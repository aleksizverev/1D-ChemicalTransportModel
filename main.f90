!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Main program
!
! - Simulate emissions and chemical reactions of gases, aerosol processes as well as 
!   transport of gases and aerosol particles within the planetary boundary layer with a
!   column model.
! - Check Fortran conventions at http://www.fortran90.org/src/best-practices.html
! - Check code conventions at
!   http://www.cesm.ucar.edu/working_groups/Software/dev_guide/dev_guide/node7.html
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

program main

use aerosol_mod
use chemistry_mod
use meteorology_mod


implicit none

!-----------------------------------------------------------------------------------------
! Control variables (can be moved to an input file in future)
!-----------------------------------------------------------------------------------------
logical :: use_emission   = .true.
logical :: use_chemistry  = .true.
logical :: use_deposition = .false.
logical :: use_aerosol    = .true.
character(len=255), parameter :: input_dir  = './input'
character(len=255), parameter :: output_dir = './output'

!-----------------------------------------------------------------------------------------
! Constants
!-----------------------------------------------------------------------------------------
! Double precision
! http://en.wikipedia.org/wiki/Double_precision_floating-point_format
integer, parameter :: dp = selected_real_kind(15, 307)

! Physics constants
real(dp), parameter :: PI     = 2*asin(1.0_dp)                  ! the constant pi
real(dp), parameter :: grav   = 9.81_dp                         ! [m s-2], gravitation
real(dp), parameter :: Rgas   = 8.3144598_dp                    ! [J mol-1 K-1], universal gas constant
real(dp), parameter :: NA     = 6.022140857e23_dp               ! [molec mol-1], Avogadro's number 
real(dp), parameter :: mm_air = 28.96e-3_dp                     ! [kg mol-1], mean molar mass of air
real(dp), parameter :: kb     = 1.38064852e-23_dp               ! [m2 kg s-2 K-1], Boltzmann constant
real(dp), parameter :: Cp     = 1012.0_dp                       ! [J kg-1 K-1], air specific heat at constant pressure,
real(dp), parameter :: p00    = 1.01325e5_dp                    ! [Pa], reference pressure at surface
real(dp), parameter :: nu_air = 1.59e-5_dp                      ! [m2 s-1], kinematic viscosity of air
real(dp), parameter :: Omega  = 2*PI/(24.0_dp*60.0_dp*60.0_dp)  ! [rad s-1], Earth angular speed
real(dp), parameter :: lambda = 300.0_dp                        ! maximum mixing length, meters
real(dp), parameter :: vonk   = 0.4_dp                          ! von Karman constant, dimensionless
real(dp), parameter :: ppb = 1e-9_dp

real(dp), parameter :: ug = 10.0d0, vg = 0.0d0  ! [m s-1], geostrophic wind

! Latitude and longitude of Hyytiala
real(dp), parameter :: latitude_deg  = 61.8455d0  ! [degN]
real(dp), parameter :: longitude_deg = 24.2833d0  ! [degE]
real(dp), parameter :: latitude      = latitude_deg  * PI/180.0d0  ! [rad]
real(dp), parameter :: longitude     = longitude_deg * PI/180.0d0  ! [rad]

real(dp), parameter :: fcor = 2*Omega*sin(latitude)  ! Coriolis parameter at Hyytiala

!-----------------------------------------------------------------------------------------
! Grid parameters
!-----------------------------------------------------------------------------------------
integer, parameter :: nz = 50  ! [-], number of height levels

! Model height levels, [m]
real(dp), parameter, dimension(nz) :: &
  hh = (/    0,   10,   20,   30,   40,   50,   60,   70,   80,   90, &
           100,  120,  140,  160,  180,  200,  230,  260,  300,  350, &
           400,  450,  500,  550,  600,  650,  700,  800,  900, 1000, &
          1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, &
          2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000 /)

real(dp), parameter :: hc = 10.0_dp  ! [m], canopy height

! Chemistry constants
real(dp), dimension(neq, nz) :: cons, cons_tmp !//TODO

! Atmospheric oxygen, N2 and H2O are kept constant:
real(dp), dimension(nz) :: Mair                                    ! Air molecules concentration [molecules/cm3]
real(dp), dimension(nz) :: O2                                      ! Oxygen concentration [molecules/cm3]
real(dp), dimension(nz) :: N2                                      ! Nitrogen concentration [molecules/cm3]
real(dp), dimension(nz) :: H2O                                     ! Water molecules [molecules/cm3]

!-----------------------------------------------------------------------------------------
! Time variables
!-----------------------------------------------------------------------------------------
integer, parameter :: one_hour = 60*60  ! [s], one hour in seconds

real(dp) :: time                  ! [s], current time
real(dp) :: time_start, time_end  ! [s], start and end time

real(dp) :: dt         ! [s], time step for main loop, usually is equal to meteorology time step
real(dp) :: dt_emis, em_t    ! [s], time step for emission calculation
real(dp) :: dt_chem    ! [s], time step for chemistry calculation
real(dp) :: dt_depo    ! [s], time step for deposition calculation
real(dp) :: dt_aero    ! [s], time step for aerosol calculation
real(dp) :: dt_output  ! [s], time step for output

real(dp) :: time_start_emission    ! [s], time to start calculating emission
real(dp) :: time_start_chemistry   ! [s], time to start calculating chemistry
real(dp) :: time_start_deposition  ! [s], time to start calculating deposition
real(dp) :: time_start_aerosol     ! [s], time to start calculating aerosol

integer :: daynumber_start  ! [day], start day of year
integer :: daynumber        ! [day], current day of year

integer :: counter  ! [-], counter of time steps

!-----------------------------------------------------------------------------------------
! Meteorology variables
!-----------------------------------------------------------------------------------------
real(dp), dimension(nz  ) :: uwind, &  ! [m s-1], u component of wind
                             vwind, &  ! [m s-1], v component of wind
                             theta     ! [K], potential temperature
real(dp), dimension(nz  ) :: temp, &   ! [K], air temperature
                             pres      ! [Pa], air pressure
real(dp), dimension(nz-1  ) :: K, K_m, K_h

integer :: i, j, layer !  used for loops


!-----------------------------------------------------------------------------------------
! Emission variables
!-----------------------------------------------------------------------------------------
real(dp) :: F_monoterpene, F_isoprene, C_L, C_T

!-----------------------------------------------------------------------------------------
! Aerosol variables
!-----------------------------------------------------------------------------------------

real(dp), dimension(nz) :: CS_H2SO4, CS_ELVOC
real(dp), dimension (nz, nr_cond) :: cond_vapour          ! Concentration of condensable vapours [molec/m^3]
REAL(dp), DIMENSION(nz) :: PN, PM, PV  ! Total particle number [# m-3] and mass concentration [kg m-3]
REAL(dp), DIMENSION(nz, nr_bins) :: particle_conc   ! number concentration in each size bin

!-----------------------------------------------------------------------------------------
! Initialization
!-----------------------------------------------------------------------------------------

call time_init()                 ! initialize time
call meteorology_init()          ! initialize meteorology

do i = 1, nz
  call Aerosol_init(diameter, particle_mass, particle_volume, particle_conc(i, :), &
                    particle_density, nucleation_coef, molecular_mass, molar_mass, &
                    molecular_volume, molecular_dia, mass_accomm)
end do
call open_files()        ! open output files
call write_files(time)   ! write initial values

cons = 0.0d0
cons_tmp = 0.0d0
CS_H2SO4 = 1d-3
CS_ELVOC = 1d-3

!-----------------------------------------------------------------------------------------
! Start main loop
!-----------------------------------------------------------------------------------------
do while (time <= time_end)
  !---------------------------------------------------------------------------------------
  ! Meteorology
  !---------------------------------------------------------------------------------------
  ! Set lower boundary condition
  call surface_values(theta(1), time+dt)  ! theta = temperature at the surface

  ! Update meteorology
  ! call updateParametersK1(uwind, vwind, theta, hh, dt, ug, vg, fcor)
  ! call updateParametersK2(uwind, vwind, theta, hh, dt, ug, vg, fcor, lambda, vonk, K)
  call updateParametersK3(uwind, vwind, theta, hh, dt, ug, vg, fcor, lambda, vonk, K_m, K_h)

  !---------------------------------------------------------------------------------------
  ! Emission
  !---------------------------------------------------------------------------------------
  ! Start to calculate emission after time_start_emission
  ! Compute emission part every dt_emis, multiplying 1000 to convert s to ms to make mod easier
  temp = theta - (grav/Cp)*hh
  if ( use_emission .and. time >= time_start_emission ) then
    if ( mod( nint((time - time_start_emission)*1000.0d0), nint(dt_emis*1000.0d0) ) == 0 ) then
      call calculate_emmisions(F_monoterpene, F_isoprene, temp(2), C_L, C_T, time)
    end if
  end if

  if ( use_emission .and. (.not. use_chemistry) ) then
    ! Add emission to the number concentrations of compounds
  end if

  !---------------------------------------------------------------------------------------
  ! Deposition
  !---------------------------------------------------------------------------------------
  ! Start to calculate gas dry deposition velocity after time_start_deposition
  ! Compute deposition part every dt_depo, multiplying 1000 to convert s to ms to make mod easier
  if ( use_deposition .and. time >= time_start_deposition ) then
    if ( mod( nint((time - time_start_deposition)*1000.0d0), nint(dt_depo*1000.0d0) ) == 0 ) then
      ! Calculate deposition velocity

      ! Remove deposited concentration at level 2 which includes canopy and soil
    end if
  end if

  !---------------------------------------------------------------------------------------
  ! Chemistry
  !---------------------------------------------------------------------------------------
  ! Start to calculate chemical reactions only after some time to save the computation time
  ! Compute chemistry part every dt_chem, multiplying 1000 to convert s to ms to make mod easier
  if ( use_chemistry .and. time >= time_start_chemistry ) then
    if ( mod( nint((time - time_start_chemistry)*1000.0d0), nint(dt_chem*1000.0d0) ) == 0 ) then

      Mair = pres*NA / (Rgas*temp) * 1d-6 
      O2   = 0.21d0*Mair                     ! Oxygen concentration [molecules/cm3]
      N2   = 0.78d0*Mair                     ! Nitrogen concentration [molecules/cm3]
      H2O  = 1.0D16                          ! Water molecules [molecules/cm3]
      ! Initial state for concentrations
      cons(1, :)  = 24.0d0   * Mair * ppb     ! O3 concentration
      cons(5, :)  = 0.2d0    * Mair * ppb     ! NO2
      cons(6, :)  = 0.07d0   * Mair * ppb     ! NO
      cons(9, :)  = 100.0d0  * Mair * ppb     ! CO
      cons(11, :) = 1759.0d0 * Mair * ppb     ! CH4
      cons(20, :) = 0.5d0    * Mair * ppb     ! SO2

      ! Solve chemical equations for each layer except boundaries //TODO 
      do layer = 2, nz - 1

        if (layer .GT. 2) then 
          F_isoprene = 0.0_dp
          F_monoterpene = 0.0_dp
        end if

        do i = 1,nz 
          do j = 1,neq
             if (cons(j,i) .LT. 0.0) then
                cons(j,i) = 0.0_dp
             endif
            enddo
        enddo

        call chemistry_step(cons(1:neq, layer), time, time+dt_chem, O2(layer), N2(layer), Mair(layer), &
                            H2O(layer), temp(layer), get_exp_coszen(time, daynumber, latitude), &
                            F_monoterpene , F_isoprene, CS_H2SO4(layer), CS_ELVOC(layer))

      end do
    end if  ! every dt_chem
  end if

  ! Update concentrations of gas phase compounds if any of these processes are considered
  ! Deposition should not be used alone because it calculates nothing in that case
  if (use_emission .or. use_chemistry) then
    ! Trick to make bottom flux zero
    cons(1:neq,1) = cons(1:neq,2)

    ! Mixing of chemical species
    do i = 2, nz-1
      do j = 1, neq
      cons_tmp(j, i) = cons(j, i)  +   dt*(K_h(i) * (cons(j, i+1) - cons(j, i))/(hh(i+1) - hh(i)) - & 
                                             K_h(i-1) * (cons(j, i) - cons(j, i-1))/(hh(i) - hh(i-1))) / & 
                                            ((hh(i+1)-hh(i-1))*0.5)
      end do
    end do 
    cons = cons_tmp

    ! Set the constraints above again for output
  end if

  !---------------------------------------------------------------------------------------
  ! Aerosol
  !---------------------------------------------------------------------------------------
  ! Start to calculate aerosol processes only after some time to save the computation time
  ! Compute aerosol part every dt_aero, multiplying 1000 to convert s to ms to make mod easier
  if ( use_aerosol .and. time >= time_start_aerosol ) then
    if ( mod( nint((time - time_start_aerosol)*1000.0d0), nint(dt_aero*1000.0d0) ) == 0 ) then
      ! Nucleation, condensation, coagulation and deposition of particles
      
      do layer = 2, nz-1
        call Nucleation(particle_conc(layer, :), particle_volume, nucleation_coef, cond_vapour(layer, :), dt_aero)

        call Coagulation(dt_aero, particle_conc(layer, :), diameter, temp(layer), pres(layer), particle_mass)

        call Condensation(dt_aero, temp(layer), pres(layer), mass_accomm, molecular_mass, &
                        molecular_volume, molar_mass, molecular_dia, particle_mass, particle_volume, &
                        particle_conc(layer, :), diameter, cond_vapour(layer, :), CS_H2SO4(layer), CS_ELVOC(layer))


        PN(layer) = sum(particle_conc(layer, :))*1D-6                 ! [# cm-3], total particle number concentration
        PM(layer) = sum(particle_conc(layer, :)*particle_mass)*1D9     ! [ug m-3], total particle mass concentration
        PV(layer) = sum(particle_conc(layer, :)*particle_volume)*1D9*1000
      end do
    end if

    ! //TODO DEFINE CS, COND_VAPOUR IN MAIN AND SPECIFY THEM FOR DIFFERENT LAYERS

    ! Trick to make bottom flux zero

    ! Concentrations can not be lower than 0 [molec m-3]

    ! Mixing of aerosol particles

    ! Set the constraints above again for output

    ! Update related values, e.g., total number concentration, total mass concentration

  end if

  !---------------------------------------------------------------------------------------
  ! Ending loop actions
  !---------------------------------------------------------------------------------------
  ! Advance to next time step
  time = time + dt

  ! Write data every dt_output [s]
  if ( mod( nint((time - time_start)*1000.0d0), nint(dt_output*1000.0d0) ) == 0 ) then
    write(*, '(a8,f8.3,a8)') 'time = ', time/one_hour, '   hours'
    call write_files(time)
  end if

  ! Count loop number
  counter = counter + 1

end do

!-----------------------------------------------------------------------------------------
! Finalization
!-----------------------------------------------------------------------------------------
! Close all the opened files
call close_files()

! Count total time steps
write(*,*) counter,'time steps'


contains


!-----------------------------------------------------------------------------------------
! subroutine open_files()
!
! Open needed files
!-----------------------------------------------------------------------------------------
subroutine open_files()
  logical :: dir_exist

  ! Create a new directory if it does not exist
  inquire(file=trim(adjustl(output_dir)), exist=dir_exist)
  if (.not. dir_exist) then
    ! This line may change for different operating systems
    call system('mkdir ' // trim(adjustl(output_dir)))
  end if

  ! Open files to write output results
  open( 8,file=trim(adjustl(output_dir))//'/time.dat' ,status='replace',action='write')
  open( 9,file=trim(adjustl(output_dir))//'/hh.dat'   ,status='replace',action='write')
  open(10,file=trim(adjustl(output_dir))//'/uwind.dat',status='replace',action='write')
  open(11,file=trim(adjustl(output_dir))//'/vwind.dat',status='replace',action='write')
  open(12,file=trim(adjustl(output_dir))//'/theta.dat',status='replace',action='write')
  open(13,file=trim(adjustl(output_dir))//'/K.dat',status='replace',action='write')
  open(14,file=trim(adjustl(output_dir))//'/F_isoprene.dat',status='replace',action='write')
  open(15,file=trim(adjustl(output_dir))//'/F_monoterpene.dat',status='replace',action='write')
  open(16,file=trim(adjustl(output_dir))//'/C_L.dat',status='replace',action='write')
  open(17,file=trim(adjustl(output_dir))//'/C_T.dat',status='replace',action='write')
  open(18,file=trim(adjustl(output_dir))//'/OH.dat',status='replace',action='write')
  open(19,file=trim(adjustl(output_dir))//'/isoprene.dat',status='replace',action='write')
  open(20,file=trim(adjustl(output_dir))//'/alphapinene.dat',status='replace',action='write')
  open(21,file=trim(adjustl(output_dir))//'/HO2.dat',status='replace',action='write')
  open(22,file=trim(adjustl(output_dir))//'/H2SO4.dat',status='replace',action='write')
  open(23,file=trim(adjustl(output_dir))//'/ELVOC.dat',status='replace',action='write')
  open(24,file=trim(adjustl(output_dir))//'/PM.dat',status='replace',action='write')
  open(25,file=trim(adjustl(output_dir))//'/PN.dat',status='replace',action='write')
  open(26,file=trim(adjustl(output_dir))//'/PV.dat',status='replace',action='write')
end subroutine open_files


!-----------------------------------------------------------------------------------------
! subroutine write_files(time)
!
! Write data to files at time
!-----------------------------------------------------------------------------------------
subroutine write_files(time)
  real(dp) :: time  ! current time
  character(255) :: outfmt_one_scalar, outfmt_two_scalar, outfmt_level, outfmt_mid_level

  ! Output real data with scientific notation with 16 decimal digits
  outfmt_one_scalar = '(es25.16e3)'                               ! for scalar
  write(outfmt_level     , '(a, i3, a)') '(', nz  , 'es25.16e3)'  ! for original levels
  write(outfmt_mid_level , '(a, i3, a)') '(', nz-1, 'es25.16e3)'  ! for middle levels
  write(outfmt_two_scalar, '(a, i3, a)') '(', 2   , 'es25.16e3)'  ! for two scalars

  ! Only output hh once
  if (time == time_start) then
    write(9, outfmt_level) hh
  end if

  ! Output every output time step
  write( 8, outfmt_one_scalar) time/(24*one_hour)  ! [day]
  write(10, outfmt_level     ) uwind
  write(11, outfmt_level     ) vwind
  write(12, outfmt_level     ) theta
  write(13, outfmt_level     ) K
  write(14, outfmt_one_scalar) F_isoprene
  write(15, outfmt_one_scalar) F_monoterpene
  write(16, outfmt_one_scalar) C_L
  write(17, outfmt_one_scalar) C_T
  write(18, outfmt_level     ) cons(3, :)
  write(19, outfmt_level     ) cons(13, :)
  write(20, outfmt_level     ) cons(23, :)
  write(21, outfmt_level     ) cons(8, :)
  write(22, outfmt_level     ) cons(21, :)
  write(23, outfmt_level     ) cons(25, :)
  write(24, outfmt_level     ) PM
  write(25, outfmt_level     ) PN
  write(26, outfmt_level     ) PV
end subroutine write_files


!-----------------------------------------------------------------------------------------
! subroutine Close_Files()
!
! Close files
!-----------------------------------------------------------------------------------------
subroutine close_files()
  close(8)
  close(9)
  close(10)
  close(11)
  close(12)
  close(13)
  close(14)
  close(15)
  close(16)
  close(17)
  close(18)
  close(19)
  close(20)
  close(21)
  close(22)
  close(23)
end subroutine close_files


!-----------------------------------------------------------------------------------------
! subroutine time_init()
!
! Time initiation
!-----------------------------------------------------------------------------------------
subroutine time_init()
  ! Basic time variables
  time_start = 0.0d0
  time_end   = 5.0d0 * 24.0d0 * one_hour
  time       = time_start

  ! Time steps
  dt        = 0.5d0
  dt_emis   = 0.5d0
  dt_chem   = 10.0d0
  dt_depo   = 10.0d0
  dt_aero   = 10.0d0
  dt_output = 3600.0d0

  ! Day number
  daynumber_start = 31+28+31+30+31+30+10  ! day is July 10th, 2011
  daynumber       = daynumber_start

  ! Start time for each process
  time_start_emission   = 3*24*one_hour
  time_start_chemistry  = 3*24*one_hour
  time_start_deposition = 3*24*one_hour
  time_start_aerosol    = 3*24*one_hour

  ! Loop number
  counter = 0
end subroutine time_init


!-----------------------------------------------------------------------------------------
! subroutine meteorology_init()
!
! Meteorology initiation
!-----------------------------------------------------------------------------------------
subroutine meteorology_init()
  ! Wind velocity
  uwind         = 0.0d0
  uwind(nz)     = ug
  uwind(2:nz-1) = uwind(nz) * hh(2:nz-1)/hh(nz)

  vwind         = 0.0d0
  vwind(nz)     = vg

  K = 0.0d0
  K_m = 0.0d0
  K_h = 0.0d0

  ! Potential temperature
  theta     = 273.15d0 + 25.0d0
  theta(nz) = 273.15d0 + 30.0d0

  ! Air temperature and pressure
  temp = theta - (grav/Cp)*hh
  pres = barometric_law(p00, temp, hh)
end subroutine meteorology_init


!-----------------------------------------------------------------------------------------
! Get the surface values from the input data file
! Now only the temperature is used.
!-----------------------------------------------------------------------------------------
subroutine surface_values(temperature, time)

  ! (Note: can also get water concentrantion, in ppt, if modify this
  ! subroutine to also use column 8)
  !
  ! Data is taken from:
  ! http://avaa.tdata.fi/web/smart

  real(dp), intent(in)            :: time ! input, in seconds
  real(dp), intent(out)           :: temperature ! output, in Kelvin
  logical, save                   :: first_time = .true.
  real(dp), dimension(8,50), save :: surface_data
  real(dp), dimension(50), save   :: temperature_data
  real(dp), parameter             :: seconds_in_day = 24*60*60
  real(dp), parameter             :: seconds_in_30min = 30*60
  integer                         :: index
  real(dp) :: time24h, time30min, time24plus15, temp1, temp2, x

  ! Only when called for the first time, read in data from file
  ! With this trick, we don't need to open the file in the main program
  IF (first_time) THEN
     open(30, file=trim(adjustl(input_dir))//'/hyytiala_20110710-t_h2o.dat', status='old')
     read(30, *) surface_data
     temperature_data(1:50) = surface_data(7,1:50) ! in Celcius
     first_time = .false.
  end IF

  time24h = modulo(time, seconds_in_day) ! time modulo 24 hours
  time24plus15 = time24h + 15*60 ! time since 23:45 previous day
  time30min = modulo(time24plus15, seconds_in_30min)
  index = 1 + floor(time24plus15/seconds_in_30min)

  temp1 = temperature_data(index)
  temp2 = temperature_data(index + 1)
  x = time30min/seconds_in_30min

  ! linear interpolation between previous and next temperature data value
  temperature = temp1 + x*(temp2 - temp1) + 273.15_dp  ! now in Kelvin

end subroutine surface_values


subroutine calculate_emmisions(F_monoterpene, F_isoprene, T, C_L, C_T, em_t)
  real(dp) :: F_monoterpene, F_isoprene, T, C_L, C_T, Q, exp_coszen, em_t
  real(dp), parameter :: betha = 0.09_dp
  real(dp), parameter :: T_s = 303.15_dp                  ! [Kelvin]
  real(dp), parameter :: T_m = 314                        ! [Kelvin]
  real(dp), parameter :: D_m = 0.0538_dp                  ! Foliar density [g/cm^2]
  real(dp), parameter :: em_factor = 100.0_dp             ! Standart emission protential [ng/g/h]
  real(dp), parameter :: alpha = 0.0027                
  real(dp), parameter :: c_l1 = 1.006_dp                  
  real(dp), parameter :: c_t1 = 1000*95.0_dp              ! [J/mol]
  real(dp), parameter :: c_t2 = 1000*230.0_dp             ! [J/mol]


  ! Calculating isoprene flux
  exp_coszen = get_exp_coszen(em_t, daynumber, latitude)
  Q = 1000.0 * exp_coszen
  C_L = alpha * c_l1 * Q/sqrt(1 + alpha**(2) * Q**(2))
  C_T = exp( (c_t1*(T-T_s))/(Rgas*T*T_s) ) /(1 + exp( (c_t2*(T-T_m))/(Rgas*T*T_s) ))

  F_isoprene = D_m * em_factor * C_L * C_T * NA/68 * 1e-9 /3600 /1000
  
  ! Calculation monoterpene flux
  F_monoterpene = D_m * em_factor * exp(betha*(T-T_s))  * NA/136 * 1e-9 /3600 /1000

end subroutine calculate_emmisions

!-----------------------------------------------------------------------------------------
! Calculate the radiation related quantities
!-----------------------------------------------------------------------------------------
real(dp) function get_exp_coszen(time,daynumber,latitude)
  real(dp), intent(in) :: time,latitude
  INTEGER, intent(in) :: daynumber
  real(dp) :: hourangle,zenith,coszen
  hourangle = get_hourangle(time)
  zenith = solar_zenith_angle(hourangle,daynumber,latitude)
  coszen = cos(zenith)
  IF (coszen > 0) THEN  ! sun is above horizon
     get_exp_coszen = exp(-0.575_dp/coszen)
  ELSE
     get_exp_coszen = 0.0_dp
  endIF
end function get_exp_coszen


real(dp) function get_hourangle(time)
  real(dp), intent(in) :: time
  real(dp), parameter :: one_day = 24*one_hour
  get_hourangle = modulo(time,one_day)/one_day * 2 * pi - pi
end function get_hourangle


real(dp) function solar_zenith_angle(hourangle,daynumber,latitude)
  ! http://en.wikipedia.org/wiki/Solar_elevation_angle
  ! http://en.wikipedia.org/wiki/Position_of_the_Sun
  INTEGER, intent(in) :: daynumber
  real(dp), intent(in) :: hourangle,latitude
  real(dp) :: declination,elevation
  real(dp), parameter :: to_rad = pi/180.0_dp

  declination = -23.44_dp * to_rad * cos(2 * pi * (daynumber + 10)/365.0_dp)
  elevation = cos(hourangle)*cos(declination)*cos(latitude) &
       + sin(declination)*sin(latitude)
  solar_zenith_angle = pi/2.0_dp - elevation
  ! Notes:
  ! - Not tested near equador or on the southern hemisphere.
  ! - solar_zenith_angle can be larger than pi/2, it just means
  !   the sun is below horizon.
  ! - solar_zenith_angle assumes time is in local solar time, which
  !   is usually not exactly true
end function solar_zenith_angle


!-----------------------------------------------------------------------------------------
! Other functions
!-----------------------------------------------------------------------------------------
function barometric_law(p00, tempK, h) result(p)
  real(dp), intent(in) :: p00, tempK(nz), h(nz)
  real(dp) :: p(nz)
  real(dp) :: dh(nz)

  dh(2:nz) = h(2:nz) - h(1:nz-1)

  p(1) = p00
  do i=2, nz
    p(i) = p(i-1)*exp(-mm_air*grav/(Rgas*(tempK(i-1)+tempK(i))/2.0d0)*dh(i))
  end do
end function barometric_law

end program main
