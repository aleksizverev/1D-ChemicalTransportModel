module meteorology_mod

private :: dp, nz
integer, parameter :: dp = selected_real_kind(15, 307)
integer, parameter :: nz = 50

contains 

    subroutine updateParametersK1(uwind, vwind, theta, hh, dt, ug, vg, fcor)
        
        real(dp), dimension(nz) :: hh
        real(dp), dimension(nz) :: uwind, uwind_tmp, vwind, vwind_tmp, theta, theta_tmp
        real(dp) :: dt 
        real(dp) :: ug, vg, fcor

        uwind_tmp = uwind
        vwind_tmp = vwind
        theta_tmp = theta
 
        do i = 2, nz-1
            uwind_tmp(i) = uwind(i) + dt*(fcor*(vwind(i)-vg) + (5 * (uwind(i+1) - uwind(i))/(hh(i+1) - hh(i)) - & 
                                                                5 * (uwind(i) - uwind(i-1))/(hh(i) - hh(i-1))) / & 
                                                                ((hh(i+1)-hh(i-1))*0.5))
        
            vwind_tmp(i) = vwind(i) + dt*(-fcor*(uwind(i)-ug) + (5 * (vwind(i+1) - vwind(i))/(hh(i+1) - hh(i)) - & 
                                                                 5 * (vwind(i) - vwind(i-1))/(hh(i) - hh(i-1))) / & 
                                                                ((hh(i+1)-hh(i-1))*0.5))
        
            theta_tmp(i) = theta(i)  +   dt*(5 * (theta(i+1) - theta(i))/(hh(i+1) - hh(i)) - & 
                                             5 * (theta(i) - theta(i-1))/(hh(i) - hh(i-1))) / & 
                                            ((hh(i+1)-hh(i-1))*0.5)
        end do

        uwind = uwind_tmp
        vwind = vwind_tmp
        theta = theta_tmp

    end subroutine updateParametersK1


    subroutine updateParametersK2(uwind, vwind, theta, hh, dt, ug, vg, fcor, lambda, vonk, K)
        
        real(dp), dimension(nz) ::  hh
        real(dp), dimension(nz) ::  uwind, uwind_tmp, vwind, vwind_tmp, &
                                    theta, theta_tmp
                                    
        real(dp), dimension(nz-1) :: K                           
        real(dp) :: dt 
        real(dp) :: ug, vg, fcor
        real(dp) :: lambda, vonk, blackadar 

        uwind_tmp = uwind
        vwind_tmp = vwind
        theta_tmp = theta

        do i = 1, nz-1
            blackadar = (vonk*(hh(i)+hh(i+1))*0.5)/(1.0 + (vonk*(hh(i)+hh(i+1))*0.5)/lambda)

            K(i) = (blackadar**2.0) * sqrt( ((uwind(i+1)-uwind(i))/(hh(i+1)-hh(i)))**2.0 + &
                                            ((vwind(i+1)-vwind(i))/(hh(i+1)-hh(i)))**2.0 ) 
        end do

        do i = 2, nz-1

            uwind_tmp(i) = uwind(i) + dt*(fcor*(vwind(i)-vg) + (K(i) * (uwind(i+1) - uwind(i))/(hh(i+1) - hh(i)) - & 
                                                                K(i-1) * (uwind(i) - uwind(i-1))/(hh(i) - hh(i-1))) / & 
                                                                ((hh(i+1)-hh(i-1))*0.5))
        
            vwind_tmp(i) = vwind(i) + dt*(-fcor*(uwind(i)-ug) + (K(i) * (vwind(i+1) - vwind(i))/(hh(i+1) - hh(i)) - & 
                                                                 K(i-1) * (vwind(i) - vwind(i-1))/(hh(i) - hh(i-1))) / & 
                                                                ((hh(i+1)-hh(i-1))*0.5))
        
            theta_tmp(i) = theta(i)  +   dt*(K(i) * (theta(i+1) - theta(i))/(hh(i+1) - hh(i)) - & 
                                             K(i-1) * (theta(i) - theta(i-1))/(hh(i) - hh(i-1))) / & 
                                            ((hh(i+1)-hh(i-1))*0.5)
        end do

        uwind = uwind_tmp
        vwind = vwind_tmp
        theta = theta_tmp

    end subroutine updateParametersK2

    subroutine updateParametersK3(uwind, vwind, theta, hh, dt, ug, vg, fcor, lambda, vonk, K_m, K_h)
        
        real(dp), dimension(nz) ::  hh
        real(dp), dimension(nz) ::  uwind, uwind_tmp, vwind, vwind_tmp, &
                                    theta, theta_tmp
                                    
        real(dp), dimension(nz-1) :: K_m, K_h                      
        real(dp) :: dt 
        real(dp) :: ug, vg, fcor
        real(dp) :: lambda, vonk, blackadar, Ri, f_m, f_h

        uwind_tmp = uwind
        vwind_tmp = vwind
        theta_tmp = theta

        do i = 1, nz-1
            blackadar = (vonk*(hh(i)+hh(i+1))*0.5)/(1.0 + (vonk*(hh(i)+hh(i+1))*0.5)/lambda)

            Ri = 9.81_dp / ((theta(i+1)+theta(i))*0.5) * &
                 (theta(i+1)-theta(i))*(hh(i+1)-hh(i)) / &
                 ((uwind(i+1)-uwind(i))**2.0 + (vwind(i+1)-vwind(i))**2.0)
            
            if (Ri < 0.0) then
                f_m = (1.0-16.0*Ri)**(0.5)
                f_h = (1.0-16.0*Ri)**(3.0/4.0)  
            else if (Ri >= 0 .and. Ri < 0.2) then
                f_m = max((1-5*Ri)**2, 0.1_dp)
                f_h = max((1-5*Ri)**2, 0.1_dp)
            else
                f_m = 0.1_dp
                f_h = 0.1_dp
            endif


            K_m(i) = (blackadar**2.0) * sqrt( ((uwind(i+1)-uwind(i))/(hh(i+1)-hh(i)))**2.0 + &
                                            ((vwind(i+1)-vwind(i))/(hh(i+1)-hh(i)))**2.0 ) * f_m
                                            

            K_h(i) = (blackadar**2.0) * sqrt( ((uwind(i+1)-uwind(i))/(hh(i+1)-hh(i)))**2.0 + &
                                            ((vwind(i+1)-vwind(i))/(hh(i+1)-hh(i)))**2.0 ) * f_h
                                                                       
        end do

        do i = 2, nz-1

            uwind_tmp(i) = uwind(i) + dt*(fcor*(vwind(i)-vg) + (K_m(i) * (uwind(i+1) - uwind(i))/(hh(i+1) - hh(i)) - & 
                                                                K_m(i-1) * (uwind(i) - uwind(i-1))/(hh(i) - hh(i-1))) / & 
                                                                ((hh(i+1)-hh(i-1))*0.5))
        
            vwind_tmp(i) = vwind(i) + dt*(-fcor*(uwind(i)-ug) + (K_m(i) * (vwind(i+1) - vwind(i))/(hh(i+1) - hh(i)) - & 
                                                                 K_m(i-1) * (vwind(i) - vwind(i-1))/(hh(i) - hh(i-1))) / & 
                                                                ((hh(i+1)-hh(i-1))*0.5))
        
            theta_tmp(i) = theta(i)  +   dt*(K_h(i) * (theta(i+1) - theta(i))/(hh(i+1) - hh(i)) - & 
                                             K_h(i-1) * (theta(i) - theta(i-1))/(hh(i) - hh(i-1))) / & 
                                            ((hh(i+1)-hh(i-1))*0.5)
        end do

        uwind = uwind_tmp
        vwind = vwind_tmp
        theta = theta_tmp

    end subroutine updateParametersK3

end module meteorology_mod