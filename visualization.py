import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Emissions -----------------------------------------------------
df_time = pd.read_csv('output/time.dat')
def_emissions_isoprene = pd.read_csv('output/F_isoprene.dat')
def_emissions_monoterpene = pd.read_csv('output/F_monoterpene.dat')
plt.plot(df_time[71:], def_emissions_isoprene[71:])
plt.plot(df_time[71:], def_emissions_monoterpene[71:])
plt.title("40")
#----------------------------------------------------------------

#
# plt.figure(1)
# df_iso = pd.read_csv('output/isoprene.dat', delimiter='\s+')
# ax1 = plt.subplot(211)
# ax1.plot(df_iso.iloc[:, 1][70:])
# ax1.plot(df_iso.iloc[:, 5][70:])
# ax1.set_title("Isoprene")
#
# df_alpha = pd.read_csv('output/alphapinene.dat', delimiter='\s+')
# ax2 = plt.subplot(212)
# ax2.plot(df_alpha.iloc[:, 1][70:])
# ax2.plot(df_alpha.iloc[:, 5][70:])
# ax2.set_title("Alphapinene")


# Chemistry -----------------------------------------------------
# plt.figure(2)
# df_OH = pd.read_csv('output/OH.dat', delimiter='\s+')
# ax3 = plt.subplot(221)
# ax3.plot(df_OH.iloc[:, 1][70:])
# ax3.plot(df_OH.iloc[:, 5][70:])
# ax3.set_title("OH")
#
# df_HO2 = pd.read_csv('output/HO2.dat', delimiter='\s+')
# ax4 = plt.subplot(222)
# ax4.plot(df_HO2.iloc[:, 1][70:])
# ax4.plot(df_HO2.iloc[:, 5][70:])
# ax4.set_title("HO2")
#
# df_H2SO4 = pd.read_csv('output/H2SO4.dat', delimiter='\s+')
# ax5 = plt.subplot(223)
# ax5.plot(df_H2SO4.iloc[:, 1][70:])
# ax5.plot(df_H2SO4.iloc[:, 5][70:])
# ax5.set_title("H2SO4")
#
# df_elvoc = pd.read_csv('output/ELVOC.dat', delimiter='\s+')
# ax6 = plt.subplot(224)
# ax6.plot(df_elvoc.iloc[:, 1][70:])
# ax6.plot(df_elvoc.iloc[:, 5][70:])
# ax6.set_title("ELVOC")
# ----------------------------------------------------------------


# Aerosol --------------------------------------------------------
# plt.figure(3)
# df_PN = pd.read_csv('output/PN.dat', delimiter='\s+')
# ax7 = plt.subplot(111)
# ax7.plot(df_PN.iloc[:, 1][72:])
# ax7.grid(True, which="both")
# ax7.minorticks_on()
# ax7.set_ylabel('Total PN (cm^-3)')
#
# plt.figure(4)
# df_PM = pd.read_csv('output/PM.dat', delimiter='\s+')
# ax8 = plt.subplot(111)
# ax8.plot(df_PM.iloc[:, 1][72:])
# ax8.grid(True, which="both")
# ax8.minorticks_on()
# ax8.set_ylabel('Total PM')
#
# plt.figure(5)
# df_particle_conc = pd.read_csv('output/particle_conc.dat', delimiter='\s+')
# df_diam = pd.read_csv('output/diameter.dat') * 1e9
# ax9 = plt.subplot(111)
# ax9.loglog(df_diam, df_particle_conc.iloc[119][0:99]*1e-6)
# ax9.grid(True, which="both")
# ax9.set_ylabel('\delta N (cm^-3)')
#
#
# fig10 = plt.figure(6)
# ax10 = plt.subplot(111)
# df_time = pd.read_csv('output/time.dat')
# X,Y = np.meshgrid(df_time[70:], df_diam)
#
# p = ax10.pcolor(X, Y, np.log10((df_particle_conc.iloc[71:, 0:98]*1e-6).T), vmin=0)
# ax10.set_yscale('log')
# cb = fig10.colorbar(p, ax=ax10)
# #----------------------------------------------------------------

plt.show()
