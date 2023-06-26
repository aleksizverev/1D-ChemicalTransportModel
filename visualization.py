import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Emissions -----------------------------------------------------
# def_emissions_isoprene = pd.read_csv('output/F_isoprene.dat')
# def_emissions_monoterpene = pd.read_csv('output/F_monoterpene.dat')
# plt.plot(def_emissions_isoprene)
# plt.plot(def_emissions_monoterpene)
#----------------------------------------------------------------


# plt.figure(1)
# df_iso = pd.read_csv('output/isoprene.dat', delimiter='\s+')
# ax1 = plt.subplot(211)
# ax1.plot(df_iso.iloc[:, 1][70:])
# ax1.plot(df_iso.iloc[:, 5][70:])
# ax1.set_title("Isoprene")

# df_alpha = pd.read_csv('output/alphapinene.dat', delimiter='\s+')
# ax2 = plt.subplot(212)
# ax2.plot(df_alpha.iloc[:, 1][70:])
# ax2.plot(df_alpha.iloc[:, 5][70:])
# ax2.set_title("Alphapinene")


# # Chemistry -----------------------------------------------------
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

df_time = pd.read_csv('output/time.dat')

# Aerosol --------------------------------------------------------
plt.figure(3)
df_PN = pd.read_csv('output/PN.dat', delimiter='\s+')
ax7 = plt.subplot(111)
ax7.plot(df_time[97:], df_PN.iloc[:, 1][97:])
ax7.grid(True, which="both")
ax7.minorticks_on()
ax7.set_ylabel('Particle number concentration (cm^-3)')
ax7.set_xlabel("Time (days)")

plt.figure(4)
df_PV = pd.read_csv('output/PV.dat', delimiter='\s+')
ax8 = plt.subplot(111)
ax8.plot(df_time[90:], df_PV.iloc[:, 1][90:])
ax8.grid(True, which="both")
ax8.minorticks_on()
ax8.set_ylabel('Particle volume concentration (μm^3/cm^-3))')
ax8.set_xlabel("Time (days)")

# df_time = pd.read_csv('output/time.dat')
# df_particle_conc = pd.read_csv('output/particle_conc.dat', delimiter='\s+')
# df_diam = pd.read_csv('output/diameter.dat')
#
#
# # size distribution plot
# fig10 = plt.figure(6)
# ax10 = plt.subplot(111)
# X,Y = np.meshgrid(df_time[70:], df_diam)
#
# df_conc_cutted = df_particle_conc.iloc[71:, :98]*1e-6 # [# cm-3], number concentration
# df_diam = df_diam[:98]
#
# # print(np.log10(df_diam.values[1]) - np.log10(df_diam.values[0]))
# # print(np.log10(df_diam.values[2]) - np.log10(df_diam.values[1]))
# # print(np.log10(df_diam.values[3]) - np.log10(df_diam.values[2]))
# # # dlogDp = []
# # # for i in range(len(df_diam)-2):
# # #     dlogDp.append(np.log10(df_diam.values[i+1]) - np.log10(df_diam.values[i]))
# # # print(dlogDp)
#
# df_dn_dlogdp = df_conc_cutted/((np.log10(df_diam.values.max()*100) - np.log10(df_diam.values.min()*100)))
#
# p = ax10.pcolor(X, Y, df_dn_dlogdp.T)
# cb = fig10.colorbar(p, ax=ax10)
# ax10.set_yscale('log')
#
#
# #particle number size distribution evolution
# plt.figure(5)
# ax9 = plt.subplot(111)
# ax9.plot(df_diam[:50], df_dn_dlogdp.iloc[:, :50].T)
# ax9.set_xscale('log')
# ax9.set_ylabel('dN/dLogDp (cm^-3)')
# ax9.grid(True, which="both")


plt.show()













# m, n = df_particle_conc.iloc[71:, 0:98].shape
# for i in range(1, m):
#     for j in range(1, n):
#         df_dn_dp.iloc[i, j] /= np.log10(df_diam.iloc[0, j])

# print(df_diam.iloc[0, 0])