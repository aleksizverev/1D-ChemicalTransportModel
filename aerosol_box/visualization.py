import pandas as pd
import matplotlib.pyplot as plt

#Emissions ------------------------------------------------------
# def_emissions_isoprene = pd.read_csv('output/F_isoprene.dat')
# def_emissions_monoterpene = pd.read_csv('output/F_monoterpene.dat')
# plt.plot(def_emissions_isoprene)
# plt.plot(def_emissions_monoterpene)
#----------------------------------------------------------------

df_time = pd.read_csv('output/time.dat')

plt.figure(1)
df_PN = pd.read_csv('output/PN.dat')
ax1 = plt.subplot(211)
ax1.plot(df_time, df_PN)
ax1.set_ylabel('Total PN (cm^-3)')


df_PV = pd.read_csv('output/PV.dat', delimiter='\s+')
ax2 = plt.subplot(212)
ax2.plot(df_time, df_PV)
ax2.set_ylabel('PV (mu*m^3/cm^3)')


# plt.figure(2)
# df_OH = pd.read_csv('output/OH.dat', delimiter='\s+')
# ax3 = plt.subplot(221)
# ax3.plot(df_OH.iloc[:, 1][70:])
# ax3.plot(df_OH.iloc[:, 5][70:])
# ax3.set_title("OH")

# df_HO2 = pd.read_csv('output/HO2.dat', delimiter='\s+')
# ax4 = plt.subplot(222)
# ax4.plot(df_HO2.iloc[:, 1][70:])
# ax4.plot(df_HO2.iloc[:, 5][70:])
# ax4.set_title("HO2")

# df_H2SO4 = pd.read_csv('output/H2SO4.dat', delimiter='\s+')
# ax5 = plt.subplot(223)
# ax5.plot(df_H2SO4.iloc[:, 1][70:])
# ax5.plot(df_H2SO4.iloc[:, 5][70:])
# ax5.set_title("H2SO4")

# df_elvoc = pd.read_csv('output/ELVOC.dat', delimiter='\s+')
# ax6 = plt.subplot(224)
# ax6.plot(df_elvoc.iloc[:, 1][70:])
# ax6.plot(df_elvoc.iloc[:, 5][70:])
# ax6.set_title("ELVOC")


plt.show()