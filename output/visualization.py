import pandas as pd
import matplotlib.pyplot as plt

df1 = pd.read_csv('output/F_monoterpene.dat')
df2 = pd.read_csv('output/F_isoprene.dat')
df3 = pd.read_csv('output/C_L.dat')
df4 = pd.read_csv('output/C_T.dat')

plt.plot(df1[70:])
plt.plot(df2[70:])
plt.grid()
plt.legend(["monoterpene", "isoprene"])
# plt.plot(df3)
# plt.plot(df4)
plt.show()