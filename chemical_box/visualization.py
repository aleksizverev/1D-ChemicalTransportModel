import pandas as pd
import matplotlib.pyplot as plt

df1 = pd.read_csv('output/concentrations.dat', delimiter='\s+')
df2 = pd.read_csv('output/concentrations_corr.dat', delimiter='\s+')
# df2 = pd.read_csv('output/F_isoprene.dat')
# df3 = pd.read_csv('output/C_L.dat')
# df4 = pd.read_csv('output/C_T.dat')

# print(df1.iloc[:, 0])
for i in range(25):
    plt.plot(df1.iloc[:, i])
    plt.plot(df2.iloc[:, i])
    plt.show()