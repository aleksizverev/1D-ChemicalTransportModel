import pandas as pd
import matplotlib.pyplot as plt


df5 = pd.read_csv('output/isoprene.dat', delimiter='\s+')


# print(df5.iloc[:, 1][75:])
plt.plot(df5.iloc[:, 1][75:])
plt.plot(df5.iloc[:, 5][75:])
plt.show()