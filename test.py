import pandas as pd 


data1 = {"a":[2.,2.,2.,2.],
         "b":[4.,4.,4.,4.],
         "c":[8.,8.,8.,8.]}
data2 = {"a":[4.],
         "b":[4.],
         "c":[4.]}

df1 = pd.DataFrame(data1)
df2 = pd.DataFrame(data2) 

print(df1)
print()
print(df2.values)

# print(df1.T.div(df2.iloc[0], axis='columns'))

print(df1/df2.values)
print(df2.shape, df1.shape)
