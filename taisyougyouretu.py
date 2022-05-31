import numpy as np

def get_symmetric(a, use_tril=True):
    if use_tril:
        a = np.tril(a)
    else:
        a = np.triu(a)
    return a + a.T - np.diag(a.diagonal())

row=4 #対称行列の行数
low=2 #対称行列の要素の最小値
high=10000 #対称行列の要素の最大値
filename='../data/f.csv' #出力先のファイル名
a=np.random.randint(low,high,(row,row))

for i in range(row):
    a[(i,i)]=1

Q=get_symmetric(a)



print(Q)
print(len(Q))

np.savetxt(filename, Q, delimiter=',')