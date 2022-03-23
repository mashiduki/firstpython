# ------ Import necessary packages ----
import networkx as nx
from collections import defaultdict
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
import matplotlib.pyplot as plt
import csv
import numpy as np
import pandas as pd
import sys


#パラメータの設定/////////////////////////////////////////////////////////////////////////////////////////////////////////////


csvdata = pd.read_csv("f.csv", header = None) #f.csvに書かれた対称行列をcsvdataに読み込む

size = len(csvdata)  #行列のサイズ
lsize = 1 #小ブロックのサイズ
bsize = 3 #大ブロックのサイズ

m = size//lsize #mは小ブロックの行数
n = size//bsize #nは大ブロックの行数

filename='f_swap.csv' #出力先のファイル名
#現在出力は小ブロック1の時のみに対応

print(csvdata)
print(size)

#ペナルティを設定
a=35000 #同じ大ブロックを選ばないための条件 細かくみていく
b=35000 #大ブロックのサイズを等しくするための条件





#パラメータここまで////////////////////////////////////////////////////////////////////////////////////////////////////////
print('行列のサイズは%d' % size )
f = np.array #fの型を定義
#Fb = np.zeros((bsize,bsize), dtype = None, order='C') いらない
f = csvdata.values #fを読み込んだ対称行列
F = [0]*m #F[i]はi行目の小ブロックの和



#F[i]を生成
for i in range(size//lsize):
    for j in range(size):
        F[i]+=f[(i, j)]

print(F)
print(f)


#制約条件のチェック用変数
smallblock_num = size/lsize
bigblock_num = size/bsize
sample_num = smallblock_num* bigblock_num
block_count = 0   
block_count2 = 0
sample_num_int = int(sample_num)
pre_block_count = -1
lbcounter = 0
bigblock_num_int = int(bigblock_num)


sum_all=0
for i in range(size):
    for j in range(size):
        sum_all=sum_all+f[i,j]*f[i,j]

sum_taikaku=0
for i in range(bigblock_num_int):
    for j in range(bsize):
        for k in range(bsize):
            sum_taikaku = sum_taikaku + f[j+(bsize*i), k+(bsize*i)]*f[j+(bsize*i), k+(bsize*i)]

print("before")
print(sum_all-sum_taikaku)
before_H=sum_all-sum_taikaku

# ------- Set up our QUBO dictionary -------
# Initialize our Q matrix
Q = defaultdict(int)

#目的関数
for i in range(0, m):
    for k in range(0, n):
        p=k*m + i
        Q[(p+1,p+1)] += 2*F[i]
        for j in range(0, m):
            q=k*m + j
            Q[(p+1,q+1)]+=-2*f[i,j]


#制約条件1
for i in range(0, m):
    for k in range(0, n):
        p=k*m + i
        for l in range(0, n):
            if(l != k):
                q=l*m + i
                Q[(p+1,q+1)] += a
        Q[(p+1,p+1)]+= -a

    

#制約条件2
for k in range(0, n):
    for i in range(0, m):
        p=k*m+i
        for j in range(0, m):
            if(j != i):
                q=k*m + j
                Q[(p+1,q+1)] += b
        Q[(p+1,p+1)]+=(1-(2*m/n))*b

print(Q)

dwavecount = 0

# ------- Run our QUBO on the QPU -------
# Set up QPU parameters
numruns = 1000 #←アニーリング回数の設定
# Run the QUBO on the solver from your config file
sampler=EmbeddingComposite(DWaveSampler()) #　　←サンプラーを設定（D-Waveの構造にマップします）
response = sampler.sample_qubo(Q, num_reads=numruns) #←Qをマップし、D-Waveで1000回アニーリング実施
for s,e,o,c in response.data(['sample', 'energy', 'num_occurrences', 'chain_break_fraction']):
    print(s, e, o, c) #←1000回のアニーリングの集計結果（頂点組み合わせと、
#そのエネルギー値、出現数及びchain break fraction値（※参照2）をプリント
    if dwavecount==0:
        tmp = s #tmpに一番エネルギーの低いサンプルの結果（辞書型）を代入
    dwavecount= dwavecount + 1
    break #breakをすると一番エネルギーの低いもののみ表示。消すと全て表示（デフォルトではbreak無し)

vsum = 0
for v in tmp.values():
    vsum += v

print('1の数の合計', vsum)

#制約条件のチェック用変数
smallblock_num = size/lsize
bigblock_num = size/bsize
sample_num = smallblock_num* bigblock_num
block_count = 0   
block_count2 = 0
sample_num_int = int(sample_num)
pre_block_count = -1
f_tmp=np.array(f)
f_swap=np.array(f_tmp)
lbcounter = 0
bigblock_num_int = int(bigblock_num)


print(sample_num_int)

for i in range(sample_num_int):       #制約条件のチェック()
    block_count += tmp[i+1]
    block_count2 += 1
    if block_count2 == smallblock_num:
        print(i+1,'までのブロックの1の数')
        print(block_count)
        if pre_block_count != -1:
            if pre_block_count != block_count:
                print('ERROR 制約条件を満たしていない')
                sys.exit()
        pre_block_count = block_count
        block_count2=0
        block_count = 0

for i in range(bigblock_num_int):
    for j in range(size):
        if tmp[j+(i*size)+1]==1:
            for k in range(size):
                f_tmp[lbcounter,k]=f[j,k]
            lbcounter=lbcounter+1

print(f_tmp)

counter=0

for i in range(size):
    for j in range(size):
        if f_tmp[i,j]<2:
            print(f_tmp[i,j],i,j)
            for k in range(size):
                f_swap[k, counter]=f_tmp[k,j]
            counter=counter+1
            

print(f_swap)
sum_taikaku=0
for i in range(bigblock_num_int):
    for j in range(bsize):
        for k in range(bsize):
            sum_taikaku = sum_taikaku + f_swap[j+(bsize*i), k+(bsize*i)]*f_swap[j+(bsize*i), k+(bsize*i)]

sum_all=0
for i in range(size):
    for j in range(size):
        sum_all=sum_all+f_swap[i,j]*f_swap[i,j]

hami=0
hami=sum_all-sum_taikaku
print("before")
print(before_H)
print("after")
print(sum_all-sum_taikaku)


np.savetxt(filename, f_swap, delimiter=',')