from re import T
from dwave.system import DWaveSampler,EmbeddingComposite
import numpy as np
from pyqubo import Array,Constraint,Placeholder

N=5
Q=np.array([[1000,20,20,50,40],
            [30,1000,10,30,20],
            [20,10,1000,30,20],
            [50,30,20,1000,10],
            [40,20,20,10,1000]])
x=Array.create('x',shape=(N,N),vartype='BINARY')
cost=0
#目的関数
for i in range(N):
    for j in range(N):
        cost+=Q[i][j]*x[i][j]
#制約1
constr1=0
for i in range(N):
    constr1+=(np.sum(x[i])-1)**2
#制約2
constr2=0
for i in range(N):
    constr2+=(np.sum(x.T[i])-1)**2

cost_func=cost+Placeholder('lambda')*Constraint(constr1,label='constr1')+Placeholder('lambda')*Constraint(constr1,label='constr2')
model=cost_func.compile()
feed_dict={'lambda':100.0}
qubo,offset=model.to_qubo(feed_dict=feed_dict)
response=DWaveSampler().sample_qubo(qubo,num_reads=100)
for s,e,o in response.data(['sample','energy','num_occurrences']):
    print(s,e,o)