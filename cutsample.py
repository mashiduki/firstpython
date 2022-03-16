import networkx as nx
from collections import defaultdict
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite

G=nx.Graph()
#G.add_edges_from([(1,2),(1,3),(2,4),(3,4),(3,5),(4,5)])
G.add_edges_from([(1,2),(2,3),(3,4),(4,5),(6,7),(7,8),(8,9),(9,10),(11,12),(12,13),(13,14),(14,15),(1,6),(6,11),(2,7),(7,12),(3,8),(8,13),(4,9),(9,14),(5,10),(10,15)])
Q=defaultdict(int)

for u,v in G.edges:
    Q[(u,u)]+=-1
    Q[(v,v)]+=-1
    Q[(u,v)]+=2
print(Q)

numruns=1000
sampler=EmbeddingComposite(DWaveSampler())
response=sampler.sample_qubo(Q,num_reads=numruns)
for s,e,o,c in response.data(['sample','energy','num_occurrences','chain_break_fraction']):
    print(s,e,o,c)
