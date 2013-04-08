
import itertools
import networkx as nx
import numpy as np




""" finite set of stations """
numstations=4
br=2

I = [ 'i%d' % k for k in range(numstations) ]
    
"""
generate vehicle states with associated locations, 
and transformation costs
"""
M = nx.Graph()
A = nx.DiGraph()

X = []
config_count = itertools.count()
for i in I :
    for t in range(br) :
        x = 'x%d' % config_count.next()
        M.add_edge( i, x )
        X.append( x )
        
for x,y in itertools.product( X, X ) :
    A.add_edge( x,y, cost=np.random.rand() )
    
APSP = nx.all_pairs_dijkstra_path_length( A )
for x,y, data in A.edges_iter( data=True ) :
    if data.get('cost') > APSP[x][y] : A.remove_edge( x, y )
    
DEMAND_TYPES = [ pair for pair in itertools.product( I, I ) ]
NUMTYPES = len( DEMAND_TYPES )
TRIALS = 100
P = np.ones( NUMTYPES ) / NUMTYPES
outcome = np.random.multinomial( TRIALS, P )
demand_hist =  dict( zip( DEMAND_TYPES, outcome ) )








