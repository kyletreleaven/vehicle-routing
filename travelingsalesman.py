
import ktutils
import itertools
import networkx as nx

import numpy as np

""" Note to self: for BASIC persistence, use primitive data types stored in a dict, and pickle the dict. """
import pickle

from matching import max_weight_matching




def tsp_mstheuristic( graph, depot ) :
    # the simple MST heurisitc, should exhibit the right rate of growth at least
    mst = nx.minimum_spanning_tree( graph, weight='weight' )
    gg = nx.MultiGraph()
    gg.add_edges_from( mst.edges_iter() )
    gg.add_edges_from( mst.edges_iter() )
    circuit = nx.eulerian_circuit( gg, depot )
    tour = []
    for i,j in nx.eulerian_circuit( gg, depot ) :
        if i not in tour :
            tour.append( i )
    tour.append( depot )
    return tour

def tourlength( tour, lut ) :
    links = ( np.linalg.norm( lut[j] - lut[i] ) for i,j in slidingtuple( tour ) )
    return sum( links )




def node_graph( coords ) :
    enum = ( (i,c) for i,c in enumerate( coords ) )
    lut = dict( enum )
    graph = nx.Graph()
    graph.add_nodes_from( lut )
    return graph, lut

def get_metricgraph( coords ) :
    graph, lut = node_graph( coords )
    
    for i,j in itertools.combinations( lut, 2 ) :
        weight = np.linalg.norm( lut[j] - lut[i] )
        graph.add_edge( i, j, weight=weight )
    return graph, lut

# tsp shortcut
def batch_tsp( batch, depot_location=np.zeros(2) ) :
    tour_batch = [ depot_location ] + batch
    graph, lut = get_metricgraph( tour_batch )
    tour = tsp_mstheuristic( graph, 0 )
    return tour, lut














class TravelingSalesmanChristophides(object) :
    """
    Input:
        -E, a metric graph (completely-connected, satisfying triangle identity)
        
    Output:
        -a cyclic ordering of the nodes of E
    """
    class node(object) :
        """ helper class, when generic graph nodes are needed """
        def __init__( self, label=None ) :
            if label : self.label = label
            
            
    @classmethod
    def TSP(cls, graph, weight_attr='weight' ) :
        assert not graph.is_multigraph()
        cls.graph = graph
            
        # Step 1: Find a minimum cost spanning tree for the nodes {ni}
        MST = nx.algorithms.minimum_spanning_tree( graph, weight=weight_attr )
        cls.MST = MST
        
        # Step 2: Identify nodes of odd degree in the spanning tree,
        # and perform a minimum cost matching on these nodes.
        ODD_NODES = [ u for ( u, deg ) in MST.degree_iter() if deg % 2 == 1 ]
        oddnode_graph = graph.subgraph( ODD_NODES )
        # want a *minimum* weight matching
        match_graph = nx.Graph()
        cls.match_graph = match_graph
        for u, v, data in oddnode_graph.edges_iter( data=True ) :
            weight = -data[weight_attr]
            match_graph.add_edge( u, v, weight=weight )
        #oddnode_graph = graph.subgraph( ODD_NODES ).copy()
        #for u, v, data in oddnode_graph.edges_iter( data=True ) : data[weight_attr] *= -1
        
        MATCH = max_weight_matching( oddnode_graph, weight='weight', maxcardinality=True )
        cls.MATCH = MATCH
        
        # Step 3: Define Ett as the multiset of spanning and matching edges.
        eulerian_graph = nx.MultiGraph()
        cls.eulerian_graph = eulerian_graph
        
        eulerian_graph.add_edges_from( MST.edges_iter() )
        eulerian_graph.add_edges_from( MATCH.edges_iter() )
        
        # Step 4: Form a Hamiltonian cycle
        circuit = nx.eulerian_circuit( eulerian_graph )
        walk = [ edge for edge in circuit ]     # flatten
        cls.walk = walk
        #cls.circuit = circuit
        
        # Step 5: Make the cycle Hamiltonian by shortcutting.
        tour = []
        for u, v in walk :
            if u not in tour :
                tour.append( u )
                
        return tour
        
        




