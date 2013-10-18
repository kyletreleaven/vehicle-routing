
import itertools

import networkx as nx


def cyclefactor( demands, transport_graph, cost_label='cost' ) :
    # Step 1: Find a minimum cost bipartite matching
    # between the multisets of heads and tails of arcs.
    # (Should this be in shortest path sense?)
    
    # transport_graph may be a Graph or DiGraph
    APSP = nx.all_pairs_dijkstra_path_length( transport_graph, weight=cost_label )
    
    # demands must be a MultiDiGraph
    ARCS = demands.edges( keys=True )
    
    bp_graph = nx.Graph()
    for arc1, arc2 in itertools.product( ARCS, ARCS ) :
        head, tail = arc1[1], arc2[0]
        cost = APSP[head][tail]
        bp_graph.add_edge( (arc1,'head'), (arc2,'tail'), weight = -cost )   # we're interested in *minimum* weight
        
    MATCH = nx.matching.max_weight_matching( bp_graph, maxcardinality=True )
    
    # Step 2: post-poning til *later*
    
    # Step 3(a): Let the m disjoint cycles be Ri, i = 1..m,each represented by a single node ni.
    cycles = []
    BAG = demands.edges( keys=True )
    while len( BAG ) > 0 :
        cycle = []
        arc = BAG[0]
        while True :
            BAG.remove( arc )
            cycle.append( arc )
            next_arc, _ = MATCH[ (arc,'head') ]
            if next_arc in cycle : break
            arc = next_arc
        cycles.append( cycle )
        
    # user is responsible for the post-processing to generate Dijkstra shortest connecting paths.
    return cycles



