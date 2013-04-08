
import itertools

import numpy as np
import networkx as nx

import cvxpy

import group_cyclefactor_approx


""" convenience """
class node : pass

def eulerian_circuit_verbose( multidigraph, source=None ) :
    registry = dict()
    for u,v,key in multidigraph.edges_iter( keys=True ) :
        pair = u,v
        id = u,v,key
        if not pair in registry : registry[pair] = []
        registry[ pair ].append( id )
        
    for pair in nx.eulerian_circuit( multidigraph, source ) :
        options = registry[ pair ]
        yield options.pop(0)
        

def cyclegraph( cycle ) :
    mdg = nx.MultiDiGraph()
    #
    for u,v, key in cycle : mdg.add_edge( u,v, key )
    #
    for edge1, edge2 in zip( cycle, cycle[1:] + cycle[0] ) :
        _,v,__ = edge1
        u,_,__ = edge2
        mdg.add_edge( v, u, CONNECT_ONLY=True )
    #
    return mdg

def cycle_edges( cycle ) :
    enroute = nx.MultiDiGraph()
    balance = nx.MultiDiGraph()
    #
    for u,v, key in cycle :
        enroute.add_edge( u,v, key )
    for edge1, edge2 in zip( cycle, cycle[1:] + [ cycle[0] ] ) :
        _,v,__ = edge1
        u,_,__ = edge2
        #
        balance.add_edge( v, u )
    #
    return enroute, balance


""" Exports """

def largearcs_connecting_heuristic( cycles, transport_graph, cost_label='cost' ) :
    APSP = nx.all_pairs_dijkstra_path( transport_graph, weight=cost_label )
    APSPL = nx.all_pairs_dijkstra_path_length( transport_graph, weight=cost_label )
    
    # Step 3(b): Form the inter-node distance from the original edge costs
    # d(ni,nj) = min { c'(u,v) | u \in Ri, v \in Rj }.
    # Associate with (ni,nj) the edge (u,v) yielding minimum cost.
    inter_node_distance = nx.Graph()
    for cycle in cycles :
        u = node()
        u.cycle = cycle     # do we need this?... eh, it's cheap
        inter_node_distance.add_node( u )
        #
        u.enroute, u.balance = cycle_edges( cycle )
        
        # also want to store all nodes visited while *EMPTY*
        u.nodes = set()
        for x,y in u.balance.edges_iter() :
            u.nodes.update( APSP[x][y] )
        #for x,y,key, data in u.graph.edges_iter( keys=True, data=True ) :
            #if not data.get( 'CONNECT_ONLY', False ) : continue
            #u.nodes.update( APSP[x][y] )
            
    NODES = inter_node_distance.nodes()
    for u, v in itertools.combinations( NODES, 2 ) :
        options = [ ( APSPL[x][y] + APSPL[y][x], (x,y) ) for x,y in itertools.product( u.nodes, v.nodes ) ]
            # round trip cost
        cost, edge = min( options )
        inter_node_distance.add_edge( u, v, cost=cost, edge=edge )
        
    # Step 4: Find a minimum cost spanning tree on inter-node distance graph
    MST = nx.algorithms.minimum_spanning_tree( inter_node_distance, weight='cost' )
    
    # Step 2: Initialize PRETOUR to be empty. For each edge in the matching,
    # associate a direction (from head to tail); insert into PRETOUR
    # (Also, insert the arcs...)
    eulerian = nx.MultiDiGraph()
    for u in NODES :
        for x,y, key in u.enroute.edges_iter( keys=True ) :
            eulerian.add_edge( x,y, key, SERVICE=True )
            
    for u in NODES :
        for x,y in u.balance.edges_iter() :
            eulerian.add_edge( x, y )
    for _,__,data in MST.edges_iter( data=True ) :
        x,y = data.get('edge')
        eulerian.add_edge( x, y )
        eulerian.add_edge( y, x )
        
    try :
        tour = []
        for edge in eulerian_circuit_verbose( eulerian ) :
            if not eulerian.get_edge_data( *edge ).get('SERVICE', False ) : continue
            tour.append( edge )
            
        return tour
    
    except :
        return eulerian





def group_stackercrane_approx( demand_hist, M, transport_graph, cost_label='cost' ) :
    cycles = group_cyclefactor_approx.group_cyclefactor_approx( demand_hist, M, transport_graph, cost_label )
    tour = largearcs_connecting_heuristic( cycles, transport_graph, cost_label )
    return tour






if __name__ == '__main__' :
    
    import testcase
    
    
    demand_hist = testcase.demand_hist
    M = testcase.M
    transport_graph = testcase.A
    
    import group_cyclefactor_approx as gcfp
    
    cycles_old = gcfp.group_cyclefactor_approx( demand_hist, M, transport_graph )
    
    sctour = group_stackercrane_approx( demand_hist, M, transport_graph )
    cycles = [ sctour ]
    
    
    station_tours = []
    for cycle in cycles :
        tour = []
        for x,y, _ in cycle :
            i = M.neighbors(x)[0]
            j = M.neighbors(y)[0]
            tour.append( (i,j) )
        station_tours.append( tour )
    
    freq_dict = dict()
    for tour in station_tours :
        for i,j in tour :
            if not freq_dict.has_key( (i,j) ) : freq_dict[ (i,j) ] = 0
            freq_dict[ (i,j) ] += 1

    
    






