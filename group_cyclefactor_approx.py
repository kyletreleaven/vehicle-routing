
import itertools

import numpy as np
import networkx as nx

import cvxpy

import cyclefactor


DEBUG = True


""" convenience utility """

def greedy_chairman_rounding( fractional_target ) :
    """
    should be a dictionary, from keys to real nums
    we assume that the values sum to an *integer*
    """
    target = fractional_target      # alias
    hist_out = dict()
    
    # initialize with lower bounds
    REMAIN_f = sum( target.values() )
    REMAIN = int( round( REMAIN_f ) )
    # this part was/is a little dicey; was "flooring", but 1.999999 became 1 instead of 2
    # as long as the precision is good enough, round() should have the desired effect
    #print REMAIN_f, REMAIN
    
    for key, x in target.iteritems() :
        min_occup = int( np.floor( x ) )
        hist_out[ key ] = min_occup
        REMAIN -= min_occup
        
    # distribute remaining demands greedily
    for k in range(REMAIN) :
        options = [ ( x - target[key], key ) for key, x in hist_out.iteritems() ]
        _, key = min( options )
        hist_out[key] += 1
        
    return hist_out




""" specialized utilities """

def construct_lp_relaxation( demand_hist, M, transport_graph, cost_label='cost' ) :
    """
    demand_hist is a dictionary from ordered station pairs (i,j) -> whole integers
    M is a membership graph with edges (i,x) for all corresponding stations-state pairs
    A is a DiGraph on X, with cost given by property 'cost'
    """
    APSP = nx.all_pairs_dijkstra_path_length( transport_graph, weight=cost_label )
    
    cost = 0.
    constr = []
    
    """ create enroute variables """
    S = nx.DiGraph()        # service "edges"
    for (i,j), NUM in demand_hist.iteritems() :
        VARS = []
        for x,y in itertools.product( M.neighbors_iter(i), M.neighbors_iter(j) ) :
            var = cvxpy.variable()
            distance = APSP[x][y]
            cost += var * distance
            constr.append( cvxpy.geq( var, 0. ) )   # needs to be non-negative
            
            S.add_edge( x, y, var=var, distance=distance )
            VARS.append( var )
            
        constr.append( cvxpy.eq( cvxpy.sum(VARS), NUM ) )   # total such trips should sum to num
            
    """ create balance variables """
    B = nx.DiGraph()
    HEADS = [ h for h, d in S.in_degree_iter() if d > 0 ]
    TAILS = [ t for t, d in S.out_degree_iter() if d > 0 ]
    for y,x in itertools.product( HEADS, TAILS ) :
        var = cvxpy.variable()
        cost += var * APSP[y][x]
        constr.append( cvxpy.geq( var, 0. ) )
        
        B.add_edge( y,x, var=var )
        
    """ generate flow conservation constraints """
    for x in S.nodes_iter() :
        sin = sum( [ data['var'] for _,__,data in S.in_edges_iter( x, data=True ) ] )
        sout = sum( [ data['var'] for _,__,data in S.out_edges_iter( x, data=True ) ] )
        bin = sum( [ data['var'] for _,__,data in B.in_edges_iter( x, data=True ) ] )
        bout = sum( [ data['var'] for _,__,data in B.out_edges_iter( x, data=True ) ] )
        
        constr.append( cvxpy.eq( sin, bout ) )
        constr.append( cvxpy.eq( sout, bin ) )
        
    #service_vars = S        # alias
    return cost, constr, S






""" Exports """

def group_cyclefactor_approx( demand_hist, M, transport_graph, cost_label='cost' ) :
    """ construct and solve the LP relaxation of the IP """
    cost, constr, SERVICE = construct_lp_relaxation( demand_hist, M, transport_graph, cost_label=cost_label )
    prog = cvxpy.program( cvxpy.minimize( cost ), constr )
    res = prog.solve( quiet=DEBUG )
    # still need to debug this garbage!!!! wtf is going on here?
    
    
    """ convert to proper format and "round" the solution """
    fractional = dict()
    
    trips_hist = dict()
    for (i,j) in demand_hist :
        target = dict()
        for edge in itertools.product( M.neighbors_iter(i), M.neighbors_iter(j) ) :
            temp = SERVICE.get_edge_data( *edge ).get('var').value
            fractional[ edge ] = temp
            target[ edge ] = temp
            
        demand_update = greedy_chairman_rounding( target )
        trips_hist.update( demand_update )
        
        
    """ construct and solve the resulting cycle factor instance """
    service_edges = nx.MultiDiGraph()
    for (x,y), NUMTRIPS in trips_hist.iteritems() :
        weight = nx.dijkstra_path_length( transport_graph, x, y, weight=cost_label )
        for k in range(NUMTRIPS) :
            service_edges.add_edge( x,y, cost=weight )
    
    cycles = cyclefactor.cyclefactor( service_edges, transport_graph )
    return cycles




""" Testing """

if __name__ == '__main__' :
    import testcase
    
    demand_hist = testcase.demand_hist
    M = testcase.M
    transport_graph = testcase.A
    
    cycles = group_cyclefactor_approx( demand_hist, M, transport_graph )
    
    
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

    
    
    
    #
