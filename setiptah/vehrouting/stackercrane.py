

import itertools
import networkx as nx

from setiptah.matching import max_weight_matching




class StackerCraneFHK(object) :
    """
    Input:
        -E, an undirected nx.Graph object, E;
        -A, an nx.DiGraph or nx.MultiDiGraph object
            (the algorithm will convert nx.DiGraph to nx.MultiDiGraph);
        -edge/arc attribute 'cost' will be used as attribute of input function c
        
        other:
        -vs, a designated starting node
        
            
    Output:
        -nothing yet
        
    Comments:
        A represents the multiset of arcs that the stacker crane tour must traverse (each once, exactly);
        E represents a set of edges that the stacker crane tour may traverse arbitrarily, at known cost
    """
    class node(object) :
        """ helper class, when generic graph nodes are needed """
        def __init__( self, label=None ) :
            if label : self.label = label
            
    @classmethod
    def bipartite_matching_algorithm(cls, bp_graph ) :
        return nx.matching.max_weight_matching( bp_graph, maxcardinality=True )
            
            
    @classmethod
    def PREPROCESS(cls, E, A, vs ) :
        """
        inputs:
            E, an undirected 'cost'-weighted graph;
            A, a directed 'cost'-weighted graph or multigraph
            vs, a starting node, or None
        """
        workE = E.copy()
        workA = nx.MultiDiGraph( A ) # A.copy()
        
        if vs == None : vs, head = workA.edges_iter().next()
        
        # Step 1: If vs is not an end point of any arc, create a terminal vertex vf
        # and an arc <vf,vs> of cost zero.
        # For all vi in V, create edges (vi,vf) of cost c(vi,vs).
        if A.in_degree( vs, weight=None ) <= 0 :
            vf = cls.node()
            workA.add_edge( vf, vs, cost=0.0 )
            for i,j, data in E.edges( vs, data=True ) :
                if i is vs :
                    t = j
                else :
                    t = i
                workE.add_edge( t, vf, data )
                
        # Step 2: Insert all vertices which are endpoints of arcs into Vt.
        # Calculate the shortest path b/w every pair of vertices vi and vj in V',
        # insert (vi,vj) into Et, and set c(vi,vj) to the cost obtained.
        paths_dict = nx.algorithms.all_pairs_shortest_path_length( workE )
        EE = nx.Graph()
        for i, temp in paths_dict.iteritems() :
            for j, cost in temp.iteritems() :
                EE.add_edge( i, j, cost=cost )
                
        return EE, workA, vf
    
    
    @classmethod
    def CYCLEFACTOR(cls, E, A ) :
        cls.E, cls.A = E, A
        
        # Step 1: Find a minimum cost bipartite matching
        # between the multisets of heads and tails of arcs.
        # (Should this be in shortest path sense?)
        ARCS = A.edges( keys=True )
        
        bp_graph = nx.Graph()
        for arc1, arc2 in itertools.product( ARCS, ARCS ) :
            head, tail = arc1[1], arc2[0]
            if E.has_edge( head, tail ) :
                cost = E.get_edge_data( head, tail )['cost']
                bp_graph.add_edge( (arc1,'head'), (arc2,'tail'), weight = -cost )    # we're interested in *minimum* weight
        #MATCH = nx.matching.max_weight_matching( bp_graph, maxcardinality=True )
        MATCH = StackerCraneFHK.bipartite_matching_algorithm( bp_graph )
        
        # debug
        cls.bp_graph = bp_graph
        cls.MATCH = MATCH
        
        # Step 2: *later*
        
        # Step 3(a): Let the m disjoint cycles be Ri, i = 1..m,each represented by a single node ni.
        cycles = []
        BAG = A.edges( keys=True )
        while len( BAG ) > 0 :
            cycle = []
            arc = BAG[0]
            while True :
                BAG.remove( arc )
                cycle.append( arc )
                next_arc, throw = MATCH[ (arc,'head') ]
                if next_arc in cycle : break
                arc = next_arc
            cycles.append( cycle )
            
        # debug
        cls.cycles = cycles
        
        return cycles
    
    
    @classmethod
    def LARGEARCS(cls, E, A ) :
        cls.E, cls.A = E, A
        ARCS = A.edges( keys=True )
        
        # Steps 1, 2, and 3(a) are the "Cycle Factor Problem"
        cycles = cls.CYCLEFACTOR( E, A )
        
        # Step 3(b): Form the inter-node distance from the original edge costs
        # d(ni,nj) = min { c'(u,v) | u \in Ri, v \in Rj }.
        # Associate with (ni,nj) the edge (u,v) yielding minimum cost.
        inter_node_distance = nx.Graph()
        for cycle in cycles :
            u = cls.node()
            u.cycle = cycle
            u.graph = nx.MultiDiGraph()
            #u.graph.add_edges_from( cycle )
            edges_iter = ( (i,j,key, {} ) for (i,j,key) in cycle )
            u.graph.add_edges_from( edges_iter )
            #
            inter_node_distance.add_node( u )
            
        # debug
        cls.inter_node_distance = inter_node_distance
            
        NODES = inter_node_distance.nodes()
        for u, v in itertools.combinations( NODES, 2 ) :
            nodes_u = u.graph.nodes()
            nodes_v = v.graph.nodes()
            opt = [ ( E.get_edge_data( i,j )['cost'], (i,j) )
                    for i,j in itertools.product( nodes_u, nodes_v )
                    if E.has_edge( i, j ) ]
            cost, edge = min( opt )
            inter_node_distance.add_edge( u, v, cost=cost, edge=edge )
            
        # Step 4: Find a minimum cost spanning tree on inter-node distance graph
        MST = nx.algorithms.minimum_spanning_tree( inter_node_distance, weight='cost' )
        # debug
        cls.MST = MST
        
        # Step 2: Initialize PRETOUR to be empty. For each edge in the matching,
        # associate a direction (from head to tail); insert into PRETOUR
        # (Also, insert the arcs...)
        PRETOUR = nx.MultiDiGraph()
        for ( arc1, side ), ( arc2, throw ) in cls.MATCH.iteritems() :
            if not side == 'head' : continue
            head, tail = arc1[1], arc2[0]
            PRETOUR.add_edge( head, tail )
            
        PRETOUR.add_edges_from( A.edges_iter() )
        
        # Step 5: Rename each spanning edge in terms of the original vertices.
        # Make two copies of each edge, one for each direction; insert into Ett.
        for u,v in MST.edges_iter() :
            i,j = MST.get_edge_data( u, v )['edge']
            PRETOUR.add_edge( i, j )
            PRETOUR.add_edge( j, i )
            
        # debug
        cls.PRETOUR = PRETOUR
        
        # Step 6: POSTPROCESS
        return PRETOUR
        #tour = cls.POSTPROCESS( PRETOUR )
        #return tour
        
        
        
        
    @classmethod
    def SMALLARCS(cls, E, A ) :
        cls.E, cls.A = E, A
        
        # Step 1(a): Represent each arc <ti,hi> by a node ni.
        # For each pair of nodes ni and ni, define c'(ni,nj) as the min of
        # c(ti,tj), c(ti,hj), c(hi,tj), and c(hi,hj).
        ARCS = A.edges( keys=True )
        NODES = [ cls.node() for i in range( len( ARCS ) ) ]
        NODE_MAP = {}
        for a, u in zip( ARCS, NODES ) :
            NODE_MAP[a] = u
            NODE_MAP[u] = a
            
        node_graph = nx.Graph()
        ct_graph = nx.Graph()
        cls.ct_graph = ct_graph
        for u, v in itertools.combinations( NODES, 2 ) :
            # produces [ (ti,tj), (ti,hj), (hi,tj), (hi,hj) ]
            arci, arcj = NODE_MAP[u], NODE_MAP[v]
            options = []
            #print ''
            for edge in itertools.product( arci[:2], arcj[:2] ) :
                #print edge[0], edge[1]
                if E.has_edge( *edge ) : 
                    cost = E.get_edge_data( *edge )['cost']
                    options.append( ( cost, edge ) )
                if len( options ) <= 0 : continue
                cost, edge = min( options )
                #
                node_graph.add_edge( u, v )
                ct_graph.add_edge( u, v, cost=cost, edge=edge )
                
        # Step 1(b): Perform an all shortest paths algorithm using c' 
        # to find the inter-node distance d(ni,nj), between ni and nj.
        # Associate with each edge (ni,nj) the edges in the shortest path.
        paths = nx.all_pairs_dijkstra_path( ct_graph, weight='cost' )
        lengths = nx.all_pairs_dijkstra_path_length( ct_graph, weight='cost' )
        
        d_graph = nx.Graph()
        cls.d_graph = d_graph
        for u, v in itertools.combinations( NODES, 2 ) :
            try :
                path = paths[u][v]
                length = lengths[u][v]
                d_graph.add_edge( u, v, cost=length, path=path )
                
            except KeyError :
                continue
            
        # Step 2: Find a minimum cost spanning tree for the nodes {ni}, using the distance function d.
        for u, v, data in node_graph.edges_iter( data=True ) :
            data['weight'] = d_graph.get_edge_data( u, v )['cost']
        MST = nx.algorithms.minimum_spanning_tree( node_graph, weight='weight' )
        cls.MST = MST
        
        # Step 3: Identify nodes of odd degree in the spanning tree,
        # and perform a minimum cost matching on these nodes using the distance function d.
        ODD_NODES = [ u for ( u, deg ) in MST.degree_iter() if deg % 2 == 1 ]
        oddnode_graph = node_graph.subgraph( ODD_NODES )
        for u, v, data in oddnode_graph.edges_iter( data=True ) :
            data['weight'] = -d_graph.get_edge_data( u, v )['cost']
        MATCH = max_weight_matching( oddnode_graph, weight='weight', maxcardinality=True )
        
        # debug
        cls.oddnode_graph = oddnode_graph
        cls.MATCH = MATCH
        
        # Step 4(a): Rename the spanning edges and matching edges in terms of V'.
        # Define Gtt = (V,Ett,A), where Ett is the multiset of spanning and matching edges.
        Ett = nx.MultiGraph()
        cls.Ett = Ett
        
        iter1 = MST.edges_iter()
        iter2 = MATCH.edges_iter()
        for u, v in itertools.chain( iter1, iter2 ) :
            edge = ct_graph.get_edge_data( u, v )['edge']
            Ett.add_edge( *edge )
            
        # Step 4(b):
        # Consider degree of nodes in Gtt.
        # For all arcs <u,v> whose endpoints have odd degree, add the edges (u,v) to Ett,
        # and associate with it the direction opposite of that arc.
        # These arcs requiring no such edge shall be called even arcs, with total cost CAt...
        def degree( i ) :
            deg = A.degree(i)
            if Ett.has_node(i) : deg += Ett.degree(i)
            return deg
        
        Ett_dir = nx.MultiDiGraph()
        cls.Ett_dir = Ett_dir
        EVEN_ARCS, ODD_ARCS = [], []
        cls.EVEN_ARCS = EVEN_ARCS
        CtA = 0.0
        for arc in ARCS :
            i,j,key = arc
            # only need to check *one* of the endpoints; each arc has endpoints of same degree parity
            if degree( i ) % 2 == 1 :
                ODD_ARCS.append( arc )
                Ett_dir.add_edge( j, i )
            else :
                EVEN_ARCS.append( arc )
                CtA += A.get_edge_data( *arc )['cost']
        
        
        # Step 5(a): Find a tour which exactly covers Ett union A, ignoring [*only*] even arc directions.
        """
        This part gave me some *serious* trouble. Originally, I ignored the direction of *all* edges,
        believing that *any* eulerian cycle would traverse the doubled odd arcs once in each direction.
        That assumption is provably false, much to my chagrin.
        My solution has been (1) to compute eulerian cycles within connected components of 
        the union graph of EVEN_ARCS and Ett,
        (2) add the doubled odd arcs, and (3) compute an eulerian cycle on the resulting MultiDiGraph.
        The new approach prevents doubled odd arcs from being traversed both the *same* direction
        (that wreaked havoc on the rest of the algorithm from time to time);
        the approach retains optimality.
        (I am not *sure* that there can be strictly more than one connected component;
        however, this way we are prepared for anything.)
        """
        eulerian_component_graph = nx.MultiGraph()
        eulerian_component_graph.add_edges_from( ( arc[:2] for arc in EVEN_ARCS ) )
        eulerian_component_graph.add_edges_from( Ett.edges_iter() )
        
        eulerian_walk_graph = nx.MultiDiGraph()
        cls.eulerian_component_graph = eulerian_component_graph
        #
        for cmp in nx.connected_component_subgraphs( eulerian_component_graph ) :
            circuit = nx.eulerian_circuit( cmp )
            eulerian_walk_graph.add_edges_from( circuit )
        #
        eulerian_walk_graph.add_edges_from( ( arc[:2] for arc in ODD_ARCS ) )
        eulerian_walk_graph.add_edges_from( Ett_dir.edges_iter() )
        
        circuit = nx.eulerian_circuit( eulerian_walk_graph )
        walk = [ edge for edge in circuit ]     # flatten
        cls.walk = walk
        #cls.circuit = circuit
        
        
        # Step 6(b): If the cost of the even arcs which are traversed backwards is more than (1/2)CAt,
        # then reverse the direction of the tour.
        """
        Most of the following code is simply to *detect* the backward-traversed arcs.
        Mis-diagnosis of issues with early implementations of the prequel mixed-multi-graph-eulerian-cycle algorithm,
        led me to focus significant attention "fixing" the following code.
        (I replaced a greedy implementation (if-Falsed-out, below), with a "bullet-proof" version based on bipartite matching;
        I have been reluctant to revisit the old version, though the greedy algorithm would be surely faster.
        An even better alternative might be to merge the eulerian cycles detected in the above code, which shouldn't be too hard.)
        """
        workA = A.copy()
        workEtt = Ett.copy()
        workEtt_dir = Ett_dir.copy()
        FRWD, BKWD = [], []
        Ettt = nx.MultiDiGraph()
        #
        if False :
            # I *thought* this greedy method would work, but I had some issues with non-eulerian cycles.
            # (The greedy algorithm *could* still work... )
            for i, j in walk :
                if workA.has_edge( i, j ) :
                    key = workA.edge[i][j].iterkeys().next()
                    arc = i,j,key
                    if arc in EVEN_ARCS : FRWD.append( arc )
                    workA.remove_edge( *arc )
                    
                elif workEtt_dir.has_edge( i, j ) :
                    workEtt_dir.remove_edge( i, j )
                    
                elif workEtt.has_edge( i, j ) :
                    workEtt.remove_edge( i, j )
                    Ettt.add_edge( i, j )       # actually, part of Step 6.
                    
                else : # graph A better have (j,i)...
                    key = workA.edge[j][i].iterkeys().next()
                    arc = j,i,key
                    #if arc in EVEN_ARCS : BKWD.append( arc )
                    assert arc in EVEN_ARCS # if arc is not even, then we're doing something wrong already.
                    BKWD.append( arc )
                    workA.remove_edge( *arc )
                    
        else :      
            # The problem might not require *this* much machinery, but for now I just need the code to work...
            #
            # The idea of this code is to create a bipartite graph,
            # with the ordered eulerian-walk edges on one side, and the "physical" arcs and edges of A, Ett, and Ett_dir on the other.
            # (Edges exist between "physical" edges and the eulerian walk edges that they "could explain".)
            # Thus, a bipartite matching gives us a mapping of the walk onto the "physical" edges.
            
            # assign indexable placeholders to each "physical" edge
            for arc in workA.edges_iter( keys=True ) :
                workA.get_edge_data( *arc )['circuit_node'] = cls.node( repr(arc) )
            for arc in workEtt_dir.edges_iter( keys=True ) :
                workEtt_dir.get_edge_data( *arc )['circuit_node'] = cls.node( repr(arc) )
            for edge in workEtt.edges_iter( keys=True ) :
                workEtt.get_edge_data( *edge )['circuit_node'] = cls.node( repr(edge) )
                
            # connect each walk edge to its candidate "physical" edges within similarity_graph (bipartite)
            similarity_graph = nx.Graph()
            cls.similarity_graph = similarity_graph
            for idx, edge in enumerate( walk ) :
                i,j = edge
                
                # connect to arcs of A in forward direction
                data_tree = workA.get_edge_data( i, j, default={} )
                for data in data_tree.itervalues() :
                    similarity_graph.add_edge( idx, data['circuit_node'] )
                        
                # connect to *odd*-arc companions in forward direction
                data_tree = workEtt_dir.get_edge_data( i, j, default={} )
                for data in data_tree.itervalues() :
                    similarity_graph.add_edge( idx, data['circuit_node'] )
                
                # connect to *even* arcs of A in the backward direction
                data_tree = workA.get_edge_data( j, i, default={} )
                for key, data in data_tree.iteritems() :
                    arc = j,i,key
                    if arc in EVEN_ARCS :
                        similarity_graph.add_edge( idx, data['circuit_node'] )
                            
                # connect to undirected edges
                data_tree = workEtt.get_edge_data( i, j, default={} )
                for data in data_tree.itervalues() :
                    similarity_graph.add_edge( idx, data['circuit_node'] )
                        
            # compute the matching between walk edges and "physical" edges
            eulerian_match = nx.max_weight_matching( similarity_graph )
            cls.eulerian_match = eulerian_match
            
            # identify the forward- vs. backward- traversed arcs
            for arc in EVEN_ARCS :
                i,j,key = arc
                data = workA.get_edge_data( *arc )
                node = data['circuit_node']
                edge_idx = eulerian_match[ node ]
                u,v = walk[edge_idx]
                
                if u==i and v==j :
                    FRWD.append( arc )
                elif u==j and v==i :
                    BKWD.append( arc )
                else :
                    raise Exception('bad matching (eulerian walk)!')    # good, this doesn't seem to ever happen!
                
            # associate direction with the undirected edges
            for edge in workEtt.edges_iter( keys=True ) :
                u = workEtt.get_edge_data( *edge )['circuit_node']
                idx = eulerian_match[u]
                Ettt.add_edge( *walk[idx] )
                        
                        
                        
        # detection is done---check for smaller direction
        backward_cost = sum([ A.get_edge_data( *arc )['cost'] for arc in BKWD ])
        if backward_cost > .5 * CtA :
            FRWD, BKWD = BKWD, FRWD
            Ettt.reverse( copy=False )
            
        #debug
        cls.FRWD, cls.BKWD = FRWD, BKWD
        
        
        # Step 6(c): Associate with each undirected edge the direction of the tour. (Done.)
        # For every arc that is incorrectly traversed, add two edges to Ett,
        # both with associated direction opposite that of the arc.
        PRETOUR = nx.MultiDiGraph()
        PRETOUR.add_edges_from( A.edges_iter() )
        PRETOUR.add_edges_from( Ett_dir.edges_iter() )
        PRETOUR.add_edges_from( Ettt.edges_iter() )
        for i,j,key in BKWD :
            for times in range(2) :
                PRETOUR.add_edge( j, i )
        
        # Step 7: POSTPROCESS
        cls.PRETOUR = PRETOUR
        return PRETOUR
        
        
        
        
    @classmethod
    def POSTPROCESS(cls, digraph_eul, A ) :
        """
        This implementation is different than what is specified in FHK.
        It simply returns an ordering of arcs in A.
        The process of moving from arc to arc by shortest paths should be straightforward.
        (Nevertheless, the steps of FHK are written below in comments.)
        """
        # Step 1: Given a set of arcs and edges with associated direction,
        # from which a tour can be constructed, find the tour.
        
        # Step 2: For any series of two or more consecutive edges in the tour,
        # replace them with one edge, thus eliminating intermediate vertices. (???)
        
        # Step 3: For any two vertices vi and vj anchoring an edge in the tour,
        # insert any vertices that would create a shorter path between vi and vj.
        # (Undo Step 2 of PREPROCESS.)
        
        # Step 4: List the tour beginning at vs, the initial vertex, and finishing at vf,
        # the terminal vertex.
        
        # Alternative Method:
        arcs = {}
        for arc in A.edges( keys=True ) :
            arcs.setdefault( arc[:2], [] ).append( arc )
            
        new_walk = []
        for edge in nx.eulerian_circuit( digraph_eul ) :
            if edge in arcs :
                longver = arcs[edge]
                new_walk.append( longver.pop(0) )
                if len( longver ) <= 0 : del arcs[edge]
        return new_walk
    
    
    @classmethod
    def CRANE(cls, E, A, vs=None ) :
        if vs == None : vs = A.nodes_iter().next() # any, irrelevant node
        Et, At, vs = cls.PREPROCESS( E, A, vs )
        ans_large = cls.LARGEARCS( Et, At )
        ans_small = cls.SMALLARCS( Et, At )
        cls.largearcs, cls.smallarcs = ans_large, ans_small
        # find out which is smaller...
        
        # for now, just return the largearcs solution
        return cls.POSTPROCESS( ans_large, At )
        
        
        
        
    """ DEBUGGING UTILITIES """
    @classmethod
    def display_Gtt(cls, Ett_dir=True ) :
        #lay = nx.layout.circular_layout( cls.A )
        lay = nx.layout.spring_layout( cls.A )
        nx.draw( cls.A, pos=lay )
        nx.draw_networkx_edges( cls.Ett, pos=lay )
        if Ett_dir :
            nx.draw_networkx_edges( cls.Ett_dir, pos=lay )
            
    @classmethod
    def draw_similarility_graph( cls ) :
        labels = {}
        layout = {}
        lhs_count = itertools.count()
        rhs_count = itertools.count()
        
        graph = cls.similarity_graph
        for node in graph.nodes() :
            if type( node ) == type(0) :
                labels[node] = cls.walk[node]
                layout[node] = np.array( (0, lhs_count.next() ) )
            else :
                labels[node] = node.label
                layout[node] = np.array( (1, rhs_count.next() ) )
        nx.draw_networkx_edges( graph, pos=layout )
        nx.draw_networkx_labels( graph, layout, labels )
        
        
        
        
        
        
        
        
        
        
        
        
        
        
if __name__ == '__main__' :
    
    import sys, time
    import numpy as np
    
    import matplotlib as mpl
    mpl.rcParams['ps.useafm'] = True
    mpl.rcParams['pdf.use14corefonts'] = True
    mpl.rcParams['text.usetex'] = True
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import axes3d
    from matplotlib.collections import PolyCollection
    
    #from dynpdp import DemandBatch, PickupDeliveryDemandBatch, StackerCraneTour
    #from dynpdp.random.demand import demands_centeredcube
    # debug
    fhk = StackerCraneFHK
    
    
    class node(object) : pass
    
    POINT_LUT = dict()
    indexer = itertools.count()
    
    def distance( u, v ) :
        return np.linalg.norm( POINT_LUT[v] - POINT_LUT[u] )
    
    
    """ instance generation """
    N = 10
    PICKS, DELVS = [], []
    for i in range(N) :
        for side in [ PICKS, DELVS ] :
            u = indexer.next()
            POINT_LUT[u] = np.random.uniform(size=2)
            side.append( u )
            
    A = nx.MultiDiGraph()
    for i, p in enumerate( PICKS ) :
        d = DELVS[i]
        cost = distance( p, d )
        A.add_edge( p,d, cost=cost )
        
    E = nx.Graph()
    for x, y in itertools.combinations_with_replacement( PICKS + DELVS, 2 ) :
        cost = distance( x, y )
        E.add_edge( x, y, cost=cost )
        
    vs = PICKS[0]

    def display_tour( tour ) :
        #lay = nx.layout.circular_layout( tour )
        labels = dict( [ (i,i) for i in POINT_LUT ] )
        nx.draw_networkx_nodes( tour, POINT_LUT, PICKS, node_color='r', alpha=.6 ) #marker='x' )
        nx.draw_networkx_nodes( tour, POINT_LUT, DELVS, node_color='b', alpha=.6 ) #marker='o' )
        nx.draw_networkx_labels( tour, POINT_LUT, labels )
        nx.draw_networkx_edges( tour, POINT_LUT )

    
    #stuff = StackerCraneFHK.PREPROCESS( E, A, vs )
    LA = StackerCraneFHK.LARGEARCS( E, A )    # LARGEARCS seems to be working!
    SA = StackerCraneFHK.SMALLARCS( E, A )
    
    
    
