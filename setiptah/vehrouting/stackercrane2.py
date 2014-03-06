
import itertools

import networkx as nx
#from setiptah.matching import max_weight_matching

"""
distance is assumed to be a metric distance function
arcs is a list of arcs
getHead and getTail can be applied to arcs to obtain the head and tail
(coordinates which can be fed to distance)

kwargs in case need to override matching algorithm
"""

# "symbols"
if False :
    ORIGINAL, CYCLE, TREE = 'O', 'C', 'T'
    TAIL, HEAD = 'tail', 'head'
    WALK, AGENT = 'walk', 'agent'
else :
    ORIGINAL, CYCLE, TREE = 0, 1, 2
    TAIL, HEAD = 0, 1
    WALK, AGENT = 'w', 'a'




def CYCLE2WALK( cycle ) :
    walk = []
    chain = cycle.copy()
    
    j = chain.iterkeys().next()
    while j in chain :
        i = j
        walk.append(i)
        j = chain.pop(i)
        
    return walk

def WALK2CYCLE( walk ) :
    return { i: j for i, j in zip( walk, walk[1:] + walk[:1] ) }


def SPLITCHAINS( chains ) :
    res = []
    
    M = chains.copy()
    while len( M ) > 0 :
        curr = {}
        j = M.iterkeys().next()
        while j not in curr :
            i = j
            j = M[i]
            curr[i] = j
            
        res.append( curr )
        for i in curr : del M[i]
        
    return res


def CYCLEFACTOR( arcs, getTail, getHead, distance, **kwargs ) :
    """
    obtain a cycle factor, expressed as a chain, or dictionary-based
    adjacency list
    """
    # need to match heads to tails
    graph = nx.Graph()
    TAIL, HEAD = 0, 1
    for i, arc1 in enumerate( arcs ) :
        for j, arc2 in enumerate( arcs ) :
            graph.add_edge( (i,HEAD), (j,TAIL),
                            weight = -distance( getHead(arc1), getTail(arc2) ) )
                            
    temp = nx.max_weight_matching( graph, maxcardinality=True )
    chains = { i : temp[(i,HEAD)][0] for i in xrange(len(arcs)) }
    return SPLITCHAINS( chains )




def LARGEARCS( arcs, getTail, getHead, distance, **kwargs ) :
    """
    obtain a stacker crane tour
    """
    graph = nx.MultiDiGraph()
    arcgen = itertools.count()
        
    # add original arcs to graph
    for i, arc in enumerate( arcs ) :
        tail = (i,TAIL)
        head = (i,HEAD)
        key = (i,ORIGINAL)
        graph.add_edge( tail, head, key )
        
    # compute cycles arc and add to graph
    cycles = CYCLEFACTOR( arcs, getTail, getHead, distance, **kwargs )
    for cycle in cycles :
        for i, j in cycle.iteritems() :
            tail = (i,HEAD)
            head = (j,TAIL)
            key = arcgen.next(), CYCLE
            graph.add_edge( tail, head, key )
            
    # construct super nodes, compute tree arcs and add to graph
    def cyclepoints( u ) :
        for i in cycles[u] :
            arc = arcs[i]
            yield getTail( arc ), (i,TAIL)
            yield getHead( arc ), (i,HEAD)
            
    def cycleconn( u, v ) :
        options = ( ( distance(x,y) + distance(y,x), tagx, tagy )
                    for x, tagx in cyclepoints(u)
                    for y, tagy in cyclepoints(v) )
        return min( options )
        
    #print len(cycles)
    cyclegraph = nx.Graph()
    for u, v in itertools.combinations( xrange( len(cycles) ), 2 ) :
        wt, tagx, tagy = cycleconn( u, v )
        #print u, v, tagx, tagy
        cyclegraph.add_edge( u, v, weight=wt, ends=[ tagx, tagy ] )
        
    MST = nx.minimum_spanning_tree( cyclegraph )
    for _,__, data in MST.edges_iter( data=True ) :
        tagx, tagy = data.get('ends')
        frwd = arcgen.next(), TREE
        bkwd = arcgen.next(), TREE
        
        graph.add_edge( tagx, tagy, frwd )
        graph.add_edge( tagy, tagx, bkwd )
        
    tour = []
    euler = nx.eulerian_circuit( graph )
    traversal = [ edge for edge in euler ]
    for edge in traversal :
        data = graph.get_edge_data( *edge )
        #print data
        arcidx, arctype = key = data.iterkeys().next()
        if arctype == ORIGINAL : tour.append( arcidx )
        u, v = edge
        graph.remove_edge( u, v, key )
        
    tour = { i : j for i, j in zip( tour, tour[1:] + tour[:1] ) }
    return tour







def CHAINCOST( tour, arcs, getTail, getHead, distance ) :
    total = 0.
    
    for i, j in tour.iteritems() :
        arci = arcs[i]
        arcj = arcs[j]
        
        # destination to destination
        total += distance( getHead(arci), getTail(arcj) )
        total += distance( getTail(arcj), getHead(arcj) )
        
    return total
    
TOURLENGTH = CHAINCOST      # the same thing, really
    
    
def WALKCOST( walk, arcs, getTail, getHead, distance ) :
    if not len( walk ) > 0 :
        return 0.
    
    def arccost( i ) :
        arc = arcs[i]
        return distance( getTail(arc), getHead(arc) )
    def edgecost( i, j ) :
        arci, arcj = arcs[i], arcs[j]
        return distance( getHead(arci), getTail(arcj) )
    
    i = walk[0]
    total = arccost( i )
    for i, j in zip( walk[:-1], walk[1:] ) :
        total += edgecost( i, j )
    total += arccost( j )
    
    return total
    
    
    
def SPLITTOUR( tour, N, arcs, getTail, getHead, distance, start=None ) :
    """
    tour should be a single cyclic chain
    agents is a dictionary with agents as keys and locations as values
    """
    walk = CYCLE2WALK( tour )
    #print walk
    
    def arccost( i ) :
        arc = arcs[i]
        return distance( getTail(arc), getHead(arc) )
    
    def edgecost( i, j ) : 
        arci, arcj = arcs[i], arcs[j]
        return distance( getHead(arci), getTail(arcj) )
    
    
    tourlen = TOURLENGTH( tour, arcs, getTail, getHead, distance )
    share = float(tourlen) / N
    #
    seqs = []
    for k in xrange(N) :
        frag = []
        fraglen = 0.
        
        j = None
        while len( walk ) > 0 :
            i = walk[0]
            
            c = 0.
            if not j is None : c += edgecost( j, i )
            if fraglen + c >= share : break
            
            walk.pop(0)
            frag.append( i )
            c += arccost( i )
            fraglen += c
            #
            j = i
            
        seqs.append( frag )
        
    return seqs




def ASSIGNFRAGS( seqs, agents, arcs, getTail, getHead, distance ) :
    # this should be easier
    graph = nx.Graph()
    
    # for each agent and its location
    for agent, x in agents.iteritems() :
        def options( seq ) :
            for startidx, arcidx in enumerate( seq ) :
                arc = arcs[arcidx]
                yield distance( x, getTail( arc ) ), startidx
                
        # for each sequence and its index
        for i, seq in enumerate( seqs ) :
            # get the optimal starting index in the sequence
            w, startidx = min( options( seq ) )
            # add a semantic edge to the graph
            graph.add_edge( (i,WALK), (agent,AGENT), weight=-w,
                        startidx=startidx )
            
    # get the min total cost matching (minimax would be better...)
    match = nx.max_weight_matching( graph, maxcardinality=True )
    
    assign = {}
    for agent in agents :
        i,_ = match[(agent,AGENT)]
        data = graph.get_edge_data( (i,WALK), (agent,AGENT) )
        startidx = data.get('startidx')
        
        seq = seqs[i]
        assign[agent] = seq[startidx:] + seq[:startidx]
    
    return assign
    
    
    
    
    
    

def kLARGEARCS( arcs, agents, getTail, getHead, distance ) :
    # compute largearcs tour
    tour = LARGEARCS( arcs, getTail, getHead, distance )
    
    # split the tour into fragments
    n = len( agents )
    seqs = SPLITTOUR( tour, n, arcs, getTail, getHead, distance )
    
    # assign the fragments
    assign = ASSIGNFRAGS( seqs, agents, arcs, getTail, getHead, distance )
    return assign


    
    
    
    
if __name__ == '__main__' :
    import sys, time
    import numpy as np
    
    import matplotlib as mpl
    if False :
        mpl.rcParams['ps.useafm'] = True
        mpl.rcParams['pdf.use14corefonts'] = True
        mpl.rcParams['text.usetex'] = True
        
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import axes3d
    from matplotlib.collections import PolyCollection
    
    plt.close('all')
    
    
    
    class Arc :
        def __init__(self, orig, dest ) :
            self.origin = orig
            self.destination = dest
            
        def getOrigin(self) : return self.origin
        def getDestination(self) : return self.destination
        
        def __repr__(self) :
            print '<%s,%s>' % ( repr(self.origin), repr(self.destination) )
    
    
    sample = lambda : np.random.rand(2)
    def distance( x, y ) : return np.linalg.norm( y - x )
    
    arcs = [ Arc( sample(), sample()) for i in xrange(100) ]
    getHead = lambda arc : arc.getOrigin()
    getTail = lambda arc : arc.getDestination()
    
    agents = { i : sample() for i in xrange(5) }
    
    #cycfact = CYCLEFACTOR( arcs, getTail, getHead, distance )
    
    # since everything's Euclidean...
    inst = nx.DiGraph()
    for k, arc in enumerate( arcs ):
        inst.add_edge( (k,TAIL), (k,HEAD) )
        
    pos = {}
    for k, arc in enumerate( arcs ) :
        pos[(k,TAIL)] = getTail( arc )
        pos[(k,HEAD)] = getHead( arc )
        
    plt.figure()
    nx.draw( inst, pos=pos )
    
    if False :
        tour = LARGEARCS( arcs, getTail, getHead, distance )
        tourgraph = nx.DiGraph()
        for i, j in tour.iteritems() :
            tourgraph.add_edge( (i,HEAD), (j,TAIL) )
            tourgraph.add_edge( (j,TAIL), (j,HEAD) )
            
        plt.figure()
        nx.draw( tourgraph, pos=pos )
    
    
    walks = kLARGEARCS( arcs, agents, getTail, getHead, distance )
    costs = { i : WALKCOST( walk, arcs, getTail, getHead, distance )
             for i, walk in walks.iteritems() } 
    
    if True :
        for agent, x in agents.iteritems() :
            pos[(agent,AGENT)] = x
        
        
        walkgraphs = {}
        for agent in agents :
            walk = walks[agent]
            if not len( walk ) > 0 : continue
            
            walkgraph = nx.DiGraph()
            i = walk[0]
            walkgraph.add_edge( (agent,AGENT), (i,TAIL) )
            for i, j in zip( walk[:-1], walk[1:] ) :
                walkgraph.add_edge( (i,TAIL), (i,HEAD) )
                walkgraph.add_edge( (i,HEAD), (j,TAIL) )
            # this should show the last one
            walkgraph.add_edge( (j,TAIL), (j,HEAD) )
            
            walkgraphs[agent] = walkgraph
            
        plt.figure()
        for agent in agents :
            nx.draw( walkgraphs[agent], pos=pos )
        
    

