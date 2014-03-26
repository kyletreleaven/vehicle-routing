
import itertools
import networkx as nx




def SHORTESTEDGE( u, v, graph, weight_attr='weight', data=False ) :
    edge_data = graph.get_edge_data(u,v)
    options = ( ( val.get(weight_attr), key )
                for key, val in edge_data.iteritems() )
    _, edgekey = min( options )
    #edgekey = min( edge_data, key=lambda k : edge_data[k].get(weight_attr) )
    
    if not data :
        return edgekey
    else :
        return edgekey, edge_data[edgekey] 
    

def PATHLEN( path, graph, weight_attr='weight' ) :
    """ just assume graph is a multi-graph """
    total = 0.
    for u, v in zip( path[:-1], path[1:] ) :
        _, data = SHORTESTEDGE(u,v, graph, weight_attr, data=True )
        total += data.get(weight_attr)
        
    return total



def ODDNODES( graph ) :
    return [ u for u, d in graph.degree_iter() if d % 2 == 1 ]





ORIGINAL = 0
EXTRA = 1

def ROUTEINSPECTION( graph, weight_attr='weight' ) :
    """
    graph should be a weighted multi-graph (undirected)
    """
    workgraph = nx.MultiGraph()
    
    # add an original copy of every edge to the workgraph
    for u, v, key in graph.edges_iter( keys=True ) :
        workgraph.add_edge( u, v, (ORIGINAL,key) )
        
    # if non-Eulerian, add joining paths between odd-degree vertices
    T = ODDNODES( graph )
    # metric graph
    met = nx.Graph()
    for s in T :
        for t in T :
            if t == s : continue
            path = nx.shortest_path( graph, s, t, weight=weight_attr )
            pathlen = PATHLEN( path, graph, weight_attr )
            met.add_edge( s, t, weight=-pathlen, path=path )    # want min cost
    match = nx.max_weight_matching( met, maxcardinality=True )
    
    extras = itertools.count()
    while len( match ) > 0 :
        s, t = match.iteritems().next()
        del match[s]
        del match[t]    # have to kill both ends
        
        path = met.get_edge_data(s,t).get('path')
        for u,v in zip( path[:-1], path[1:] ) :
            #edgekey = SHORTESTEDGE(u,v, graph, weight_attr )
            #idx = len( extras )
            #extras.append(edgekey)
            idx = extras.next()
            workgraph.add_edge(u,v, (EXTRA,idx) )
            
    #return workgraph

    # traverse
    walk = []
    
    for u, v in nx.eulerian_circuit( workgraph ) :
        edge_data = workgraph.get_edge_data(u,v)
        which, key = datakey = edge_data.iterkeys().next()
        workgraph.remove_edge(u,v, datakey )
        
        if which == ORIGINAL :
            edge_key = key
        elif which == EXTRA :
            edge_key = SHORTESTEDGE(u,v, graph, weight_attr )
        
        if not len( walk ) > 0 :
            walk.append(u)
        walk.extend([edge_key,v])
        
    return walk






if __name__ == '__main__' :
    
    
    g = nx.MultiGraph()
    g.add_edge(0,1,'A',weight=1.)
    g.add_edge(0,1,'B',weight=1.)
    g.add_edge(1,2,'C',weight=1.)
    
    #g.add_path([0,1,2,3], weight=1. )
    #g.add_path([1,2])
    #g.add_path([0,1], weight=1. )
    
    walk = ROUTEINSPECTION(g)




