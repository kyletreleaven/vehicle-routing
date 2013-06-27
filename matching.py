
import numpy as np
import networkx as nx




def max_weight_matching( graph, weight=None, maxcardinality=False ) :
    """
    a more featureful max_weight_matching algorithm (just simple pre/post-processing);
    accepts graph, or multigraph; returns the same type!
    to get a vertex u's match, use neighbors(u)[0] instead of the dictionary;
    comments: the algorithm makes a copy of the input
    """
    assert isinstance( graph, nx.Graph )
    MULTIGRAPH = isinstance( graph, nx.MultiGraph )
    multigraph = nx.MultiGraph( graph )
    
    # PREPROCESS
    workgraph = nx.Graph()
    for u,v, key, new_data in multigraph.edges_iter( data=True, keys=True ) :
        if weight == None :
            if workgraph.has_edge( u, v ) : continue
            workgraph.add_edge( u, v, { 'edge' : (u,v,key) } )
            
        else :
            new_value = new_data[weight]
            if workgraph.has_edge( u, v ) :
                data = workgraph.get_edge_data( u, v )
                if new_value < data[weight] :
                    data['weight'] = new_value
                    data['edge'] = (u,v,key)
            else :
                data = { 'weight' : new_value, 'edge' : (u,v,key) }
                workgraph.add_edge( u, v, data )
                
    # RUN networkx ALGORITHM
    match = nx.max_weight_matching( workgraph, maxcardinality )
    
    # POSTPROCESS
    if MULTIGRAPH :
        res = nx.MultiGraph()
    else :
        res = nx.Graph()
        
    while len( match ) > 0 :
        u, v = match.popitem()
        del match[v]
        
        i,j,key = workgraph.get_edge_data( u, v )['edge']
        data = multigraph.get_edge_data( i,j,key )
        if MULTIGRAPH :
            res.add_edge( i,j,key, data )
        else :
            res.add_edge( i,j, data )
            
    return res
    
def match_to_dict( match ) :
    """
    in case you simply *must* to have a match in dictionary form;
    assumes a proper input, e.g. from max_weight_matching above;
    thus, no input checking; use at own risk!
    """
    res = {}
    for i,j in match.edges_iter() :
        res[i] = j
        res[j] = i
    return res




