"""
Network constructor and information propagator (ncip) contains functions for constructing the information network and entropy for the Reactome TopLevelPathways
"""

import os,sys
import numpy as np
import math
import networkx as nx
import pandas as pd
from scipy import sparse
import scipy.stats as st
from CoRe.cypher_commands import command_set

def construct_edges(inputs,outputs,property,neglect_ids=[],neglect_class=[]):
    """
    Converts reactions into a local graph topology.
    If a reaction is: A + B -> C + D + E,
    then A and B are two input nodes to three separate output nodes C, D, and E.
    If the reaction is chemically balanced then the edges are bidirectional,
    else edges are unidirectional.

    Parameters
    ----------
    inputs: neo4j result
        Result object from neo4j containing Reactome physical entities nodes that are inputs to an event.

    outputs: neo4j result
        object from neo4j containing Reactome physical entities nodes that are outputs to an event.

    property: string
        Property name of the Reactome physical entity that is used to identify the nodes.

    neglect_ids: list
        stId of Reactome physical entities that should be neglected from the network.

    neglect_class: list
        schemaClass of Reactome physical entities that should be neglected from the network.

    Returns
    -------
    edge_nodes: dict
        Dict containing the fields ['inputs','outputs','input_names','output_names'].
        'inputs' and 'outputs' for the input and output node property values of an edge.
        'input names' and 'output_names' are tuples (displayName,schemaClass) with values from
        the Reactome Data Schema, https://reactome.org/content/schema/PhysicalEntity.
    """

    edge_nodes = {}
    edge_nodes['inputs'] = []            # Property value as input node id
    edge_nodes['outputs'] = []           # Property value as output node id
    edge_nodes['input_names'] = []       # Physical entity name of input nodes
    edge_nodes['output_names'] = []      # Physical entity name of output nodes

    node_in_names = []
    node_out_names = []

    for res_in in inputs:
        node_in = res_in[property]

        if node_in not in neglect_ids and res_in['schemaClass'] not in neglect_class:

            for res_out in outputs:
                node_out = res_out[property]

                if node_out not in neglect_ids  and res_out['schemaClass'] not in neglect_class:
                    edge_nodes['inputs'].append(node_in)
                    edge_nodes['outputs'].append(node_out)
                    edge_nodes['input_names'].append((res_in['displayName'],res_in['schemaClass']))
                    edge_nodes['output_names'].append((res_out['displayName'],res_out['schemaClass']))

    return edge_nodes

def construct_exo_edges(inputs,property,neglect_ids=[],neglect_class=[]):
    """
    Converts reactions into a local graph topology.
    If a reaction is: A + B -> C + D + E,
    then A and B are two input nodes to three separate output nodes C, D, and E.
    If the reaction is chemically balanced then the edges are bidirectional,
    else edges are unidirectional.

    Parameters
    ----------
    inputs: neo4j result
        Result object from neo4j containing Reactome physical entities nodes that are inputs to an event.

    property: string
        Property name of the Reactome physical entity that is used to identify the nodes.

    neglect_ids: list
        stId of Reactome physical entities that should be neglected from the network.

    neglect_class: list
        schemaClass of Reactome physical entities that should be neglected from the network.

    Returns
    -------
    edge_nodes: dict
        Dict containing the fields ['inputs','outputs','input_names','output_names'].
        'inputs' and 'outputs' for the input and output node property values of an edge.
        'input names' and 'output_names' are tuples (displayName,schemaClass) with values from
        the Reactome Data Schema, https://reactome.org/content/schema/PhysicalEntity.
    """

    edge_nodes = {}
    edge_nodes['inputs'] = []            # Property value as input node id
    edge_nodes['outputs'] = []           # Property value as output node id
    edge_nodes['input_names'] = []       # Physical entity name of input nodes
    edge_nodes['output_names'] = []      # Physical entity name of output nodes

    node_in_names = []
    node_out_names = []

    for res_in in inputs:
        node_in = res_in[property]

        if node_in not in neglect_ids and res_in['schemaClass'] not in neglect_class:
            edge_nodes['inputs'].append(node_in)
            edge_nodes['outputs'].append(node_in)
            edge_nodes['input_names'].append((res_in['displayName'],res_in['schemaClass']))
            edge_nodes['output_names'].append((res_in['displayName'],res_in['schemaClass']))

    return edge_nodes

def state_entropy(state,sources):
    """
    Computes the entropy of the network state with respect to the maximum entropy state.
    For binary valued network nodes, the state is: {P(v=1),P(v=0)}, where v=1 means the physical
    entity at high concentration (or copy number). Similary, v=0 means the physical entity is
    at low concentration (or copy number). The total number of nodes in the network is n and the
    number of states for each node is m.

    Nodes that act as sources are neglected from the entropy calculation, to determine only
    the effect of the information propagation through the network.

    Parameters
    ----------
    state: array_like
        A [mn,1] numpy array.

    sources: list
        Node numbers that are sources in the network.

    Returns
    -------
    H: float
        Entropy of the network state in bits.
    """

    H = 0.0

    maxEnt_state = np.array([0.5,0.5])

    for i in range(0,int(state.shape[0]/2)):
        if i not in sources:
            H += st.entropy(state[2*i:2*i+2,0],maxEnt_state,base=2)

    return H

def get_final_state(C,X,sources):
    """
    Computes the final state of the network using the initial state X and the information transfer
    matrix Tmat. The state of the sources are fixed during the information transfer. The number of
    nodes in the network is n and the number of states for each node is m.

    Parameters
    ----------
    C: numpy matrix
        A [mn,mn] numpy array.

    X: array_like
        A [mn,1] numpy array.

    sources: list
        Node numbers that are sources in the network, whose probabilistic state, e.g. {P[v=1],P[v=0]}
        has to be fixed during information propagation.

    Returns
    -------
    X: array_like
        Stationary state of the network.

    steps: int
        Number of steps it took to reach the stationary state of the network.
    """
    tol = 1.0
    steps = 0

    while tol>1e-6:
        X_old = np.array(X, copy=True)
        X = C.dot(X)

        for s in sources:
            X[2*s:2*(s+1),:] = X_old[2*s:2*(s+1),:]

        tol = np.linalg.norm(X-X_old,2)

        steps += 1

    return X, steps

def make_graph(node_data,edge_data):
    """
    Creates a networkx graph using the connecitivity data and assigns a color value to each edge
    depending the subpathway they belong to.

    Parameters
    ----------
    node_data: DataFrame
        Contains the node ids (as Reactome stId) and attributes (as Reactome displayName) for the network.

    edge_data: DataFrame
        Contains the input and the outputs that make each edge of the graph.

    Returns
    -------
    G_d: networkx DiGraph
        NetworkX directed graph.
    """
    G_d = nx.DiGraph()
    edge_category = {}
    edge_names = {}
    edge_model = {}

    #for n1,n2,prop,sc,name in zip(edge_data['input'],edge_data['output'],edge_data['category'],edge_data['schemaClass'],edge_data['name']):
    for i in range(0,len(edge_data['input'])):
        if G_d.has_edge(edge_data['input'][i],edge_data['output'][i]):
            pass
        else:
            G_d.add_edge(edge_data['input'][i],edge_data['output'][i])

            edge_names[(edge_data['input'][i],edge_data['output'][i])] = edge_data['name'][i]

            if 'binding' in edge_data['category'][i]:
                category_name = 'binding'
            elif 'Negative' in edge_data['schemaClass'][i]:
                category_name = 'negative regulation'
            elif 'Positive' in edge_data['schemaClass'][i]:
                category_name = 'positive regulation'
            elif 'phorylat' in edge_data['name'][i]:
                category_name = 'phosphorylation'
            elif 'olymeris' in edge_data['schemaClass'][i]:
                category_name = 'polymerisation'
            else:
                category_name = 'other'

            edge_category[(edge_data['input'][i],edge_data['output'][i])] = category_name

            if 'Inhibit' in edge_data['name'][i] or 'Negative' in edge_data['schemaClass'][i]:
            #if 'Exocytosis' in edge_data['name'][i] or 'Inhibit' in edge_data['name'][i] or 'Negative' in edge_data['schemaClass'][i]:
                edge_model[(edge_data['input'][i],edge_data['output'][i])] = 'down'
            else:
                edge_model[(edge_data['input'][i],edge_data['output'][i])] = 'up'

            #if edge_data['category'][i]=='binding' or edge_data['category'][i]=='dissociation' or edge_data['category'][i]=='PPI':
            #if prop!='PositiveRegulation' and prop!='NegativeRegulation':
            if edge_data['schemaClass'][i]=='Reaction' or edge_data['category'][i]=='PPI':# and 'Exocytosis' not in edge_data['name'][i]:
                G_d.add_edge(edge_data['output'][i],edge_data['input'][i])

                edge_names[(edge_data['output'][i],edge_data['input'][i])] = edge_data['name'][i]
                edge_category[(edge_data['output'][i],edge_data['input'][i])] = category_name
                edge_model[(edge_data['output'][i],edge_data['input'][i])] = edge_model[(edge_data['input'][i],edge_data['output'][i])]

    nx.set_edge_attributes(G_d, edge_category, "category")
    nx.set_edge_attributes(G_d, edge_names, "name")
    nx.set_edge_attributes(G_d, edge_model, "channel")

    node_names = {}
    node_genes = {}
    node_w_seq = {}

    for id,n,g,s in zip(node_data['node'],node_data['name'],node_data['reference gene'],node_data['genome-encoded']):
        node_names[id] = n
        node_genes[id] = g
        node_w_seq[id] = s

    nx.set_node_attributes(G_d, node_names, "name")
    nx.set_node_attributes(G_d, node_genes, "gene")
    nx.set_node_attributes(G_d, node_w_seq, "sequenced")

    return G_d

def construct_C(G_d,rho=0,h=1):
    """
    Creates the information transition matrix for the graph using the adjacency matrix and
    interference due to multiple information sources.

    Parameters
    ----------
    G_d: networkx DiGraph
        Graph of a Reactome TopLevelPathway consisting of ReactionLikeEvents and regulations.

    rho: float
        Communication error in each edge.

    Returns
    -------
    C: numpy matrix
        Information transition matrix of the network.
    """

    code_length = int(2**h)

    up_regulation, down_regulation = np.eye((code_length)), np.rot90(np.eye((code_length)))

    up_regulation *= (1 - (code_length-1)*rho)*up_regulation
    down_regulation *= (1 - (code_length-1)*rho)*down_regulation

    up_regulation += rho*(np.ones((code_length,code_length)) - np.eye(code_length))
    down_regulation += rho*(np.ones((code_length,code_length)) - np.rot90(np.eye(code_length)))

    n = len(G_d.nodes())
    C = np.zeros(shape=(code_length*n,code_length*n))
    A = np.transpose(nx.to_numpy_matrix(G_d))

    node_list = list(G_d.nodes)

    for n1 in range(0,A.shape[0]):
        in_degree = G_d.in_degree(node_list[n1])
        for n2 in range(0,A.shape[0]):
            if A[n1,n2]>0:
                if G_d.get_edge_data(node_list[n2],node_list[n1])['channel']=='down':
                    C[code_length*n1:code_length*(n1+1),code_length*n2:code_length*(n2+1)] = (1.0/in_degree)*down_regulation
                else:
                    C[code_length*n1:code_length*(n1+1),code_length*n2:code_length*(n2+1)] = (1.0/in_degree)*up_regulation

                # if G_d.get_edge_data(node_list[n2],node_list[n1])['category']=='NegativeRegulation':
                #     C[2*n1:2*n1+2,2*n2:2*n2+2] = (1.0/in_degree)*down_regulation
                # else:
                #     C[2*n1:2*n1+2,2*n2:2*n2+2] = (1.0/in_degree)*up_regulation

    # Apply bounday condition
    for n in range(0,len(node_list)):
        if G_d.in_degree(node_list[n])==0:
            C[code_length*n:code_length*(n+1),code_length*n:code_length*(n+1)] = np.identity(code_length)

    C_sparse = sparse.csr_matrix(C)

    return C_sparse

def add_edge_and_nodes(edge_nodes,edge_data,node_data,session,command_set,edge_info,edge_log):
    """
    Adds information about an edge to the edge_data and node_data.

    Parameters
    ----------
    edge_nodes: dict
        Contains the node ids (as Reactome stId) and attributes (as Reactome displayName) for the network.

    edge_data: DataFrame
        Contains the input and the outputs that make each edge of the graph.

    node_data: DataFrame
        Contains the node ids (as Reactome stId) and attributes (as Reactome displayName) for the network.

    session: neo4j session
        To query the neo4j reactome database

    command_set: dict
        Set of commands to query the neo4j reactome database.

    edge_info: list
        Contains Reactome stId, superpathways, category, and schemaClass of the ReactionLikeEvent.
    """

    r = edge_info[0]
    superpathways = edge_info[1]
    c = edge_info[2]
    s = edge_info[3]

    for i in range(0,len(edge_nodes['inputs'])):
        if (edge_nodes['inputs'][i]+','+edge_nodes['outputs'][i]) in edge_log:
            pass
        else:
            edge_log.append(edge_nodes['inputs'][i]+','+edge_nodes['outputs'][i])

            edge_data['input'].append(edge_nodes['inputs'][i])
            edge_data['output'].append(edge_nodes['outputs'][i])
            edge_data['reaction'].append(r)
            edge_data['superpathways'].append(superpathways)
            edge_data['category'].append(c)
            edge_data['schemaClass'].append(s)

            if edge_nodes['inputs'][i] not in node_data['node']:
                node_data['node'].append(edge_nodes['inputs'][i])
                node_data['name'].append(edge_nodes['input_names'][i][0])
                node_data['schemaClass'].append(edge_nodes['input_names'][i][1])

                if edge_nodes['input_names'][i][1]=='EntityWithAccessionedSequence':
                    node_data['genome-encoded'].append(True)

                    ref_gene_command = command_set['reference_gene'].replace('#',edge_nodes['inputs'][i])
                    ref_gene_result = session.run(ref_gene_command)

                    for ref_res in ref_gene_result:
                        gene_dbId = str(ref_res[0]['dbId'])

                        try:
                            gene_name = ref_res[0]['geneName'][0]
                        except TypeError:
                            gene_name = ref_res[0]['name'][0]

                    node_data['reference gene'].append(gene_name)
                else:
                    node_data['genome-encoded'].append(False)
                    node_data['reference gene'].append('None')

            if edge_nodes['outputs'][i] not in node_data['node']:
                node_data['node'].append(edge_nodes['outputs'][i])
                node_data['name'].append(edge_nodes['output_names'][i][0])
                node_data['schemaClass'].append(edge_nodes['output_names'][i][1])

                if edge_nodes['output_names'][i][1]=='EntityWithAccessionedSequence':
                    node_data['genome-encoded'].append(True)

                    ref_gene_command = command_set['reference_gene'].replace('#',edge_nodes['outputs'][i])
                    ref_gene_result = session.run(ref_gene_command)

                    for ref_res in ref_gene_result:
                        gene_dbId = str(ref_res[0]['dbId'])

                        try:
                            gene_name = ref_res[0]['geneName'][0]
                        except TypeError:
                            gene_name = ref_res[0]['name'][0]

                    node_data['reference gene'].append(gene_name)
                else:
                    node_data['genome-encoded'].append(False)
                    node_data['reference gene'].append('None')

def add_edges(edge_nodes,edge_data,node_data,session,reaction_info,edge_log):
    """
    Adds information about an edge to the edge_data and node_data.

    Parameters
    ----------
    edge_nodes: dict
        Contains the node ids (as Reactome stId) and attributes (as Reactome displayName) for the network.

    edge_data: DataFrame
        Contains the input and the outputs that make each edge of the graph.

    node_data: DataFrame
        Contains the node ids (as Reactome stId) and attributes (as Reactome displayName) for the network.

    session: neo4j session
        To query the neo4j reactome database

    command_set: dict
        Set of commands to query the neo4j reactome database.

    reaction_info: list
        Row entry of the ReactionLikeEvents dataframe containing stId, category, schemaClass, and name.
    """

    ref_gene_info = {}
    ref_gene_keys = []

    for i in range(0,len(edge_nodes['inputs'])):
        if edge_nodes['input_names'][i][1]=='EntityWithAccessionedSequence':
            if edge_nodes['inputs'][i] not in ref_gene_keys:
                ref_gene_command = command_set['reference_gene'].replace('#',edge_nodes['inputs'][i])
                result = session.run(ref_gene_command)
                refgene_inputs = [res[0] for res in result]
                #edge_log_name_in = str(refgene_inputs[0]['dbId'])

                try:
                    gene_name = refgene_inputs[0]['geneName'][0]
                except TypeError:
                    gene_name = refgene_inputs[0]['name'][0]

                edge_log_name_in = gene_name

                ref_gene_keys.append(edge_nodes['inputs'][i])
                ref_gene_info[edge_nodes['inputs'][i]] = (edge_log_name_in,gene_name)

        if edge_nodes['output_names'][i][1]=='EntityWithAccessionedSequence':
            if edge_nodes['outputs'][i] not in ref_gene_keys:
                ref_gene_command = command_set['reference_gene'].replace('#',edge_nodes['outputs'][i])
                result = session.run(ref_gene_command)
                refgene_outputs = [res[0] for res in result]
                #edge_log_name_out = str(refgene_outputs[0]['dbId'])

                try:
                    gene_name = refgene_outputs[0]['geneName'][0]
                except TypeError:
                    gene_name = refgene_outputs[0]['name'][0]

                edge_log_name_out = gene_name

                ref_gene_keys.append(edge_nodes['outputs'][i])
                ref_gene_info[edge_nodes['outputs'][i]] = (edge_log_name_out,gene_name)

    #print('Reference gene search completed')

    for i in range(0,len(edge_nodes['inputs'])):
        edge_log_name_in = edge_nodes['inputs'][i]
        edge_log_name_out = edge_nodes['outputs'][i]

        if edge_nodes['input_names'][i][1]=='EntityWithAccessionedSequence':
            edge_log_name_in = ref_gene_info[edge_nodes['inputs'][i]][0]

        if edge_nodes['output_names'][i][1]=='EntityWithAccessionedSequence':
            edge_log_name_out = ref_gene_info[edge_nodes['outputs'][i]][0]

        edge_tag = edge_log_name_in+','+edge_log_name_out

        if edge_tag in edge_log:
            pass
        else:
            edge_log.append(edge_tag)

            edge_data['input'].append(edge_log_name_in)
            edge_data['output'].append(edge_log_name_out)
            edge_data['reaction'].append(reaction_info['reaction'])
            edge_data['category'].append(reaction_info['category'])
            edge_data['schemaClass'].append(reaction_info['schemaClass'])
            edge_data['name'].append(reaction_info['name'])

            if edge_data['input'][-1] not in node_data['node']:
                node_data['node'].append(edge_data['input'][-1])
                node_data['name'].append(edge_nodes['input_names'][i][0])
                node_data['schemaClass'].append(edge_nodes['input_names'][i][1])

                if edge_nodes['input_names'][i][1]=='EntityWithAccessionedSequence':
                    node_data['genome-encoded'].append(True)

                    node_data['reference gene'].append(ref_gene_info[edge_nodes['inputs'][i]][1])
                else:
                    node_data['genome-encoded'].append(False)
                    node_data['reference gene'].append('None')
                    #node_data['reference stId'].append('None')

            if edge_data['output'][-1] not in node_data['node']:
                node_data['node'].append(edge_data['output'][-1])
                node_data['name'].append(edge_nodes['output_names'][i][0])
                node_data['schemaClass'].append(edge_nodes['output_names'][i][1])

                if edge_nodes['output_names'][i][1]=='EntityWithAccessionedSequence':
                    node_data['genome-encoded'].append(True)

                    node_data['reference gene'].append(ref_gene_info[edge_nodes['outputs'][i]][1])
                else:
                    node_data['genome-encoded'].append(False)
                    node_data['reference gene'].append('None')
                    #node_data['reference stId'].append('None')

def save_network(node_data,edge_data,pathway_nametag,network_type):
    """
    Adds information about an edge to the edge_data and node_data.

    Parameters
    ----------
    node_data: dict
        Contains the node ids (as Reactome stId) and attributes (as Reactome displayName) for the network.

    edge_data: DataFrame
        Contains the input and the outputs that make each edge of the graph.

    pathway_name: string
        Name of the reactome top level pathway.

    network_type: string
        To query the neo4j reactome database
    """
    df_e = pd.DataFrame(edge_data)

    df_e.to_pickle(pathway_nametag+'_'+network_type+'-edges.pkl')
    df_e.to_csv(pathway_nametag+'_'+network_type+'-edges.csv',index=None)

    df_n = pd.DataFrame(node_data)
    df_n.to_pickle(pathway_nametag+'_'+network_type+'-nodes.pkl')
    df_n.to_csv(pathway_nametag+'_'+network_type+'-nodes.csv',index=None)

    G_d = make_graph(node_data,edge_data)

    print("Communication network "+pathway_nametag+"-"+network_type+".gml created.")

    nx.write_gml(G_d,pathway_nametag+"-"+network_type+".gml")
