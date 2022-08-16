"""
Network constructor and information propagator (ncip) contains functions for constructing the information network and entropy for the Reactome TopLevelPathways.
"""

import os,sys
import numpy as np
import math
import networkx as nx
import pandas as pd
from scipy import sparse
import scipy.stats as st
from CoRe.cypher_commands import command_set
from tqdm import tqdm
#import BA_C.BA as BA

class ncip:
    def __init__(self,_property=None):
        self.edge_data = {}
        self.edge_data['input'] = []
        self.edge_data['output'] = []
        self.edge_data['reaction'] = []
        self.edge_data['name'] = []
        #self.edge_data['category'] = []
        self.edge_data['schemaClass'] = []
        self.edge_data['module'] = []

        self.property = _property

        self.node_data = {}
        self.node_data['node'] = []
        self.node_data['name'] = []
        self.node_data['genome-encoded'] = []
        self.node_data['schemaClass'] = []
        self.node_data['reference gene'] = []
        self.edge_log = []

    def load_data(self,_edge_data,_node_data):
        self.edge_data = _edge_data
        self.node_data = _node_data

    def load_graph(self,filename):
        self.G_d = nx.read_gml(filename)

    def construct_edges(self,inputs,outputs,neglect_ids=[],neglect_class=[]):
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

        self.edge_nodes = {}
        self.edge_nodes['inputs'] = []            # Property value as input node id
        self.edge_nodes['outputs'] = []           # Property value as output node id
        self.edge_nodes['input_names'] = []       # Physical entity name of input nodes
        self.edge_nodes['output_names'] = []      # Physical entity name of output nodes

        node_in_names = []
        node_out_names = []

        for res_in in inputs:
            node_in = res_in[self.property]

            if node_in not in neglect_ids and res_in['schemaClass'] not in neglect_class:

                for res_out in outputs:
                    node_out = res_out[self.property]

                    if node_out not in neglect_ids  and res_out['schemaClass'] not in neglect_class:
                        self.edge_nodes['inputs'].append(node_in)
                        self.edge_nodes['outputs'].append(node_out)
                        self.edge_nodes['input_names'].append((res_in['displayName'],res_in['schemaClass']))
                        self.edge_nodes['output_names'].append((res_out['displayName'],res_out['schemaClass']))

    def construct_exo_edges(self,inputs,neglect_ids=[],neglect_class=[]):
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

        self.edge_nodes = {}
        self.edge_nodes['inputs'] = []            # Property value as input node id
        self.edge_nodes['outputs'] = []           # Property value as output node id
        self.edge_nodes['input_names'] = []       # Physical entity name of input nodes
        self.edge_nodes['output_names'] = []      # Physical entity name of output nodes

        node_in_names = []
        node_out_names = []

        for res_in in inputs:
            node_in = res_in[self.property]

            if node_in not in neglect_ids and res_in['schemaClass'] not in neglect_class:
                self.edge_nodes['inputs'].append(node_in)
                self.edge_nodes['outputs'].append(node_in)
                self.edge_nodes['input_names'].append((res_in['displayName'],res_in['schemaClass']))
                self.edge_nodes['output_names'].append((res_in['displayName'],res_in['schemaClass']))

    def add_edges(self,session,reaction_info):
        """
        Adds information about an edge to the edge_data and node_data.

        Parameters
        ----------
        session: neo4j session
            To query the neo4j reactome database

        command_set: dict
            Set of commands to query the neo4j reactome database.

        reaction_info: list
            Row entry of the ReactionLikeEvents dataframe containing stId, category, schemaClass, and name.
        """

        ref_gene_info = {}
        ref_gene_keys = []
        species_name = {}

        for i in range(0,len(self.edge_nodes['inputs'])):
            if self.edge_nodes['input_names'][i][1]=='EntityWithAccessionedSequence':
                if self.edge_nodes['inputs'][i] not in ref_gene_keys:
                    ref_gene_command = command_set['reference_gene'].replace('#',self.edge_nodes['inputs'][i])
                    result = session.run(ref_gene_command)
                    refgene_inputs = [res[0] for res in result]

                    spec_command = "MATCH(n:EntityWithAccessionedSequence {stId:'"+ self.edge_nodes['inputs'][i] +"'}) RETURN n.speciesName"
                    result = session.run(spec_command)
                    spec_name = [res[0] for res in result][0]

                    species_name[self.edge_nodes['inputs'][i]] = spec_name

                    try:
                        gene_name = refgene_inputs[0]['geneName'][0]
                    except TypeError:
                        gene_name = refgene_inputs[0]['name'][0]

                    edge_log_name_in = gene_name

                    ref_gene_keys.append(self.edge_nodes['inputs'][i])
                    ref_gene_info[self.edge_nodes['inputs'][i]] = (edge_log_name_in,gene_name)

            if self.edge_nodes['output_names'][i][1]=='EntityWithAccessionedSequence':
                if self.edge_nodes['outputs'][i] not in ref_gene_keys:
                    ref_gene_command = command_set['reference_gene'].replace('#',self.edge_nodes['outputs'][i])
                    result = session.run(ref_gene_command)
                    refgene_outputs = [res[0] for res in result]

                    spec_command = "MATCH(n:EntityWithAccessionedSequence {stId:'"+ self.edge_nodes['outputs'][i] +"'}) RETURN n.speciesName"
                    result = session.run(spec_command)
                    spec_name = [res[0] for res in result][0]

                    species_name[self.edge_nodes['outputs'][i]] = spec_name

                    try:
                        gene_name = refgene_outputs[0]['geneName'][0]
                    except TypeError:
                        gene_name = refgene_outputs[0]['name'][0]

                    edge_log_name_out = gene_name

                    ref_gene_keys.append(self.edge_nodes['outputs'][i])
                    ref_gene_info[self.edge_nodes['outputs'][i]] = (edge_log_name_out,gene_name)

        #print('Reference gene search completed')

        for i in range(0,len(self.edge_nodes['inputs'])):
            edge_log_name_in = self.edge_nodes['inputs'][i]
            edge_log_name_out = self.edge_nodes['outputs'][i]

            if self.edge_nodes['input_names'][i][1]=='EntityWithAccessionedSequence':
                try:
                    edge_log_name_in = ref_gene_info[self.edge_nodes['inputs'][i]][0]
                except KeyError:
                    pass

            if self.edge_nodes['output_names'][i][1]=='EntityWithAccessionedSequence':
                try:
                    edge_log_name_out = ref_gene_info[self.edge_nodes['outputs'][i]][0]
                except KeyError:
                    pass

            edge_tag = edge_log_name_in+','+edge_log_name_out

            if edge_tag in self.edge_log:
                pass
            else:
                self.edge_log.append(edge_tag)

                self.edge_data['input'].append(edge_log_name_in)
                self.edge_data['output'].append(edge_log_name_out)
                self.edge_data['reaction'].append(reaction_info['reaction'])
                #self.edge_data['category'].append(reaction_info['category'])
                self.edge_data['schemaClass'].append(reaction_info['schemaClass'])
                self.edge_data['name'].append(reaction_info['name'])
                self.edge_data['module'].append(reaction_info['module'])

                if self.edge_data['input'][-1] not in self.node_data['node']:
                    self.node_data['node'].append(self.edge_data['input'][-1])
                    self.node_data['name'].append(self.edge_nodes['input_names'][i][0])
                    self.node_data['schemaClass'].append(self.edge_nodes['input_names'][i][1])

                    if self.edge_nodes['input_names'][i][1]=='EntityWithAccessionedSequence':
                        if species_name[self.edge_nodes['inputs'][i]]=='Homo sapiens':
                            self.node_data['reference gene'].append(ref_gene_info[self.edge_nodes['inputs'][i]][1])
                            self.node_data['genome-encoded'].append(True)
                        else:
                            self.node_data['genome-encoded'].append(False)
                            self.node_data['reference gene'].append('None')
                    else:
                        self.node_data['genome-encoded'].append(False)
                        self.node_data['reference gene'].append('None')
                        #node_data['reference stId'].append('None')

                if self.edge_data['output'][-1] not in self.node_data['node']:
                    self.node_data['node'].append(self.edge_data['output'][-1])
                    self.node_data['name'].append(self.edge_nodes['output_names'][i][0])
                    self.node_data['schemaClass'].append(self.edge_nodes['output_names'][i][1])

                    if self.edge_nodes['output_names'][i][1]=='EntityWithAccessionedSequence':
                        if species_name[self.edge_nodes['outputs'][i]]=='Homo sapiens':
                            self.node_data['reference gene'].append(ref_gene_info[self.edge_nodes['outputs'][i]][1])
                            self.node_data['genome-encoded'].append(True)
                        else:
                            self.node_data['genome-encoded'].append(False)
                            self.node_data['reference gene'].append('None')
                    else:
                        self.node_data['genome-encoded'].append(False)
                        self.node_data['reference gene'].append('None')

    def make_graph(self):
        """
        Creates a networkx graph using the connecitivity data and assigns a color value to each edge
        depending the subpathway they belong to.
        """
        self.G_d = nx.DiGraph()
        edge_category = {}
        edge_names = {}
        edge_model = {}
        edge_module = {}
        edge_reaction = {}

        check1, check2 = 0, 0

        #for n1,n2,prop,sc,name in zip(edge_data['input'],edge_data['output'],edge_data['category'],edge_data['schemaClass'],edge_data['name']):
        for i in tqdm(range(0,len(self.edge_data['input']))):
            #if self.edge_data['input'][i]!=self.edge_data['output'][i]:
            if self.G_d.has_edge(self.edge_data['input'][i],self.edge_data['output'][i]):
                pass
            else:
                #if self.edge_data['input'][i]!=self.edge_data['output'][i]:
                #if (self.edge_data['input'][i]==self.edge_data['output'][i]):
                #    pass
                #else:
                self.G_d.add_edge(self.edge_data['input'][i],self.edge_data['output'][i])

                edge_names[(self.edge_data['input'][i],self.edge_data['output'][i])] = self.edge_data['name'][i]

                if 'binding' in self.edge_data['name'][i] or 'Binding' in self.edge_data['name'][i]:
                    category_name = 'binding'
                elif 'Negative' in self.edge_data['schemaClass'][i] or 'nhibit' in self.edge_data['name'][i]:
                    category_name = 'negative regulation'
                elif 'Positive' in self.edge_data['schemaClass'][i] or 'egulation' in self.edge_data['name'][i]:
                    category_name = 'positive regulation'
                elif 'phorylat' in self.edge_data['name'][i]:
                    category_name = 'phosphorylation'
                elif 'Expression'  in self.edge_data['name'][i]:
                    category_name = 'expression'
                elif self.edge_data['name'][i]=='PPI':
                    category_name = 'PPI'
                else:
                    category_name = 'other'

                edge_category[(self.edge_data['input'][i],self.edge_data['output'][i])] = category_name
                edge_module[(self.edge_data['input'][i],self.edge_data['output'][i])] = self.edge_data['module'][i]
                edge_reaction[(self.edge_data['input'][i],self.edge_data['output'][i])] = self.edge_data['reaction'][i]


                if 'nhibit' in self.edge_data['name'][i] or 'Negative' in self.edge_data['schemaClass'][i]:
                #if 'Exocytosis' in edge_data['name'][i] or 'Inhibit' in edge_data['name'][i] or 'Negative' in edge_data['schemaClass'][i]:
                    edge_model[(self.edge_data['input'][i],self.edge_data['output'][i])] = 'down'
                else:
                    edge_model[(self.edge_data['input'][i],self.edge_data['output'][i])] = 'up'

                #if True==False:
                #if edge_data['category'][i]=='binding' or edge_data['category'][i]=='dissociation' or edge_data['category'][i]=='PPI':
                #if prop!='PositiveRegulation' and prop!='NegativeRegulation':
                if self.edge_data['schemaClass'][i]=='Reaction' or self.edge_data['schemaClass'][i]=='PPI':# and 'Exocytosis' not in edge_data['name'][i]:
                    self.G_d.add_edge(self.edge_data['output'][i],self.edge_data['input'][i])

                    edge_names[(self.edge_data['output'][i],self.edge_data['input'][i])] = self.edge_data['name'][i]
                    edge_category[(self.edge_data['output'][i],self.edge_data['input'][i])] = category_name
                    edge_model[(self.edge_data['output'][i],self.edge_data['input'][i])] = edge_model[(self.edge_data['input'][i],self.edge_data['output'][i])]
                    edge_module[(self.edge_data['output'][i],self.edge_data['input'][i])] = self.edge_data['module'][i]
                    edge_reaction[(self.edge_data['output'][i],self.edge_data['input'][i])] = self.edge_data['reaction'][i]

        nx.set_edge_attributes(self.G_d, edge_category, "category")
        nx.set_edge_attributes(self.G_d, edge_names, "name")
        nx.set_edge_attributes(self.G_d, edge_model, "channel")
        nx.set_edge_attributes(self.G_d, edge_module, "module")
        nx.set_edge_attributes(self.G_d, edge_reaction, "reaction")

        node_names = {}
        node_genes = {}
        node_w_seq = {}
        node_types = {}

        for index,row in self.node_data.iterrows():
            id = row['node']

            node_names[id] = row['name']
            node_genes[id] = row['reference gene']
            node_w_seq[id] = row['genome-encoded']
            node_types[id] = row['schemaClass']

        nx.set_node_attributes(self.G_d, node_names, "name")
        nx.set_node_attributes(self.G_d, node_genes, "gene")
        nx.set_node_attributes(self.G_d, node_w_seq, "sequenced")
        nx.set_node_attributes(self.G_d, node_types, "class")

    def state_entropy(self,state,sources,reference_state=[]):
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

        maxEnt_state = (1.0/float(self.code_length))*np.ones(shape=(self.code_length,))

        if len(reference_state)==0:
            for i in range(0,int(state.shape[0]/self.code_length)):
                if i not in sources:
                    H += st.entropy(state[self.code_length*i:self.code_length*(i+1),0],maxEnt_state,base=2)
                    #print(state[self.code_length*i:self.code_length*(i+1),0])
        else:
            for i in range(0,int(state.shape[0]/self.code_length)):
                if i not in sources:
                    this_entropy = st.entropy(state[self.code_length*i:self.code_length*(i+1),0],reference_state[self.code_length*i:self.code_length*(i+1),0],base=2)
                    if this_entropy!=np.inf:
                        H += this_entropy
                    #print(state[self.code_length*i:self.code_length*(i+1),0])

        return H

    def entropy_drop_and_gain(self,state_1,state_2,reference_state,sources):
        """
        Computes the difference between two Kullback-Leibeler divergences, D_KL(state_1||reference_state), and
        D_KL(state_2||reference_state), for each unconstrained node in the network.

        Parameters
        ----------
        state: array_like
            A [mn,1] numpy array.

        sources: list
            Node numbers that are sources in the network.

        Returns
        -------
        drop: float
            Sum of all difference in the relative entropy values for all the nodes in the network where, state_2 > state_1.

        gain: float
            Sum of all difference in the relative entropy values for all the nodes in the network where, state_1 <= state_2.
        """

        H_drop = 0.0
        H_gain = 0.0

        maxEnt_state = (1.0/float(self.code_length))*np.ones(shape=(self.code_length,))

        for i in range(0,int(state_1.shape[0]/self.code_length)):
            if i not in sources:
                this_entropy_1 = st.entropy(state_1[self.code_length*i:self.code_length*(i+1),0],reference_state[self.code_length*i:self.code_length*(i+1),0],base=2)

                this_entropy_2 = st.entropy(state_2[self.code_length*i:self.code_length*(i+1),0],reference_state[self.code_length*i:self.code_length*(i+1),0],base=2)

                if this_entropy_1!=np.inf and this_entropy_2!=np.inf:
                    if this_entropy_1>=this_entropy_2:
                        H_gain += this_entropy_1 - this_entropy_2
                    else:
                        H_drop += this_entropy_1 - this_entropy_2
                    #print(state[self.code_length*i:self.code_length*(i+1),0])

        return H_drop, H_gain

    def get_final_state(self,X,sources):
        """
        Computes the final state of the network using the initial state X and the information transfer
        matrix Tmat. The state of the sources are fixed during the information transfer. The number of
        nodes in the network is n and the number of states for each node is m.

        Parameters
        ----------
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

        while tol>1e-4:
            X_old = np.array(X, copy=True)
            X = self.C_sparse.dot(X)

            for s in sources:
                X[2*s:2*(s+1),:] = X_old[2*s:2*(s+1),:]

            tol = np.linalg.norm(X-X_old,1)

            steps += 1

        return X, steps

    def construct_C(self,rho=0,h=1,neglect_modules=[]):
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

        self.code_length = int(2**h)
        m_rho = rho/(self.code_length - 1)

        up_regulation, down_regulation = np.eye((self.code_length)), np.rot90(np.eye((self.code_length)))

        up_regulation *= (1 - rho)*up_regulation
        down_regulation *= (1 - rho)*down_regulation

        up_regulation += m_rho*(np.ones((self.code_length,self.code_length)) - np.eye(self.code_length))
        down_regulation += m_rho*(np.ones((self.code_length,self.code_length)) - np.rot90(np.eye(self.code_length)))

        null_communicator = (1.0/float(self.code_length))*np.ones((self.code_length,self.code_length))

        # bao = BA()
        # bao.set_response(up_regulation)
        # print(bao.get_CC())
        #
        # bao.set_response(down_regulation)
        # print(bao.get_CC())

        n = len(self.G_d.nodes())
        C = np.zeros(shape=(self.code_length*n,self.code_length*n))
        A = np.transpose(nx.to_numpy_matrix(self.G_d))

        node_list = list(self.G_d.nodes)

        node_idx = {}

        for i in range(0,len(node_list)):
            node_idx[node_list[i]] = i

        # for n1 in range(0,A.shape[0]):
        #     in_degree = self.G_d.in_degree(node_list[n1])
        #     for n2 in range(0,A.shape[0]):
        #         if A[n1,n2]>0:
        #             if self.G_d.get_edge_data(node_list[n2],node_list[n1])['channel']=='down':
        #                 C[self.code_length*n1:self.code_length*(n1+1),self.code_length*n2:self.code_length*(n2+1)] = (1.0/in_degree)*down_regulation
        #             else:
        #                 C[self.code_length*n1:self.code_length*(n1+1),self.code_length*n2:self.code_length*(n2+1)] = (1.0/in_degree)*up_regulation

        for e in self.G_d.edges():
            i, j = node_idx[e[0]], node_idx[e[1]]

            if 'PPI' in self.G_d.get_edge_data(e[0],e[1])['category']:
                w = 1.0
            else:
                w = 1.0

            if self.G_d.get_edge_data(e[0],e[1])['channel']=='down':
                C[self.code_length*j:self.code_length*(j+1),self.code_length*i:self.code_length*(i+1)] = w*down_regulation
            else:
                C[self.code_length*j:self.code_length*(j+1),self.code_length*i:self.code_length*(i+1)] = w*up_regulation

            if self.G_d.get_edge_data(e[0],e[1])['module'] in neglect_modules:
                C[self.code_length*j:self.code_length*(j+1),self.code_length*i:self.code_length*(i+1)] = w*null_communicator

        for i in range(0,C.shape[0]):
            total_row = np.sum(C[i,:])

            if total_row>0.0:
                C[i,:] *= 1.0/total_row

        # Apply bounday condition
        for n in range(0,len(node_list)):
            if self.G_d.in_degree(node_list[n])==0:
                C[self.code_length*n:self.code_length*(n+1),self.code_length*n:self.code_length*(n+1)] = np.identity(self.code_length)

        self.C_sparse = sparse.csr_matrix(C)
        self.up_regulation = up_regulation
        self.down_regulation = down_regulation

    def disconnect_nodes(self,classname,sources):
        """
        Removes nodes of the specified class from the network communication (or transition) matrix.

        Parameters
        ----------
        classname: string
            Class of the nodes to be disconnected from the network communication matrix.

        sources: list
            List of source nodes which remains to be connected to the communication matrix.
        """

        nodes = list(self.G_d.nodes)

        C_dense = self.C_sparse.todense()

        for n in self.G_d.nodes(data=True):
            if classname==n[1]['class'] and n[0] not in sources:
                i = nodes.index(n[0])

                C_dense[:,self.code_length*i:self.code_length*(i+1)] = 0
                C_dense[self.code_length*i:self.code_length*(i+1),self.code_length*i:self.code_length*(i+1)] = np.identity(self.code_length)

        for i in range(0,C_dense.shape[0]):
            s = np.sum(C_dense[i,:])
            if s>0.0:
                C_dense[i,:] *= 1.0/np.sum(C_dense[i,:])
            else:
                C_dense[i,i] = 1.0

        self.C_sparse = sparse.csr_matrix(C_dense)

    def disconnect_drug_nodes(self,skip='None'):
        """
        Removes nodes of the specified class from the network communication (or transition) matrix.

        Parameters
        ----------
        classname: string
            Class of the nodes to be disconnected from the network communication matrix.

        sources: list
            List of source nodes which remains to be connected to the communication matrix.
        """

        drugs = []
        drug_ids = []

        nodes = list(self.G_d.nodes)

        all_nodes = self.G_d.nodes(data=True)

        for n in self.G_d.nodes(data=True):
            if 'Drug' in n[1]['class'] and skip not in n[1]['name']:
                drugs.append(n[1]['name'].split(' [')[0])
                drug_ids.append(n[0])

        total_edges_list = []
        total_edges_to_delete = 0
        total_deletions = []

        for e in self.G_d.edges:
            if e[0] in drug_ids:
                total_edges_to_delete += 1
                total_edges_list.append(self.G_d[e[0]][e[1]]['reaction'])
                total_deletions.append((nodes.index(e[0]),nodes.index(e[1])))

        C_dense = self.C_sparse.todense()

        for deletion in total_deletions:
            i,j = deletion[0], deletion[1]

            C_dense[self.code_length*j:self.code_length*(j+1),self.code_length*i:self.code_length*(i+1)] = np.zeros(shape=(self.code_length,self.code_length,))

        for i in range(0,C_dense.shape[0]):
            s = np.sum(C_dense[i,:])
            if s>0.0:
                C_dense[i,:] *= 1.0/np.sum(C_dense[i,:])
            else:
                C_dense[i,i] = 1.0

        self.C_sparse = sparse.csr_matrix(C_dense)

    def save_edges_and_nodes(self,pathway_nametag,network_type):
        """
        Saves the network nodes and edges as csv and pickled files.

        Parameters
        ----------
        pathway_name: string
            Name of the Reactome top level pathway.

        network_type: string
            Additional description about the network.
        """
        df_e = pd.DataFrame(self.edge_data)

        df_e.to_pickle(pathway_nametag+'_'+network_type+'-edges.pkl')
        df_e.to_csv(pathway_nametag+'_'+network_type+'-edges.csv',index=None)

        df_n = pd.DataFrame(self.node_data)
        df_n.to_pickle(pathway_nametag+'_'+network_type+'-nodes.pkl')
        df_n.to_csv(pathway_nametag+'_'+network_type+'-nodes.csv',index=None)

    def save_network(self,pathway_nametag,network_type):
        """
        Saves the network as a gml file.

        Parameters
        ----------
        pathway_name: string
            Name of the Reactome top level pathway.

        network_type: string
            Additional description about the network.
        """

        print("Communication network "+pathway_nametag+"-"+network_type+".gml created.")

        nx.write_gml(self.G_d,pathway_nametag+"-"+network_type+".gml")

    def initialize_network_state(self,additional_sources,source_state):
        """
        Saves the network as a gml file.

        Parameters
        ----------
        pathway_name: string
            Name of the Reactome top level pathway.

        network_type: string
            Additional description about the network.
        """

        network_state = (1.0/float(self.code_length))*np.ones(shape=(C.shape[0],1))
        network_sources = []

        node_list = self.G_d.nodes

        for n in additional_sources:
            try:
                i = node_list.index(n)

                network_state[self.code_length*i:self.code_length*(i+1),0] = source_state

                network_sources.append(i)

            except ValueError:
                pass

    def write_nodes_and_edges(self,filetag=''):
        ef = open(filetag+'-edges.tsv','w')
        print('Input\tOutput\tReactionName\tChannel\tReactionId',file=ef)

        for e in self.G_d.edges(data=True):
            outstring = e[0]+'\t'+e[1]+'\t'+e[2]['name']+'\t'+e[2]['channel']+'\t'+e[2]['reaction']
            print(outstring,file=ef)

        ef.close()

        nf = open(filetag+'-nodes.tsv','w')
        print('Node\tNodeName\tClass\tSequenced',file=nf)

        for n in self.G_d.nodes(data=True):
            if n[1]['sequenced']==0:
                outstring = n[0]+'\t'+n[1]['name']+'\t'+n[1]['class']+'\t'+str(n[1]['sequenced'])
            else:
                outstring = n[0]+'\t'+n[0]+'\t'+n[1]['class']+'\t'+str(n[1]['sequenced'])

            print(outstring,file=nf)

        nf.close()
