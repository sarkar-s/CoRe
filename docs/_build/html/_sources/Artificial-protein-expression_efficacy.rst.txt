Artificial protein expression screening
=======================================

**Computes the change in network relative entropy from viral PPIs after
setting immune system proteins, one at a time, to high abundance**

.. code:: ipython3

    import os, sys
    import numpy as np
    import scipy as sp
    import pandas as pd
    import copy as copy
    from tqdm.notebook import tqdm
    import math
    import scipy.stats as st
    
    from CoRe import reader
    from CoRe.ncip import ncip
    from CoRe.BA_C import BA
    
    import importlib
    
    import networkx as nx
    import matplotlib.pyplot as plt
    import json

.. code:: ipython3

    data_directory = "./Examples/Immune_System"
    os.chdir(data_directory)
    
    edge_data = pd.read_pickle('Immune_System_medium-PPI-edges.pkl')
    node_data = pd.read_pickle('Immune_System_medium-PPI-nodes.pkl')

.. code:: ipython3

    remake_graph = False
    
    if remake_graph==False:
        netObj = ncip()
        netObj.load_graph('Immune_System-medium-PPI.gml')
    else:
        netObj = ncip()
        netObj.load_data(edge_data,node_data)
        netObj.make_graph()
        netObj.save_network(pathway_nametag,network_type)

**All immune system communication network proteins that have PPI with
SARS-CoV-2 proteins**

.. code:: ipython3

    f = open('SARS_CoV2-Immune_System_interactions.json')
    SARS_nodes = json.load(f)
    f.close()
    
    all_sars_nodes = []
    
    for s in SARS_nodes.keys():
        all_sars_nodes += SARS_nodes[s]
        
    all_sars_nodes = list(set(all_sars_nodes))
    
    print(all_sars_nodes)


.. parsed-literal::

    ['GOLGA7', 'RIPK1', 'NEU1', 'NLRX1', 'PVR', 'RAB18', 'HECTD1', 'RAB7A', 'ITGB1', 'IMPDH2', 'ERP44', 'RAB14', 'ELOC', 'CSNK2B', 'RALA', 'AP2A2', 'RNF41', 'CYB5R3', 'ELOB', 'ECSIT', 'ANO6', 'PTGES2', 'STOM', 'HMOX1', 'GGH', 'RAB5C', 'NPC2', 'GLA', 'RHOA', 'TBK1', 'SLC27A2', 'IL17RA', 'SLC44A2', 'TOMM70', 'RAB10', 'EIF4E2']


**Specifying the reference state and construction of the global
transition matrix**

.. code:: ipython3

    initial_state_type = 'maxEnt'
    
    errorname = '0.0'
    rho = float(errorname)
    
    input_bits = 1
    code_length = int(2**input_bits)
    
    max_entropy_state = (1.0/float(code_length))*np.ones(shape=(code_length,))
    
    low_state = np.zeros(shape=(code_length,))
    low_state[-1] = 1.0
    
    high_state = np.zeros(shape=(code_length,))
    high_state[0] = 1.0
    
    if initial_state_type=='high':
        initial_state = high_state
    elif initial_state_type=='low':
        initial_state = low_state
    else:
        initial_state = max_entropy_state
    
    print(high_state,low_state)
    
    netObj.construct_C(rho,h=input_bits,neglect_modules=[])
    node_list = list(netObj.G_d.nodes)


.. parsed-literal::

    [1. 0.] [0. 1.]


**Disconnect all drugs from the network initially**

.. code:: ipython3

    netObj.disconnect_drug_nodes()

**Compute the reference stationary state of the network**

.. code:: ipython3

    initial_network_state = np.zeros(shape=(netObj.C_sparse.shape[0],1))
    network_sources = {}
    
    for n in range(0,len(node_list)):
        initial_network_state[code_length*n:code_length*(n+1),0] = initial_state
        
    network_sources = []
    
    reference_final_state, steps = netObj.get_final_state(initial_network_state,[])
    reference_final_entropy = netObj.state_entropy(reference_final_state,[])
    print('Reference state relative entropy: ',reference_final_entropy)


.. parsed-literal::

    Reference state relative entropy:  0.0


**Set the SARS-CoV-2 nodes in the network to low abundance**

.. code:: ipython3

    network_state = np.zeros(shape=(netObj.C_sparse.shape[0],1))
    network_sources = []
    
    for n in range(0,len(node_list)):
        network_state[code_length*n:code_length*(n+1),0] = initial_state
    
    for k in tqdm(SARS_nodes.keys()):
        for n in SARS_nodes[k]:
            try:
                i = node_list.index(n)
    
                network_state[netObj.code_length*i:netObj.code_length*(i+1),0] = low_state
    
                if i not in network_sources:
                    network_sources.append(i)
            except ValueError:
                pass



.. parsed-literal::

      0%|          | 0/17 [00:00<?, ?it/s]


**Relative entropy of the total network and number of steps to
stationary state.**

.. code:: ipython3

    final_state, steps = netObj.get_final_state(network_state,network_sources)
    final_entropy = netObj.state_entropy(final_state,network_sources)
    print(final_entropy)


.. parsed-literal::

    65.4157289254401


**Compute stationary state of the network due to SARS-CoV-2 PPIs and
proteins**

The proteins in the Reactome database were set to the state {1,0} to
compute the stationary state, and the subsequent change in the network
relative entropy.

.. code:: ipython3

    node_class = nx.get_node_attributes(netObj.G_d,"class")
    node_n = list(netObj.G_d.nodes())
    
    c = 0
    
    for i in range(0,len(node_n)):
        nn = node_n[i]
        if node_class[nn]=='EntityWithAccessionedSequence':
            relH = st.entropy(final_state[netObj.code_length*i:netObj.code_length*(i+1),0],max_entropy_state,base=2)
            
            if relH>0.01:
                c += 1

.. code:: ipython3

    all_sources = []
    
    for n in netObj.G_d.nodes(data=True):
        if n[1]['class']=='EntityWithAccessionedSequence' and n[0] not in all_sars_nodes:
            all_sources.append((n[0],netObj.G_d.in_degree(n[0])))

.. code:: ipython3

    entropy_shifts = {}
    H_drops = {}
    H_gains = {}
    
    for s_pair in tqdm(all_sources):
        s = s_pair[0]
        additional_source_nodes = [s]
        
        netObj.construct_C(rho,h=input_bits)
        netObj.disconnect_nodes('ChemicalDrug',additional_source_nodes)
        netObj.disconnect_nodes('ProteinDrug',additional_source_nodes)
        
        network_state = np.zeros(shape=(netObj.C_sparse.shape[0],1))
        network_sources = []
        
        for n in range(0,len(node_list)):
            network_state[code_length*n:code_length*(n+1),0] = initial_state
    
        for k in SARS_nodes.keys():
            for n in SARS_nodes[k]:
                try:
                    i = node_list.index(n)
    
                    network_state[netObj.code_length*i:netObj.code_length*(i+1),0] = low_state
    
                    network_sources.append(i)
                except ValueError:
                    pass
    
            for n in additional_source_nodes:
                try:
                    i = node_list.index(n)
    
                    network_state[netObj.code_length*i:netObj.code_length*(i+1),0] = high_state
    
                    network_sources.append(i)
                except ValueError:
                    pass
            
        entropy_shifts[s] = 0.0
    
        this_state, steps = netObj.get_final_state(network_state,network_sources)
        this_entropy = netObj.state_entropy(this_state,network_sources)
        H_drop, H_gain = netObj.entropy_drop_and_rise(this_state,final_state,reference_final_state,network_sources)
                    
        entropy_shifts[s] = this_entropy
            
        H_drops[s] = H_drop
        H_gains[s] = H_gain



.. parsed-literal::

      0%|          | 0/1122 [00:00<?, ?it/s]


.. code:: ipython3

    try:
        os.chdir('./counter_entropic_shift')
    except OSError:
        os.mkdir('./counter_entropic_shift')
        os.chdir('./counter_entropic_shift')

.. code:: ipython3

    node_data = nx.get_node_attributes(netObj.G_d,"name")
    node_class = nx.get_node_attributes(netObj.G_d,"class")
    
    of = open('high_all_protein_shifts.csv','w')
    
    print('Gene,Relative Entropy,In Degree',file=of)
    
    print('Ref,'+str(final_entropy)+',0',file=of)
    
    for s in all_sources:
        if node_class[s[0]]=="Complex":
            this_name = node_data[s[0]]
            this_name = this_name.replace(',',';')
        else:
            this_name = s[0]
            
        print(this_name+','+str(entropy_shifts[s[0]])+','+str(s[1]),file=of)
        
    of.close()
    
    of = open('split_all_high_protein_shifts-'+initial_state_type+'.csv','w')
    
    print('Protein,Drop,Gain',file=of)
    
    for s in all_sources:
        if node_class[s[0]]=="Complex":
            this_name = node_data[s[0]]
            this_name = this_name.replace(',',';')
        else:
            this_name = s[0]
            
        print(this_name+','+str(H_drops[s[0]])+','+str(H_gains[s[0]]),file=of)
        
    of.close()

