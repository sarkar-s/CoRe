Scanning drugs to counter signaling from viral PPIs
===================================================

**Computes the change in network relative entropy from viral PPIs after
setting drug state to high abundance**

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

**Drugs that are part of the immune system in the Reactome database**

.. code:: ipython3

    drugs = []
    drug_names = []
    
    all_nodes = netObj.G_d.nodes(data=True)
    
    for n in netObj.G_d.nodes(data=True):
        if 'Drug' in n[1]['class'] and 'lat' not in n[1]['name'] and n[0] not in drug_names:
            print(n[0],n[1]['name'],n[1]['class'])
            drugs.append((n[0],n[1]['class'],n[1]['name']))
            drug_names.append(n[0])


.. parsed-literal::

    R-ALL-9679737 natalizumab [extracellular region] ProteinDrug
    R-ALL-9681784 anakinra [extracellular region] ProteinDrug
    R-ALL-9681764 isunakinra [extracellular region] ProteinDrug
    R-ALL-9678829 baricitinib [cytosol] ChemicalDrug
    R-ALL-9678779 tofacitinib [cytosol] ChemicalDrug
    R-ALL-9678901 ibrutinib [cytosol] ChemicalDrug
    R-ALL-9678960 acalabrutinib [cytosol] ChemicalDrug
    R-ALL-9678786 ruxolitinib [cytosol] ChemicalDrug
    R-ALL-9715049 sarilumab [extracellular region] ProteinDrug
    R-ALL-9715070 satralizumab [extracellular region] ProteinDrug
    R-ALL-9681301 tocilizumab [extracellular region] ProteinDrug
    R-ALL-9679801 amlexanox [cytosol] ChemicalDrug
    R-ALL-9724689 omalizumab [extracellular region] ProteinDrug
    R-ALL-9678991 tacrolimus [cytosol] ChemicalDrug
    R-ALL-9679471 camrelizumab [extracellular region] ProteinDrug
    R-ALL-9679434 cemiplimab [extracellular region] ProteinDrug
    R-ALL-9679411 nivolumab [extracellular region] ProteinDrug
    R-ALL-9678628 HCQ [cytosol] ChemicalDrug
    R-ALL-9717004 delgocitinib [cytosol] ChemicalDrug


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

    ['RAB7A', 'ELOC', 'ERP44', 'ECSIT', 'ANO6', 'AP2A2', 'RHOA', 'PTGES2', 'IL17RA', 'SLC44A2', 'CYB5R3', 'TOMM70', 'RAB10', 'TBK1', 'HECTD1', 'RIPK1', 'SLC27A2', 'ELOB', 'GLA', 'ITGB1', 'GGH', 'NLRX1', 'EIF4E2', 'HMOX1', 'RAB14', 'IMPDH2', 'NEU1', 'RALA', 'CSNK2B', 'STOM', 'RNF41', 'PVR', 'NPC2', 'GOLGA7', 'RAB5C', 'RAB18']


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


**Compute stationary state of the network due to SARS-CoV-2 PPIs**

.. code:: ipython3

    final_state, steps = netObj.get_final_state(network_state,network_sources)
    final_entropy = netObj.state_entropy(final_state,network_sources,reference_final_state)
    print(final_entropy)


.. parsed-literal::

    65.4157289254401


**Compute stationary state of the network due to SARS-CoV-2 PPIs and
drugs**

The drugs in the Reactome database were set to the state {1,0} to
compute the stationary state, and the subsequent change in the network
relative entropy.

.. code:: ipython3

    entropy_shifts = {}
    H_drops = {}
    H_gains = {}
    
    for s_pair in tqdm(drugs):
        s = s_pair[0]
        additional_source_nodes = [s]
        
        netObj.load_graph('Immune_System-medium-PPI.gml')
        netObj.disconnect_drug_nodes(skip=s)
        
        netObj.construct_C(rho,h=input_bits)
        
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
        this_entropy = netObj.state_entropy(this_state,network_sources,reference_final_state)
        H_drop, H_gain = netObj.entropy_drop_and_rise(this_state,final_state,reference_final_state,network_sources)
                    
        entropy_shifts[s] = this_entropy
        
        H_drops[s] = H_drop
        H_gains[s] = H_gain



.. parsed-literal::

      0%|          | 0/19 [00:00<?, ?it/s]


.. code:: ipython3

    try:
        os.chdir('./counter_entropic_shift')
    except OSError:
        os.mkdir('./counter_entropic_shift')
        os.chdir('./counter_entropic_shift')

.. code:: ipython3

    node_data = nx.get_node_attributes(netObj.G_d,"name")
    node_class = nx.get_node_attributes(netObj.G_d,"class")
    
    of = open('high_all_drug_shifts-'+initial_state_type+'.csv','w')
    
    print('Drug,Relative Entropy,Drug Type',file=of)
    
    print('Ref,'+str(final_entropy)+',0',file=of)
    
    for s in drugs:
        if node_class[s[0]]=="Complex":
            this_name = node_data[s[0]]
            this_name = this_name.replace(',',';')
        else:
            this_name = s[2].split(' [')[0]
            
        print(this_name+','+str(entropy_shifts[s[0]])+','+str(s[1]),file=of)
        
    of.close()
    
    of = open('split_all_drug_shifts-'+initial_state_type+'.csv','w')
    
    print('Drug,Drop,Gain,Drug Type',file=of)
    
    for s in drugs:
        if node_class[s[0]]=="Complex":
            this_name = node_data[s[0]]
            this_name = this_name.replace(',',';')
        else:
            this_name = s[2].split(' [')[0]
            
        print(this_name+','+str(H_drops[s[0]])+','+str(H_gains[s[0]])+','+str(s[1]),file=of)
        
    of.close()

