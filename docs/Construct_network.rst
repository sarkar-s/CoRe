Construct network
=================

**We use the data in the Reactome database for a Top Level Pathway to
construct a network where the nodes represent physical entities
connected by events, which include reactions, regulation,
(de)polymerisation, synthesis, and degradations.**

.. code:: ipython3

    import os, sys
    from json import dumps
    import logging
    import pandas as pd
    import numpy as np
    import scipy as sp
    import math
    from tqdm.notebook import tqdm
    
    from CoRe.ncip import ncip
    from CoRe.cypher_commands import command_set
    
    from neo4j import GraphDatabase, basic_auth
    
    import networkx as nx
    import matplotlib.pyplot as plt

**Create reactome graph database object using neo4j**

.. code:: ipython3

    reactome = GraphDatabase.driver("bolt://127.0.0.1:7687", auth=("neo4j", "test"))
    
    session = reactome.session()

**All top level pathways in Reactome**

.. code:: ipython3

    result = session.run("MATCH(tp:TopLevelPathway {speciesName:'Homo sapiens'}) RETURN DISTINCT tp.displayName")
    
    all_TopLevelPathways = [res[0] for res in result]
        
    all_TopLevelPathways




.. parsed-literal::

    ['Autophagy',
     'Cell Cycle',
     'Cell-Cell communication',
     'Cellular responses to stimuli',
     'Chromatin organization',
     'Circadian Clock',
     'Developmental Biology',
     'Digestion and absorption',
     'Disease',
     'DNA Repair',
     'DNA Replication',
     'Drug ADME',
     'Extracellular matrix organization',
     'Gene expression (Transcription)',
     'Hemostasis',
     'Immune System',
     'Metabolism',
     'Metabolism of proteins',
     'Metabolism of RNA',
     'Muscle contraction',
     'Neuronal System',
     'Organelle biogenesis and maintenance',
     'Programmed Cell Death',
     'Protein localization',
     'Reproduction',
     'Sensory Perception',
     'Signal Transduction',
     'Transport of small molecules',
     'Vesicle-mediated transport']



**Select a Top Level Pathway**

.. code:: ipython3

    selected_TopLevelPathway = 'Immune System'
    pathway_nametag = selected_TopLevelPathway.replace(' ','_')
    
    output_directory = '/Users/swarnavo/Research/Reactome-Graph-Database/HumanData/'+ pathway_nametag
    
    try:
        os.chdir(output_directory)
    except FileNotFoundError:
        os.mkdir(output_directory)
        os.chdir(output_directory)

**Find the subpathways of the top level pathway**

.. code:: ipython3

    command = command_set['subpathways'].replace('#',selected_TopLevelPathway)
    result = session.run(command)
    
    subpathways = [res[0] for res in result]
        
    subpathways




.. parsed-literal::

    ['Innate Immune System',
     'Cytokine Signaling in Immune system',
     'Adaptive Immune System']



**Collect all ReactionLikeEvents in the selected TopLevelPathway, these
constitute the edges, or information channels, of the network.**

.. code:: ipython3

    command = command_set['events'].replace('#',selected_TopLevelPathway)
    
    result = session.run(command)
        
    all_results = [res[0] for res in result]
    
    data = {}
    data['reaction'] = [res['stId'] for res in all_results if res['isInDisease']==False]
    #data['category'] = [res['category'] for res in all_results if res['isInDisease']==False]
    data['name'] = [res['displayName'] for res in all_results if res['isInDisease']==False]
    data['schemaClass'] = [res['schemaClass'] for res in all_results if res['isInDisease']==False]
    data['module'] = [selected_TopLevelPathway for res in all_results if res['isInDisease']==False]
    
    df = pd.DataFrame(data)
    df.to_pickle(pathway_nametag+'-ReactionLikeEvents.pkl')
    df.to_csv(pathway_nametag+'-ReactionLikeEvents.csv',index=None)
    
    print('Total ReactionLikeEvents in',selected_TopLevelPathway,': ',len(data['reaction']))


.. parsed-literal::

    Total ReactionLikeEvents in Immune System :  1622


.. code:: ipython3

    c = 0
    
    for ii in range(0,len(df['reaction'])):
        command_reg = command_set['regulation'].replace('#',df['reaction'][ii])
        result_reg = session.run(command_reg)
    
        for res in result_reg:
            schemaClass = res[0]['schemaClass']
            stId = res[0]['stId']
            n = res[0]['displayName']
    
            if stId!=None:
                c += 1
                
    print(c)


.. parsed-literal::

    176


**Identify inputs and outputs to each ReactionLikeEvent in the selected
TopLevelPathway**

The inputs and the outputs are physical entities that form the nodes of
the network. There are 3 options for querying the inputs and ouputs
(*network_type*): **coarse** - does not break down physical entities
complexes, defined set, and candidate set. **medium** - break downs
defined and candidate sets into individual compoments. **fine** - break
downs complexes, defined and candidate sets into individual compoments.

**neglect_class** - list of Reactome schemaClass of physial entities to
be neglected from the graph. If network_type=medium, then we DefinedSet
and CandidateSet has to be neglected because their components are being
included individually. Additionally, we neglect the SimpleEntities,
*e.g.* ATP, ADP, H2O, from the network.

.. code:: ipython3

    network_type = 'medium'
    neglect_class = ['DefinedSet','CandidateSet']
    field_value = 'stId' # Reactome database object attribute to use as a nodename
    
    # Create network object
    netmaker = ncip(field_value)
    
    for ii in tqdm(range(0,len(df['reaction']))):
        command_in = command_set['input'][network_type].replace('#',df['reaction'][ii])
        result_in = session.run(command_in)
        input_results = [res[0] for res in result_in]
    
        command_out = command_set['output'][network_type].replace('#',df['reaction'][ii])
        result_out = session.run(command_out)
        output_results = [res[0] for res in result_out]
    
        command_sp = command_set['superpathways'].replace('#',df['reaction'][ii])
        result_sp = session.run(command_sp)
        
        if 'Exocytosis' in df['name'][ii]:
            netmaker.construct_exo_edges(input_results,neglect_class=neglect_class)
        else:
            netmaker.construct_edges(input_results,output_results,neglect_class=neglect_class)
    
        netmaker.add_edges(session,df.loc[ii])
    
        command_reg = command_set['regulation'].replace('#',df['reaction'][ii])
        result_reg = session.run(command_reg)
    
        for res in result_reg:
            schemaClass = res[0]['schemaClass']
            stId = res[0]['stId']
            n = res[0]['displayName']
    
            if stId!=None:
                command_1 = command_set['regulator'].replace('#',schemaClass)
                command = command_1.replace('%',stId)
    
                result_reg = session.run(command)
                reg_results = [result[0] for result in result_reg]
    
                reg_edge_nodes = netmaker.construct_edges(reg_results,output_results,neglect_class=neglect_class)
    
                edge_info = [stId,schemaClass,'Regulation',n]
                edge_info = {'reaction':stId, 'schemaClass':schemaClass, 'name':n,'module':df.loc[ii]['module']}
    
                netmaker.add_edges(session,edge_info)
    
    print('Total nodes in',selected_TopLevelPathway,': ',len(netmaker.node_data['node']))
    print('Total edges in',selected_TopLevelPathway,': ',len(netmaker.edge_data['input']))



.. parsed-literal::

      0%|          | 0/1622 [00:00<?, ?it/s]


.. parsed-literal::

    Total nodes in Immune System :  3367
    Total edges in Immune System :  16557


**Create and save networkx graph of the top level pathway**

.. code:: ipython3

    netmaker.save_edges_and_nodes(pathway_nametag,network_type)

