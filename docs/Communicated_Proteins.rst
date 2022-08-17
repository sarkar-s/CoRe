Communicated proteins
=====================

.. code:: ipython3

    import os
    from json import dumps
    import logging
    import pandas as pd
    import numpy as np
    import copy

    import matplotlib.pyplot as plt
    from matplotlib import cm

    from CoRe.cypher_commands import command_set

    import networkx as nx
    import json

.. code:: ipython3

    current_directory = os.getcwd()

    selected_TopLevelPathway = 'Immune System'
    pathway_nametag = selected_TopLevelPathway.replace(' ','_')

    network_type = 'medium-PPI'

    state_type = 'maxEnt'

    data_directory = current_directory + "/Examples/"+pathway_nametag
    os.chdir(data_directory)

    f = open('SARS_CoV2-'+pathway_nametag+'_interactions.json')
    SARS_nodes = json.load(f)
    f.close()

    communicated_proteins = pd.read_csv(state_type+'-SARS_CoV2_'+pathway_nametag+'_'+network_type+'_affected_genes0.0.csv')
    print(list(communicated_proteins))
    print(communicated_proteins['node_ids'].count())


.. parsed-literal::

    ['node_ids', 'node_index', 'SARS-CoV2 Nsp2', 'SARS-CoV2 Nsp7', 'SARS-CoV2 Nsp12', 'SARS-CoV2 Nsp13', 'SARS-CoV2 Nsp14', 'SARS-CoV2 ORF3a', 'SARS-CoV2 M', 'SARS-CoV2 ORF8', 'SARS-CoV2 ORF9b', 'SARS-CoV2 ORF9c', 'SARS-CoV2 ORF10']
    91


.. code:: ipython3

    all_ref_gene_names = {}

    i = 0
    for gen_name in communicated_proteins['node_ids']:
        try:
            all_ref_gene_names[gen_name].append(i)
        except KeyError:
            all_ref_gene_names[gen_name] = [i]

        i += 1

.. code:: ipython3

    # Condense selected protein list
    first_indices = []
    indices_to_drop = []

    for k in all_ref_gene_names.keys():
        first_indices.append(all_ref_gene_names[k][0])

        for d in all_ref_gene_names[k][1:]:
            indices_to_drop.append(d)

.. code:: ipython3

    communicated_proteins = communicated_proteins.drop(indices_to_drop)
    all_ref_gene_names = communicated_proteins['node_ids'].to_list()

.. code:: ipython3

    sars_proteins = list(communicated_proteins)[2:]
    print(sars_proteins)


.. parsed-literal::

    ['SARS-CoV2 Nsp2', 'SARS-CoV2 Nsp7', 'SARS-CoV2 Nsp12', 'SARS-CoV2 Nsp13', 'SARS-CoV2 Nsp14', 'SARS-CoV2 ORF3a', 'SARS-CoV2 M', 'SARS-CoV2 ORF8', 'SARS-CoV2 ORF9b', 'SARS-CoV2 ORF9c', 'SARS-CoV2 ORF10']


.. code:: ipython3

    SARS_indirect_nodes = {}
    SARS_indirect_nodes_wts = {}

    for k in SARS_nodes.keys():
        SARS_indirect_nodes[k] = []
        SARS_indirect_nodes_wts[k] = []

    print(SARS_nodes.keys())

    for s in sars_proteins:
        d = communicated_proteins[s].to_numpy()

        for i in range(0,d.shape[0]):
            if d[i]>0.0 and all_ref_gene_names[i] not in SARS_indirect_nodes[s]:
                SARS_indirect_nodes[s].append(all_ref_gene_names[i])
                SARS_indirect_nodes_wts[s].append(d[i])


.. parsed-literal::

    dict_keys(['SARS-CoV2 Nsp2', 'SARS-CoV2 Nsp7', 'SARS-CoV2 Nsp8', 'SARS-CoV2 Nsp10', 'SARS-CoV2 Nsp12', 'SARS-CoV2 Nsp13', 'SARS-CoV2 Nsp14', 'SARS-CoV2 Nsp15', 'SARS-CoV2 Spike', 'SARS-CoV2 ORF3a', 'SARS-CoV2 E', 'SARS-CoV2 M', 'SARS-CoV2 ORF8', 'SARS-CoV2 ORF9b', 'SARS-CoV2 ORF9c', 'SARS-CoV2 N', 'SARS-CoV2 ORF10'])


.. code:: ipython3

    SARS_affected_refgenes = copy.deepcopy(SARS_nodes)

    for s in sars_proteins:
        d = communicated_proteins[s].to_numpy()

        for i in range(0,d.shape[0]):
            if d[i]>0.0 and all_ref_gene_names[i] not in SARS_affected_refgenes[s]:
                SARS_affected_refgenes[s].append(all_ref_gene_names[i])

        print(s,'\t',len(SARS_affected_refgenes[s]))

    json_obj = json.dumps(SARS_affected_refgenes)

    f = open(state_type+'-SARS_CoV2_total_'+pathway_nametag+'_'+network_type+'_interactions.json','w')
    f.write(json_obj)
    f.close()


.. parsed-literal::

    SARS-CoV2 Nsp2 	 3
    SARS-CoV2 Nsp7 	 28
    SARS-CoV2 Nsp12 	 7
    SARS-CoV2 Nsp13 	 3
    SARS-CoV2 Nsp14 	 7
    SARS-CoV2 ORF3a 	 2
    SARS-CoV2 M 	 33
    SARS-CoV2 ORF8 	 46
    SARS-CoV2 ORF9b 	 3
    SARS-CoV2 ORF9c 	 5
    SARS-CoV2 ORF10 	 4
