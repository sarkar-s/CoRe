SARS-CoV2 interactions
======================

.. code:: ipython3

    import os, sys
    import numpy as np
    import scipy as sp
    import pandas as pd

    import math

    import networkx as nx
    import matplotlib.pyplot as plt

    from neo4j import GraphDatabase, basic_auth
    from CoRe.cypher_commands import command_set
    from CoRe.ncip import ncip

    import json

**Read SARS-CoV2 proteins and human proteins that have protein-protein
interaction**

.. code:: ipython3

    data_directory = "./Examples"
    os.chdir(data_directory)

    data = pd.read_csv('SARS-Cov2_proteins.csv')

    Cov2_proteins = data['Bait'].to_list()

    interaction_proteins = data['PreyGene'].to_list()

    SARS_CoV2_pnames = pd.read_csv('ordered_SARS-CoV-2_proteins.csv')['Bait'].to_list()
    print(SARS_CoV2_pnames)

    CoV2_interactions = {}

    for CoV2_protein in SARS_CoV2_pnames:
        CoV2_interactions[CoV2_protein] = []

    for a,b in zip(Cov2_proteins,interaction_proteins):
        CoV2_interactions[a].append(b)

    # create json object from dictionary
    json_obj = json.dumps(CoV2_interactions)

    # open file for writing, "w"
    f = open('SARS_CoV2-interactions.json','w')

    # write json object to file
    f.write(json_obj)

    # close file
    f.close()


.. parsed-literal::

    ['SARS-CoV2 Nsp1', 'SARS-CoV2 Nsp2', 'SARS-CoV2 Nsp3', 'SARS-CoV2 Nsp4', 'SARS-CoV2 Nsp5', 'SARS-CoV2 Nsp5_C145A', 'SARS-CoV2 Nsp6', 'SARS-CoV2 Nsp7', 'SARS-CoV2 Nsp8', 'SARS-CoV2 Nsp9', 'SARS-CoV2 Nsp10', 'SARS-CoV2 Nsp11', 'SARS-CoV2 Nsp12', 'SARS-CoV2 Nsp13', 'SARS-CoV2 Nsp14', 'SARS-CoV2 Nsp15', 'SARS-CoV2 Spike', 'SARS-CoV2 ORF3a', 'SARS-CoV2 ORF3b', 'SARS-CoV2 E', 'SARS-CoV2 M', 'SARS-CoV2 ORF6', 'SARS-CoV2 ORF7a', 'SARS-CoV2 ORF8', 'SARS-CoV2 ORF9b', 'SARS-CoV2 ORF9c', 'SARS-CoV2 N', 'SARS-CoV2 ORF10']


Read the NetworkX graph of the chosen top level pathway.

.. code:: ipython3

    selected_pathway = 'Immune System'
    pathway_nametag = selected_pathway.replace(' ','_')

    network_type = 'medium'

    os.chdir(pathway_nametag)

    netObj = ncip()
    netObj.load_graph(pathway_nametag+"-"+network_type+".gml")

.. code:: ipython3

    node_list = netObj.G_d.nodes.data()

    ref_genes = [x[0] for x in node_list if x[1]['sequenced']!=0]

    print(len(ref_genes))


.. parsed-literal::

    1468


**Identify graph nodes that match with the protein-protein interaction
data.**

.. code:: ipython3

    CoV2_immune_interactions_node_id = {}

    print('Directly interacting reference gene products:\n')

    for CoV2_protein in SARS_CoV2_pnames:
        CoV2_immune_interactions_node_id[CoV2_protein] = []

        #for hp in CoV2_interactions[CoV2_protein]:
        #    if hp in ref_genes and hp not in CoV2_immune_interactions_node_id[CoV2_protein]:
        #        CoV2_immune_interactions_node_id[CoV2_protein].append(hp)

        CoV2_immune_interactions_node_id[CoV2_protein] = [hp for hp in CoV2_interactions[CoV2_protein] if hp in ref_genes and hp not in CoV2_immune_interactions_node_id[CoV2_protein]]

        if len(CoV2_immune_interactions_node_id[CoV2_protein])==0:
            del CoV2_immune_interactions_node_id[CoV2_protein]
        else:
            print(CoV2_protein,CoV2_immune_interactions_node_id[CoV2_protein])

    all_interacting_nodes = []
    for v in CoV2_immune_interactions_node_id.values():
        all_interacting_nodes += v


.. parsed-literal::

    Directly interacting reference gene products:

    SARS-CoV2 Nsp2 ['SLC27A2', 'EIF4E2']
    SARS-CoV2 Nsp7 ['CYB5R3', 'RALA', 'RAB5C', 'RAB7A', 'RAB10', 'RAB14', 'RHOA', 'PTGES2', 'RAB18']
    SARS-CoV2 Nsp8 ['HECTD1']
    SARS-CoV2 Nsp10 ['AP2A2']
    SARS-CoV2 Nsp12 ['RIPK1']
    SARS-CoV2 Nsp13 ['TBK1']
    SARS-CoV2 Nsp14 ['GLA', 'IMPDH2']
    SARS-CoV2 Nsp15 ['RNF41']
    SARS-CoV2 Spike ['GOLGA7']
    SARS-CoV2 ORF3a ['HMOX1']
    SARS-CoV2 E ['SLC44A2']
    SARS-CoV2 M ['STOM', 'ANO6']
    SARS-CoV2 ORF8 ['ITGB1', 'PVR', 'NPC2', 'GGH', 'IL17RA', 'NEU1', 'ERP44']
    SARS-CoV2 ORF9b ['TOMM70']
    SARS-CoV2 ORF9c ['NLRX1', 'ECSIT']
    SARS-CoV2 N ['CSNK2B']
    SARS-CoV2 ORF10 ['ELOC', 'ELOB']


.. code:: ipython3

    # create json object from dictionary
    json = json.dumps(CoV2_immune_interactions_node_id)

    f = open('SARS_CoV2-'+pathway_nametag+'_interactions.json','w')

    f.write(json)

    f.close()
