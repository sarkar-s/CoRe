Add human PPI to the network
============================

This notebook identfies the proteins in the communication network built
from Reactome and adds edges due protein-protein interactions (PPI). The
human PPI data is taken from `The Human Reference Protein Interactome
Mapping Project <http://www.interactome-atlas.org/download>`__.

.. code:: ipython3

    import os, sys
    import numpy as np
    import scipy as sp
    import pandas as pd
    import copy as copy
    from tqdm.notebook import tqdm
    import math
    import scipy.stats as st

    from CoRe.ncip import ncip
    from CoRe import reader

    import matplotlib.pyplot as plt

    import networkx as nx

**Read protein-protein interaction data.**

.. code:: ipython3

    data_directory = "./Examples"
    os.chdir(data_directory)

    huPPI = pd.read_csv('ppi-genes.csv',header=None)

**Read the nodes and the edges of the reactome information network.**

.. code:: ipython3

    selected_pathway = 'Immune System'
    pathway_nametag = selected_pathway.replace(' ','_')

    network_type = 'medium'

    data_directory = "./" + pathway_nametag
    os.chdir(data_directory)

    edge_data = pd.read_pickle(pathway_nametag+'_'+network_type+'-edges.pkl')
    node_data = pd.read_pickle(pathway_nametag+'_'+network_type+'-nodes.pkl')
    all_nodes = node_data['node']

**Identify the nodes in the network that have PPI, and add them as
additional edges to the network.**

.. code:: ipython3

    c = 0

    for i,row in tqdm(huPPI.iterrows()):
        hu1, hu2 = row[0], row[1]

        x = all_nodes[all_nodes == hu1]
        y = all_nodes[all_nodes == hu2]

        if len(x)>0 and len(y)>0:
            df_temp = pd.DataFrame([[hu1,hu2,'PPI','PPI','PPI','PPI']],columns=['input','output','reaction','name','schemaClass','module'])
            edge_data = pd.concat([edge_data,df_temp],sort=False,ignore_index=True)

            c += 1

    print('Total protein-protein interactions to add:',c)



.. parsed-literal::

    0it [00:00, ?it/s]


.. parsed-literal::

    Total protein-protein interactions to add: 532


**Create the communication network with the combined reactome and
interactome edges.**

.. code:: ipython3

    netmaker = ncip()
    netmaker.load_data(edge_data,node_data)

    netmaker.save_edges_and_nodes(pathway_nametag,network_type)

    netmaker.make_graph()
    netmaker.save_network(pathway_nametag,network_type)


.. parsed-literal::

    100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 20813/20813 [00:01<00:00, 14931.90it/s]


.. parsed-literal::

    Communication network Immune_System-medium.gml created.


**Remove self-loops that are not regulatory reactions.**

.. code:: ipython3

    sl = list(nx.selfloop_edges(netmaker.G_d))

    for e in sl:
        if 'regulation' not in netmaker.G_d[e[0]][e[1]]['category']:
            netmaker.G_d.remove_edge(e[0],e[1])

**Remove isolates, or nodes that are note connected to any other node in
the network.**

.. code:: ipython3

    ni = list(nx.isolates(netmaker.G_d))

    for n in ni:
        netmaker.G_d.remove_node(n)

**Save the network with the PPI edges.**

.. code:: ipython3

    netmaker.save_network(pathway_nametag,network_type+'-PPI')
    netmaker.save_edges_and_nodes(pathway_nametag,network_type+'-PPI')


.. parsed-literal::

    Communication network Immune_System-medium-PPI.gml created.
