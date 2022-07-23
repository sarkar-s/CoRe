.. CoReI documentation master file, created by
   sphinx-quickstart on Sat Jul  2 17:19:44 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to CoReI's documentation!
=================================

.. image:: ./ImmuneSystem-category.png
  :width: 140
  :align: right
  :alt: Immune System network

**Co**\ mmunication through the **Re**\ actome and Interactome (CoRe) is an open-source python-based software tools to model human biological pathways
as information transferring networks. CoReI uses the `Reactome <https://reactome.org>`_ and the `Interactome <http://interactome.dfci.harvard.edu/H_sapiens/>`_ databases
to construct communication network models of human biology. The current version computes the information transferred due protein-protein interactions (PPIs) with viral
proteins, demonstrated using the SARS-CoV-2 PPIs as a case study.

CoReI uses the gene sets from the c5 collection in the `MSigDB <https://www.gsea-msigdb.org/gsea/msigdb/>`_ database to identify the significant biological processes using
Gene Ontology over-representation analysis.

:Author:
  `Swarnavo Sarkar <https://swarnavosarkar.com>`_

:Download:
 `github.com/sarkar-s/CoReI <https://github.com/sarkar-s/CoRe.git>`_

.. toctree::
   :maxdepth: 3
   :caption: Contents:

   Construction_of_Top_Level_Pathway_Network.ipynb
   Add-human_PPI.ipynb
   SARS-CoV-2_human_PPI.ipynb
   Information_transfer_in_the_network.ipynb
   Communicated_Proteins.ipynb
   communicated_ORA.ipynb
   ncip
   enGO
   reader
   dependencies



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
