"""
Entropic enrichment of Gene Ontology(GO) gene sets
"""

import numpy as np
import math
import networkx as nx
import pandas as pd
from glob import glob
import os, sys
from scipy.stats import hypergeom as hg
#from multipy.fdr import qvalue
from scipy.interpolate import UnivariateSpline
from statsmodels.stats.multitest import fdrcorrection
from itertools import compress

def qvalue(pvals, threshold=0.05, verbose=True):
    """Function for estimating q-values from p-values using the Storey-
    Tibshirani q-value method (2003).
    Input arguments:
    ================
    pvals       - P-values corresponding to a family of hypotheses.
    threshold   - Threshold for deciding which q-values are significant.
    Output arguments:
    =================
    significant - An array of flags indicating which p-values are significant.
    qvals       - Q-values corresponding to the p-values.
    """

    """Count the p-values. Find indices for sorting the p-values into
    ascending order and for reversing the order back to original."""
    m, pvals = len(pvals), np.asarray(pvals)
    ind = np.argsort(pvals)
    rev_ind = np.argsort(ind)
    pvals = pvals[ind]

    # Estimate proportion of features that are truly null.
    kappa = np.arange(0, 0.96, 0.01)
    pik = [sum(pvals > k) / (m*(1-k)) for k in kappa]
    cs = UnivariateSpline(kappa, pik, k=3, s=None, ext=0)
    pi0 = float(cs(1.))
    if (verbose):
        print('The estimated proportion of truly null features is %.3f' % pi0)

    """The smoothing step can sometimes converge outside the interval [0, 1].
    This was noted in the published literature at least by Reiss and
    colleagues [4]. There are at least two approaches one could use to
    attempt to fix the issue:
    (1) Set the estimate to 1 if it is outside the interval, which is the
        assumption in the classic FDR method.
    (2) Assume that if pi0 > 1, it was overestimated, and if pi0 < 0, it
        was underestimated. Set to 0 or 1 depending on which case occurs.
    Here we have chosen the first option, since it is the more conservative
    one of the two.
    """
    if (pi0 < 0 or pi0 > 1):
        pi0 = 1
        print('Smoothing estimator did not converge in [0, 1]')

    # Compute the q-values.
    qvals = np.zeros(np.shape(pvals))
    qvals[-1] = pi0*pvals[-1]
    for i in np.arange(m-2, -1, -1):
        qvals[i] = min(pi0*m*pvals[i]/float(i+1), qvals[i+1])

    # Test which p-values are significant.
    significant = np.zeros(np.shape(pvals), dtype='bool')
    significant[ind] = qvals<threshold

    """Order the q-values according to the original order of the p-values."""
    qvals = qvals[rev_ind]
    return significant, qvals

def readGOsets(GO_file,GO_category):
    """
    Reads the gene ontology data set as a python dictionary.

    Parameters
    ----------
    GO_file: string
        Name of the file containing Gene Ontology gene sets.

    GO_category: string
        Name of the gene ontology category, one among the three, 'GOBP', 'GOCC', or 'GOMF'.

    Returns
    -------
    GO_BPs: dict
        Dict with gene ontology gene set as keys and the list of associated genes as dictionary entries.

    total_unique_genes: list
        Names of unique genes in the total gene ontology data set.
    """

    rf = open(GO_file,'r')
    all_lines = rf.readlines()
    rf.close()

    GO_BPs = {}

    total_unique_genes = []

    for l in all_lines:
        this_set = l.rstrip('\r\n').split('\t')

        if GO_category in this_set[0]:
            #GO_BPs[this_set[0]] = this_set[2:]
            GO_BPs[this_set[0]] = pd.DataFrame(this_set[2:])

            GO_BPs_dict[this_set[0]] = this_set[2:]

            for g in this_set[2:]:
                if g not in total_unique_genes:
                    total_unique_genes.append(g)

    total_gene_size = len(total_unique_genes)

    return GO_BPs, total_unique_genes

def findGOembedding(GO_sets,outputfile):
    """
    Identifies that gene sets that are embedded within other gene sets. The gene sets that are not embedded
    in any other gene set are returned as dictionary keys with an empty list as the entry.

    Parameters
    ----------
    GO_sets: dict
        Dictionary with gene set name as keys and the list of associated genes as entries.

    outputfile: string
        Name of the file to store the output.

    Returns
    -------
    GO_embedding: dict
        Dict with the index of gene ontology set name as keys and the list of indices of the gene
        sets that contains the key gene set.
    """

    GO_embedding = {}

    GO_BP_names = list(GO_sets.keys())

    lf = open('GOBP_list.csv','w')
    lf.close()

    for i in range(0,len(GO_BP_names)):
        bp = GO_BP_names[i]
        GO_embedding[i] = []

    for i_1 in range(0,len(GO_BP_names)):
        bp_1 = GO_BP_names[i_1]
        this_gene_set_size = len(GO_sets[bp_1])

        outline = str(i_1)

        for i_2 in range(0,len(GO_BP_names)):
            bp_2 = GO_BP_names[i_2]

            if len(GO_BPs[bp_2])>this_gene_set_size:
                check = 0

                for gene in GO_sets[bp_1]:
                    if gene in GO_sets[bp_2]:
                        check += 1

                if check==this_gene_set_size:
                    GO_embedding[i_1].append(i_2)

                    outline += ','+str(i_2)

        print(bp_1,file=open('GOBP_list.csv','a'))

    wf = open(outputfile,'w')

    for bp in GO_embedding.keys():
        outline = bp

        for bp_in in GO_embedding[bp]:
            outline += ','+bp_in

        print(outline,file=wf)

    wf.close()

    return GO_embedding

def findGOcontainer(GO_embedding,outputfile):
    """
    Identifies that gene sets that contains other gene sets. The gene sets that do not contain other
    gene sets are returned as dictionary keys with an empty list as the entry.

    Parameters
    ----------
    GO_sets: dict
        Dictionary with gene set index as keys and the list of gene set that contains it as entries.

    outputfile: string
        Name of the file to store the output.

    Returns
    -------
    GO_container: dict
        Dict with the index of gene ontology set name as keys and the list of indices of the gene sets
        that are contained in the key gene set.
    """

    GO_container = {}

    for bp in GO_embedding.keys():
        GO_container[bp] = []

    for bp_1 in GO_container.keys():

        for bp_2 in GO_embedding.keys():
            if bp_1 in GO_embedding[bp_2]:
                GO_container[bp_1].append(bp_2)

    wf = open(outputfile,'w')

    for bp in GO_container.keys():
        outline = str(bp)

        for bp_in in GO_container[bp]:
            outline += ','+str(bp_in)

        print(outline,file=wf)

    wf.close()

    return GO_container

def MinMaxGOsets(GO_embedding,GO_container,GO_BP_names):
    """
    Identifies the minimum gene sets that do not contain any other gene sets and the maximum gene sets that
    are not contained in any other gene sets.

    Parameters
    ----------
    GO_embedding: dict
        Dict with the index of gene ontology set name as keys and the list of indices of the gene
        sets that contains the key gene set.

    GO_container: dict
        Dict with the index of gene ontology set name as keys and the list of indices of the gene sets
        that are contained in the key gene set.

    GO_BP_names: list
        Names of gene ontology gene sets.
    """

    wf = open('minimum_BP_set.csv','w')

    for bp in GO_container.keys():
        if len(GO_container[bp])==0:
            print(GO_BP_names[bp],file=wf)

    wf.close()

    wf = open('maximum_BP_set.csv','w')

    for bp in GO_embedding.keys():
        if len(GO_embedding[bp])==0:
            print(GO_BP_names[bp],file=wf)

    wf.close()

def readGOBPs(GO_directory):
    """
    Reads the gene sets associated with Gene Ontology Biological Processes.

    Parameters
    ----------
    GO_directory:
        Directory containing the GO data files from `Moleculary Signatures Database <https://www.gsea-msigdb.org/gsea/msigdb/>`_.
        This directory contains a set of .csv files with GOBP names as filenames containing the list of associated genes.

    Returns
    -------
    GO_BPs: dict
        Dict with GOBP names as keys and the list of associated genes as the dictionary entry.
    """
    os.chdir(GO_directory)

    all_GO_files = glob('*.csv')

    GO_BPs = {}
    GO_BPs_dict = {}

    for f in all_GO_files:
        k = f.strip('.csv')

        GO_BPs[k] = pd.read_csv(f,header=None)

        GO_BPs_dict[k] = pd.read_csv(f,header=None)[0].to_list()

    return GO_BPs, GO_BPs_dict

def total_genes(GOBPs):
    """
    Determines the list of unique genes across all the GO gene sets.

    Parameters
    ----------
    GO_BPs: dict
        Dict with GOBP names as keys and the list of associated genes as the dictionary entry.

    Returns
    -------
    unique_genes: list
        Name of genes present in the gene set database.
    """

    gsets = list(GOBPs.keys())

    union = pd.DataFrame(GOBPs[gsets[0]])

    for g in gsets[1:]:
        union = union.merge(GOBPs[g],how='outer').drop_duplicates([0])

    return union[0].to_list()

def p_adjust_bh(p):
    """
    Benjamini-Hochberg p-value correction for multiple hypothesis testing.

    Parameters
    ----------
    p: array_like
        A list or array of p-values, from Fisher's exact text, for multiple hypothesis.

    Returns
    -------
    q: array_like
        Adjusted positive False Discovery Rate (pFDR).
    """
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]

    by_orig = by_descend.argsort()

    steps = float(len(p)) / np.arange(float(len(p)), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))

    return q[by_orig]

def compute_q_values(p_values,go_names,go_tags,alpha=0.01,return_all=False):
    """
    Benjamini-Hochberg p-value correction for multiple hypothesis testing.

    Parameters
    ----------
    p_values: array_like
        A list or array of p-values, from Fisher's exact text, for multiple hypothesis.

    Returns
    -------
    q_values: array_like
        Sorted positive False Discovery Rate corrected for multiple hypothesis testing.

    p_values_go: array_like
        Names of the associated gene ontology biological processes.
    """

    tmp = sorted(zip(go_tags, go_names, p_values),key=lambda x: x[2])
    p_values_go_tags = [x[0] for x in tmp]
    p_values_go = [x[1] for x in tmp]
    p_values = [x[2] for x in tmp]

    p_values_arr = np.array(p_values)
    q_values = p_adjust_bh(p_values)
    rej, qs = fdrcorrection(p_values,alpha)

    q_tags = list(compress(p_values_go_tags,rej))
    q_gos = list(compress(p_values_go,rej))
    q_vs = qs[rej]

    #_, q_values = qvalue(p_values,threshold=0.01)

    if return_all==False:
        return q_tags, q_gos, q_vs
    else:
        return p_values_go_tags, p_values_go, q_values

def compute_p_values(sources,GO_BPs,interaction_set,total_genes,minimum_GOBP=False,size_threshold=np.inf,full=False):
    """
    Benjamini-Hochberg p-value correction for multiple hypothesis testing.

    Parameters
    ----------
    sources: list
        Names of factors that are causing the information transfer in the network.

    minimum_GOBP: list
        Names of gene sets for biological processes at the lowest level, i.e. these sets do not contain other gene sets.

    GO_BPs: array_like
        Gene sets for Gene Ontology Biological Processes.

    interaction_set: dict
        Set of genes receiving information from the sources.

    total_genes: int
        Total number of unique genes across all gene sets.

    Returns
    -------
    go_names: dict
        Gene Ontology Biological Processes that are over-represented by the sources.

    p_values: dict
        Fisher's exact test p-value for Gene Ontology over-represenation analysis.
    """

    go_tags, go_names, p_values = {}, {}, {}

    if minimum_GOBP==False:
        all_gene_sets = list(GO_BPs.keys())
    else:
        all_gene_sets = minimum_GOBP

    for s_g in sources:
        len_interaction = int(interaction_set[s_g][0].count())

        go_tags[s_g] = []
        go_names[s_g] = []
        p_values[s_g] = []

        #for go in minimum_GOBP:
        for go in all_gene_sets:
            if int(GO_BPs[go].count())<=size_threshold:
                intersection = pd.merge(GO_BPs[go], interaction_set[s_g], how='inner').drop_duplicates([0])

                len_intersection = int(intersection[0].count())

                if len_intersection>0:
                    len_gene_BP = int(GO_BPs[go].count())

                    pvalue = hg.sf(len_intersection-1, total_genes, len_gene_BP, len_interaction)

                    p_values[s_g].append(pvalue)
                    go_names[s_g].append(go.replace('_',' ').replace('GOBP ',''))
                    go_tags[s_g].append(go)
                    #if s_g=='SARS-CoV2 nsp7':
                    #    print(s_g,go,len_intersection,pvalue,intersection[0])

    return go_tags, go_names, p_values

    def read_embedding(filename):
        wf = open(filename)
        all_lines = wf.readlines()
        wf.close()

        embed_idx = {}

        for l in all_lines:
            all_values = l.rstrip('\r\n').split(',')

            embed_idx[int(all_values[0])] = []

            for k in all_values[1:]:
                embed_idx[int(all_values[0])].append(int(k))

        return embed_idx
