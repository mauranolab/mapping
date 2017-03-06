import numpy as np
import collections
import pandas as pd 
import os
#This takes in the barcode file and counts the number of occurences of each barcode.
#It then uses the directional adjacency https://cgatoxford.files.wordpress.com/2015/08/schematic_25-e1443714121688.png
#Most of this code is taken from https://www.cgat.org/downloads/public/training/notebooks/umi_simulation/Simulating_umi_deduping.html In [18]: from the Jupyter notebook

#I have only tested on a subset of the barcodes 
#head -10000 BS00067A-5xBGlo_K562d4_2hDpn_iPCR.barcodes.txt > first100K.txt
with open('first10K_BS74A') as f:
    firstColumn = [line.split('\t')[0] for line in f]
		
	

myUMI = collections.Counter(firstColumn)

#This command counts the number of each unique adapter and sorts accordingly  

#From a subset of 100000
#Counter({'CTAGTTGTGGGATC': 1045, 'TTGAAGTCGGAAGA': 293, 'TCGGAAGAGATCGG': 165, 'ACAGTTGAAGTCGG': 47, 'GTTTCATACCACTG': 47, 'AGAACGCTTCAGAA': 46, 'TTAGTTTATTTTCC': 45, 'AGTGTATCGTATTG': 44, 'TCAACATTTACAGT': 43, 'AAAATAAACAGTTG': 43, 'CGAAGCGTGGTACA': 42, 'TCTGCTTTGATGTA': 42, 'TTACCCTTAAGCTA': 42, 'ATTGGTCTGTGTTT': 41, 'GTACAGTTGAAGTC': 4.......

def edit_dist(first, second):
    ''' returns the edit distance/hamming distances between
    its two arguements '''
    dist = sum([not a == b for a, b in zip(first, second)])
    return dist
		


def breadth_first_search(node, adj_list):
    searched = set()
    found = set()
    queue = set()
    queue.update((node,))
    found.update((node,))
    
    while len(queue)>0:
        node=(list(queue))[0]
        found.update(adj_list[node])
        queue.update(adj_list[node])
        searched.update((node,))
        queue.difference_update(searched)
            
    return found
		


def dedup_dir_adj(Counter):
    def get_adj_list_directional_adjacency(umis, counts):
        return {umi: [umi2 for umi2 in umis if edit_dist(umi, umi2) == 1 and
                      counts[umi] >= (counts[umi2]*2)-1] for umi in umis}  
    def get_connected_components_adjacency(graph, Counter):
        found = list()
        components = list()
        for node in sorted(graph, key=lambda x: Counter[x], reverse=True):
            if node not in found:
                component = breadth_first_search(node, graph)
                found.extend(component)
                components.append(component)
        return components
    def remove_umis(adj_list, cluster, nodes):
        '''removes the specified nodes from the cluster and returns
        the remaining nodes '''
        # list incomprehension: for x in nodes: for node in adj_list[x]: yield node
        nodes_to_remove = set([node
                               for x in nodes
                               for node in adj_list[x]] + nodes)
        return cluster - nodes_to_remove
    def reduce_clusters_directional_adjacency(adj_list, clusters, counts):
        n = 0
        for cluster in clusters:
            n+=1
        return n
    adj_list = get_adj_list_directional_adjacency(Counter.keys(), Counter)
    clusters = get_connected_components_adjacency(adj_list, Counter)
    count = reduce_clusters_directional_adjacency(adj_list, clusters, Counter)
    return clusters
    return count
    return adj_list
    
	

dedup_dir_adj(myUMI)