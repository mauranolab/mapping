#!/bin/env python3
from sys import argv
import sys
import re
import argparse
import csv
#import numpy as np
import collections
#import pandas as pd 
import os
import cython
import leven 
from leven import levenshtein
from collections import defaultdict
from umi_tools._dedup_umi import edit_distance

#This takes in the barcode file and counts the number of occurences of each barcode.
#It then uses the directional adjacency https://cgatoxford.files.wordpress.com/2015/08/schematic_25-e1443714121688.png


#TODO Make sure to ignore input lines with "" barcode
# edit_distance('abcd'.encode(),'aaad'.encode())

###
#this code is taken from https://www.cgat.org/downloads/public/training/notebooks/umi_simulation/Simulating_umi_deduping.html In [18]: from the Jupyter notebook
def edit_dist(first, second):
       ''' returns the edit distance/hamming distances between
       its two arguements '''
       dist = levenshtein(first, second)
       return dist
		
#def edit_dist(first, second):
#       '''Using Numpy instead to speed up''' 
#              np.dstack((np.array(firstColumn),np.array(firstColumn)))
#       dist = sum([not a == b for a, b in np.dstack((np/arrafirst,second))])
#       return dist


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
       
#get_adj_list_directional_adjacency
#        ''' identify all umis within the hamming distance threshold
#        and where the counts of the first umi is > (2 * second umi counts)-1'''
def dedup_dir_adj(Counter):
       def get_adj_list_directional_adjacency(umis, counts):
#              return {umi: [umi2 for umi2 in umis if counts[umi] >= (counts[umi2]*2)-1 and edit_distance(umi.encode(), umi2.encode()) == 1]
#                                      for umi in umis}

#Create list for Umi1 and append Umi2 to each match             
             Umi_adj_Dict = {} #Create directory
             for umi in umis:
             	umiCounts = counts[umi]
             	umiCode=umi.encode()
             	nearbyUmis = [] #add Umi as key in list
             	[nearbyUmis.append(umi2) for umi2 in umis if umi != umi2 and umiCounts >= (counts[umi2]*2)-1 and edit_distance(umiCode, umi2.encode()) == 1]
             	Umi_adj_Dict[umi] = nearbyUmis
             return Umi_adj_Dict

##Create list for Umi1 and append Umi2 to each match             
#             #print(umis) 
#             #print(counts)
#             Umi_adj_Dict = {} #Create directory
#             for umi in umis:
#             	#print(umi)
#             	umiCounts = counts[umi]
#             	#print(umiCounts)
#             	umiCode=umi.encode()
#             	#print(umiCounts)
#             	nearbyUmis = [] #add Umi as key in list
#             	for umi2 in umis: #List comprehension here
#             		if umi != umi2: 
#             			if umiCounts >= (counts[umi2]*2)-1:
#             				if edit_distance(umiCode, umi2.encode()) == 1:
#             					nearbyUmis.append(umi2) #If this take time, preallocate
#             	Umi_adj_Dict[umi] = nearbyUmis
#             return Umi_adj_Dict
#             	
#
             	

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


###
def replace_dedup(data, bcColNum, myUMIcounts, deduped_UMI, wr):
       #Create index of all deduplicated barcodes
       dedupMulti = [s for s in deduped_UMI if len(s) > 1]
       demul_index = defaultdict(list)
       for i,x in enumerate(dedupMulti):
              for y in x :
                     demul_index[y].append(i)
       
       for line in data:
              oldBC = line[bcColNum]
              results = demul_index[oldBC]
              matches = [dedupMulti[x] for x in results]
              if len(matches)>0 and oldBC!="":
                     if len(matches)==1:
                            barCount = dict((k, myUMIcounts[k]) for k in matches[0])
                            newBC = max(barCount, key=barCount.get)
                            #print(line)
                     if len(matches)>1:
                            newBC=""
                            #print('Ambigious') 
                     line[bcColNum] = line[bcColNum].replace(oldBC, newBC)
              wr.writerows([line])


def process_lines(input_data, wr):
       #Checks which barcodes are the most common before deduplication
       bcColNum = col
       myUMIcounts = collections.Counter([line[bcColNum] for line in input_data])
       
       if args.verbose:
              print("\nNumber of unique barcodes:", len(myUMIcounts), file=sys.stderr)
              print("The 10 most common UMIs before directional adjacency correction:",myUMIcounts.most_common(10), file=sys.stderr)
              print('\nDeduping barcodes', file=sys.stderr)

       #Deduplicate the barcodes using directional adjacency
       deduped_UMI = dedup_dir_adj(myUMIcounts)
       
       if args.verbose:
              print("\nNumber of unique barcodes after deduplication:", len(deduped_UMI), file=sys.stderr)
              print("\nAssociated barcodes with the top 10 UMIs directional adjacency correction:",'\n', deduped_UMI[:10], file=sys.stderr)
              print("\nReplacing deduplicated barcodes", file=sys.stderr)
       
       replace_dedup(input_data, bcColNum, myUMIcounts, deduped_UMI, wr)


def process_lines_byGroup(lastGroup, startRow, index, wr):
       if args.verbose:
              print("Processing ", lastGroup, " on lines ", startRow, "-", index-1, file=sys.stderr, sep="")
       if startRow==index-1 or lastGroup=="":
              #Only one BC in this group (or it's "") so just print
              if args.verbose:
                     print("Only one BC in this group so just print", file=sys.stderr, sep="")
              wr.writerows([input_data[startRow:(index-1)]])
       else:
              process_lines(input_data[startRow:(index-1)], wr)


###Main
parser = argparse.ArgumentParser(prog = "DeDup_adjacencyOutput", description = "Deduplicate and correct barcodes", add_help=True)
parser.add_argument('inputfilename', action='store', help='input filename. Format: tab-delimited with barcode sequences (other columns are passed through)')
parser.add_argument("--col", action='store', type=int, default=1, help = "Column with barcodes in it [%(default)s]")
parser.add_argument("--groupcol", action='store', type=int, help = "BC will only be collapsed for sequences having the same value in this column. NB: must be sorted on this column. [%(default)s]")
parser.add_argument("--verbose", action='store_true', default=False, help = "Verbose mode")
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
parser.add_argument('-o','--output', action='store', dest='output',help='Deduplicated barcode file')


try:
          args = parser.parse_args()
except argparse.ArgumentError as exc:
          print(exc.message, '\n', exc.argument, file=sys.stderr)
          sys.exit(2)

col = args.col-1
if args.groupcol is None:
       groupcol = None
else:
       groupcol = args.groupcol-1


print("Deduplicating and correcting barcodes.\n", file=sys.stderr)
print(args, file=sys.stderr)

if args.inputfilename=="-":
          inputfile = sys.stdin
else:
          inputfile = open(args.inputfilename, 'r') 

input_data = inputfile.readlines()
input_data = [line.rstrip().split('\t') for line in input_data]
## open file handle 

       ## open file handle 
outfile = open(args.output, 'w')
wr = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep, skipinitialspace=True)



if groupcol is None:
       process_lines(input_data, wr)
else:
       startRow = 0
       lastGroup = None
       for index in range(len(input_data)):
              line = input_data[index]
              curGroup = line[groupcol]
              if curGroup != lastGroup:
                     if index>0:
                            process_lines_byGroup(lastGroup, startRow, index, wr)
                     lastGroup = curGroup
                     startRow = index
       #Do last BC (index+1 because this would be done one past the last iteration)
       process_lines_byGroup(lastGroup, startRow, index+1, wr)


print("\nAll barcodes have been deduplicated", file=sys.stderr)
