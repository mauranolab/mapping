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




inputfile = open('SortMerged100K', 'r')###TEST


input_data = inputfile.readlines()
input_data = [line.rstrip().split('\t') for line in input_data]
## open file handle 


bcColNum = 0
myUMIcounts = collections.Counter([line[bcColNum] for line in input_data])


##def get_adj_list_directional_adjacency(umis, counts):
##adj_list = get_adj_list_directional_adjacency(Counter.keys(), Counter)
#
#UmiEncoded = [(umi, umi.encode()) for umi in myUMIcounts.keys()] #TODO encode all umis form start in list
##print(UmiEncoded)
##print(type(UmiEncoded))
##print(umis)
#Umi_adj_Dict = {} #Create directory TODO encode all umis form start in list
#umisKeys = list(myUMIcounts.keys())
#for i in range(1,len(myUMIcounts.keys())): #TODO len of dictionary  
#	umi = umisKeys[i]#TODO get keys from dict_keys
#	#umi=umis[i]
#	#print(umi)
#	umiCounts = myUMIcounts[umi]
#	#print(umiCounts)
#	#umiCode=umi.encode()
#	umiCode = UmiEncoded[i]
#	nearbyUmis = [] #add Umi as key in list
#	for j in range(1, len(myUMIcounts.keys())): #List comprehension here
#		umi2 = umisKeys[j]
#		#print(umi2)
#		umi2Code = UmiEncoded[j]
#		if umi != umi2 and umiCounts >= (myUMIcounts[umi2]*2)-1 and edit_distance(umiCode, umi2Code) == 1:
#			nearbyUmis.append(umi2) #If this take time, preallocate
#	Umi_adj_Dict[umi] = nearbyUmis
#	print(Umi_adj_Dict)
#
#return Umi_adj_Dict
#
#
#

umis = myUMIcounts.keys()
counts = myUMIcounts

#
#UmisEncoded = [(umi, umi.encode(), counts[umi],(counts[umi]*2)-1) for umi in umis]  #NEW
#Umi_adj_Dict = {} #Create directory 
##umisKeys = list(umis) #Convert the dynamic dict.keys object into indexed list 
#for i in range(0,len(umis)):
#	umi = UmisEncoded[i][0]#get keys from dict_keys
#	umiCounts = UmisEncoded[i][2]
#	umiEncoded = UmisEncoded[i][1]
#	nearbyUmis =[umi2[0]  for umi2 in UmisEncoded if edit_distance(umiEncoded, umi2[1]) == 1 and umiCounts >= (umi2[3] #NEW
#	Umi_adj_Dict[umi] = nearbyUmis
#	
#print(Umi_adj_Dict)
#
#
UmisEncoded = [(umi, umi.encode(), counts[umi],(counts[umi]*2)-1) for umi in umis]  #NEW
Umi_adj_Dict = {} #Create directory 
#umisKeys = list(umis) #Convert the dynamic dict.keys object into indexed list 
for i in range(0,len(umis)):
	umi = UmisEncoded[i][0]#get keys from dict_keys
	umiCounts = UmisEncoded[i][2]
	umiEncoded = UmisEncoded[i][1]
	nearbyUmis =[umi2[0]  for umi2 in UmisEncoded if edit_distance(umiEncoded, umi2[1]) == 1 and umiCounts >= (umi2[3]] #NEW
	Umi_adj_Dict[umi] = nearbyUmis
	
print(Umi_adj_Dict)


apa ='AACG', 'GTAC', 'ATCC', 'ATCT', 'GCAC', 'AACC', 'TACC', 'GAAC', 'ATAC'
for j in range(len(apa)):
	for i in range(len(apa)):
		print(apa[j], apa[i],levenshtein(apa[j],apa[i]))





def get_connected_components_adjacency(graph, Counter):
       found = list()
       components = list()
       for node in sorted(graph, key=lambda x: Counter[x], reverse=True):
              print(node)
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
              process_lines(input_data[startRow:(index-1)], wr)#Inputs multiple lines into process_lines


#TODO --groupcol only works if it's 0
if groupcol is None:
       process_lines(input_data, wr)

groupcol=0
startRow = 0
lastGroup = None
for index in range(len(input_data)):
       line = input_data[index]#Gives the line number from input_data
       #print(line)#TEST
       curGroup = line[groupcol]#Assigns the barcode to curGroup for each line 
       #print(curGroup)#TEST if groupcol 1 it prints the FASTQ header fir the read
       if curGroup != lastGroup:#Checks if current group is the same as last group
              if index>0:# if index is over the first line
                     print(lastGroup, startRow, index)
                     #process_lines_byGroup(lastGroup, startRow, index, wr)
              lastGroup = curGroup
              startRow = index
       #Do last BC (index+1 because this would be done one past the last iteration)
       process_lines_byGroup(lastGroup, startRow, index+1, wr)


#Do last BC (index+1 because this would be done one past the last iteration)
process_lines_byGroup(lastGroup, startRow, index+1, wr)


print("\nAll barcodes have been deduplicated", file=sys.stderr)




####
#test for myUMIcounts = collections.Counter([line[2] for line in input_data[99770:(99836-1)]])
####
myUMIcounts = collections.Counter([line[2] for line in input_data[99770:(99836-1)]])




found = list()
components = list()
for node in sorted(Umi_adj_Dict, key=lambda x: myUMIcounts[x], reverse=True):
       print(node)
       if node not in found:
              component = breadth_first_search(node, myUMIcounts)
              found.extend(component)
              components.append(component)


def breadth_first_search(node, adj_list):
       searched = set()
       found = set()
       queue = set()
       queue.update((node,))
       found.update((node,))
       
       while len(queue)>0:
              node=(list(queue))[0]
              found.update(Umi_adj_Dict[node])
              queue.update(Umi_adj_Dict[node])
              searched.update((node,))
              queue.difference_update(searched)
                     
       return found
       
       
       
       





#TEST 

apa ='AACG', 'GTAC', 'ATCC', 'ATCT', 'GCAC', 'AACC', 'TACC', 'GAAC', 'ATAC'

apa = [('TCGGAAGAGCACATGT', 'TCGGAAGAGCAAACGT', 'TCGGAAGAGCACACGT', 'TCGGAAGAGCACACAT', 'TCGGAAGCGCACAAGT', 'TCGGAAGAGCACACGC', 'TCGGAAGAGCACAAGT', 'TCGGAAGAGCACACGA', 'TTGGAAGAGCACACGT'), ('AGACCCTGTTTCTCTA', 'AGACCCTGTTTCTCCA', 'AGACCCTGTTTCCCCA', 'AGACCCTGTTTCTTCA', 'AGACCCTGTTTCCCTA', 'AGACTCTGTTTCTCTA', 'AGACACTGTTTCTCCA', 'ATACCCTGTTTCTCTA'), ('AAGAGCACACGTCTGA', 'AAGAGCTCACGTCTGA', 'AGGAGCACACGTCTGA', 'AAGAGTACACGTCTGA', 'AAGAGCACAAGTCTGA', 'AAGAGCACACGTCTGG', 'GAGAGCACACGTCTGA'), ('ATCGGAAGAGCACACA', 'ATAGGAAGAGCACACG', 'ATCGGAAGAGCACACG', 'ATCGGAAGAGAACACG', 'AACGGAAGAGCACACG'), ('ACTAGTTGTGGGATCT', 'ACTAATTGTGGGATCT', 'GCTAGTTGTGGGATCT', 'ACCAGTTGTGGGATCT', 'CCTAGTTGTGGGATCT', 'ACTCGTTGTGGGATCT', 'TCTAGTTGTGGGATCT'), ('AGACTGTTTCTCTAGC', 'AGACTATTTCTCTAGC')]

for j in range(len(apa)):
	for i in range(len(apa[j])):
		for k in range(len(apa[j])):
			print(apa[j][i], apa[j][k],levenshtein(apa[j][i],apa[j][k]))




