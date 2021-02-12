#!/bin/env python3
from sys import argv
import sys
import re
import argparse
import csv
import collections
import os
import cython
import leven 
from leven import levenshtein
from collections import defaultdict
from umi_tools._dedup_umi import edit_distance

#This takes in the barcode file and counts the number of occurences of each barcode.
#It then uses the directional adjacency https://cgatoxford.files.wordpress.com/2015/08/schematic_25-e1443714121688.png
#Based on https://github.com/CGATOxford/UMI-tools


#TODO Make sure to ignore input lines with "" barcode

###
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


def breadth_first_search2(node, adj_list):
    searched = set()
    queue = set([node])
    
    while len(queue) > 0:
        node = queue.pop()
        neighbors = adj_list.get(node, [])
        searched.add(node)
        queue.update(neighbors)
        queue.difference_update(searched)
    return searched


def dedup_dir_adj(Counter):
        def get_adj_list_directional_adjacency(umis, counts):
            UmisEncoded = [(umi, umi.encode(), counts[umi],(counts[umi]*2)-1) for umi in umis] #Get all barcode information here 
            Umi_adj_Dict = {} #Create directory 
            for i in range(0,len(umis)): #TODO changed 1 to 0, otherwise it skips the first entry
                umi = UmisEncoded[i][0]#get keys from dict_keys
                umiEncoded = UmisEncoded[i][1]
                umiCounts = UmisEncoded[i][2]
                nearbyUmis =[umi2[0] for umi2 in UmisEncoded if umiCounts >= umi2[3] and edit_distance(umiEncoded, umi2[1]) == 1] #
                Umi_adj_Dict[umi] = nearbyUmis
                #print(Umi_adj_Dict)
            return Umi_adj_Dict


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
        
        adj_list = get_adj_list_directional_adjacency(Counter.keys(), Counter)
        clusters = get_connected_components_adjacency(adj_list, Counter)
        return clusters

def dedup_dir_adj2(count):
    ## _get_adj_list_directional
    ## TODO: implement UMI-tools build_substr_idx optimization?
    umis = sorted(count.keys(), key = lambda x: count[x], reverse = True)
    umis_encoded = [umi.encode() for umi in umis]
    adj_list = { umi: [] for umi in umis }
    for i in range(len(umis)):
        for j in range(i + 1, len(umis)):
            if edit_distance(umis_encoded[i], umis_encoded[j]) > 1:
                continue
            if count[umis[i]] >= (2 * count[umis[j]] - 1):
                adj_list[umis[i]].append(umis[j])
            if count[umis[j]] >= (2 * count[umis[i]] - 1):
                adj_list[umis[j]].append(umis[i])
    ## _get_connected_components
    found = set()
    components = list()
    for umi in umis:
        if umi in found:
            continue
        component = breadth_first_search2(umi, adj_list)
        found.update(component)
        components.append(component)
    return components


def iter_indexed(umis, umi_length):
    ## define slices
    cs, r = divmod(umi_length, 2)
    slices = [(0, cs + r), (cs + r, umi_length)]
    ## build index
    substr_idx = collections.defaultdict(lambda: collections.defaultdict(set))
    for idx in slices:
        for i, u in enumerate(umis):
            u_sub = u[slice(*idx)]
            substr_idx[idx][u_sub].add(i)
    ## iterate over index
    for i, u in enumerate(umis):
        neighbors = set()
        for idx, sub_map in substr_idx.items():
            u_sub = u[slice(*idx)]
            neighbors = neighbors.union(substr_idx[idx][u_sub])
        neighbors.difference_update(list(range(i + 1)))
        for nbr in neighbors:
            yield i, nbr


def iter_pairwise(umis):
    for i in range(len(umis)):
        for j in range(i+1, len(umis)):
            yield i, j


def dedup_dir_adj3(count):
    ## _get_adj_list_directional
    umis = sorted(count.keys(), key = lambda x: count[x], reverse = True)
    umis_encoded = [umi.encode() for umi in umis]
    adj_list = { umi: [] for umi in umis }
    if len(umis) > 25:
        it = iter_indexed(umis, len(umis[0]))
    else:
        it = iter_pairwise(umis)
    for i, j in it:
        if edit_distance(umis_encoded[i], umis_encoded[j]) > 1:
            continue
        if count[umis[i]] >= (2 * count[umis[j]] - 1):
            adj_list[umis[i]].append(umis[j])
        if count[umis[j]] >= (2 * count[umis[i]] - 1):
            adj_list[umis[j]].append(umis[i])
    ## _get_connected_components
    found = set()
    components = list()
    for umi in umis:
        if umi in found:
            continue
        component = breadth_first_search2(umi, adj_list)
        found.update(component)
        components.append(component)
    return components

###
def replace_dedup(data, bcColNum, myUMIcounts, deduped_UMI, wr):
    global numLinesProcessed
    global numNonEmptyBClinesProcessed
    global numAmbiguousLines
    
    #Create index of all deduplicated barcodes
    dedupMulti = [s for s in deduped_UMI if len(s) > 1]
    demul_index = defaultdict(list)
    for i,x in enumerate(dedupMulti):
        for y in x :
            demul_index[y].append(i)
    
    for line in data:
        numLinesProcessed += 1
        oldBC = line[bcColNum]
        if oldBC != "":
            numNonEmptyBClinesProcessed += 1
            results = demul_index[oldBC]
            matches = [dedupMulti[x] for x in results]
            if len(matches)>0:
                if len(matches)==1:
                    barCount = dict((k, myUMIcounts[k]) for k in matches[0])
                    if max(barCount, key=barCount.get) == min(barCount, key=barCount.get):
                        newBC = sorted(barCount)[0]
                    else:
                        newBC = max(barCount, key=barCount.get) ##TODO Barcodes with same number replace
                    #print(line)
                if len(matches)>1:
                    newBC = ""
                    numAmbiguousLines += 1
                    #print('Ambigious')
                line[bcColNum] = newBC
        wr.writerows([line])


def process_lines(input_data, wr):
    global numUniqueBCs
    global numUniqueBCsAfterDedup
    #BUGBUG I think numUniqueBCs/numUniqueBCsAfterDedup double count the same BC that shows up in multiple groups
    
    #Checks which barcodes are the most common before deduplication
    bcColNum = col
    myUMIcounts = collections.Counter([line[bcColNum] for line in input_data])
    #if args.verbose:
        #print(myUMIcounts, file=sys.stderr)
    numUniqueBCs += len(myUMIcounts)
    if args.verbose:
        print("\nNumber of unique barcodes:", len(myUMIcounts), file=sys.stderr)
        print("The 10 most common UMIs before directional adjacency correction:", myUMIcounts.most_common(10), file=sys.stderr)
        print('\nDeduping barcodes', file=sys.stderr)
    #Deduplicate the barcodes using directional adjacency
    deduped_UMI = dedup_dir_adj(myUMIcounts)
    numUniqueBCsAfterDedup += len(deduped_UMI)
    if args.verbose:
        print("\nNumber of unique barcodes after deduplication:", len(deduped_UMI), file=sys.stderr)
        print("\nAssociated barcodes with the top 10 UMIs directional adjacency correction:",'\n', deduped_UMI[:10], file=sys.stderr)
        print("\nReplacing deduplicated barcodes", file=sys.stderr)
    replace_dedup(input_data, bcColNum, myUMIcounts, deduped_UMI, wr)


def process_lines_byGroup(group, startRow, index, wr):
    global numLinesProcessed
    global numNonEmptyBClinesProcessed
    
    bcColNum = col
    
    if args.verbose:
        print("Processing ", group, " on lines ", startRow, "-", index-1, file=sys.stderr, sep="")
    if startRow==index-1 or group=="":
        #Only one BC in this group or it's "" so just print
        
        #Not pretty but since we are not calling process_lines, we do need to update the statistics
        for i in range(startRow, index):
            numLinesProcessed += 1
            if input_data[i][bcColNum] != "":
                numNonEmptyBClinesProcessed +=1
        #BUGBUG numUniqueBCs and numUniqueBCsAfterDedup would need to get updated too
        
        #BUGBUG why aren't the subsequent two ranges startRow:(index-1)
        
        if args.verbose:
            print("Only one BC in this group so just print", file=sys.stderr, sep="")
        wr.writerows(input_data[startRow:index]) ###Writes only the unique barcodes
    else:
        #print(input_data[startRow:index], startRow,index)
        process_lines(input_data[startRow:index], wr) #FIX Removed index-1 now it prints the last entry per group


###Argument parsing
version="1.1"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog = "DeDup_adjacencyOutput", description = "Deduplicate and correct barcodes", allow_abbrev=False)
    parser.add_argument('inputfilename', action='store', help='input filename. Format: tab-delimited with barcode sequences (other columns are passed through)')
    parser.add_argument("--col", action='store', type=int, default=1, help = "Column with barcodes in it [%(default)s]")
    parser.add_argument("--groupcols", action='store', type=str, help = "Comma-separated list of column numbers. BCs will only be collapsed for sequences having the same values in these columns. NB: must be sorted on these columns, with priority given to lower-numbered columns. [%(default)s]")
    parser.add_argument("--verbose", action='store_true', default=False, help = "Verbose mode")
    parser.add_argument('--version', action='version', version='%(prog)s ' + version)
    parser.add_argument('-o', '--output', action='store', dest='output',help='Deduplicated barcode file')


    try:
        args = parser.parse_args()
    except argparse.ArgumentError as exc:
        print(exc.message, '\n', exc.argument, file=sys.stderr)
        sys.exit(2)

    col = args.col-1
    if args.groupcols is None:
        groupcols = None
    else:
        groupcols = [int(i)-1 for i in args.groupcols.split(",")]
        groupcols.sort()


    ###Main
    print("[AdjacencyDeDup] Deduplicating and correcting barcodes.", file=sys.stderr)
    print("[AdjacencyDeDup] " + str(args), file=sys.stderr)


    ## open file handles
    if args.inputfilename=="-":
        inputfile = sys.stdin
    else:
        inputfile = open(args.inputfilename, 'r') 

    input_data = inputfile.readlines()
    input_data = [line.rstrip("\n").split('\t') for line in input_data]

    if args.output=="-":
        outfile = sys.stdout
    else:
        outfile = open(args.output, 'w')
    wr = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep, skipinitialspace=True)


    numLinesProcessed = 0
    numNonEmptyBClinesProcessed = 0
    numAmbiguousLines = 0
    numGroupsRead = 0
    numUniqueBCs = 0
    numUniqueBCsAfterDedup = 0

    if groupcols is None:
        process_lines(input_data, wr)
    else:
        startRow = 0
        lastGroup = None
        minCols = max(groupcols)
        for index in range(len(input_data)):
            line = input_data[index]
            if len(line) < minCols: continue # Skip entries too short
            #Make key by concatenating contents from all groupcols specified
            curGroup = "+".join([ str(line[groupcol]) for groupcol in groupcols ])
            #index points to the first record in input_data of the next group (not the one we are processing)
            if curGroup != lastGroup:
                if index > 0:
                    numGroupsRead += 1
                    process_lines_byGroup(lastGroup, startRow, index, wr)
                lastGroup = curGroup
                startRow = index
        if lastGroup is not None:
            #Do last BC (index+1 because this would be done one past the last iteration)
            process_lines_byGroup(lastGroup, startRow, index+1, wr)

    print("[AdjacencyDeDup] All barcodes have been deduplicated (" + str(numLinesProcessed) + " lines processed, " + str(numNonEmptyBClinesProcessed) + " with non-empty BCs)", file=sys.stderr)
    print("[AdjacencyDeDup] Number of unique BCs pre-deduplication:", numUniqueBCs, file=sys.stderr)
    print("[AdjacencyDeDup] Number of unique BCs post-deduplication:", numUniqueBCsAfterDedup, file=sys.stderr)
    if groupcols is not None:
        print("[AdjacencyDeDup] " + str(numGroupsRead) + " unique groups were found", file=sys.stderr)
    print("[AdjacencyDeDup] " + str(numAmbiguousLines) + " reads with ambiguous BCs were masked out", file=sys.stderr)
    #TODO report number unique BCs before/after


    inputfile.close()
    outfile.close()
