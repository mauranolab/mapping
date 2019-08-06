#!/bin/env python3.5
from sys import argv
import sys
import argparse
import os
import csv
import networkx as nx
#import matplotlib.pyplot as plt
#from networkx import drawing
#from networkx.drawing.nx_pylab import draw_networkx



###Command line arguments
version="1.0"

parser = argparse.ArgumentParser(prog = "genotypeClones.py", description = "Graph approach to generate unified list of which BCs are present in which cells", allow_abbrev=False)
parser.add_argument('--inputfilename', action='store', help='input filename. Format: tab-delimited with barcode sequences, ')
parser.add_argument('--outputfilename', action='store', help='Tab-delimited list of clones and the cells/BCs they include')

parser.add_argument("--verbose", action='store_true', default=False, help = "Verbose mode")
parser.add_argument('--version', action='version', version='%(prog)s ' + version)


try:
    args = parser.parse_args()
except argparse.ArgumentError as exc:
    print(exc.message, '\n', exc.argument, file=sys.stderr)
    sys.exit(1)


inputfilename = args.inputfilename
if inputfilename=="-":
    inputfile = sys.stdin
else:
    inputfile = open(inputfilename, 'r') 

outputfilename = args.outputfilename
if outputfilename=="-":
    outfile = sys.stdout
else:
    outfile = open(outputfilename, 'w')

wr = csv.DictWriter(outfile, delimiter='\t', lineterminator=os.linesep, skipinitialspace=True, fieldnames=['clone', 'cellBC', 'BC', 'count', 'ncells', 'nBCs'])
wr.writeheader()


input_data = inputfile.readlines()
input_data = [line.rstrip("\n").split('\t') for line in input_data]

#Bipartite graph between BCs and cells
G = nx.Graph()
for line in input_data:
    bc = line[0]
    cellbc = line[1]
    count = line[2]
    
    G.add_node(bc, type="BC")
    G.add_node(cellbc, type="cell")
    G.add_edge(bc, cellbc, weight=count)


print("[genotypeClones] Read " + str(len([x for x in G.nodes() if G.nodes[x]['type'] == 'cell'])) + " unique cells", file=sys.stderr)
print("[genotypeClones] Read " + str(len([x for x in G.nodes() if G.nodes[x]['type'] == 'BC'])) + " unique BCs", file=sys.stderr)


#Iterate over connected components (clones) and print out all neighbors
cloneid = 0
totalCells = 0
totalBCs = 0
totalCount = 0
while(len(G.nodes()) > 0):
    cloneid += 1
    nodes = nx.algorithms.components.node_connected_component(G, list(G.nodes)[0])
    #Get all edges between these nodes
    edges = [x for x in G.edges if x[0] in nodes and x[1] in nodes]
    
    cells = [x for x in nodes if G.nodes[x]['type'] == 'cell']
    bcs = [x for x in nodes if G.nodes[x]['type'] == 'BC']
    #Num UMIs supporting total clone/BCs
    count = sum([int(G.edges[x]['weight']) for x in edges])
    G.remove_nodes_from(nodes)
    
    totalCells += len(cells)
    totalBCs += len(bcs)
    totalCount += count
    
    wr.writerow({ 'clone': 'clone-' + str(cloneid), 'cellBC': ",".join(cells), 'BC': ",".join(bcs), 'count': count, 'ncells': len(cells), 'nBCs': len(bcs) })


print("[genotypeClones] Identified " + str(cloneid) + " unique clones, with an average of ", round(totalCells/cloneid, 2), " cells, ", round(totalBCs/cloneid,2), " BCs, ", round(totalCount/cloneid, 2), " reads", sep="", file=sys.stderr)


inputfile.close()
outfile.close()
