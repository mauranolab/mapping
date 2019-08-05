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

G = nx.Graph()

#Bipartite graph between BCs and cells
for line in input_data:
    bc = line[0]
    cellbc = line[1]
    count = line[2]
    
    G.add_node(bc, type="BC")
    G.add_node(cellbc, type="cell")
    G.add_edge(bc, cellbc, weight=count)

#Iterate over cell nodes and print out all neighbors
cloneid = 0
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
    wr.writerow({ 'clone': 'clone-' + str(cloneid), 'cellBC': ",".join(cells), 'BC': ",".join(bcs), 'count': count, 'ncells': len(cells), 'nBCs': len(bcs) })
    #BUGBUG printing some weird output


print("[genotypeClones] Identified " + str(cloneid) + " unique clones", file=sys.stderr)


inputfile.close()
outfile.close()
