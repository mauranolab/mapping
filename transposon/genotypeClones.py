#!/bin/env python3.5
from sys import argv
import sys
import argparse
import os
import csv
import networkx as nx
import statistics

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import matplotlib.cm as cmx
from networkx import drawing
from networkx.drawing.nx_pylab import draw_networkx
from networkx import edge_betweenness_centrality as betweenness


#BUGBUG don't think connected_component_subgraphs presents subgraphs in deterministic order


def printGraph(G, filename, edge_color='weight'):
    #print("[genotypeClones] Printing graph ", filename, sep="", file=sys.stderr)
    nodeColorDict = {'BC': 'darkblue', 'cell': 'darkred'}
#    node_sizes = [node[1]*25000 for node in G.nodes.data('weight')]
    node_colors = [mcolors.to_rgba(nodeColorDict[node[1]]) for node in G.nodes.data('type')]
    edge_weights = [edge[2] for edge in G.edges.data('weight')]
    
    kwds = {'edgelist': G.edges,
        'font_size': 10,
        'node_color': node_colors,
        'node_size': 25,
        'width': 2.8,
#        'edge_vmin': 0, #min/max edge weight for color scale
#        'edge_vmax': 0.5
        }
    
    edge_colormap = plt.get_cmap('Blues')
    if edge_color=="color":
        edge_colors = [mcolors.to_rgba(edge[2]) for edge in G.edges.data('color')]
        kwds['edge_color'] = edge_colors
    elif edge_color=="weight":
        kwds['edge_color'] = edge_weights
        kwds['edge_cmap'] = edge_colormap
    else:
        print("ERROR", edge_color)
    
    fig = plt.figure()
    fig.set_size_inches(7, 7)
    
    #James originally had kamada_kawai_layout but seems to perform badly
    #https://stackoverflow.com/questions/14283341/how-to-increase-node-spacing-for-networkx-spring-layout shows how to space out the components a bit
    drawing.nx_pylab.draw_networkx(G, pos=nx.spring_layout(G,k=0.25,iterations=50, weight='weight'), **kwds)
    
    sm = plt.cm.ScalarMappable(cmap=edge_colormap, norm=mcolors.NoNorm(vmin=0, vmax=100))
    sm._A = []
    plt.colorbar(sm, shrink=0.7)
    
    plt.title(filename + " ({} nodes)".format(len(G.nodes)), fontsize=14, x=0.5, y=1.02)
    
    plt.savefig(filename)
    
    matplotlib.pyplot.close(fig)


#Remove nodes without edges
def pruneOrphanNodes(G):
    #BUGBUG removes everything?
    isolates = [node for node in nx.isolates(G)]
    print("[genotypeClones] Dropped ", len([node for node in isolates if G.nodes[node]['type']=='cell']), " unconnected cells and ", len([node for node in isolates if G.nodes[node]['type']=='BC']), " unconnected BCs from graph", sep="", file=sys.stderr)
    G.remove_nodes_from(isolates)


def pruneEdgesLowPropOfReads(G, minPropOfReads, type='BC'):
    countsremoved = 0
    edgesremoved = 0
    for node in G.nodes:
        if G.nodes[node]['type'] == type:
            nodeweight = G.nodes[node]['weight']
            edges = [edge for edge in G.edges([node]) if G.edges[edge]['weight'] / nodeweight < minPropOfReads ]
            filteredEdgeWeight = sum([G.edges[x]['weight'] for x in edges])
            countsremoved += filteredEdgeWeight
            edgesremoved += len(edges)
            
            #Update node weights
            for edge in edges:
                G.nodes[edge[0]]['weight'] -= G.edges[edge]['weight']
                G.nodes[edge[1]]['weight'] -= G.edges[edge]['weight']
            
            G.remove_edges_from(edges)
    
    print("[genotypeClones] Pruned ", edgesremoved, " edges (", len(G.edges), " left) representing <", minPropOfReads, " of all UMIs for a given ", type, "; ", countsremoved, " UMIs removed", sep="", file=sys.stderr)
    print("[genotypeClones] Average UMIs per edge: ", str(statistics.mean([ G.edges[x]['weight'] for x in G.edges])), sep="", file=sys.stderr)
    
    pruneOrphanNodes(G)


def breakUpWeaklyConnectedCommunities(G, minCentrality, doGraph=False):
    #Break up weakly connected communities
    precloneid = 0
    edgesToDrop = []
    countsremoved = 0
    nCommunities = 0
    #Downside of doing filtering in a separate pass is that it is harder to debug why some clusters aren't broken up
    for subG in nx.connected_component_subgraphs(G):
        bcs = [x for x in subG.nodes if subG.nodes[x]['type'] == 'BC']
        cells = [x for x in subG.nodes if subG.nodes[x]['type'] == 'cell']
        
        nodesToPrune = []
        
        if len(bcs) > 0 and len(cells) > 0:
            precloneid += 1
            
            centrality = betweenness(subG, weight='weight')
            #Start with edges with lowest UMIs
            for edge in sorted(nx.bridges(subG), key=lambda e: subG.edges[e]['weight'], reverse=False):
                subGedgeless = subG.copy()
                subGedgeless.remove_edge(edge[0], edge[1])
                leftNodes = nx.algorithms.components.node_connected_component(subGedgeless, edge[0])
                rightNodes = nx.algorithms.components.node_connected_component(subGedgeless, edge[1])
                
                leftCells = [ node for node in leftNodes if subG.nodes[node]['type'] == 'cell']
                rightCells = [ node for node in rightNodes if subG.nodes[node]['type'] == 'cell']
                
                leftBCs = [ node for node in leftNodes if subG.nodes[node]['type'] == 'BC']
                rightBCs = [ node for node in rightNodes if subG.nodes[node]['type'] == 'BC']
                
                leftReads = sum([ subG.nodes[node]['weight'] for node in leftCells])
                rightReads = sum([ subG.nodes[node]['weight'] for node in rightCells])
                
                #Don't create orphan components with no cells or BCs
                #Make sure we don't remove both edges from a BC node by not pruning >1 edge from any given node
                if edge[0] not in nodesToPrune and edge[1] not in nodesToPrune and len(leftCells) >= 1 and len(rightCells) >= 1 and len(leftBCs) >= 1 and len(rightBCs) >= 1:
                    #Separate if for tunable filters to facilitate tuning
                    #Identify communities by removing bridge edges based centrality metric.
                    #TODO arbitrary read cutoff for edges. Should min BC/cells in each component be a proportion rather than 1?
                    if centrality[edge] > minCentrality and G.edges[edge]['weight'] / min(leftReads, rightReads) <= args.maxpropreads :
                        nCommunities += 1
                        countsremoved += subG.edges[edge]['weight']
                        edgesToDrop.append(edge)
                        subG.nodes[edge[0]]['weight'] -= subG.edges[edge]['weight']
                        subG.nodes[edge[1]]['weight'] -= subG.edges[edge]['weight']
                        subG.edges[edge]['color'] = 'darkred'
                        nodesToPrune.append(edge[0])
                        nodesToPrune.append(edge[1])
                    
                        print("[genotypeClones] ", 'preclone-' + str(precloneid), ' ', edge, " weight:", subG.edges[edge]['weight'], ", L:", str(len(leftNodes)), " (", str(len(leftCells)), " cells, ", leftReads, " reads), R:", str(len(rightNodes)), " (", str(len(rightCells)), " cells, ", rightReads, " reads), centrality:", centrality[edge], sep="", file=sys.stderr)
            
            ###Never implemented Louvain communities (did ok but tended to split up smaller graphs)
            #pip install --upgrade --user python-louvain
            #from community import community_louvain
            #comp = community_louvain.best_partition(subG, weight='weight', resolution=1000)
            #nCommunities = len(set(comp.values()))
            #comp = community.girvan_newman(subG#, most_valuable_edge=most_central_edge)
            #nCommunities = len([c for c in comp])
            
            ###Print graph
            if doGraph and args.printGraph is not None:
                printGraph(subG, args.printGraph + '/preclone-' + str(precloneid) + '.png', edge_color='color')
    
    G.remove_edges_from(edgesToDrop)
    print("[genotypeClones] Created ", nCommunities, " new clones by pruning ", len(edgesToDrop), " edges (", len(G.edges), " left) ", countsremoved, " UMIs removed", sep="", file=sys.stderr)
    
    return nCommunities


###Command line arguments
version="1.0"

parser = argparse.ArgumentParser(prog = "genotypeClones.py", description = "Graph approach to generate unified list of which BCs are present in which cells", allow_abbrev=False)
parser.add_argument('--inputfilename', action='store', help='input filename. Format: tab-delimited with barcode sequences, . First line must be header (unused)')
parser.add_argument('--outputfilename', action='store', help='Tab-delimited list of clones and the cells/BCs they include')

parser.add_argument('--minreads', action='store', type=int, default=2, help='Min UMI filter for input file')
parser.add_argument('--minPropOfReads', action='store', type=float, default=0.15, help='Each BC-cell edge must represent at least this proportion of UMIs for BC')
parser.add_argument('--minCentrality', action='store', type=float, default=0.2, help='Each BC-cell edge must represent at least this proportion of UMIs for BC')
parser.add_argument('--maxpropreads', action='store', type=int, default=0.1, help='Edges joining communities must have fewer than this number of UMIs as proportion of the smaller community they bridge')

parser.add_argument('--printGraph', action='store', type=str, help='Plot a graph for each clone into this directory')

parser.add_argument("--verbose", action='store_true', default=False, help = "Verbose mode")
parser.add_argument('--version', action='version', version='%(prog)s ' + version)


try:
    args = parser.parse_args()
except argparse.ArgumentError as exc:
    print(exc.message, '\n', exc.argument, file=sys.stderr)
    sys.exit(1)

print("[genotypeClones] " + str(args), file=sys.stderr)

if args.printGraph is not None:
    os.makedirs(args.printGraph, exist_ok=True)

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

###Initialize undirected bipartite graph between BCs and cells
#I think python3.6 might obviate need for OrderedGraph but we are trying to make results deterministic.
G = nx.OrderedGraph()
totalReads = 0
#Expects header line
for line in input_data[1:]:
    bc = line[0]
    cellbc = line[1]
    count = int(line[2])
    
    if count >= args.minreads:
        totalReads += count
        
        #Node weight is sum of edge weights
        if bc in G:
            if G.nodes[bc]['type'] != "BC":
                print("WARNING collision between BC and cellBC")
            G.nodes[bc]['weight'] += count
        else:
            G.add_node(bc, type="BC", weight=count)
        
        if cellbc in G:
            if G.nodes[cellbc]['type'] != "cell":
                print("WARNING collision between cellBC and BC")
            G.nodes[cellbc]['weight'] += count
        else:
            G.add_node(cellbc, type="cell", weight=count)
        
        G.add_edge(bc, cellbc, weight=count, color="lightgray")


print("[genotypeClones] Read ", str(len([x for x in G.nodes if G.nodes[x]['type'] == 'cell'])), " unique cells", sep="", file=sys.stderr)
print("[genotypeClones] Read ", str(len([x for x in G.nodes if G.nodes[x]['type'] == 'BC'])), " unique BCs", sep="", file=sys.stderr)
print("[genotypeClones] Read ", str(totalReads), " total UMIs", sep="", file=sys.stderr)
print("[genotypeClones] Initialized ", len(G.edges), " total edges", sep="", file=sys.stderr)

print("[genotypeClones] Average UMIs per BC: ", str(statistics.mean([ int(G.nodes[x]['weight']) for x in G.nodes if G.nodes[x]['type'] == 'BC'])), sep="", file=sys.stderr)
print("[genotypeClones] Average UMIs per BC-cell edge: ", str(statistics.mean([ int(G.edges[x]['weight']) for x in G.edges])), sep="", file=sys.stderr)


###Filtering
#This filter seems more stringent on the individual libraries than the aggregate one
pruneEdgesLowPropOfReads(G, args.minPropOfReads, type='BC')
#TODO this drops a lot of otherwise unconnected BCs, maybe keep?
pruneEdgesLowPropOfReads(G, 0.01, type='cell')
breakUpWeaklyConnectedCommunities(G, minCentrality=args.minCentrality, doGraph=True)
#Do twice to break up some of the bigger graphs since we don't iterate internally, 3x doesn't do anything else
breakUpWeaklyConnectedCommunities(G, minCentrality=args.minCentrality, doGraph=False)

#TODO doublet detection by finding cell nodes that join strongly connected components? Maybe need higher cell density


###Iterate over connected components (clones) and print out all neighbors
cloneid = 0
totalCells = 0
totalCellsSkipped = 0
totalBCs = 0
totalBCsSkipped = 0
totalCount = 0


for subG in nx.connected_component_subgraphs(G):
    bcs = [x for x in subG.nodes if subG.nodes[x]['type'] == 'BC']
    cells = [x for x in subG.nodes if subG.nodes[x]['type'] == 'cell']
    #Sort lists so that print is deterministic
    bcs.sort()
    cells.sort()
    
    cloneid += 1
    
    #Don't output any BCs without cells and vice versa
    #Num UMIs supporting total clone/BCs
    count = sum([subG.edges[x]['weight'] for x in subG.edges])
    totalCells += len(cells)
    totalBCs += len(bcs)
    totalCount += count
    
    if args.printGraph:
        printGraph(subG, args.printGraph + '/clone-' + str(cloneid) + '.png', edge_color='weight')
    
    ret = { 'clone': 'clone-' + str(cloneid), 'cellBC': ",".join(cells), 'BC': ",".join(bcs), 'count': count, 'ncells': len(cells), 'nBCs': len(bcs)}
#    print(ret, file=sys.stderr)
    wr.writerow(ret)


print("[genotypeClones] Identified ", str(cloneid), " unique clones, ", totalCells, " cells, and ", totalBCs, " BCs. Average of ", round(totalCells/cloneid, 2), " cells, ", round(totalBCs/cloneid, 2), " BCs, ", round(totalCount/cloneid, 2), " UMIs", sep="", file=sys.stderr)


inputfile.close()
outfile.close()
