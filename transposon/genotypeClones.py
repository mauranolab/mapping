#!/bin/env python3.5
from sys import argv
import sys
import argparse
import os
import csv
import statistics
import collections

import networkx as nx
from networkx import drawing
from networkx.drawing.nx_pylab import draw_networkx
from networkx import edge_betweenness_centrality as betweenness

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import matplotlib.cm as cmx
#Make pyplot export text as truetype 2 for illustrator
#http://jonathansoma.com/lede/data-studio/matplotlib/exporting-from-matplotlib-to-open-in-adobe-illustrator/
matplotlib.rcParams['pdf.fonttype'] = 42
#Otherwise default is DejaVuSans
#BUGBUG seems to have no effect
matplotlib.rc('font', family='sans-serif')
matplotlib.rc('font', serif='Helvetica')

version="1.0"

#BUGBUG don't think connected_component_subgraphs presents subgraphs in deterministic order
#TODO doublet detection by finding cell nodes that join strongly connected components? Maybe need higher cell density
#TODO sort output


###Initialize undirected bipartite graph between BCs and cells
#1) Node attributes
#    type - BC / cell
#    weight - number of UMIs involving this BC/cell
#2) edge attributes
#    weight - number of UMIs supporting this link
#    color - used to distinguish trimmed edges in plots
def initializeGraphFromInput(inputfilename, minCount):
    if inputfilename=="-":
        inputfile = sys.stdin
    else:
        inputfile = open(inputfilename, 'r') 
    
    try:
        input_data = inputfile.readlines()
        input_data = [line.rstrip("\n").split('\t') for line in input_data]
        
        #We are trying to make results deterministic. I think python3.6 might obviate need for OrderedGraph
        G = nx.OrderedGraph()
        totalReads = 0
        #Expects header line
        #TODO actually use header
        for line in input_data[1:]:
            bc = line[0]
            cellbc = line[1]
            count = int(line[2])
            
            if count >= minCount:
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
    finally:
        inputfile.close()
    
    print("[genotypeClones] Read ", str(len([x for x in G.nodes if G.nodes[x]['type'] == 'cell'])), " unique cells", sep="", file=sys.stderr)
    print("[genotypeClones] Read ", str(len([x for x in G.nodes if G.nodes[x]['type'] == 'BC'])), " unique BCs", sep="", file=sys.stderr)
    print("[genotypeClones] Read ", str(totalReads), " total UMIs", sep="", file=sys.stderr)
    print("[genotypeClones] Initialized ", len(G.edges), " total edges", sep="", file=sys.stderr)
    print("[genotypeClones] Average UMIs per BC: ", str(statistics.mean([ int(G.nodes[x]['weight']) for x in G.nodes if G.nodes[x]['type'] == 'BC'])), sep="", file=sys.stderr)
    print("[genotypeClones] Average UMIs per BC-cell edge: ", str(statistics.mean([ int(G.edges[x]['weight']) for x in G.edges])), sep="", file=sys.stderr)
    
    return G


#Remove nodes without edges (BCs without cells and vice versa)
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


def breakUpWeaklyConnectedCommunities(G, minCentrality, maxPropReads, doGraph=False, verbose=True, graphOutput=None):
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
                    if centrality[edge] > minCentrality and G.edges[edge]['weight'] / min(leftReads, rightReads) <= maxPropReads:
                        nCommunities += 1
                        countsremoved += subG.edges[edge]['weight']
                        edgesToDrop.append(edge)
                        subG.nodes[edge[0]]['weight'] -= subG.edges[edge]['weight']
                        subG.nodes[edge[1]]['weight'] -= subG.edges[edge]['weight']
                        subG.edges[edge]['color'] = 'darkred'
                        nodesToPrune.append(edge[0])
                        nodesToPrune.append(edge[1])
                        
                        if verbose:
                            print("[genotypeClones] ", 'preclone-' + str(precloneid), ' ', edge, " weight:", subG.edges[edge]['weight'], ", L:", str(len(leftNodes)), " (", str(len(leftCells)), " cells, ", leftReads, " reads), R:", str(len(rightNodes)), " (", str(len(rightCells)), " cells, ", rightReads, " reads), centrality:", centrality[edge], sep="", file=sys.stderr)
            
            ###Never implemented Louvain communities (did ok but tended to split up smaller graphs)
            #pip install --upgrade --user python-louvain
            #from community import community_louvain
            #comp = community_louvain.best_partition(subG, weight='weight', resolution=1000)
            #nCommunities = len(set(comp.values()))
            #comp = community.girvan_newman(subG#, most_valuable_edge=most_central_edge)
            #nCommunities = len([c for c in comp])
            
            ###Print graph
            if graphOutput is not None:
                printGraph(subG, filename=printGraph + '/preclone-' + str(precloneid), edge_color='color')
    
    G.remove_edges_from(edgesToDrop)
    print("[genotypeClones] Created ", nCommunities, " new clones by pruning ", len(edgesToDrop), " edges (", len(G.edges), " left) ", countsremoved, " UMIs removed", sep="", file=sys.stderr)
    
    return nCommunities


def printGraph(G, filename=None, fig=None, node_color='type', node_color_dict={'BC': 'darkblue', 'cell': 'darkred'}, edge_color='weight', edge_color_cmap="Blues", with_labels=True):
    #print("[genotypeClones] Printing graph ", filename, sep="", file=sys.stderr)
    # nodeColorDict = 
#    node_sizes = [node[1]*25000 for node in G.nodes.data('weight')]
    node_colors = [mcolors.to_rgba(node_color_dict.get(node[1], "lightgray"))
                   for node in G.nodes.data(node_color)]
    edge_weights = [edge[2] for edge in G.edges.data('weight')]
    
    kwds = {
        'edgelist': G.edges,
        'font_size': 8,
        'nodelist': G.nodes,
        'node_color': node_colors,
        'node_size': 25,
        'width': 2,
        'with_labels': with_labels
#        'edge_vmin': 0, #min/max edge weight for color scale
#        'edge_vmax': 0.5
        }
    
    edge_colormap = plt.get_cmap(edge_color_cmap)
    if edge_color=="color":
        edge_colors = [mcolors.to_rgba(edge[2]) for edge in G.edges.data('color')]
        kwds['edge_color'] = edge_colors
    elif edge_color=="weight":
        kwds['edge_color'] = edge_weights
        kwds['edge_cmap'] = edge_colormap
    else:
        print("ERROR", edge_color)
    
    if fig is None:
        fig = plt.figure()
        fig.set_size_inches(7, 7)
    
    #James originally had kamada_kawai_layout but seems to perform badly
    #https://stackoverflow.com/questions/14283341/how-to-increase-node-spacing-for-networkx-spring-layout shows how to space out the components a bit
    pos = nx.spring_layout(G, k=0.25, iterations=25, weight='weight')
    nx.draw_networkx(G, pos=pos, **kwds)
    
    sm = plt.cm.ScalarMappable(cmap=edge_colormap, norm=mcolors.NoNorm(vmin=0, vmax=100))
    sm._A = []
    plt.colorbar(sm, shrink=0.7)
    ## add node legend
    ax = fig.gca()
    for label in node_color_dict:
        ax.scatter([0], [0], 
            color=mcolors.to_rgba(node_color_dict[label]),
            label=label)
    plt.legend(loc="upper right")
    
    nCells = len([ node for node in G.nodes if G.nodes[node]['type'] == 'cell'])
    nBCs = len([node for node in G.nodes if G.nodes[node]['type'] == 'BC'])
    
    plt.title("({} cells and {} BCs)".format(nCells, nBCs), fontsize=14, x=0.5, y=1.02)
    
    if filename is not None:
        plt.title(os.path.basename(filename) +
            " ({} cells and {} BCs)".format(nCells, nBCs), fontsize=14, x=0.5, y=1.02)
        plt.savefig(filename + '.png')
        plt.savefig(filename + '.pdf')
        plt.close(fig)


###Iterate over connected components to define clones; returns a dict of dicts. The key to the top-level dict is the clone name; the value includes the clone subgraph, and summary info like UMI count, bcs, cells
def identifyClones(G):
    cloneid = 0
    totalCells = 0
    totalCellsSkipped = 0
    totalBCs = 0
    totalBCsSkipped = 0
    totalUMIs = 0
    
    clones = collections.OrderedDict()
    
    for subG in nx.connected_component_subgraphs(G):
        cloneid += 1
        clonename = 'clone-' + str(cloneid).zfill(4)
        
        bcs = [x for x in subG.nodes if subG.nodes[x]['type'] == 'BC']
        cells = [x for x in subG.nodes if subG.nodes[x]['type'] == 'cell']
        #Sort lists so that print is deterministic
        bcs.sort()
        cells.sort()
        
        #Num UMIs supporting total clone/BCs
        umi_count = sum([subG.edges[x]['weight'] for x in subG.edges])
        totalUMIs += umi_count
        totalCells += len(cells)
        totalBCs += len(bcs)
        
        clones[clonename] = { 'clones': subG, 'umi_count': umi_count , 'bcs': bcs, 'cells': cells }
        
    print("[genotypeClones] Identified ", str(cloneid), " unique clones, ", totalCells, " cells, and ", totalBCs, " BCs. Average of ", round(totalCells/cloneid, 2), " cells, ", round(totalBCs/cloneid, 2), " BCs, ", round(totalUMIs/cloneid, 2), " UMIs per clone", sep="", file=sys.stderr)
    
    return clones


#Print three output files and optionally graphs for each clone
#TODO make each output file optional
def writeOutputFiles(clones, output, outputlong, outputwide, graphOutput=None):
    try:
        #Prepare file IO
        longoutfile = open(outputlong, 'w')
        longwr = csv.DictWriter(longoutfile, delimiter='\t', lineterminator=os.linesep, skipinitialspace=True, fieldnames=['BC', 'clone', 'count', 'nCells'])
        longwr.writeheader()
            
        outputfilename = output
        if outputfilename=="-":
            outfile = sys.stdout
        else:
            outfile = open(outputfilename, 'w')
        outwr = csv.DictWriter(outfile, delimiter='\t', lineterminator=os.linesep, skipinitialspace=True, fieldnames=['BC', 'cellBC', 'count', 'clone'])
        outwr.writeheader()
        
        wideoutfile = open(outputwide, 'w')
        widewr = csv.DictWriter(wideoutfile, delimiter='\t', lineterminator=os.linesep, skipinitialspace=True, fieldnames=['BCs', 'cellBCs', 'clone', 'count', 'nedges', 'nBCs', 'ncells'])
        widewr.writeheader()
        
        for clonename, clones in clones.items():
            subG = clones['clones']
            count = clones['umi_count']
            bcs = clones['bcs']
            cells = clones['cells']
            
            if graphOutput:
                printGraph(subG, filename=graphOutput + '/' + clonename, edge_color='weight')
            
            for bc in bcs:
                longwr.writerow({ 'BC': bc, 'clone': clonename, 'count': sum([subG.edges[x]['weight'] for x in subG.edges([bc])]), 'nCells': len(subG.edges([bc]))})
                
                #Start with edges with highest UMIs
                for edge in sorted(subG.edges([bc]), key=lambda e: subG.edges[e]['weight'], reverse=True):
                    if subG.nodes[edge[0]]['type'] == 'BC':
                        cellBC = edge[1]
                    else:
                        cellBC = edge[0]
                    outwr.writerow({ 'BC': bc, 'cellBC': cellBC, 'clone': clonename, 'count': subG.edges[edge]['weight'] })
            
            widewr.writerow({ 'BCs': ",".join(bcs), 'cellBCs': ",".join(cells), 'clone': clonename, 'count': umi_count, 'nedges': len(subG.edges), 'nBCs': len(bcs), 'ncells': len(cells) })
    
    finally:
        outfile.close()
        longoutfile.close()
        wideoutfile.close()


###Utility functions not currently used in command line operation
## Given a set of nodes, it expands their neighborhood to up a certain degree of distance, which allows exploring specific node proximities
def expandNeighborhood(G, seednodes, degree = 1, to_skip = set()):
    #Initialize neighborhood with seed nodes
    neighborhood = set(seednodes)
    to_skip = set(to_skip)
    for _ in range(degree):
        tmp = neighborhood.union(
            set(n for x in neighborhood for n in nx.neighbors(G, x)))
        neighborhood = tmp - to_skip
    return neighborhood


## Assign assigns values to nodes based on a dictionary of node -> value, it can be used to assign nodes transfections for example
def assignToNodes(G, key, annotationdict, default=None):
    for n in G.nodes:
        val = annotationdict.get(n, default)
        G.nodes[n][key] = val


## Convert bipartite graph into a Jaccard index graph of the cells
def toJaccard(G):
    def jaccard(G, u, v):
        unbrs = set(G[u])
        vnbrs = set(G[v])
        return float(len(unbrs & vnbrs)) / len(unbrs | vnbrs)
    
    cells = [n for n in G.nodes if G.nodes[n]['type'] == 'cell']
    return nx.bipartite.generic_weighted_projected_graph(G, cells, weight_function=jaccard)


###Command line operation
if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog = "genotypeClones.py", description = "Graph approach to generate unified list of which BCs are present in which cells", allow_abbrev=False)
    parser.add_argument('--inputfilename', action='store', help='input filename. Format: tab-delimited with barcode sequences. First line must be header (unused)')
    parser.add_argument('--outputlong', action='store', help='Tab-delimited list of BC counts totalled per clone')
    parser.add_argument('--outputwide', action='store', help='Tab-delimited list of clones and the cells/BCs they include')
    parser.add_argument('--output', action='store', help='Tab-delimited list of clone, cell, BC links - filtered version of barcode.counts.byCell file. Can be - for stdout.')
    
    parser.add_argument('--minreads', action='store', type=int, default=2, help='Min UMI filter for input file')
    parser.add_argument('--minPropOfBCReads', action='store', type=float, default=0.15, help='Each BC-cell edge must represent at least this proportion of UMIs for BC')
    parser.add_argument('--minPropOfCellReads', action='store', type=float, default=0.01, help='Each BC-cell edge must represent at least this proportion of UMIs for cell')
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
    
    
    ###Initialize graph, filter, write output
    G = initializeGraphFromInput(inputfilename=args.inputfilename, minCount=args.minreads)
    
    #This filter seems more stringent on the individual libraries than the aggregate one
    pruneEdgesLowPropOfReads(G, args.minPropOfBCReads, type='BC')
    #TODO this drops a lot of otherwise unconnected BCs, maybe keep?
    pruneEdgesLowPropOfReads(G, args.minPropOfCellReads, type='cell')
    breakUpWeaklyConnectedCommunities(G, minCentrality=args.minCentrality, maxPropReads=args.maxpropreads, verbose=args.verbose, graphOutput=None)
    #Do twice to break up some of the bigger graphs since we don't iterate internally, 3x doesn't do anything else
    breakUpWeaklyConnectedCommunities(G, minCentrality=args.minCentrality, maxPropReads=args.maxpropreads, verbose=args.verbose, graphOutput=args.printGraph)
    
    clones = identifyClones(G)
    writeOutputFiles(clones, args.output, args.outputlong, args.outputwide, args.printGraph)
