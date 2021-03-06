#!/bin/env python
from sys import argv
import sys
import re
import argparse
import os
import csv
import statistics
import collections
import gzip
import pickle

import networkx as nx
from networkx import drawing
from networkx.drawing.nx_pylab import draw_networkx

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

version="1.1"

#TODO Sort output to stdout so I don't need mlr

###Initialize undirected bipartite graph between BCs and cells
#1) Node attributes
#    type - BC / cell
#    weight - number of UMIs involving this BC/cell
#2) edge attributes
#    weight - number of UMIs supporting this link
#    color - used to distinguish trimmed edges in plots
def initializeGraphFromInput(inputfilename, minCount):
    print("[genotypeClones] Initializing graph", sep="", file=sys.stderr)
    
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
                        print("[genotypeClones] WARNING collision between BC and cellBC", sep="", file=sys.stderr)
                    G.nodes[bc]['weight'] += count
                else:
                    G.add_node(bc, type="BC", weight=count)
                
                if cellbc in G:
                    if G.nodes[cellbc]['type'] != "cell":
                        print("[genotypeClones] WARNING collision between cellBC and BC", sep="", file=sys.stderr)
                    G.nodes[cellbc]['weight'] += count
                else:
                    G.add_node(cellbc, type="cell", weight=count)
                
                G.add_edge(bc, cellbc, weight=count, color="lightgray")
    finally:
        inputfile.close()
    
    summarizeGraph(G)
    
    return G


def summarizeGraph(G):
    bcs = [x for x in G.nodes if G.nodes[x]['type'] == 'BC']
    cells = [x for x in G.nodes if G.nodes[x]['type'] == 'cell']
    
    totalUMIs = sum([G.edges[e]['weight'] for e in G.edges])
    
    #BUGBUG fails with divide by 0 error if graph is empty
    print("[genotypeClones] summarizeGraph - ", len(G.nodes), " nodes (", len(cells), " cells, ", len(bcs), " BCs), ", len(G.edges), " edges with ", totalUMIs, " UMIs (avg. ", round(totalUMIs/len(bcs), 1), " per BC, and ", round(totalUMIs/len(G.edges), 1), " per edge).", sep="", file=sys.stderr)


def filterNodesFromFile(G, filename, keep=True):
    try:
        if re.search('\.gz$', filename) is not None:
            maskfile = gzip.open(filename, 'rt')
        else:
            maskfile = open(filename, 'r')
        mask_data = maskfile.readlines()
        mask_data = [line.rstrip("\n") for line in mask_data]
        mask_data = set(mask_data) #for O(1) lookup performance below
        maskfile.close()
        
        if keep:
            #Apply whitelist only to BCs
            nodesToRemove = [ node for node in G.nodes if node not in mask_data and G.nodes[node]['type']=='BC' ]
            nodesPresentInGraph =  len([ node for node in G.nodes if node in mask_data ])
        else:
            nodesToRemove = [ node for node in G.nodes if node in mask_data ]
            nodesPresentInGraph = len(nodesToRemove)
        
        ncells = len([x for x in nodesToRemove if G.nodes[x]['type'] == 'cell'])
        nbcs = len([x for x in nodesToRemove if G.nodes[x]['type'] == 'BC'])
        
        print("[genotypeClones] ", sep="", file=sys.stderr)
        print("[genotypeClones] filterNodesFromFile - ", "whitelist " if keep else "blacklist ", filename, " with ", len(mask_data), " entries, matching ", nodesPresentInGraph, " nodes. Removing ", len(nodesToRemove), " nodes (", nbcs, " BCs and ", ncells, " cells).", sep="", file=sys.stderr)
        
        remove_edges(G, G.edges(nodesToRemove))
        G.remove_nodes_from(nodesToRemove)
        
        summarizeGraph(G)
    except FileNotFoundError as e:
        print("[genotypeClones] WARNING Problem opening filter file ", filename, sep="", file=sys.stderr)



#Remove nodes without edges (BCs without cells and vice versa)
def pruneOrphanNodes(G):
    #BUGBUG removes everything?
    isolates = [node for node in nx.isolates(G)]
    print("[genotypeClones] pruneOrphanNodes - Dropped ", len([node for node in isolates if G.nodes[node]['type']=='cell']), " unconnected cells and ", len([node for node in isolates if G.nodes[node]['type']=='BC']), " unconnected BCs from graph", sep="", file=sys.stderr)
    #By definition there are no edges to update
    G.remove_nodes_from(isolates)
    summarizeGraph(G)


#Remove list of edges, keeping node weights in sync
def remove_edges(G, edgesToRemove):
    #Convert to set to speed-up inclusion check
    edgesToRemove = set(edgesToRemove)
    #Also include the edges in the reverse order [(x, y) => (y, x)] to make sure they match the edges emited by G.edges
    edgesToRemove.update([(y, x) for (x, y) in edgesToRemove])
    #Update node weights first
    #Make sure we only subtract umis once per edge, regardless of representation (i.e. node order) in edgesToRemove
    for edge in [e for e in G.edges if e in edgesToRemove]:
        edgeweight = G.edges[edge]['weight']
        G.nodes[edge[0]]['weight'] -= edgeweight
        G.nodes[edge[1]]['weight'] -= edgeweight
    G.remove_edges_from(edgesToRemove)


def identifyEdgesLowPropOfReads(G, minPropOfReads, type='BC'):
    edgesToRemove = set()
    for node in G.nodes:
        if G.nodes[node]['type'] == type:
            nodeweight = G.nodes[node]['weight']
            curEdgesToRemove = [ edge for edge in G.edges([node]) if G.edges[edge]['weight'] / nodeweight < minPropOfReads ]
            for edge in curEdgesToRemove:
                edgesToRemove.add(edge)
    
    countsremoved = sum([G.edges[x]['weight'] for x in edgesToRemove])
    
    print("[genotypeClones] identifyEdgesLowPropOfReads - Identified ", len(edgesToRemove), " edges representing <", minPropOfReads, " of all UMIs for a given ", type, "; ", countsremoved, " UMIs removed", sep="", file=sys.stderr)
    
    return edgesToRemove


def pruneEdgesLowPropOfReads(G, minPropOfBCReads, minPropOfCellReads):
    print("[genotypeClones] ", sep="", file=sys.stderr)
    
    #Make a first pass to identify edges for removal so each edge is assessed against the incoming counts
    edgesToRemove = set()
    edgesToRemove.update(identifyEdgesLowPropOfReads(G, minPropOfReads=minPropOfBCReads, type="BC"))
    edgesToRemove.update(identifyEdgesLowPropOfReads(G, minPropOfReads=minPropOfCellReads, type="cell"))
    
    #Remove edges all at once
    remove_edges(G, edgesToRemove)
    print("[genotypeClones] pruneEdgesLowPropOfReads - Removed ", len(edgesToRemove), " total edges", sep="", file=sys.stderr)
    pruneOrphanNodes(G)

# Within a set of destiny nodes and a source node, return a dictionary mapping key attribute values to the proportion of total UMIs each attribute value represent
def computeKeyRate(G, source, dest, key):
    total = 0
    rate = dict()
    ## Count UMI per key value
    for x in dest:
        value = G.nodes[x][key]
        weight = G.edges[(source, x)]["weight"]
        total += weight
        rate[value] = rate.get(value, 0) + weight
    ## Compute rate by key value
    for value in rate.keys():
        rate[value] /= total
    return rate


def pruneConflictingEdges(G, transfectionKey, maxConflict):
    edges_to_remove = set()
    ## TODO: remove hard-coding of labels to skip
    skip = set(["conflicting", "uninformative", "None"])
    
    print("[genotypeClones] ", sep="", file=sys.stderr)
    
    for cell in [x for x in G.nodes if G.nodes[x]['type'] == "cell"]:
        bcs = [x for (_, x) in G.edges(cell) if G.nodes[x][transfectionKey] not in skip]
        transfection_rate = computeKeyRate(G, cell, bcs, transfectionKey)
        ## Skip cells with BCs from single transfection or no transfection and those with too high conflict rate
        if len(transfection_rate) > 1:
            transfection_majority = max(transfection_rate, key = transfection_rate.get)
            if (1 - transfection_rate[transfection_majority]) > maxConflict:
                continue
            keep = set([transfection_majority]) | skip
            #Remove all edges connecting to BCs not in the majority
            edges_to_remove.update([(cell, x) for x in bcs if G.nodes[x][transfectionKey] not in keep])
    remove_edges(G, edges_to_remove)
    print("[genotypeClones] pruneConflictingEdges - Removed ", len(edges_to_remove), " total edges", sep="", file=sys.stderr)
    pruneOrphanNodes(G)


def pruneConflictingNodes(G, transfectionKey, type = "BC"):
    nodes_to_remove = set()
    ## TODO: remove hard-coding of labels to skip
    skip = set(["conflicting", "uninformative", "None"])
    
    print("[genotypeClones] ", sep="", file=sys.stderr)
    
    for node in [x for x in G.nodes if G.nodes[x]["type"] == type]:
        neighboring_bcs = set()
        if type == "BC":
            #If we are starting from BCs, then we consider BCs connected to neighbor cells
            for cell in G.neighbors(node):
                neighboring_bcs.update(G.neighbors(cell))
        else:
            #If we are starting from cells, then we consider directly connected BCs
            neighboring_bcs.update(G.neighbors(node))
        neighboring_bcs_transfection = set([G.nodes[x][transfectionKey] for x in neighboring_bcs])
        neighboring_bcs_transfection = neighboring_bcs_transfection - skip
        if len(neighboring_bcs_transfection) > 1:
            nodes_to_remove.add(node)
    edges_to_remove = list(G.edges(nodes_to_remove))
    remove_edges(G, edges_to_remove)
    G.remove_nodes_from(nodes_to_remove)
    print("[genotypeClones] pruneConflictingNodes - Removed ", len(edges_to_remove), " edges and ", len(nodes_to_remove), " nodes of type ", type, ".", sep="", file=sys.stderr)
    pruneOrphanNodes(G)


# Compute centrality of bridge edges given the number of nodes on left and right side
def bridgeCentrality(left, right):
    n = left + right
    return (left * right) / (n * (n-1) / 2)
    

#Break up weakly connected communities
def breakUpWeaklyConnectedCommunities(G, minCentrality, maxPropReads, doGraph=False, verbose=True, graphOutput=None):
    precloneid = 0
    edgesToDrop = []
    countsremoved = 0
    nCommunities = 0
    
    print("[genotypeClones] ", sep="", file=sys.stderr)
    
    #Downside of doing filtering in a separate pass is that it is harder to debug why some clusters aren't broken up
    for subG in [G.subgraph(c) for c in sorted(nx.connected_components(G), key=len, reverse=True)]:
        bcs = [x for x in subG.nodes if subG.nodes[x]['type'] == 'BC']
        cells = [x for x in subG.nodes if subG.nodes[x]['type'] == 'cell']
        
        nodesToPrune = set()
        if len(bcs) > 0 and len(cells) > 0:
            precloneid += 1
            #Identify bridges and sort nodes to ensure consistency
            bridges = [tuple(sorted(e)) for e in nx.bridges(subG)]
            #Start with edges with lowest UMIs
            #Include edge vertices in sort key, to ensure consistent order
            for edge in sorted(bridges, key=lambda e: (subG.edges[e]['weight'], e[0], e[1]), reverse=False):
                #Skip border bridges where either node has no other neighbor
                if subG.degree[edge[0]] == 1 or subG.degree[edge[1]] == 1:
                    continue
                #Produce a view of subG, hiding the bridge edge
                subGedgeless = nx.restricted_view(subG, [], [edge])
                leftNodes = nx.algorithms.components.node_connected_component(subGedgeless, edge[0])
                rightNodes = nx.algorithms.components.node_connected_component(subGedgeless, edge[1])
                
                leftCells = [ node for node in leftNodes if subG.nodes[node]['type'] == 'cell' ]
                rightCells = [ node for node in rightNodes if subG.nodes[node]['type'] == 'cell' ]
                
                leftBCs = [ node for node in leftNodes if subG.nodes[node]['type'] == 'BC' ]
                rightBCs = [ node for node in rightNodes if subG.nodes[node]['type'] == 'BC' ]
                
                leftReads = sum([ subG.nodes[node]['weight'] for node in leftCells ])
                rightReads = sum([ subG.nodes[node]['weight'] for node in rightCells ])
                
                #Don't create orphan components with no cells or BCs
                #Make sure we don't remove both edges from a BC node by not pruning >1 edge from any given node
                if edge[0] not in nodesToPrune and edge[1] not in nodesToPrune and len(leftCells) >= 1 and len(rightCells) >= 1 and len(leftBCs) >= 1 and len(rightBCs) >= 1:
                    centrality = bridgeCentrality(len(leftNodes), len(rightNodes))
                    #Separate if for tunable filters to facilitate tuning
                    #Identify communities by removing bridge edges based centrality metric.
                    if centrality > minCentrality and G.edges[edge]['weight'] / min(leftReads, rightReads) <= maxPropReads:
                        nCommunities += 1
                        countsremoved += subG.edges[edge]['weight']
                        edgesToDrop.append(edge)
                        subG.edges[edge]['color'] = 'darkred'
                        nodesToPrune.add(edge[0])
                        nodesToPrune.add(edge[1])
                        
                        if verbose:
                            print("[genotypeClones] ", 'preclone-' + str(precloneid), ' ', edge, " weight:", subG.edges[edge]['weight'], ", L:", str(len(leftNodes)), " (", str(len(leftCells)), " cells, ", leftReads, " reads), R:", str(len(rightNodes)), " (", str(len(rightCells)), " cells, ", rightReads, " reads), centrality:", centrality, sep="", file=sys.stderr)
            
            ###Never implemented Louvain communities (did ok but tended to split up smaller graphs)
            #pip install --upgrade --user python-louvain
            #from community import community_louvain
            #comp = community_louvain.best_partition(subG, weight='weight', resolution=1000)
            #nCommunities = len(set(comp.values()))
            #comp = community.girvan_newman(subG#, most_valuable_edge=most_central_edge)
            #nCommunities = len([c for c in comp])
            
            ###Print graph
            if graphOutput is not None:
                printGraph(subG, filename=graphOutput + '/preclone-' + str(precloneid), edge_color='color', **printGraph_kwds)
    
    remove_edges(G, edgesToDrop)
    print("[genotypeClones] breakUpWeaklyConnectedCommunities Created ", nCommunities, " new clones by pruning ", len(edgesToDrop), " edges (", len(G.edges), " left) ", countsremoved, " UMIs removed", sep="", file=sys.stderr)
    
    return nCommunities


def printGraph(G, filename=None, fig=None, node_color='type', node_color_dict={'BC': 'darkblue', 'cell': 'darkred'}, edge_color='weight', edge_color_cmap="Blues", with_labels=True, with_legend=True):
    #print("[genotypeClones] Printing graph ", filename, sep="", file=sys.stderr)
    # nodeColorDict = 
#    node_sizes = [ node[1]*25000 for node in G.nodes.data('weight') ]
    node_colors = [ mcolors.to_rgba(node_color_dict.get(node[1], "lightgray")) for node in G.nodes.data(node_color) ]
    edge_weights = [ edge[2] for edge in G.edges.data('weight') ]
    
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
        print("[genotypeClones] ERROR", edge_color, sep="", file=sys.stderr)
    
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
    if with_legend:
        ## add node legend
        ax = fig.gca()
        for label in node_color_dict:
            ax.scatter([], [], color=mcolors.to_rgba(
                node_color_dict[label]), label=label)
        plt.legend(loc="upper right")
    
    nCells = len([ node for node in G.nodes if G.nodes[node]['type'] == 'cell'])
    nBCs = len([node for node in G.nodes if G.nodes[node]['type'] == 'BC'])
    
    plt.title("({} cells and {} BCs)".format(nCells, nBCs), fontsize=14, x=0.5, y=1.02)
    
    if filename is not None:
        plt.title(os.path.basename(filename) + " ({} cells and {} BCs)".format(nCells, nBCs), fontsize=14, x=0.5, y=1.02)
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
    maxCloneSize = 0
    
    clones = collections.OrderedDict()
    
    #Go through connected components in descending order of total #BCs + # cells
    for subG in [G.subgraph(c) for c in sorted(nx.connected_components(G), key=len, reverse=True)]:
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
        if len(cells) > maxCloneSize:
            maxCloneSize = len(cells)
        
        clones[clonename] = { 'clone': subG, 'umi_count': umi_count , 'bcs': bcs, 'cells': cells }
    
    print("[genotypeClones] ", sep="", file=sys.stderr)
    print("[genotypeClones] Identified ", str(cloneid), " unique clones, max clone had ", maxCloneSize, " cells. Average of ", round(totalCells/cloneid, 2), " cells, ", round(totalBCs/cloneid, 2), " BCs, ", round(totalUMIs/cloneid, 2), " UMIs per clone", sep="", file=sys.stderr)
    
    return clones


#Print three output files and optionally graphs for each clone:
#    outwr, specified by second positional argument
#    widwr, specified by --outputwide
#    longwr, specified by --outputlong

#TODO make each output file optional
def writeOutputFiles(G, clones, output, outputlong, outputwide, cloneobj, graphOutput=None, transfectionKey=None):
    try:
        #Prepare file IO
        longoutfile = open(outputlong, 'w')
        longwr = csv.DictWriter(longoutfile, delimiter='\t', lineterminator=os.linesep, skipinitialspace=True, fieldnames=['BC', 'BC_transfection', 'clone', 'transfection', 'count', 'nCells'])
        longwr.writeheader()
        
        outputfilename = output
        if outputfilename=="-":
            outfile = sys.stdout
        else:
            outfile = open(outputfilename, 'w')
        outwr = csv.DictWriter(outfile, delimiter='\t', lineterminator=os.linesep, skipinitialspace=True, fieldnames=['BC', 'cellBC', 'count', 'clone'])
        outwr.writeheader()
        
        wideoutfile = open(outputwide, 'w')
        widewr = csv.DictWriter(wideoutfile, delimiter='\t', lineterminator=os.linesep, skipinitialspace=True, fieldnames=['BCs', 'cellBCs', 'clone', 'transfection', 'count', 'nedges', 'nBCs', 'ncells'])
        widewr.writeheader()
        
        if cloneobj is not None:
            pickle.dump(dict(graph = G, clones = clones), open(cloneobj, "wb"))
        
        for clonename, clone in clones.items():
            subG = clone['clone']
            umi_count = clone['umi_count']
            bcs = clone['bcs']
            cells = clone['cells']
            # Gather BC's transfection assignment set when transfectionKey is defined, otherwise defaults to None.
            # If --removeConflictingCells and --removeConflictingBCs are specified, then there will be a single transfection per clone
            transfection = set()
            if transfectionKey is not None:
                #Only record unique transfections, not any other labels
                skip = set(["conflicting", "uninformative", "None"])
                transfection = set(G.nodes[bc][transfectionKey] for bc in bcs if G.nodes[bc][transfectionKey] not in skip)
            if len(transfection) == 0:
                transfection.add("None")
            
            if graphOutput:
                printGraph(subG, filename=graphOutput + '/' + clonename, edge_color='weight', **printGraph_kwds)
            
            widewr.writerow({ 'BCs': ",".join(bcs), 'cellBCs': ",".join(cells), 'clone': clonename, 'count': umi_count, 'nedges': len(subG.edges), 'nBCs': len(bcs), 'ncells': len(cells), 'transfection': ",".join(transfection) })
            
            for bc in bcs:
                # Default BC_Transfection to None when transfectionKey is not defined
                bc_transfection = "None"
                if transfectionKey is not None:
                    bc_transfection = G.nodes[bc][transfectionKey]
                longwr.writerow({ 'BC': bc, 'BC_transfection': bc_transfection, 'clone': clonename, 'transfection': ",".join(transfection), 'count': sum([subG.edges[x]['weight'] for x in subG.edges([bc])]), 'nCells': len(subG.edges([bc]))})
                
                #Get all neighbors for this BC (which must be cellBCs)
                for cellBC in subG[bc]:
                    outwr.writerow({ 'BC': bc, 'cellBC': cellBC, 'clone': clonename, 'count': subG.edges[bc, cellBC]['weight'] })
    
    finally:
        outfile.close()
        longoutfile.close()
        wideoutfile.close()


###Utility functions not currently used in command line operation
## Given a list of nodes, it expands their neighborhood to up a certain degree of distance, which allows exploring specific node proximities
def expandNeighborhood(G, seednodes, degree = 1, to_skip = set()):
    #Initialize neighborhood with seed nodes
    neighborhood = set(seednodes)
    to_skip = set(to_skip)
    for _ in range(degree):
        tmp = neighborhood.union(set(n for x in neighborhood for n in nx.neighbors(G, x)))
        neighborhood = tmp - to_skip
    return neighborhood


## Assign assigns values to nodes based on a dictionary of node -> value, values for key must match the name of a cell or BC node.
def assignToNodes(G, key, annotationdict, default=None):
    hadAnnotation = 0
    missingAnnotation = 0
    for n in G.nodes:
        if n in annotationdict:
            hadAnnotation += 1
            G.nodes[n][key] = annotationdict[n]
        else:
            missingAnnotation += 1
            G.nodes[n][key] = default
    print("[genotypeClones] added annotation to " + str(hadAnnotation) + " nodes; " + str(missingAnnotation) + " missing were assigned \"" + str(default) + "\".", file=sys.stderr)


## Convert bipartite graph into a weighted Jaccard index graph of the cells
def toJaccard(G):
    def jaccard(G, u, v):
        u_bcs = set(G[u])
        v_bcs = set(G[v])
        common_bcs = u_bcs & v_bcs
        common_u_weight = sum([G.edges[(u, bc)]["weight"] for bc in common_bcs])
        common_v_weight = sum([G.edges[(v, bc)]["weight"] for bc in common_bcs])
        return float((common_u_weight + common_v_weight) / (G.nodes[u]["weight"] + G.nodes[v]["weight"]))
    cells = [n for n in G.nodes if G.nodes[n]['type'] == 'cell']
    return nx.bipartite.generic_weighted_projected_graph(G, cells, weight_function=jaccard)


## Unlink cells with low similarity (Jaccard-index)
def pruneLowJaccard(G, min_jaccard):
    J = toJaccard(G)
    low_jaccard_edges = [e for e in J.edges if J.edges[e]["weight"] < min_jaccard]
    edges_to_remove = set()
    ## iterates over cell-pairs with low Jaccard, identifies their shared BCs and for each bc select the edge to u or v with lowest weight to remove
    for (u, v) in low_jaccard_edges:
        bcs = set(G[u]) & set(G[v])
        edges_to_remove.update([(v, bc) if (G.edges[(u, bc)]["weight"] > G.edges[(v, bc)]["weight"]) else (u, bc) for bc in bcs])
    remove_edges(G, edges_to_remove)
    print("[genotypeClones] pruneLowJaccard - Removed ", len(edges_to_remove), " edges.", sep="", file=sys.stderr)
    pruneOrphanNodes(G)


###Command line operation
if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog = "genotypeClones.py", description = "Graph approach to generate unified list of which BCs are present in which cells", allow_abbrev=False)
    parser.add_argument('--inputfilename', action='store', required=True, help='input filename. Format: tab-delimited with barcode sequences. First line must be header (unused). Columns must be BC, cellBC, count.')
    parser.add_argument('--outputlong', action='store', required=True, help='output filename for tab-delimited list of BC counts totalled per clone')
    parser.add_argument('--outputwide', action='store', required=True, help='output filename for tab-delimited list of clones and the cells/BCs they include')
    parser.add_argument('--output', action='store', required=True, help='output filename for tab-delimited list of clone, cell, BC links - filtered version of barcode.counts.byCell file. Can be - for stdout.')
    parser.add_argument('--cloneobj', action='store', default=None, help='output .pickle filename to serialize native clones and G objects')
    parser.add_argument('--whitelist', action='store', default=None, help='Comma separated list filenames, each containing BCs to be retained. Format: one per line, no other columns.')
    parser.add_argument('--blacklist', action='store', default=None, help='Comma separated list filenames, each containing cellBCs or BCs to be dropped. Format: one per line, no other columns')
    parser.add_argument('--annotateclones', action='store', default=None, help='Comma separated list filenames, each containing annotation to be added to cells or cellBCs. Format: gzipped tab-delimited file with header, one header column must be "bc", another must be "transfection"')
    
    parser.add_argument('--minreads', action='store', type=int, default=2, help='Min UMI filter for input file')
    parser.add_argument('--minPropOfBCReads', action='store', type=float, default=0.15, help='Each BC-cell edge must represent at least this proportion of UMIs for BC')
    parser.add_argument('--minPropOfCellReads', action='store', type=float, default=0.02, help='Each BC-cell edge must represent at least this proportion of UMIs for cell')
    parser.add_argument('--minCentrality', action='store', type=float, default=0.2, help='Edges joining communities must have higher centrality than this to be retained')
    parser.add_argument('--maxpropreads', action='store', type=float, default=0.1, help='Edges joining communities must have fewer than this number of UMIs as proportion of the smaller community they bridge')
    
    key_arg = parser.add_argument("--transfectionKey", action="store", type=str, default=None, help="Attribute key used to identify BCs transfections")
    parser.add_argument("--removeMinorityBCsFromConflictingCells", action="store", type=float, default=None, help="For cells with BCs from more than one transfection, removes BCs from minority transfections that together represent at most this proportion of UMIs for that cell. Requires `transfectionKey`")
    parser.add_argument("--removeConflictingCells", action='store_true', default=False, help="Removes cells linked to BCs from 2+ transfections. It requires `transfectionKey`")
    parser.add_argument("--removeConflictingBCs", action='store_true', default=False, help="Removes BCs linked to cells from 2+ transfections. It requires `transfectionKey`")
    
    parser.add_argument("--minjaccard", action="store", type=float, default=0, help="Remove edges connecting cells whose BC overlap is below this Jaccard-index threshold")
    
    parser.add_argument('--printGraph', action='store', type=str, help='Plot a graph for each clone into this directory')
    
    parser.add_argument("--verbose", action='store_true', default=False, help = "Verbose mode")
    parser.add_argument('--version', action='version', version='%(prog)s ' + version)
    
    try:
        args = parser.parse_args()
        if args.removeMinorityBCsFromConflictingCells is not None and args.transfectionKey is None:
            raise ArgumentError(key_arg, "`--transfectionKey` is required to use `--removeMinorityBCsFromConflictingCells`")
        if args.transfectionKey is not None and args.removeConflictingCells and args.transfectionKey is None:
            raise ArgumentError(key_arg, "`--transfectionKey` is required to use `--removeConflictingCells`")
    except argparse.ArgumentError as exc:
        print("[genotypeClones] ", exc.message, '\n', exc.argument, sep="", file=sys.stderr)
        sys.exit(1)
    
    print("[genotypeClones] " + str(args), file=sys.stderr)
    
    if args.printGraph is not None:
        os.makedirs(args.printGraph, exist_ok=True)
    
    
    ###Initialize graph, filter, write output
    G = initializeGraphFromInput(inputfilename=args.inputfilename, minCount=args.minreads)
    
    #Add annotation if any
    #TODO generalize
    if args.annotateclones is not None:
        for f in args.annotateclones.split(","):
            print("[genotypeClones] Loading annotation from ", f, sep="", file=sys.stderr)
            inputfile_reader = csv.DictReader(gzip.open(f, 'rt'), delimiter='\t')
            bcToTransfection = {}
            for line in inputfile_reader:
                bcToTransfection[line['bc']] = line['transfection']
            
            assignToNodes(G, args.transfectionKey, bcToTransfection, default="None")
            ## Loops over all cells and assign their `transfectionKey` field to "cellBC" that facilitate coloring them in printGraph
            for n in G.nodes:
                if G.nodes[n]['type'] == 'cell':
                    G.nodes[n][args.transfectionKey] = 'cellBC'
        ## TODO hard-coded colorscheme for T0215toT0218 and T0219toT0222
        printGraph_kwds = { 'node_color': args.transfectionKey, 'node_color_dict': {'cellBC': 'black', 'conflicting': 'yellow', 'T0215A': 'orange', 'T0216B': 'purple', 'T0217B': 'green', 'T0219A': 'orange', 'T0220B': 'purple', 'T0221B': 'green', 'T0222B': 'red'} }
    else:
        printGraph_kwds = {}
    
    ## Unlink cells with low similarity
    if (args.minjaccard > 0):
        pruneLowJaccard(G, args.minjaccard)
    
    if args.removeMinorityBCsFromConflictingCells is not None:
        pruneConflictingEdges(G, args.transfectionKey, args.removeMinorityBCsFromConflictingCells)
    if args.transfectionKey is not None and args.removeConflictingCells:
        pruneConflictingNodes(G, args.transfectionKey, type = "cell")
    if args.transfectionKey is not None and args.removeConflictingBCs:
        pruneConflictingNodes(G, args.transfectionKey, type = "BC")
    
    pruneEdgesLowPropOfReads(G, minPropOfBCReads=args.minPropOfBCReads, minPropOfCellReads=args.minPropOfCellReads)
    
    if args.whitelist is not None:
        for f in args.whitelist.split(","):
            filterNodesFromFile(G, filename=f, keep=True)
    if args.blacklist is not None:
        for f in args.blacklist.split(","):
            filterNodesFromFile(G, filename=f, keep=False)
    if args.whitelist is not None or args.blacklist is not None:
        pruneOrphanNodes(G)
    
    breakUpWeaklyConnectedCommunities(G, minCentrality=args.minCentrality, maxPropReads=args.maxpropreads, verbose=args.verbose, graphOutput=None)
    #run summarizeGraph since we don't need pruneOrphanNodes
    summarizeGraph(G)
    #Do twice to break up some of the bigger graphs since we don't iterate internally, 3x doesn't do anything else
    breakUpWeaklyConnectedCommunities(G, minCentrality=args.minCentrality, maxPropReads=args.maxpropreads, verbose=args.verbose, graphOutput=args.printGraph)
    #run summarizeGraph since we don't need pruneOrphanNodes
    summarizeGraph(G)
    
    clones = identifyClones(G)
    writeOutputFiles(G, clones, args.output, args.outputlong, args.outputwide, args.cloneobj, args.printGraph, transfectionKey=args.transfectionKey)
