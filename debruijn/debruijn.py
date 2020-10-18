#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

from random import randint
import statistics
import argparse
import os
import sys
from operator import itemgetter
import random
import networkx as nx
import matplotlib.pyplot as plt
random.seed(9001)

__author__ = "Estelle Mariaux"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Estelle Mariaux"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Estelle Mariaux"
__email__ = "estelle.mariaux@hotmail.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()

# Identification des kmer unique ------------------------------

def read_fastq(fastq_file):
    """ Creates a generator of sequences
    """
    with open(fastq_file, "r") as file:
        for _ in file:
            yield next(file).strip()
            next(file)
            next(file)

def cut_kmer(sequence, kmer_size):
    """Creates a generator of kmer)
    """
    k = 0
    while kmer_size+k < len(sequence)+1:
        yield sequence[k:kmer_size+k]
        k+=1

def build_kmer_dict(fastq_file, kmer_size):
    """Creates dictionnary of kmer
    """
    kmer_dict = {}
    for seq in read_fastq(fastq_file):
        for kmer in cut_kmer(seq, kmer_size):
            if kmer in kmer_dict.keys():
                kmer_dict[kmer] += 1
            else:
                kmer_dict[kmer] = 1
    return kmer_dict

# Construction de l'arbre de De Bruijn ------------------------

def build_graph(dict_kmer):
    """Creates graph from dictionnary
    """
    graph = nx.DiGraph()
    for kmer, poids in dict_kmer.items():
        graph.add_edge(kmer[:-1],kmer[1:], weight = poids)
    #affichage du graphe, seulement sur petite donnée
    #plt.subplot(111)
    #nx.draw(graph, with_labels = True)
    #plt.savefig("graphe")
    return graph

# Parcours du graphe de De Bruijn -----------------------------

def get_starting_nodes(graph):
    """Creates list of nodes
    """
    noeud_entree = []
    for node in graph:
        if len(list(graph.predecessors(node))) == 0:
            noeud_entree.append(node)
    return noeud_entree

def get_sink_nodes(graph):
    """Creates list of nodes
    """
    noeud_sortie = []
    for node in graph:
        if len(list(graph.successors(node))) == 0:
            noeud_sortie.append(node)
    return noeud_sortie

def get_contigs(graph,list_starting, list_sink):
    """Creates list of tuple(contig, length contig)
    """
    contigs = []
    for start in list_starting:
        for sink in list_sink:
            path = list(nx.all_simple_paths(graph, source = start, target = sink))
            if path:
                contig = path[0][0]
                for i in range(1, len(path[0])):
                    contig += path[0][i][-1]
                contigs.append((contig, len(contig)))
    return contigs

def fill(text, width = 80):
    """Split text with a line return to respect fasta format
    """
    return os.linesep.join(text[i:i+width] for i in range(0,len(text), width))

def save_contigs(list_tuple, name_file):
    """Creates a txt file for the contigs
    """
    with open(name_file,"w") as file:
        for i in range(0,len(list_tuple)):
            file.write(">contig_{} len={}\n".format(i, list_tuple[i][1]))
            file.write("{}\n".format(fill(list_tuple[i][0], 80)))

# Resolution des bulles ---------------------------------------

def std(liste):
    """Computes the standard deviation of a list
    """
    return statistics.stdev(liste)

def path_average_weight(graph, path):
    """Computes the average weight of one path
    """
    weight_full = 0
    for i in range(len(path)-1):
        weight_full += graph[path[i]][path[i+1]]["weight"]
    return weight_full / (len(path)-1)

def remove_paths(graph, list_path, delete_entry_node, delete_sink_node):
    """Creates a cleaned graph
    """
    for nodes in list_path:
        if delete_entry_node:
            graph.remove_node(nodes[0])
        if delete_sink_node:
            graph.remove_node(nodes[-1])
        for node in nodes[1:-1]:
            graph.remove_node(node)
    return graph

def select_best_path(graph, list_path, list_length, list_average_weight,
                    delete_entry_node = False, delete_sink_node = False):
    """Creates the cleaned graph of unwanted paths
    """
    #test weight
    if std(list_average_weight) > 0:
        index_worst_path = list_average_weight.index(min(list_average_weight))
        path = tuple(list_path[index_worst_path])
        graph = remove_paths(graph, [path], delete_entry_node, delete_sink_node)
    #test length
    elif std(list_length) > 0:
        index_worst_path = list_length.index(min(list_length))
        path = tuple(list_path[index_worst_path])
        graph = remove_paths(graph, [path], delete_entry_node, delete_sink_node)
    #test random
    else :
        index_worst_path = randint(0,1)
        path = tuple(list_path[index_worst_path])
        graph = remove_paths(graph, [path], delete_entry_node, delete_sink_node)
    return graph

def solve_bubble(graph, node_pred, node_succ):
    """Remove multiple paths from graph between two nodes
    """
    bubble_paths = list(nx.all_simple_paths(graph, node_pred, node_succ))
    while len(bubble_paths) >= 2 :
        list_paths = [bubble_paths[0], bubble_paths[1]]
        list_length = [len(bubble_paths[0]), len(bubble_paths[1])]
        list_weight = [path_average_weight(graph, bubble_paths[0]),
                        path_average_weight(graph, bubble_paths[1])]
        graph = select_best_path(graph, list_paths, list_length, list_weight)
        bubble_paths = list(nx.all_simple_paths(graph, node_pred, node_succ))
    return graph


def simplify_bubbles(graph):
    """Gets rid of unnecessary bubbles
    """
    start_nodes = get_starting_nodes(graph)
    sink_nodes = get_sink_nodes(graph)
    for start_node in start_nodes:
        for sink_node in sink_nodes:
            paths = list(nx.all_simple_paths(graph, start_node, sink_node))
            while len(paths) > 1 :
                flag_bubble = 0
                for k, val in enumerate(paths[0]):
                    if flag_bubble == 0:
                        if val not in paths[1]:
                            node_1 = paths[0][k-1]
                            flag_bubble = 1
                    else:
                        if val in paths[1]:
                            node_2 = paths[1][k]
                            break
                graph = solve_bubble(graph, node_1, node_2)
                paths = list(nx.all_simple_paths(graph, start_node, sink_node))
    return graph

# Detection des pointes ---------------------------------------

def solve_entry_tips(graph, list_start_node):
    """Cleares the graph of unwanted starting path
    """
    last_common_node = list_start_node[0]
    node = list(graph.successors(list_start_node[0]))[0]
    while len(list(graph.successors(node))) != 0:
        if len(list(graph.predecessors(node))) > 1:
            last_common_node = node
        node = list(graph.successors(node))[0]
    if last_common_node == list_start_node[0]:
        return graph
    list_paths, list_length, list_weight = [], [], []
    for start_node in list_start_node:
        for path in nx.all_simple_paths(graph, start_node, last_common_node):
            list_paths.append(path)
            list_length.append(len(list(path)))
            list_weight.append(path_average_weight(graph, path))
    graph = select_best_path(graph, list_paths, list_length, list_weight,
            delete_entry_node = True, delete_sink_node = False)
    return graph

def solve_out_tips(graph, list_sink_node):
    """Cleares the graph of unwanted sinking path
    """
    first_common_node = list_sink_node[0]
    node = list(graph.predecessors(list_sink_node[0]))[0]
    while len(list(graph.predecessors(node))) != 0:
        if len(list(graph.successors(node))) > 1:
            first_common_node = node
        node = list(graph.predecessors(node))[0]
    if first_common_node == list_sink_node[0]:
        return graph
    list_paths, list_length, list_weight = [], [], []
    for sink_node in list_sink_node:
        for path in nx.all_simple_paths(graph, first_common_node, sink_node):
            list_paths.append(path)
            list_length.append(len(list(path)))
            list_weight.append(path_average_weight(graph, path))
    graph = select_best_path(graph, list_paths, list_length, list_weight,
            delete_entry_node = False, delete_sink_node = True)
    return graph

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    # Lecture du fichier et de construction du graphe
    dico_kmer = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(dico_kmer)

    # Resolution des bulles
    graph = simplify_bubbles(graph)

    # Resolution des pointes d'entree et de sortie
    entry_nodes = get_starting_nodes(graph)
    graph = solve_entry_tips(graph, entry_nodes)

    out_nodes = get_sink_nodes(graph)
    graph = solve_out_tips(graph, out_nodes)

    # Ecriture du/des contigs
    entry_nodes = get_starting_nodes(graph)
    out_nodes = get_sink_nodes(graph)
    contigs = get_contigs(graph, entry_nodes, out_nodes)
    save_contigs(contigs, args.output_file)

if __name__ == '__main__':
    main()
