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

import argparse
import os
import sys
import networkx as nx
import matplotlib.pyplot as plt
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics

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

#1 Creation du graphe de De Bruijn ----------------------------

def read_fastq(fastq_file):
    """ Creates a generator of sequences
    """
    with open(fastq_file, "r") as file:
        for seq in file:
            yield(next(file))
            next(file)
            next(file)

def cut_kmer(sequence, kmer_size):
    """Creates a generator of kmer)
    """
    k = 0
    while kmer_size+k < len(sequence):
        yield(sequence[k:kmer_size+k])
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
    return(kmer_dict)

def build_graph(dict_kmer):
    """Creates graph from dictionnary
    """
    G = nx.DiGraph()
    for kmer, poids in dict_kmer.items():
        G.add_edge(kmer[:-1],kmer[1:], weight = poids)
    """ #affichage du graphe, seulement sur petite donnée
    plt.subplot(111)
    nx.draw(G, with_labels = True)
    plt.savefig("test_graphe")
    """
    return(G)

#2 Parcours du graphe de De Bruijn ----------------------------
def get_starting_nodes(graph):
    """Creates list of nodes
    """
    noeud_entree = []
    for node in graph:
        if len(list(graph.predecessors(node))) == 0:
            noeud_entree.append(node)
    return(noeud_entree)
    
def get_sink_nodes(graph):
    """Creates list of nodes
    """
    noeud_sortie = []
    for node in graph:
        if len(list(graph.successors(node))) == 0:
            noeud_sortie.append(node)
    return(noeud_sortie)

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
    return(contigs)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    
    # Identification des k-mer unique
    dico_kmer = build_kmer_dict(args.fastq_file, args.kmer_size)
    #print(dico_kmer)
    
    # Construction de l'arbre de De Bruijn
    G = build_graph(dico_kmer)
    
    # Graphe de De Bruijn
    liste_entree = get_starting_nodes(G)
    liste_sortie = get_sink_nodes(G)
    #print(liste_entree)
    #print(liste_sortie)
    print(get_contigs(G, liste_entree, liste_sortie))
    
if __name__ == '__main__':
    main()
