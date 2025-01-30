#!/usr/bin/env python3
import os
import sys
import argparse
from typing import Tuple, List, Dict


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


mammal_biotypes = [
    "protein_coding", 
    "protein_coding_LoF", 
    "lncRNA", 
    "IG_C_gene", 
    "IG_D_gene", 
    "IG_J_gene", 
    "IG_LV_gene", 
    "IG_V_gene", 
    "IG_V_pseudogene", 
    "IG_J_pseudogene", 
    "IG_C_pseudogene", 
    "TR_C_gene", 
    "TR_D_gene", 
    "TR_J_gene", 
    "TR_V_gene", 
    "TR_V_pseudogene", 
    "TR_J_pseudogene"
]

bird_biotypes = [
    "protein_coding",
    "protein_coding_LoF",
    "lncRNA",
    "IG_V_gene",
    "IG_J_gene",
    "IG_V_pseudogene",
    "IG_J_pseudogene",
    "IG_C_gene",
    "IG_C_pseudogene",
    "TR_V_gene",
    "TR_J_gene",
    "TR_V_pseudogene",
    "TR_J_pseudogene",
    "TR_C_gene",
    "TR_D_gene",
    "TR_C_pseudogene"
]

amphibian_biotypes = [
   "protein_coding",
   "lncRNA",
   "IG_C_gene",
   "IG_D_gene",
   "IG_J_gene", 
   "IG_V_gene",
   "IG_V_pseudogene",
   "IG_J_pseudogene",
   "IG_C_pseudogene",
   "TR_C_gene",
   "TR_D_gene",
   "TR_J_gene",
   "TR_V_gene", 
   "TR_V_pseudogene",
   "TR_J_pseudogene"
]

fish_biotypes = [
   "protein_coding",
   "lncRNA",
   "IG_C_gene",
   "IG_D_gene", 
   "IG_J_gene",
   "IG_V_gene",
   "IG_V_pseudogene",
   "IG_J_pseudogene",
   "IG_C_pseudogene",
   "TR_C_gene",
   "TR_D_gene",
   "TR_J_gene", 
   "TR_V_gene",
   "TR_V_pseudogene",
   "TR_J_pseudogene",
   "IG_gene",
   "TR_gene"
]

invertibrate_biotypes = [
    "protein_coding",
    "lncRNA"
]

plant_biotypes = [
    "protein_coding", 
    "lncRNA",
    "lincRNA", 
    "pseudogene"
]

fungi_biotypes = [
    "protein_coding", "ncRNA", "pseudogene"
]

biotype_index = {
    # animals
    ## mammals
    "Rattus norvegicus" : mammal_biotypes,
    "Macaca mulatta" : mammal_biotypes,
    "Callithrix jacchus" : mammal_biotypes,
    "Troglodytes gorilla" : mammal_biotypes,
    "Equus caballus" : mammal_biotypes,
    "Canis lupus familiaris" : mammal_biotypes,
    "Bos taurus" : mammal_biotypes,
    "Ovis aries" : mammal_biotypes,
    "Sus scrofa" : mammal_biotypes,
    "Heterocephalus glaber" : mammal_biotypes,
    "Oryctolagus cuniculus" : mammal_biotypes,
    ## birds
    "Gallus gallus" : bird_biotypes,
    ## amphibians
    "Xenopus tropicalis" : amphibian_biotypes,
    ## fish
    "Danio rerio" : fish_biotypes,
    ## invertibrates
    "Drosophila melanogaster" : invertibrate_biotypes,
    "Caenorhabditis elegans" : invertibrate_biotypes,
    "Schistosoma mansoni" : invertibrate_biotypes,
    "Anopheles gambiae" : invertibrate_biotypes,
    # plants
    "Arabidopsis thaliana" : plant_biotypes,
    "Oryza sativa" : plant_biotypes,
    "Solanum lycopersicum" : plant_biotypes,
    "Zea mays" : plant_biotypes,
    # fungi
    "Saccharomyces cerevisiae" : fungi_biotypes
}



def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.
    Returns:
        argparse.Namespace containing arguments.
    """
    desc = 'Create STAR reference genome index.'
    epi = """DESCRIPTION:

    """
    parser = argparse.ArgumentParser(description=desc, epilog=epi, formatter_class=CustomFormatter)
    parser.add_argument(
        'fasta_url', type=str, help='URL to genome FASTA file'
    )
    parser.add_argument(
        'gtf_url', type=str, help='URL to genome GTF file'
    )
    parser.add_argument(
        'organism', type=str, help='Organism name'
        choices=['Homo sapiens', 'Mus musculus'
    )
    return parser.parse_args()