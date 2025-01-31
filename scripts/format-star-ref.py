#!/usr/bin/env python3
import os
import sys
import gzip
import argparse
from typing import Tuple, List, Dict, Set
from collections import Counter


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass

# global vars
mammal_biotypes = {
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
}

bird_biotypes = {
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
}

amphibian_biotypes = {
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
}

fish_biotypes = {
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
}

invertibrate_biotypes = {
    "protein_coding",
    "lncRNA"
}

plant_biotypes = {
    "protein_coding", 
    "lncRNA",
    "lincRNA"
}

fungi_biotypes = {
    "protein_coding", 
    "ncRNA"
}

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

# functions
def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.
    Returns:
        argparse.Namespace containing arguments.
    """
    desc = 'Create STAR reference genome index.'
    epi = """DESCRIPTION:
    # example
    ./scripts/format-star-ref.py \
      --organism "Macaca mulatta" \
      --fasta /home/nickyoungblut/tmp/genomes/reference_sources/Macaca_mulatta.Mmul_10.dna.toplevel.fa \
      /home/nickyoungblut/tmp/genomes/reference_sources/Macaca_mulatta.Mmul_10.113.gtf
    """
    parser = argparse.ArgumentParser(description=desc, epilog=epi, formatter_class=CustomFormatter)
    parser.add_argument(
        'gtf', type=str, help='Path to genome GTF file'
    )
    parser.add_argument(
        '--fasta', type=str, default=None,
        help='Path to genome FASTA file'
    )
    parser.add_argument(
        '--output-dir', type=str, default='star_ref', help='Output base directory', 
    )
    parser.add_argument(
        '--organism', type=str, choices=biotype_index.keys(), required=True,
        help='Organism name',
    )
    parser.add_argument(
        '--exclude-tags', type=str, nargs='+',
        default=["readthrough_transcript", "PAR"],
        help='Filter records containing this tag',
    )
    return parser.parse_args()

def process_gtf_line(
        line: str, 
        outF: 'TextIO', 
        biotypes: Set[str], 
        exclude_tags: List[str], 
        seq_names: Set[str], 
        status: Dict[str, int]
    ):
    """
    Process a single gtf line.
    Args:
        line: GTF line.
        outF: Output file handle.
        biotypes: Set of biotypes to keep.
        exclude_tags: List of tags to exclude.
        seq_names: Set of sequence names to keep.
        status: Dictionary to store status.
    """
    # simply write out header line
    if line.startswith("#"):
        outF.write(line)
        return None

    # status
    status["total_raw"] += 1

    # parse body line
    fields = line.strip().split("\t")
    attributes = {}
    for x in fields[8].split(";"):
        x = x.strip()
        if x:
            try:
                key, value = x.split(" ", 1)
                attributes[key.strip()] = value.strip('"')
            except ValueError:
                continue
    fields = fields[:8]

    # store seq names
    seq_names.add(fields[0])

    # convert gene_id
    if attributes.get("gene_id") and not attributes.get("gene_version"):
        try:
            gene_id,gene_version = str(attributes["gene_id"]).split(".")
            attributes["gene_id"] = gene_id
            attributes["gene_version"] = gene_version
        except ValueError:
            pass
    
    # filter by biotype
    biotype_fields = ['gene_biotype', 'gene_type', 'transcript_type', 'transcript_biotype']
    for field in biotype_fields:
        if attributes.get(field) and attributes[field] not in biotypes:
            status["biotype"] += 1
            try:
                status[attributes[field]] += 1
            except KeyError:
                status[attributes[field]] = 1
            return None
    
    # filter by tags
    if attributes.get("tag") and attributes.get("tag") in exclude_tags:
        status["tag"] += 1
        return None

    # write out
    attributes = "; ".join([f"{k} \"{v}\"" for k, v in attributes.items()])
    outF.write("\t".join(fields + [attributes]) + "\n")

def process_fasta(fasta: str, output_fasta: str, seq_names: Set[str]):
    """
    Process a fasta file. Check to make sure that the sequence names are in the set.
    Args:
        fasta: Path to input fasta file.
        output_fasta: Path to output fasta file.
        seq_names: Set of sequence names to keep.
    """
    print(f"Processing fasta: {os.path.basename(fasta)}", file=sys.stderr)

    if fasta.endswith(".gz"):
        _open = gzip.open
    else:
        _open = open

    with _open(fasta) as inF, gzip.open(output_fasta, 'wb') as outF:
        write = False
        for line in inF:
            try:
                line = line.decode()   
            except AttributeError:
                pass
            if line.startswith(">"):
                if line.lstrip(">").split(" ")[0] in seq_names:
                    write = True
                else:
                    write = False
                    print(f"Sequence not in GTF: {line.strip()}", file=sys.stderr)
            if write:
                outF.write(line.encode())

def main():
    # parse cli arguments
    args = parse_args()

    # output
    args.output_dir = os.path.join(args.output_dir, args.organism.replace(" ", "_"))
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir, exist_ok=True)

    # check input
    if args.fasta and not os.path.exists(args.fasta):
        sys.exit(f"Error: {args.fasta} not found")
    if not os.path.exists(args.gtf):
        sys.exit(f"Error: {args.gtf} not found")

    # iterate over gtf
    print(f"Processing GTF: {os.path.basename(args.gtf)}", file=sys.stderr)
    seq_names = set()
    status = {"total_raw": 0, "biotype": 0, "tag": 0}
    output_gtf = os.path.join(args.output_dir, f"{args.organism}.gtf")
    with open(args.gtf) as inF, open(output_gtf, 'w') as outF:
        for line in inF:
            process_gtf_line(
                line, outF, 
                biotypes=biotype_index[args.organism],
                exclude_tags=args.exclude_tags,
                seq_names=seq_names,
                status=status,
            )
        
    ## GTF processing status
    print(f"Total records in GTF: {status['total_raw']}", file=sys.stderr)
    for key in ["biotype", "tag"]:
        print(f"Filtered {status[key]} records by {key}", file=sys.stderr)
    print("-- Count of biotypes filtered --", file=sys.stderr)
    for key in status:
        if key not in ["total_raw", "biotype", "tag"]:
            print(f"{key}: {status[key]}", file=sys.stderr)

    # fasta
    if args.fasta:
        output_fasta = os.path.join(args.output_dir, f"{args.organism}.fna.gz")
        process_fasta(args.fasta, output_fasta, seq_names)


# main
if __name__ == "__main__":
    main()