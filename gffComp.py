#!/usr/bin/python3

import getopt
import os
import sys
import gzip
import shutil

# Import Biopython
try:
    from Bio import SeqIO
except ModuleNotFoundError:
    sys.exit("\nBiopython not installed\n\"pip3 install biopython\" to install it in your local environment.")

# Import GFF parser
try:
    from BCBio import GFF
except ModuleNotFoundError:
    sys.exit("\nbcbio-gff not installed\n\"pip3 install bcbio-gff\" to install it in your local environment.")


""" Extract and compare annotations products from 2 bacterial assembly GFF files. Print results in STDOUT and 
produce a binary table. """


def help_msg():
    """ Display help message """
    print('./gffComp.py [option]\n'
          'Compare GFF files products of two assemblies.\n'
          'Use either -g or -i\n'
          '-h\t\tDisplay this help\n'
          '-g, --gff\tGFF files (separated by comma, no space)\n'
          '-i, --ids\tAssembly IDs (separated by comma, no space)')


def gunzip_shutil(source_filepath, dest_filepath, block_size=65536):
    """ Unzipping function if gff downloaded via dl_gff() """
    with gzip.open(source_filepath, 'rb') as s_file, \
            open(dest_filepath, 'wb') as d_file:
        shutil.copyfileobj(s_file, d_file, block_size)


def check_arg(argv):
    """ Verify arguments prior to using methods """
    try:
        opts, args = getopt.getopt(argv, "hg:i:", ["gff=", "ids="])
    except getopt.GetoptError:
        help_msg()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            help_msg()
            sys.exit(2)
        elif opt in ("-g", "--gff"):
            if len(arg.split(",")) == 2:
                return "gff", arg
            else:
                help_msg()
                sys.exit(2)
        elif opt in ("-i", "--ids"):
            if len(arg.split(",")) == 2:
                return "ids", arg
            else:
                help_msg()
                sys.exit(2)


def dl_gff(ids):
    """Download of gff files"""

    # Import ncbi downloader
    try:
        import ncbi_genome_download as ngd
    except ModuleNotFoundError:
        sys.exit("\nncbi_genome_download not installed.\n\"pip install ncbi_genome_download\" "
                 "to install it in your local environment.")

    print("Downloading files ...")
    ngd.download(section='genbank', file_formats='all', assembly_accessions=ids, output='dl')
    print("Done!")


def prep_gff(ids):
    """ Extract zipped gff files and move them into another folder """
    try:
        os.mkdir("genDist_gff")
    except FileExistsError:
        shutil.rmtree("genDist_gff")
        os.mkdir("genDist_gff")

    ids = ids.split(",")  # Split the 2 ids in a list
    for id in ids:
        try:
            gz_file = os.listdir(f"dl/genbank/bacteria/{id}")
        except FileNotFoundError:  # If download failed
            sys.exit(2)
        gz_file.remove("MD5SUMS")
        gunzip_shutil(f"dl/genbank/bacteria/{id}/{gz_file[0]}", f"genDist_gff/{id}.gff")

    shutil.rmtree("dl")  # Delete download folder


def parsing_gff(mode, names):
    """ Parsing of gff files """

    # Variables initialization
    gff1 = None  # Variable for 1st gff file
    gff2 = None  # Variable for 2nd gff file
    genes1 = []  # Genes present in genome 1
    genes2 = []  # Genes present in genome 2
    names = names.split(",")  # split files names

    # Open files depending on mode selected
    if mode == "gff":
        gff1 = GFF.parse(open(names[0]))
        gff2 = GFF.parse(open(names[1]))
    elif mode == "ids":
        gff1 = GFF.parse(open(f"genDist_gff/{names[0]}.gff"))
        gff2 = GFF.parse(open(f"genDist_gff/{names[1]}.gff"))

    # Parsing files and extracting genes product names
    for sequence in gff1:
        for feature in sequence.features:
            if feature.type == "gene" or feature.type == "pseudogene":
                if feature.sub_features[0].qualifiers["gbkey"][0] == "CDS":
                    product = feature.sub_features[0].qualifiers["product"][0]
                    if "uncharacterised" not in product.lower() \
                            and "putative" not in product.lower() \
                            and "hypothetical" not in product.lower():
                        genes1.append(product)

    for sequence in gff2:
        for feature in sequence.features:
            if feature.type == "gene" or feature.type == "pseudogene":
                if feature.sub_features[0].qualifiers["gbkey"][0] == "CDS":
                    product = feature.sub_features[0].qualifiers["product"][0]
                    if "uncharacterised" not in product.lower() \
                            and "putative" not in product.lower() \
                            and "hypothetical" not in product.lower():
                        genes2.append(product)

    # Set common and unique genes
    common = [gene for gene in genes1 if gene in genes2]
    uniques1 = [gene for gene in genes1 if gene not in genes2]
    uniques2 = [gene for gene in genes2 if gene not in genes1]
    all_genes = genes1 + uniques2

    # Print results
    print(f"\nResults\n"
          f"Common genes: {len(common)}\n"
          f"Genes unique to {names[0]}: {len(uniques1)}\n"
          f"Genes unique to {names[1]}: {len(uniques2)}\n")

    # Binary matrix
    with open("binary_matrix.tsv", "w") as output:
        header = "\t".join(all_genes)
        output.write(f"{header}\n")

        # Reading genes for genome 1
        binary = []
        for gene in all_genes:
            if gene in uniques1 + common:
                binary.append("1")
            else:
                binary.append("0")
        binary_text = "\t".join(binary)
        output.write(f"{names[0]}\t{binary_text}\n")

        # Reading genes for genome 2
        binary = []
        for gene in all_genes:
            if gene in uniques2 + common:
                binary.append("1")
            else:
                binary.append("0")
        binary_text = "\t".join(binary)
        output.write(f"{names[1]}\t{binary_text}\n")


if __name__ == "__main__":
    try:
        mode, values = check_arg(sys.argv[1:])
    except TypeError:
        help_msg()
        sys.exit(2)

    if mode == "ids":
        dl_gff(values)
        prep_gff(values)
        parsing_gff(mode, values)
    elif mode == "gff":
        parsing_gff(mode, values)
