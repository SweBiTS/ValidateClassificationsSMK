#!/usr/bin/env python3

import argparse
import logging
import gzip
from stringmeup.taxonomy import TaxonomyTree
from os import path

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d [%H:%M:%S]')
log = logging.getLogger(path.basename(__file__))

# Arguments
parser = argparse.ArgumentParser(
    prog='FASTQ taxonomic subsetter (FastqTaxSub)',
    usage='subset_reads --mode {"single", "clade"} --output FILE --nodes FILE --names FILE --input FILE [--tax_id INT | --tax_id_file FILE]',
    description='Extract read IDs from Kraken2 output file that are classified to a user specified taxonomic ID (or its clade).')

# Input and output files
parser.add_argument(
    '--input',
    metavar='Kraken2 output',
    type=str,
    required=True,
    help='The kraken2 read-by-read output file.')
parser.add_argument(
    '--output',
    metavar='FILE',
    type=str,
    required=True,
    help='The output file name.')

# Supply relevant taxonomic ID on command line, or one or multiple taxonomic IDs
# through a text file.
tax_id_group = parser.add_mutually_exclusive_group(required=True)
tax_id_group.add_argument(
    '--tax_id',
    metavar='Taxonomic ID',
    type=int,
    help='A taxonomic ID for a single taxon or clade that you wish to get read IDs for.')
tax_id_group.add_argument(
    '--tax_id_file',
    metavar='FILE',
    type=str,
    help='Supply multiple taxonomic IDs at once. A textfile with one taxonomic ID per line.')

# Work only with the supplied taxonomic ID ('single'), or all taxonomic IDs
# present in the clade rooted at the supplied taxonomic ID.
parser.add_argument(
    '--mode',
    metavar='{single, clade}',
    type=str,
    required=True,
    choices=['single', 'clade'],
    help='Put "single" to extract reads hitting exactly the supplied taxonomic ID. Put "clade" to extract reads hitting any taxonomic ID in the clade rooted at the supplied taxonomic ID.')

# The taxonomy
parser.add_argument(
    '--names',
    metavar='FILE',
    required=True,
    help='taxonomy names dump file (names.dmp)')
parser.add_argument(
    '--nodes',
    metavar='FILE',
    required=True,
    help='taxonomy nodes dump file (nodes.dmp)')

# Parse the arguments
args = parser.parse_args()


def read_file(filename):
    """
    Wrapper to read either gzipped or ordinary text file input.
    """
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'r')


def get_tax_read_ids(line):
    """
    Takes as input a row from a kraken2 classifications output file.
    Returns the read ID and taxonomic ID of the input line.
    """
    if line.startswith('U'):
        return (0, 0)

    split_line = line.strip().split('\t')
    read_id = split_line[1]
    tax_id = int(split_line[2])

    return (read_id, tax_id)


def get_tax_ids(args):
    """
    Returns a set containing ints of taxonomic IDs. Reads either from a file
    with potentially multiple taxonomic IDs, or a single taxonomic ID from
    the command line.
    """
    tax_id_set = set()

    # Loop over each line in the supplied file and add the taxonomic IDs to
    # the list
    if args.tax_id_file:
        with read_file(args.tax_id_file) as f:
            for l in f:
                l = int(l.strip())
                tax_id_set.add(l)

    # If no file, taxonomic ID must come from the command line
    else:
        tax_id_set.add(args.tax_id)

    return tax_id_set


if __name__ == "__main__":

    # Initiate the taxonomy tree
    taxonomy_tree = TaxonomyTree(names_filename=args.names, nodes_filename=args.nodes)

    # Get which tax_ids that we are interested in
    supplied_tax_ids = get_tax_ids(args)

    # Get the sought-after tax_id(s)
    child2root_map = {}  # A dictionary that will hold <root_tax_id>: <child_tax_id> key/value pairs.
    if args.mode == 'single':
        log.info('Mode is "single"; will look for reads classified only to the supplied taxonomic ID(s).')

        # The translation between actual classified tax_id, and the root tax_id (the supplied tax_id)
        # will in the 'single' mode just be tax_id[N]: tax_id[N].
        child2root_map = {tax_id: tax_id for tax_id in supplied_tax_ids}
    else:
        # Find out which tax_ids are in the clade rooted at the supplied tax_id
        log.info('Mode is "clade"; will look for reads classified to any taxon in the clade(s) rooted at the supplied taxonomic ID(s).')

        for tax_id in supplied_tax_ids:
            clade_tax_ids = taxonomy_tree.get_clade([tax_id])[tax_id]

            # Add each tax_id in the clade to the child2root_map dictionary
            for taxon in clade_tax_ids:
                child2root_map[int(taxon)] = tax_id

    # We are interested in all tax_ids that are keys in the child2root_map dict
    tax_ids = set(child2root_map.keys())

    # Log some stuff
    log.info('Read IDs for the following taxonomic IDs will be output:')
    log.info('tax_id\t\tscientific_name')
    for tax_id in tax_ids:
        scientific_name = taxonomy_tree.get_name([tax_id])[tax_id]
        log.info('{}\t\t{}'.format(tax_id, scientific_name))
    log.info('Saving read IDs in file {}'.format(args.output))
    log.info('Started parsing file {} ...'.format(args.input))

    # Open the input file and parse it
    with read_file(args.input) as f, open(args.output, 'w') as o:
        i = 0
        c = 0

        # Parse each row in the input file
        for line in f:

            # Get the tax ID and read ID for the current row
            read_id, tax_id = get_tax_read_ids(line)

            # If the tax ID is in the clade, write out the read ID and tax_ids.
            if tax_id in tax_ids:

                # Get the root tax_id for the current row
                root_tax_id = child2root_map[tax_id]

                # Output
                output_line = '\t'.join([str(tax_id), str(root_tax_id), read_id])
                o.write(output_line + '\n')

                # Increment count
                c += 1

            # Keep track of progress
            i += 1
            if i % 10000000 == 0:
                log.info("Parsed {} lines...".format(i))

    # Some final logging
    log.info("Done parsing file.")
    log.info("Found {} reads matching the specifications (total reads: {}).".format(c, i))
