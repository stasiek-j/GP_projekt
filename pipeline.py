import logging
import subprocess
import argparse
import shutil
from download import SequencesContainer

from Bio import SeqIO

def download_proteomes():
    """Downloads proteome sequences from NCBI's Proteome database"""
    pass


if __name__ == '__main__':
    parser = argparse.ArgumentParser("Supertrees project pipeline")
    parser.add_argument('input_file', help="Input file containing names of proteomes to run analysis on", type=str)
    parser.add_argument("-n", '--num', help="Maximum number of species", type=int, default=1000)
    parser.add_argument('--output_root', help="Output root directory.", type=str, default='./data')
    args = parser.parse_args()

    logging.basicConfig()
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)


    sequences = SequencesContainer(args.input_file, args.output_root)
    sequences.download()


