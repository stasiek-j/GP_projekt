"""
Utilities for downloading data from NCBI's proteome database.
"""
import logging
import os
from urllib.parse import quote
import requests as r
import pandas as pd
from requests import HTTPError


class SequencesContainer:
    def __init__(self, path, output_path='./data'):
        self.species = pd.read_csv(path, sep='\t')
        self.output_path = output_path

        self.__preprocess()
        self.__get_taxids()
        self.__get_uniprotID()

    def __get_taxids(self):
        # uniprot = []
        tax = []
        for organism in self.species['organism']:
            logging.info(f'Processing organism {organism}')
            organism = organism.replace(' ', '_').replace(".", "")
            # organism = quote(organism)
            query = f'https://rest.uniprot.org/proteomes/stream?fields=upid%2Corganism%2Corganism_id&format=tsv&query={organism}'
            try:
                taxids = pd.read_csv(query, sep='\t')
                # uniprot.append(taxids['Proteome Id'][0])
                tax.append(taxids['Organism Id'][0])
            except KeyError or ValueError:
                tax.append('')
                # uniprot.append("")
                raise ValueError(f'No taxids found for organism {organism}')
        self.species['taxids'] = tax
        # self.species['uniprot'] = uniprot

    def __preprocess(self):
        if 'organism' not in self.species.columns:
            raise Exception('No organism column in species file')
        for column in ['taxids', 'seq', 'uniprot', 'file']:
            if column not in self.species.columns:
                self.species[column] = ""

    def __get_uniprotID(self):
        uniprot = []
        for taxid in self.species['taxids']:
            logging.debug(f'checking Uniprot ID for {taxid}')
            query = f'https://rest.uniprot.org/proteomes/stream?fields=upid%2Corganism_id&format=tsv&query={taxid}+AND+(proteome_type%3A1)'
            try:
                taxids = pd.read_csv(query, sep='\t')
                uniprot.append(taxids['Proteome Id'][0])
            except KeyError or ValueError:
                uniprot.append("")
                logging.error(f'No uniprot ID found for')
        self.species['uniprot'] = uniprot




    def download(self):
        def _download(row: pd.Series):
            logging.info(f"\n{row['uniprot']} ")
            rest_link_prot = f"https://rest.uniprot.org/uniprotkb/search?format=fasta&query=%28%28proteome%3A{str(row['uniprot'])}%29%29" #f"https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28proteome%3{row['uniprot']}%29%29"
            with r.get(rest_link_prot, stream=True) as req:
                try:
                    req.raise_for_status()
                    path = self.output_path
                    path_file = os.path.join(self.output_path, row['organism'].replace('.', '_').replace(" ", "_") + '.fasta')
                    try:
                        with open(path_file, 'wb') as f:
                            for chunk in req.iter_content(chunk_size=2 ** 20):
                                logging.debug(chunk.decode())
                                f.write(chunk)
                    except FileNotFoundError:
                        os.mkdir(path)
                        with open(path_file, 'wb') as f:
                            for chunk in req.iter_content(chunk_size=2 ** 20):
                                f.write(chunk)
                except HTTPError:
                    path_file = None
                    logging.error(f"No Proteome fasta found for organism {row['organism']}")
            return path_file

        files = []
        for row in self.species.iterrows():
            files.append(_download(row[1]))
        self.species['file'] = files
        logging.info(f'Downloaded {len(self.species)} species: examples:\n{self.species.head()}')

    def save_species(self):
        self.species.to_csv(f"{self.output_path}/species.csv", index=False)
        logging.info(f'Saved species to file {self.output_path}/species.csv')


if __name__ == '__main__':
    logging.basicConfig()
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    sequences = SequencesContainer('test')
    sequences.download()
    sequences.save_species()
