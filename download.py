"""
Utilities for downloading data from NCBI's proteome database.
"""
import logging
import os
import subprocess
import urllib
from urllib import request
from urllib.parse import quote
import requests as r
import pandas as pd
from requests import HTTPError
from tqdm import tqdm


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
            rest_link_prot = f"https://rest.uniprot.org/uniprotkb/search?format=fasta&query=%28%28proteome%3A{str(row['uniprot'])}%29%29"  # f"https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28proteome%3{row['uniprot']}%29%29"
            with r.get(rest_link_prot, stream=True) as req:
                try:
                    req.raise_for_status()
                    path = self.output_path
                    path_file = os.path.join(self.output_path,
                                             row['organism'].replace('.', '_').replace(" ", "_") + '.fasta')
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


class ProteomeDownloader:
    def __init__(self, path, output_path='./data/'):
        self.ftps = None
        self.ids = pd.read_csv(path)
        self.output_path = output_path
        pass

    def get_ftps(self, ids=None, key="assembly_accession"):
        if ids is None:
            ids = self.ids
        # organisms = pd.read_csv("lcr_organisms.csv")
        summary = pd.read_csv('assembly_summary', sep="\t")[["assembly_accession", "ftp_path"]]

        ftps = []
        for acc in ids:  # organisms["assembly_id"].values:
            list_ftps = summary[summary[key] == acc]["ftp_path"].values.tolist()
            ftps += list_ftps if list_ftps != [] else [""]

        self.ftps = ftps

        # organisms["ftp"] = ftps
        # organisms.to_csv("lcr_organisms.csv")

    def download_ftps(self, output_path='', unzip=True, n=-1):
        if self.output_path:
           output_path = self.output_path
        def list_url(iterab):
            return [x for x in iterab if isURL(x)]

        def isURL(string):
            """
            Bardzo uproszczona funkcja czy url, ale na te potrzeby wystarcza
            """
            if string == string:
                return string[:4] == 'http'
            return False

        def download(url, filename, pbar=None):
            remote = '/' + url.split("/")[-1] + '_protein.faa.gz'
            if pbar and logging.DEBUG >= logging.root.level:
                pbar.set_description(url + remote)
            else:
                logging.debug(url + remote)
            try:
                request.urlretrieve(url + remote, filename + '_protein.faa.gz')
            except (urllib.error.HTTPError, urllib.error.URLError):
                logging.debug(f"HTTP error {url}")
                pass

        try:
            assert self.ftps is not None, f"No ftps found"
            for i, ftp in enumerate(list_url(self.ftps)):
                download(ftp, output_path + str(i))
        except AssertionError:
            ftps = self.ids['ftp'][:n]
            pbar = tqdm(list_url(ftps))
            for i, ftp in enumerate(pbar):
                filename = "_".join(self.ids["name"][i].split(" ")[:2])
                download(ftp, output_path + filename, pbar)
                if unzip:
                    subprocess.run(['gunzip', os.path.join(output_path, filename + '_protein.faa.gz')])


if __name__ == '__main__':
    logging.basicConfig()
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    sequences = ProteomeDownloader('organisms.csv', 'data/EDF/')
    sequences.download_ftps()
    # sequences.save_species()
