"""
Created on 8.01.24
@author: Stanisław Janik

Scripts for running bioinformatics tools using python interface.
"""
import logging
import os.path
from abc import ABC, abstractmethod
import subprocess
from io import StringIO
from typing import List

import dendropy
from Bio import AlignIO
from Bio.Phylo.Consensus import majority_consensus
from Bio.Align import Applications
from tqdm import tqdm


def flatten(xss):
    return [x for xs in xss for x in xs]


def merge_files(path, output=None):
    if output is None:
        output = os.path.join(path, 'merged.fasta')
    with open(output, 'w+') as f:
        for file in os.listdir(path):
            if os.path.isfile(os.path.join(path, file)):
                with open(os.path.join(path, file), 'r') as fr:
                    f.write(fr.read())
    return output


def make_valid(s: str) -> str:
    for character in [' ', '.', '-', '/', '(', ')']:
        s = s.replace(character, "_")
    return s

class Tool(ABC):

    def __init__(self, path: str, parameters: dict):
        self.path = path
        self.parameters = parameters

    @abstractmethod
    def __call__(self, *args, **kwargs):
        pass


class ClusterTool(Tool):

    def __init__(self, path: str, parameters: dict):
        super().__init__(path, parameters)
        # assert '--n_clusters' in self.parameters

    @abstractmethod
    def __call__(self, *args, **kwargs):
        raise NotImplementedError


class Mmseq2(ClusterTool):
    def __init__(self, path: str, parameters: dict):
        super().__init__(path, parameters)

    def __call__(self, input_file: str | List[str], output_file_prefix: str):
        def _merge(f_in, f_out):
            with open(f_out, "w+") as fo:
                for file_in in f_in:
                    with open(file_in, "r") as fi:
                        fo.write(fi.read())


        mmseq2_path = self.path
        if isinstance(input_file, list):
            _merge(input_file, 'tmp_input')
            input_file = ['tmp_input']
        else:
            input_file = [input_file]
        mmseq2_cmd = [mmseq2_path, 'easy-cluster'] + input_file + [output_file_prefix, 'tmp/',
                                                                   *flatten([(k, v) for k, v
                                                                             in self.parameters.items()])]
        proc = subprocess.run(mmseq2_cmd, check=True)
        return proc.returncode


class MultiAlignmentTool(Tool):
    def __init__(self, path: str, parameters: dict):
        super().__init__(path, parameters)

    def __call__(self, *args, **kwargs):
        raise NotImplementedError


class MuscleTool(MultiAlignmentTool):
    def __init__(self, path: str, parameters: dict):
        super().__init__(path, parameters)

    def __call__(self, input_file: str, output_file: str):
        cline = Applications.MuscleCommandline(self.parameters)
        # subprocess.run(str(cline).split(), check=True)
        cline()
        return output_file


class MAFFTool(MultiAlignmentTool):
    def __init__(self, path: str, parameters: dict):
        super().__init__(path, parameters)

    def __call__(self, input_file: str, output_file: str):
        cline = Applications.MafftCommandline(input=input_file, **self.parameters)  # *flatten([(k, v) for k, v in self.parameters.items()]))
        ret = cline()[0]
        msa_ob = AlignIO.read(StringIO(ret), 'fasta')
        if output_file:
            AlignIO.write(msa_ob, output_file, 'fasta')
        return msa_ob


class TreeTool(Tool):
    def __init__(self, path: str, parameters: dict):
        super().__init__(path, parameters)

    def __call__(self, *args, **kwargs):
        raise NotImplementedError


class RAxMLTool(TreeTool):
    def __init__(self, path: str, parameters: dict):
        super().__init__(path, parameters)

    def __call__(self, aligned_file: str, output_path: str):
        proc = subprocess.run([self.path, '-s', aligned_file, '-w', output_path, '-m', 'PROTGAMMAGTR', '-p', '12345'],
                              stdout=subprocess.PIPE, check=True)
        return proc


class FasttreeTool(TreeTool):
    def __init__(self, path: str, parameters: dict):
        super().__init__(path, parameters)

    def __call__(self, align_path: str, output_path: str):
        if os.path.isfile(align_path):
            proc = subprocess.run([self.path, '-quiet', '-out', output_path, align_path, ], stdout=subprocess.PIPE, text=True,
                                  check=True)
            ret = proc.returncode
        elif os.path.isdir(align_path):
            assert not os.path.isfile(output_path), 'Alignment path is a directory and output path is not'
            if not os.path.exists(output_path):
                os.makedirs(output_path)
            ret = 0
            for file in tqdm(os.listdir(align_path)):
                proc = subprocess.run([self.path, '-quiet',  '-out',
                                       os.path.join(output_path, file), os.path.join(align_path, file)], stdout=subprocess.PIPE, text=True, check=True)
                ret += proc.returncode
        else:
            ret = 1
        return ret == 0


class MajorityConsensusTool(TreeTool):
    def __init__(self, path: str, parameters: dict):
        super().__init__(path, parameters)
        assert 'schema' in self.parameters.keys(), 'schema must be one of parameters.'

    def __call__(self, path):
        if os.path.isfile(path):
            trees = dendropy.TreeList.get(path=path, schema=self.parameters['schema'])
        elif os.path.isdir(path):
            merged = merge_files(path, os.path.join(path, 'merged.newick'))
            # trees = dendropy.TreeList()
            # for filename in os.listdir(path):
            #     logging.debug(filename)
            #     if os.path.isfile(os.path.join(path, filename)):
            #         trees.read_from_path(src=os.path.join(path, filename), schema=self.parameters['schema'],)
            trees = dendropy.TreeList.get(path=merged, schema=self.parameters['schema'])
        else:
            raise FileNotFoundError("`path` must be file or directory containing trees to build consensus for.")
        logging.debug(trees.taxon_namespace)
        return trees.consensus()


class SuperTreeTool(TreeTool):
    def __init__(self, path: str, parameters: dict):
        super().__init__(path, parameters)

    def __call__(self, tree_file):
        print(tree_file)
        proc = subprocess.run([self.path, "-i", tree_file, "--nogenetree"], stdout=subprocess.PIPE, text=True,
                              check=True)

        lines = proc.stdout.splitlines()
        return lines[-1]



