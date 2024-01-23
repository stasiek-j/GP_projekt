"""
Created on 8.01.24
@author: Stanisław Janik

Scripts for running bioinformatics tools using python interface.
"""
import logging
import os.path
import re
import tempfile
from abc import ABC, abstractmethod
import subprocess
from io import StringIO
from typing import List

import Bio.Application
import dendropy
from Bio import AlignIO
from Bio.Phylo.Consensus import majority_consensus
from Bio.Align import Applications
from Bio.Seq import Seq
from tqdm import tqdm


def flatten(xss):
    return [x for xs in xss for x in xs]


def merge_files(path, output=None, pat=''):
    if output is None:
        output = os.path.join(path, 'merged.fasta')
    with open(output, 'w+') as f:
        for file in os.listdir(path):
            if os.path.isfile(os.path.join(path, file)) and pat in file:
                with open(os.path.join(path, file), 'r') as fr:
                    f.write(fr.read())
    return output


def make_valid(s: str) -> str:
    for character in [' ', '.', '-', '/', '(', ')', ":", ",", ")", "(", ";", "]", "[", "'"]:
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
        tmpdir = tempfile.mkdtemp()
        mmseq2_cmd = [mmseq2_path, 'easy-cluster'] + input_file + [output_file_prefix, tmpdir,
                                                                   *flatten([(k, v) for k, v
                                                                             in self.parameters.items()])]
        proc = subprocess.run(mmseq2_cmd, check=True)
        subprocess.run(['rm', '-rf', tmpdir])
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
        cline = Applications.MafftCommandline(input=input_file,
                                              **self.parameters)
        try:
            ret = cline()[0]
        except Bio.Application.ApplicationError:
            return
        msa_ob = AlignIO.read(StringIO(ret), 'fasta')
        m = []
        for record in msa_ob:
            record.seq = Seq(str(record.seq).replace("J", "X"))
            m.append(record)

        msa_ob = AlignIO.MultipleSeqAlignment(m)
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
        output_dir = os.path.dirname(output_path)
        output_file = os.path.basename(output_path)[:-len('.fasta')]
        proc = subprocess.run([self.path, '-s', aligned_file, '-w', output_dir, '-n', output_file, '-m',
                               'PROTGAMMAGTR', '-p', '12345', *flatten([(k, v) for k, v
                                                                        in self.parameters.items()])],
                              stdout=subprocess.PIPE, check=True)
        return proc


class FasttreeTool(TreeTool):
    def __init__(self, path: str, parameters: dict):
        super().__init__(path, parameters)

    def __call__(self, align_path: str, output_path: str):
        if os.path.isfile(align_path):
            proc = subprocess.run([self.path, '-quiet', '-out', output_path, align_path, ], stdout=subprocess.PIPE,
                                  text=True,
                                  check=True)
            ret = proc.returncode
        elif os.path.isdir(align_path):
            assert not os.path.isfile(output_path), 'Alignment path is a directory and output path is not'
            if not os.path.exists(output_path):
                os.makedirs(output_path)
            ret = 0
            for file in tqdm(os.listdir(align_path)):
                proc = subprocess.run([self.path, '-quiet', '-out',
                                       os.path.join(output_path, file), os.path.join(align_path, file)],
                                      stdout=subprocess.PIPE, text=True, check=True)
                ret += proc.returncode
        else:
            ret = 1
        return ret == 0


class RapidNJTool(TreeTool):
    def __init__(self, path: str, parameters: dict):
        super().__init__(path, parameters)

    def __call__(self, align_path: str, output_path: str):
        def _bootstrap(path: str, b: int):
            with open(path, 'r') as f:
                data = f.read()
            confs = [float(x[1:-1])/b for x in re.findall(r'\)\d+:', data)]
            if sum(confs)/len(confs) < 0.75:
                os.remove(path)
                return 1
            return 0

        disc = 0
        if os.path.isfile(align_path):
            proc = subprocess.run([self.path, align_path, '-i', 'fa', '-n', '-x', output_path, *flatten([(k, v) for k, v
                                                                                    in self.parameters.items()])])
            # logging.debug("RapidNJTool finished, return code: {}".format(proc))
            if int(self.parameters['-b']) > 1:
                disc += _bootstrap(output_path, int(self.parameters['-b']))
            ret = proc.returncode
        elif os.path.isdir(align_path):
            assert not os.path.isfile(output_path), 'Alignment path is a directory and output path is not'
            if not os.path.exists(output_path):
                os.makedirs(output_path)
            ret = 0
            for file in tqdm(os.listdir(align_path)):
                proc = subprocess.run([self.path,
                    os.path.join(align_path, file), '-i', 'fa', '-n', '-x', os.path.join(output_path, file),
                    *flatten([(k, v) for k, v in self.parameters.items()])]
                )
                # logging.debug("RapidNJTool finished, return code: {}".format(proc))
                ret += proc.returncode
                if int(self.parameters['-b']) > 1:
                    disc += _bootstrap(os.path.join(output_path, file), int(self.parameters['-b']))

        else:
            ret = 1

        # logging.info(f'Discarded {disc} files with bootstrap')
        return ret == 0, disc


class MajorityConsensusTool(TreeTool):
    def __init__(self, path: str, parameters: dict):
        super().__init__(path, parameters)
        assert 'schema' in self.parameters.keys(), 'schema must be one of parameters.'

    def __call__(self, path, pat=''):
        try:
            if os.path.isfile(path):
                trees = dendropy.TreeList.get(path=path, schema=self.parameters['schema'])
            elif os.path.isdir(path):
                merged = merge_files(path, os.path.join(path, 'merged.newick'), pat=pat)
                # trees = dendropy.TreeList()
                # for filename in os.listdir(path):
                #     logging.debug(filename)
                #     if os.path.isfile(os.path.join(path, filename)):
                #         trees.read_from_path(src=os.path.join(path, filename), schema=self.parameters['schema'],)
                trees = dendropy.TreeList.get(path=merged, schema=self.parameters['schema'])
            else:
                raise FileNotFoundError("`path` must be file or directory containing trees to build consensus for.")
        except FileNotFoundError:
            return dendropy.Tree()
        logging.debug(trees.taxon_namespace)
        return trees.consensus()


class SuperTreeTool(TreeTool):
    def __init__(self, path: str, parameters: dict):
        super().__init__(path, parameters)

    def __call__(self, tree_file, pat=''):
        logging.debug(tree_file)
        with open(tree_file, 'r') as f:
            data = f.read()
            # logging.debug(data)
            data_re = re.sub(r'__\d+', "", data)
        with open(tree_file, 'w') as f:
            f.write(data_re)

        proc = subprocess.run([self.path, "-i", tree_file, "--nogenetree", ], stdout=subprocess.PIPE, text=True,
                              check=True)

        lines = proc.stdout.splitlines()
        return lines[-1]


class FasturecTreeTool(TreeTool):
    def __init__(self, path: str, parameters: dict):
        super().__init__(path, parameters)

    def __call__(self, tree_file, pat=''):
        with open(tree_file, 'r') as f:
            data = f.read()
        data = re.sub(r'[0-9]\.[0-9]+e-[0-9]+', lambda x: f"{float(x.group()): .9f}", data)
        # print(data)
        with open(tree_file, 'w') as f:
            f.write(data)
        sed = subprocess.run(["sed", "-i", "s/'//g", tree_file])
        sed = subprocess.run(["sed", "-i", 's/__[0-9]\+//g', tree_file])
        proc = subprocess.run([self.path, "-G", tree_file, "-Z", '-e', 'a', ], stdout=subprocess.PIPE, text=True,
                              check=True)

        with open('fu.txt') as f:
            line = f.readline()
        line = line.split(" ")[-1]

        # pat = r'\(\w+, *\w+\)'
        pat = r',\w+__set__'
        # line = re.sub(pat, lambda x: f"{x.group().split(',')[0]}", line)
        line = re.sub(pat, '', line)
        pat = r'\w+__set__,'
        line = re.sub(pat, '', line)
        subprocess.run(['rm', 'fu.txt'])
        return line + ';'


if __name__ == '__main__':
    test = FasturecTreeTool('../../fasturec/fasturec', {})
    print(test('../../merged.newick'))
