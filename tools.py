"""
Created on 8.01.24
@author: Stanisław Janik

Scripts for running bioinformatics tools using python interface.
"""
import os.path
from abc import ABC, abstractmethod
import subprocess
from Bio.Phylo.Consensus import majority_consensus
from Bio.Align import Applications


def flatten(xss):
    return [x for xs in xss for x in xs]


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
        assert '--n_clusters' in self.parameters

    @abstractmethod
    def __call__(self, *args, **kwargs):
        raise NotImplementedError


class Mmseq2(ClusterTool):
    def __init__(self, path: str, parameters: dict):
        super().__init__(path, parameters)

    def __call__(self, input_file: str, output_file: str):
        mmseq2_path = os.path.join(self.path, 'mmseqs2')
        mmseq2_cmd = [mmseq2_path, '--input', input_file, '--output', output_file,
                      *flatten([(k, v) for k, v in self.parameters.items()])]
        subprocess.run(mmseq2_cmd, check=True)


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
        subprocess.run(str(cline).split(), check=True)
        return output_file


class MAFFTool(MultiAlignmentTool):
    def __init__(self, path: str, parameters: dict):
        super().__init__(path, parameters)

    def __call__(self, input_file: str, output_file: str):
        cline = Applications.MafftCommandline(self.parameters)
        ret = cline(input_file)
        return ret


class TreeTool(Tool):
    def __init__(self, path: str, parameters: dict):
        super().__init__(path, parameters)

    def __call__(self, *args, **kwargs):
        raise NotImplementedError


class RAxMLTool(TreeTool):
    def __init__(self, path: str, parameters: dict):
        super().__init__(path, parameters)

    def __call__(self):
        # TODO
        raise NotImplementedError


class MLTool(TreeTool):
    # TODO
    def __init__(self, path: str, parameters: dict):
        super().__init__(path, parameters)

    def __call__(self):
        # TODO
        raise NotImplementedError


class MajorityConsensusTool(TreeTool):
    # TODO
    def __init__(self, path: str, parameters: dict):
        super().__init__(path, parameters)

    def __call__(self):
        raise NotImplementedError


class SuperTreeTool(TreeTool):
    # TODO
    def __init__(self, path: str, parameters: dict):
        super().__init__(path, parameters)

    def __call__(self):
        raise NotImplementedError
