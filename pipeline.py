import gzip
import logging
import subprocess
import argparse
import shutil
import tempfile
from collections import defaultdict
from typing import Dict, Tuple

from download import ProteomeDownloader
from tools import *

from Bio import SeqIO

TMP = tempfile.mkdtemp()
PATH_TO_DUPTREE = '/home/stasiek/Pulpit/BIBS/GP/lab9/linux-i386/duptree'




def get_clusters(merged: str, output: str = '', min_seq: int = 5) -> Tuple[Dict[str, List], Dict[str, int]]:
    clusters_file = SeqIO.parse(merged, 'fasta')
    sizes = defaultdict(int)
    clusters = {}
    curr = ''
    for seq in clusters_file:
        if not seq.seq:
            curr = seq.id
            clusters[curr] = []
        elif make_valid(seq.id) not in [x.id for x in clusters[curr]]:
            seq.description = make_valid(seq.description)
            seq.id = make_valid(seq.id)
            clusters[curr].append(seq)
            sizes[curr] += 1

    clusters = {k: v for k, v in clusters.items() if sizes[k] >= min_seq}
    sizes = {k: v if k in clusters else 0 for k, v in sizes.items()}
    if output:
        os.makedirs(os.path.dirname(output), exist_ok=True)
        for k, v in clusters.items():
            SeqIO.write(v, os.path.join(output, k.split(' ')[0].split('.')[0]), 'fasta')
    return clusters, sizes


def get_1_1_clusters(clusters: Dict[str, List], output: str = '', min_seq: int = 3) \
        -> Tuple[Dict[str, List], Dict[str, int]]:
    ret_cl, ret_sizes = defaultdict(list), defaultdict(int)
    logging.debug(len(clusters))
    for cluster in clusters:
        organisms = set()
        for seq in clusters[cluster]:
            left = seq.description.rfind('[')
            right = seq.description.rfind(']')
            if seq.description[left + 1:right] in organisms or left == -1 or right == -1:
                continue
            else:
                organisms.add(seq.description[left + 1:right])
                seq.id = seq.description[left + 1:right]
                seq.description = make_valid(seq.description)
                ret_cl[cluster].append(seq)
                ret_sizes[cluster] += 1
    logging.debug(f"Size of clusters before filtering: {len(ret_cl)}")
    ret_cl = {k: v for k, v in ret_cl.items() if len(v) >= min_seq}
    ret_sizes = {k: v for k, v in ret_sizes.items() if k in ret_cl}
    if output:
        os.makedirs(os.path.dirname(output), exist_ok=True)
        for k, v in ret_cl.items():
            SeqIO.write(v, os.path.join(output, k.split(' ')[0].split('.')[0]), 'fasta')

    logging.debug(f"Size of clusters: {len(ret_cl)}")
    return ret_cl, ret_sizes


def align_family(clusters: List | str, tool: MultiAlignmentTool, cluster_name: str =''):
    if isinstance(clusters, List):
        cluste = os.path.join(TMP, 'clusters/', f'{cluster_name if cluster_name else "family"}.fasta')
        os.makedirs(os.path.dirname(cluste), exist_ok=True)
        SeqIO.write(clusters, cluste, 'fasta')
        clusters = cluste
    msas = os.path.join(TMP, 'msas/', f'{cluster_name if cluster_name else "family"}_msa.fasta')
    os.makedirs(os.path.dirname(msas), exist_ok=True)

    return tool(clusters, msas), msas


def fasttree_family(msa: str, tool: FasttreeTool, cluster_name: str = ''):
    output = os.path.join(TMP, 'fasttrees/', f'{cluster_name if cluster_name else "family"}_fasttree.fasta')
    os.makedirs(os.path.dirname(output), exist_ok=True)
    try:
        return tool(msa, output)
    except subprocess.CalledProcessError as e:
        logging.error(f"command not working problem with {msa}")
        raise e


def reconcile(trees: str, tool: SuperTreeTool | MajorityConsensusTool):
    return tool(trees)


def align_and_tree(cls, _mafft, _fasttree, split='1_1/'):
    for cluster in tqdm(cls):
        obj, path = align_family(cls[cluster], _mafft, split + cluster)
        if not fasttree_family(path, _fasttree, split + cluster):
            raise Exception("Could not create fasttree for cluster {}".format(cluster))


if __name__ == '__main__':
    parser = argparse.ArgumentParser("Supertrees project pipeline")
    parser.add_argument('input_file', help="Input file containing names of proteomes to run analysis on", type=str)
    parser.add_argument('--output_root', help="Output root directory.", type=str, default='./data')
    args = parser.parse_args()

    logging.basicConfig()
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    logger.info("Downloading proteomes:")
    debugging_path = 'data/EDF_debug/'
    os.makedirs(debugging_path, exist_ok=True)
    downloader = ProteomeDownloader('organisms.csv', debugging_path)
    logger.debug("Created downloader")
    downloader.download_ftps(n=5, unzip=True)
    logger.info("Downloaded proteomes")

    logger.info("Starting clustering:")

    cluster_params = {'-v': '0'}
    cluster_tool = Mmseq2('mmseqs', cluster_params)
    logger.debug("Created clustering tool with params: {}".format(cluster_params))

    merge_files(debugging_path)
    os.makedirs(os.path.join(debugging_path, 'clusters/'), exist_ok=True)
    logger.debug("Merged fastas. \nStart clustering.")

    cluster_tool(os.path.join(debugging_path, 'merged.fasta'), os.path.join(debugging_path, 'clusters/'))
    logger.debug("Clustered.")

    logger.debug("Preprocessing clusters:")

    cls, siz = get_clusters(os.path.join(debugging_path, 'clusters/', '_all_seqs.fasta'),
                            os.path.join(debugging_path, 'clusters/', 'families_paralogs/'))
    logger.debug("Created clusters with paralogs and saved.")

    cls_1_1, siz_1_1 = get_1_1_clusters(cls, os.path.join(debugging_path, 'clusters/', 'families_unique/'))
    logger.debug("Created clusters without paralogs and saved.")

    logger.info("Starting aligning")

    mafft = MAFFTool('', {})
    logger.debug("Created mafft tool.")

    fasttree = FasttreeTool('fasttree', {})
    logger.debug("Created fasttree tool.")

    align_and_tree(cls, mafft, fasttree, 'paralogs/')
    logger.debug("Created alignments and trees for clusters with paralogs")

    align_and_tree(cls_1_1, mafft, fasttree)
    logger.debug("Created alignments and trees for clusters without paralogs")
    logger.info("Finished clustering and created trees using FastTree tool.")

    logger.info("Starting creating trees of genomes.")

    merge_files(os.path.join(TMP, 'fasttrees/paralogs/'), output=os.path.join(TMP, 'merged.nwick'))
    logger.debug(f"Merged files in {os.path.join(TMP, 'fasttrees/paralogs/')}")

    supertreetool = SuperTreeTool(PATH_TO_DUPTREE, {})
    logger.debug("Created supertree tool.")

    consensus = MajorityConsensusTool('', {'schema': 'newick'})
    logger.debug("Created consensus tool.")

    consensus_tree = reconcile(os.path.join(TMP, 'fasttrees/1_1/'), consensus)
    logger.debug(consensus_tree)
    logger.info(f"Created consensus tree.\n{dendropy.Tree.as_ascii_plot(consensus_tree)}")

    logger.info("Starting working on supertree")
    supertree = reconcile(os.path.join(TMP, 'merged.nwick'), supertreetool)
    logging.debug(supertree)
    logger.info(f"Created supertree.\n{dendropy.Tree.as_ascii_plot(supertree)}")

    # os.remove(TMP)
    logger.info("Finished.")



    # # sequences = SequencesContainer(args.input_file, args.output_root)
    # # sequences.download()
    #
    # # logger.info("Downloaded {} sequences".format(len(sequences)))
    #
    # cluster = Mmseq2('mmseqs', {'-v': '0'})
    # merge_fastas('data/EDF')
    # logger.debug('Done merging fastas')
    # # file_list = [os.path.join('data/EDF/', x) for x in os.listdir('data/EDF') if x.endswith('.faa')]
    # # print(file_list)
    # cluster('data/EDF/merged.fasta', 'out/')
    # mafft = MAFFTool('', {'--thread': -1})
    # # for file in file_list:
    # #     # with gzip.open(file, 'rb') as f_in:
    # #     #     with open('tmp_input', 'wb') as f_out:
    # #     #         shutil.copyfileobj(f_in, f_out)
    # #     print(f"Starting mafft for {file}")
    # #     print(mafft(f'{file}', ""))
    # #
    # # fasttree = FasttreeTool('fasttree', {})
    # # print(os.listdir('../../lab7/data/alignments'))
    # # fasttree('../../lab7/data/alignments', 'fasttree_out/')


