import gzip
import logging
import os
import subprocess
import argparse
import shutil
import tempfile
from collections import defaultdict
from typing import Dict, Tuple

from download import ProteomeDownloader
from utils import *

from Bio import SeqIO

TMP = tempfile.mkdtemp()
PATH_TO_FASTUREC = '/home/stasiek/fasturec/fasturec'
PATH_TO_RAPID = '/home/stasiek/rapidNJ-master/bin/rapidnj'


def get_clusters(merged: str, output: str = '', min_seq: int = 4, min_org: int = 1) -> Tuple[
    Dict[str, List], Dict[str, List]]:
    clusters_file = SeqIO.parse(merged, 'fasta')
    sizes = defaultdict(int)
    clusters = {}
    organisms = defaultdict(set)
    org_sizes = defaultdict(int)
    curr = ''
    for seq in clusters_file:
        if not seq.seq:
            curr = seq.id
            clusters[curr] = []
        elif make_valid(seq.id) not in [x.id for x in clusters[curr]]:
            left = seq.description.rfind('[')
            right = seq.description.rfind(']')
            organisms[curr].add(seq.description[left: right + 1])
            org_sizes[seq.description[left:right + 1]] += 1
            seq.id = seq.description[left: right + 1] + "_" + str(org_sizes[seq.description[left:right + 1]])
            if org_sizes[seq.description[left:right + 1]] < 10:
                print(seq.id)

            seq.description = ""
            seq.id = make_valid(seq.id)
            clusters[curr].append(seq)
            sizes[curr] += 1

    clusters = {k: v for k, v in clusters.items() if sizes[k] >= min_seq and len(organisms[k]) >= min_org}
    sizes = {k: v if k in clusters.items() else 0 for k, v in sizes.items()}
    clusters_ret = {k: v for k, v in clusters.items() if len(organisms[k]) >= min_org}
    # clusters = {k: clusters[k] for k in sorted(sizes.keys(), key=lambda x: sizes[x], reverse=True)[:100] if k in
    # clusters}
    if output:
        os.makedirs(os.path.dirname(output), exist_ok=True)
        for k, v in clusters_ret.items():
            SeqIO.write(v, os.path.join(output, k.split(' ')[0].split('.')[0]), 'fasta')
    return clusters, clusters_ret


def get_1_1_clusters(clusters: Dict[str, List], output: str = '', min_seq: int = 4) \
        -> Tuple[Dict[str, List], Dict[str, int]]:
    ret_cl, ret_sizes = defaultdict(list), defaultdict(int)
    logging.debug(len(clusters))
    for cluster in clusters:
        organisms = set()
        for seq in clusters[cluster]:
            # left = seq.description.rfind('[')
            # right = seq.description.rfind(']')
            if seq.id.rsplit("__")[0] in organisms:
                continue
            else:
                organisms.add(seq.id.rsplit("__")[0])
                seq.id = seq.id.rsplit("__")[0]
                seq.description = ""
                ret_cl[cluster].append(seq)
                ret_sizes[cluster] += 1
    logging.debug(f"Size of clusters before filtering: {len(ret_cl)}")
    ret_cl = {k: v for k, v in ret_cl.items() if len(v) == min_seq}
    ret_sizes = {k: v for k, v in ret_sizes.items() if k in ret_cl}
    if output:
        os.makedirs(os.path.dirname(output), exist_ok=True)
        for k, v in ret_cl.items():
            SeqIO.write(v, os.path.join(output, k.split(' ')[0].split('.')[0]), 'fasta')

    logging.debug(f"Size of clusters: {len(ret_cl)}")
    return ret_cl, ret_sizes


def align_family(clusters: List | str, tool: MultiAlignmentTool, cluster_name: str = ''):
    if isinstance(clusters, List):
        cluste = os.path.join(TMP, 'clusters/', f'{cluster_name if cluster_name else "family"}.fasta')
        os.makedirs(os.path.dirname(cluste), exist_ok=True)
        SeqIO.write(clusters, cluste, 'fasta')
        clusters = cluste
    msas = os.path.join(TMP, 'msas/', f'{cluster_name if cluster_name else "family"}_msa.fasta')
    os.makedirs(os.path.dirname(msas), exist_ok=True)

    return tool(clusters, msas), msas


def fasttree_family(msa: str, tool: TreeTool, cluster_name: str = '', name='fasttree'):
    output = os.path.join(TMP, f'{name}s/', f'{cluster_name if cluster_name else "family"}_{name}.fasta')
    os.makedirs(os.path.dirname(output), exist_ok=True)
    try:
        return tool(msa, output)
    except subprocess.CalledProcessError as e:
        logging.error(f"command not working problem with {msa}")
        raise e


def reconcile(trees: str, tool: SuperTreeTool | MajorityConsensusTool, pat=''):
    return tool(trees, pat=pat)


def align_and_tree(cls, _mafft, _fasttree, split='1_1/'):
    pbar = tqdm(cls)
    d = 0
    for cluster in pbar:
        obj, path = align_family(cls[cluster], _mafft, split + cluster)
        if obj:
            tre, disc = fasttree_family(path, _fasttree, split + cluster)
            d += disc
            pbar.set_description(f"Discarded: {d} files with bootstrap")
            if not tre:
                raise Exception("Could not create fasttree for cluster {}".format(cluster))


if __name__ == '__main__':
    parser = argparse.ArgumentParser("Supertrees project pipeline")
    parser.add_argument('input_file', help="Input file containing names of proteomes to run analysis on "
                                           "as well as urls to ftp containing proteomes.", type=str)
    parser.add_argument('--output_root', help="Output root directory.", type=str, default='./data/')
    parser.add_argument('-d', '--debug', help="Print lots of debugging statements",
                        action="store_const", dest="loglevel", const=logging.DEBUG, default=logging.WARNING)
    parser.add_argument(
        '-v', '--verbose', help="Be verbose", action="store_const", dest="loglevel", const=logging.INFO)
    parser.add_argument('-n', type=int, default=-1, help="Number of proteomes to analyse")
    # parser.add_argument('-T', type=int, default=6, help="Number of threads to use.")  # Can be used with RaxML
    parser.add_argument('--min_seq', type=int, default=4, help="Minimum number of sequences in cluster")
    # parser.add_argument("--num_seq", type=int, default=20)
    parser.add_argument('--min_org', type=int, default=2, help="Minimum number of organisms in a cluster")
    parser.add_argument('--min_seq_id', default=0.5, type=float, help="Minimum identity of sequence in clusters")
    parser.add_argument('--bootstrap', '-b', default=1, type=int, help="Number of bootstrap replicates "
                                                                       "to generate the trees.")
    parser.add_argument('--logfile', type=str, default=None, help="Logfile in which to write")
    args = parser.parse_args()

    if args.logfile is not None:
        logging.basicConfig(filename=args.logfile)
    else:
        logging.basicConfig()
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    logger.info("Downloading proteomes:")
    debugging_path = args.output_root
    os.makedirs(debugging_path, exist_ok=True)
    downloader = ProteomeDownloader(f'{args.input_file}', debugging_path)
    logger.debug("Created downloader")
    downloader.download_ftps(unzip=True, n=args.n)
    logger.info(f"Downloaded {downloader.count} proteomes")

    logger.info("Starting clustering:")

    cluster_params = {'-v': '0', '--min-seq-id': f"{args.min_seq_id}"}
    cluster_tool = Mmseq2('mmseqs', cluster_params)
    logger.debug("Created clustering tool with params: {}".format(cluster_params))

    merge_files(debugging_path)
    os.makedirs(os.path.join(debugging_path, 'clusters/'), exist_ok=True)
    logger.debug("Merged fastas. \nStart clustering.")

    cluster_tool(os.path.join(debugging_path, 'merged.fasta'), os.path.join(debugging_path, 'clusters/'))
    logger.debug("Clustered.")

    logger.debug("Preprocessing clusters:")

    cls, _ = get_clusters(os.path.join(debugging_path, 'clusters/', '_all_seqs.fasta'),
                          os.path.join(debugging_path, 'clusters/', 'families_paralogs/'),
                          min_seq=args.min_seq, min_org=args.min_org)
    logger.debug("Created clusters with paralogs and saved.")

    # cluster_tool_no_para = Mmseq2('mmseqs', {'-v': '0'})
    # os.makedirs(os.path.join(debugging_path, 'clusters_1/'), exist_ok=True)
    # cluster_tool_no_para(os.path.join(debugging_path, 'merged.fasta'), os.path.join(debugging_path, 'clusters_1/'))
    #
    # cls_1, _ = get_clusters(os.path.join(debugging_path, 'clusters_1/', '_all_seqs.fasta'),
    #                         os.path.join(debugging_path, 'clusters_1/', 'families_paralogs/'),
    #                         min_seq=args.min_seq, min_org=args.min_org)

    cls_1_1, siz_1_1 = get_1_1_clusters(cls, os.path.join(debugging_path, 'clusters_1/', 'families_unique/'),
                                        downloader.count) # if downloader.count < args.num_seq else args.num_seq)
    logger.debug("Created clusters without paralogs and saved.")

    logger.info("Starting aligning")

    mafft = MAFFTool('', {})
    logger.debug("Created mafft tool.")
    #
    fasttree = FasttreeTool('fasttreeMP', {})
    # raxml = RAxMLTool('/home/stasiek/Pulpit/BIBS/GP/standard-RAxML/raxmlHPC-PTHREADS-AVX',
    #                   {'-T': f'{args.T}', '-x': '123', '-#': '1'})
    # raxml = fasttree
    raxml = RapidNJTool(PATH_TO_RAPID, {'-b': f'{args.bootstrap}'})
    logger.debug("Created tree tool.")
    #
    align_and_tree(cls, mafft, raxml, 'paralogs/')
    logger.debug("Created alignments and trees for clusters with paralogs")
    #
    align_and_tree(cls_1_1, mafft, raxml)
    logger.debug("Created alignments and trees for clusters without paralogs")
    logger.info("Finished clustering and created trees using FastTree tool.")

    logger.info("Starting creating trees of genomes.")

    # TMP = '/tmp/tmp0m_acmic'
    merge_files(os.path.join(TMP, 'fasttrees/paralogs/'), output=os.path.join(TMP, 'merged.newick'), pat='')
    logger.debug(f"Merged files in {os.path.join(TMP, 'fasttrees/paralogs/')}")

    supertreetool = FasturecTreeTool(PATH_TO_FASTUREC, {})
    logger.debug("Created supertree tool.")

    consensus = MajorityConsensusTool('', {'schema': 'newick'})
    logger.debug("Created consensus tool.")

    consensus_tree = reconcile(os.path.join(TMP, 'fasttrees/1_1/'), consensus, '')
    logger.debug(consensus_tree)
    logger.info(f"Created consensus tree.\n{dendropy.Tree.as_ascii_plot(consensus_tree)}")
    consensus_tree.write_to_path(os.path.join(debugging_path, 'consensus_tree.nwk'), schema='newick')

    logger.info("Starting working on supertree")


    supertree = reconcile(os.path.join(TMP, 'merged.newick'), supertreetool)
    logging.debug(supertree)
    supertree = dendropy.Tree.get_from_string(supertree, schema='newick')
    logger.info(f"Created supertree.\n{dendropy.Tree.as_ascii_plot(supertree)}")
    supertree.write_to_path(os.path.join(debugging_path, 'supertree.nwk'), schema='newick')

    shutil.rmtree(TMP)
    logger.info("Finished.")

