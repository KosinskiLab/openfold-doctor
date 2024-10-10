import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
from openfold.data.parsers import *
import logging

logger = logging.getLogger(__file__)
logger.setLevel(level=logging.DEBUG)

MAX_SEQ_LEN = 100 

class SequenceCoveragePlotter:
    def __init__(self, out_dir):
        self.out_dir = out_dir               

    def _compute_sequence_coverage(self, msa):
        msa_array = np.array([list(seq) for seq in msa])
        coverage = np.sum(msa_array != '-', axis=0)
        
        return coverage

    def _plot_sequence_coverage(self, coverage, query_sequence, out_file_name):
        plt.figure(figsize=(10, 4))
        
        plt.plot(coverage, color='black', linewidth=1)
        plt.xlabel('Position in sequence')
        plt.ylabel('Sequence coverage')
        plt.title('Sequence coverage from MSA')
        
        # query sequence as x-ticks (if not too huge)
        if query_sequence:
            plt.xticks(ticks=np.arange(len(query_sequence)), labels=list(query_sequence), fontsize=8, rotation=90)
        
        plt.grid(True)
        plt.tight_layout()
        # plt.show()
        plt.savefig(os.path.join(self.out_dir, out_file_name))
        plt.close()
        logger.info(f"Sequence coverage plot saved as {os.path.join(self.out_dir, out_file_name)}")
    
    # taken from https://github.com/sokrypton/ColabFold/blob/main/colabfold/plot.py
    def _plot_msa_v2(tag, feature_dict, sort_lines=True, dpi=100):
        seq = feature_dict["msa"][0]
        if "asym_id" in feature_dict:
          Ls = [0]
          k = feature_dict["asym_id"][0]
          for i in feature_dict["asym_id"]:
            if i == k: Ls[-1] += 1
            else: Ls.append(1)
            k = i
        else:
          Ls = [len(seq)]    
        Ln = np.cumsum([0] + Ls)

        try:
            N = feature_dict["num_alignments"][0]
        except:
            N = feature_dict["num_alignments"] 
        
        msa = feature_dict["msa"][:N]
        gap = msa != 21
        qid = msa == seq
        gapid = np.stack([gap[:,Ln[i]:Ln[i+1]].max(-1) for i in range(len(Ls))],-1)
        lines = []
        Nn = []
        for g in np.unique(gapid, axis=0):
            i = np.where((gapid == g).all(axis=-1))
            qid_ = qid[i]
            gap_ = gap[i]
            seqid = np.stack([qid_[:,Ln[i]:Ln[i+1]].mean(-1) for i in range(len(Ls))],-1).sum(-1) / (g.sum(-1) + 1e-8)
            non_gaps = gap_.astype(float)
            non_gaps[non_gaps == 0] = np.nan
            if sort_lines:
                lines_ = non_gaps[seqid.argsort()] * seqid[seqid.argsort(),None]
            else:
                lines_ = non_gaps[::-1] * seqid[::-1,None]
            Nn.append(len(lines_))
            lines.append(lines_)
        
        Nn = np.cumsum(np.append(0,Nn))
        lines = np.concatenate(lines,0)
        plt.figure(figsize=(8,5), dpi=dpi)
        plt.title("Sequence coverage")
        plt.imshow(lines,
                  interpolation='nearest', aspect='auto',
                  cmap="rainbow_r", vmin=0, vmax=1, origin='lower',
                  extent=(0, lines.shape[1], 0, lines.shape[0]))
        for i in Ln[1:-1]:
            plt.plot([i,i],[0,lines.shape[0]],color="black")
        for j in Nn[1:-1]:
            plt.plot([0,lines.shape[1]],[j,j],color="black")
        
        plt.plot((np.isnan(lines) == False).sum(0), color='black')
        plt.xlim(0,lines.shape[1])
        plt.ylim(0,lines.shape[0])
        plt.colorbar(label="Sequence identity to query")
        plt.xlabel("Positions")
        plt.ylabel("Sequences")
        
        coverage_png = os.path.join(self.out_dir, f"{tag}_coverage.png")
        plt.savefig(str(coverage_png), bbox_inches='tight')
        plt.close()
        logger.info(f"Sequence coverage plot saved as {coverage_png}")


    def plot(self, a3m_file, out_file_name):
        with open(a3m_file, "r") as fp:
            msa = parse_a3m(fp.read())    
        # print(msa.sequences)

        query_sequence = msa.sequences[0] if len(msa.sequences[0]) < MAX_SEQ_LEN else None   
        coverage = self._compute_sequence_coverage(msa.sequences)    
        self._plot_sequence_coverage(coverage, query_sequence, out_file_name)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate sequence coverage plot from a3m file")
    parser.add_argument("a3m_file", help="Path to .a3m file")
    parser.add_argument("out_file_name", help="Plot file name")

    args = parser.parse_args()
    sc = SequenceCoveragePlotter(".")
    sc.plot(args.a3m_file, args.out_file_name)

