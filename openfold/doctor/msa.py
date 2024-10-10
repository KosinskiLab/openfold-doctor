import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
from openfold.data.parsers import *
import logging

logger = logging.getLogger(__file__)
logger.setLevel(level=logging.DEBUG)

MAX_SEQ_LEN = 100 

class SequenceCoverage:
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
        logger.info(f"Sequence coverage plot saved as {os.path.join(self.out_dir, out_file_name)}")

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
    sc = SequenceCoverage(".")
    sc.plot(args.a3m_file, args.out_file_name)

