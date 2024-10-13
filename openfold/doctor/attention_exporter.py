import matplotlib.pyplot as plt
import numpy as np
import os
import logging
logger = logging.getLogger(__file__)
logger.setLevel(level=logging.DEBUG)

class AttnExporter:
    def __init__(self, model, args, output_dir, avg_only=False):
        self.model = model
        self.args = args
        self.output_dir = output_dir
        self.avg_only = avg_only
        self.col_calls = 0
        self.row_calls = 0
        self.col_dir = os.path.join(output_dir, "col")
        self.row_dir = os.path.join(output_dir, "row")
        os.makedirs(self.col_dir, exist_ok=True)
        os.makedirs(self.row_dir, exist_ok=True)

        # set callbacks
        for block in self.model.evoformer.blocks:
            block.msa_att_row.mha.save_attn_callback = self._attn_row_callback
        for block in self.model.evoformer.blocks:
            if not block.no_column_attention:
                block.msa_att_col._msa_att.mha.save_attn_callback = self._attn_col_callback 

    def _attn_col_hook(self, attn_block, input, output):
        logger.debug("col hook")
    

    def _attn_row_hook(self, attn_block, input, output):
        logger.debug("row hook")
        _m = output.detach().cpu().numpy()
        logger.debug(f"m: {_m.shape}")

    def _attn_col_callback(self, a):
        _a = a.detach().cpu().numpy()
        logger.debug(f"a_col: {_a.shape}")
        num_channels = 32
        num_residues, num_heads, x_dim, y_dim = _a.shape
        for head in range(num_heads):
            avg_data = np.mean(_a[:, head, :num_channels, :num_channels], axis=0)
            filename = os.path.join(self.col_dir, f"col_attn_call_{self.col_calls}_head_{head}_res_avg.png")
            title = f"Head {head} mean over residue indices"
            self._save_heatmap(avg_data, title, filename)
            
            if not self.avg_only:
                for res_id in range(num_residues):
                    if res_id in [40, 60, 80]:  # temp. to replicate alphafold suppl
                        filename = os.path.join(self.col_dir, f"col_attn_call_{self.col_calls}_head_{head}_res_{res_id}.png")
                        title = f"Head {head} residue index {res_id}"
                        self._save_heatmap(_a[res_id, head, :num_channels, :num_channels], title, filename)
        self.col_calls += 1

    def _attn_row_callback(self, a):
        _a = a.detach().cpu().numpy()
        logger.debug(f"a_row: {_a.shape}")
        slice_dim, num_heads, x_dim, y_dim = _a.shape
        num_residues = x_dim  # == y_dim
        for head in range(num_heads):
            avg_data = np.mean(_a[:num_residues, head, :, :], axis=0)
            filename = os.path.join(self.row_dir, f"row_attn_call_{self.row_calls}_head_{head}_res_avg.png")
            title = f"Head {head} mean over residue indices"
            self._save_heatmap(avg_data, title, filename)

            if not self.avg_only:
                for res_id in range(num_residues):
                    if res_id in [14, 28, 20]:  # temp. to replicate alphafold suppl
                        filename = os.path.join(self.row_dir, f"row_attn_call_{self.row_calls}_head_{head}_res_{res_id}.png")
                        title = f"Head {head} residue index {res_id}"
                        self._save_heatmap(_a[res_id, head, :, :], title, filename)
        self.row_calls += 1

    def _save_heatmap(self, data, title, filename):
            plt.figure()
            plt.imshow(data, cmap='hot', interpolation='nearest')
            plt.title(title)
            plt.colorbar()
            plt.savefig(filename)
            plt.close()
            logger.info(f"Attention heatmap saved as {filename}")
