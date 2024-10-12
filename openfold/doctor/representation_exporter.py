import seaborn as sns
import matplotlib.pyplot as plt
import os
import torch
import logging
import subprocess
import numpy as np

logger = logging.getLogger(__file__)
logger.setLevel(level=logging.DEBUG)

class RepresentationExporter:
    def __init__(self, model, output_dir="heatmaps", export_msa=True, export_pair=True):
        self.model = model
        self.output_dir = output_dir
        self.export_msa = export_msa
        self.export_pair = export_pair
        self.iteration = 0
        os.makedirs(self.output_dir, exist_ok=True)
        #TODO these should be created only when needed...
        os.makedirs(os.path.join(self.output_dir, "msa"), exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "pair"), exist_ok=True)

        # set callback
        self.model.representation_hook = self._representation_hook

    def _representation_hook(self, msa_representation, pair_representation, stage, iteration):
        if self.export_msa:
            self._heatmap(msa_representation, stage, iteration, "msa")

        if self.export_pair:
            self._heatmap(pair_representation, stage, iteration, "pair")

    def _heatmap(self, data, stage, iteration, which):
        if data is None:
            logger.warning(f"Warning: No data for {title}.")
            return

        if isinstance(data, torch.Tensor):
            data = data.detach().cpu().numpy()
        
        stage_num = 0 if stage == "before" else 1
        frame_number = iteration * 2 + (0 if stage == "before" else 1)
        # "squash" the third dimension using average along third axis
        # N.B. msa shape: (N_alignments, seq_length, msa_dim), pair shape: (seq_length, seq_length, pair_dim)
        #TODO here we could experiment with PCA, e.g. (just jotted, not tested): 
        # from sklearn.decomposition import PCA
        # pca = PCA(n_components=1)
        # reduced_representation = pca.fit_transform(data.reshape(-1, data.shape[-1])).reshape(data.shape[:-1])
        
        # export data *averaged* along z axis
        filename = f"{which}_avg_{frame_number:02d}"
        title = f"{which} representation {stage} recycling {iteration}"
        reduced_representation = data.mean(axis=-1)
        self._save_heatmap(reduced_representation, title, which)         

        # export data, *median* along z axis
        filename = f"{which}_median_{frame_number:02d}"
        title = f"{which} representation {stage} recycling {iteration}"
        reduced_representation = np.median(data, axis=-1)
        self._save_heatmap(reduced_representation, title, which)         

        # export data, *max* along z axis
        filename = f"{which}_max_{frame_number:02d}"
        title = f"{which} representation {stage} recycling {iteration}"
        reduced_representation = np.max(data, axis=-1)
        self._save_heatmap(reduced_representation, title, which)         
   
    def _save_heatmap(self, reduced_representation, title, which):
        plt.figure(figsize=(12, 8))
        sns.heatmap(reduced_representation, cmap='viridis', cbar=True)
        plt.title(title)
        plt.xlabel("seq length")
        plt.ylabel("N alignments") if which == "msa" else plt.ylabel("seq length")
        filepath = os.path.join(os.path.join(self.output_dir, which), f"{filename}.png")
        plt.savefig(filepath)
        plt.close()
        logger.debug(f"Heatmap saved: {filepath}")

    def pngs_to_mpg(self, framerate=1):
        if self.export_msa:
            self._pngs_to_mpg("msa", framerate)
        if self.export_pair:
            self._pngs_to_mpg("pair", framerate)

    def _pngs_to_mpg(self, which, framerate=1):
        pngs_output_dir = os.path.join(self.output_dir, which)
        frame_pattern = os.path.join(pngs_output_dir, f"{which}_%02d.png")
        output_path = os.path.join(pngs_output_dir, f"{which}_movie.mp4")

        ffmpeg_cmd = [
            "ffmpeg", "-framerate", str(framerate), "-i", frame_pattern,
            "-c:v", "libx264", "-pix_fmt", "yuv420p", output_path
        ]

        try:
            subprocess.run(ffmpeg_cmd, check=True)
            logger.info(f"{which} heatmap movie saved: {output_path}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error generating movie: {e}")

