import seaborn as sns
import matplotlib.pyplot as plt
import os
import torch
import logging

logger = logging.getLogger(__file__)
logger.setLevel(level=logging.DEBUG)

class PDBExporter:
    def __init__(self, model, output_dir):
        self.model = model
        self.output_dir = output_dir
        self.total_block_calls = 0
        os.makedirs(self.output_dir, exist_ok=True)

        # set callback
        [block.register_forward_hook(self._structure_hook) for block in self.model.evoformer.blocks]

    def _structure_hook(self, module, input, output):
        logger.debug(f"input: {input}")
        logger.debug(f"output: {output}")
        logger.debug(f"ciao {self.total_block_calls}")
        self.total_block_calls += 1

