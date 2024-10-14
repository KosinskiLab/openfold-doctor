import os
import logging
import numpy as np

logger = logging.getLogger(__file__)
logger.setLevel(level=logging.DEBUG)

class MSAExporter:
    def __init__(self, model, args, output_dir):
        self.model = model
        self.args = args
        self.output_dir = output_dir

        # set callback
        [block.register_forward_hook(self._export_msa) for block in self.model.evoformer.blocks]

    def _export_msa(self, module, input, output):
        logger.debug(f"type input {type(input)}, type ouput {type(output)}")
        logger.debug(f"input {input}")
        logger.debug(f"output {output}")

