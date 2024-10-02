import seaborn as sns
import matplotlib.pyplot as plt
import os
import torch
import logging
import numpy as np
from openfold.utils.tensor_utils import tensor_tree_map
from openfold.utils.feats import atom14_to_atom37

logger = logging.getLogger(__file__)
logger.setLevel(level=logging.DEBUG)

class PDBExporter:
    def __init__(self, model, output_dir):
        self.model = model
        self.batch = None
        self.output_dir = output_dir
        self.total_block_calls = 0
        self.no_blocks = 48
        os.makedirs(self.output_dir, exist_ok=True)

        # set callback
        self.model.register_forward_pre_hook(self._batch_hook)
        [block.register_forward_hook(self._structure_hook) for block in self.model.evoformer.blocks]

    def _batch_hook(self, module, input):
        self.batch = input

    def _structure_hook(self, module, input, output):
        m, z = output
        m = m.detach()
        z = z.detach()
        s = self.model.evoformer.linear(m[..., 0, :, :]).detach()

        cycle_no = self.total_block_calls // self.no_blocks
        fetch_cur_batch = lambda t: t[..., cycle_no]
        feats = tensor_tree_map(fetch_cur_batch, self.batch)[0]  # altrimenti è una tupla; bah...

        logger.debug(f"feats: {feats}")

        # dtype = next(self.model.parameters()).dtype
        # for k in feats:
        #     if feats[k].dtype == torch.float32:
        #         feats[k] = feats[k].to(dtype=dtype)

        n_seq = feats["msa_feat"].shape[-3]
        # device = feats["target_feat"].device

        outputs = {}
        outputs["msa"] = m[..., :n_seq, :, :]
        outputs["pair"] = z
        outputs["single"] = s
        
        # Predict 3D structure
        outputs["sm"] = self.model.structure_module(
            outputs,
            feats["aatype"],
            mask=feats["seq_mask"],  # .to(dtype=s.dtype),
            inplace_safe=False,
            _offload_inference=self.model.globals.offload_inference,
        )
        outputs["final_atom_positions"] = atom14_to_atom37(
            outputs["sm"]["positions"][-1], feats
        )
        outputs["final_atom_mask"] = feats["atom37_atom_exists"]
        outputs["final_affine_tensor"] = outputs["sm"]["frames"][-1]

        self.total_block_calls += 1

