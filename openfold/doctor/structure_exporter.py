import seaborn as sns
import matplotlib.pyplot as plt
import os
import torch
import logging
import numpy as np
from openfold.utils.tensor_utils import tensor_tree_map
from openfold.utils.feats import atom14_to_atom37
from openfold.utils.script_utils import prep_output
from openfold.np import protein
from openfold.doctor.movie import ProteinMovieMaker

logger = logging.getLogger(__file__)
logger.setLevel(level=logging.DEBUG)

class PDBExporter:
    def __init__(self, model, feature_dict, feature_processor, args, output_dir):
        self.model = model
        self.feature_dict = feature_dict
        self.feature_processor = feature_processor
        self.args = args
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

        # cycle_no = self.total_block_calls // self.no_blocks
        cycle_no = self.model._cycle_no
        # iter_num = self.model.iter_num
        try:
            fetch_cur_batch = lambda t: t[..., cycle_no]
            feats = tensor_tree_map(fetch_cur_batch, self.batch)[0]  # altrimenti è una tupla; bah...
        except:
            logger.error(f"There is something fishy here... total block calls: {self.total_block_calls}")
            fetch_cur_batch = lambda t: t[..., -1]
            feats = tensor_tree_map(fetch_cur_batch, self.batch)[0]  # altrimenti è una tupla; bah...

        # logger.debug(f"feats: {feats}")

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

        # model.py: 606-610 
        try:
            outputs["asym_id"] = feats["asym_id"]
        except:
            # logger.debug("asym_id not in feats")
            pass

        # Run auxiliary heads
        outputs.update(self.model.aux_heads(outputs))

        # run_pretrained_openfold.py: ~378-393
        # Toss out the recycling dimensions --- we don't need them anymore
        processed_feature_dict = tensor_tree_map(
            lambda x: np.array(x[..., -1].cpu()),
            self.batch
        )[0]  # altrimenti è una tupla; bah...
        # logger.debug(f"processed feature dict: {processed_feature_dict}")
        outputs = tensor_tree_map(lambda x: np.array(x.cpu()), outputs)

        the_protein = prep_output(
            outputs,
            processed_feature_dict,
            self.feature_dict,
            self.feature_processor,
            self.args.config_preset,
            self.args.multimer_ri_gap,
            self.args.subtract_plddt
        )

        self._save_structure(the_protein)

        self.total_block_calls += 1


    def _save_structure(self, the_protein):
        file_suffix = "evoformer.pdb"
        if self.args.cif_output:
            file_suffix = "evoformer.cif"
        # recycle = self.total_block_calls // self.no_blocks
        # block = self.total_block_calls % self.no_blocks
        output_path = os.path.join(
            self.output_dir, f'{self.total_block_calls:03d}_{file_suffix}'
        )

        with open(output_path, 'w') as fp:
            if self.args.cif_output:
                fp.write(protein.to_modelcif(the_protein))
            else:
                fp.write(protein.to_pdb(the_protein))
        logger.info(f"Output written to {output_path}...")

    def make_movie(self):
        mmaker = ProteinMovieMaker(
            input_directory=self.output_dir,
            frame_duration_seconds=self.args.frame_duration_seconds,
            low_res=self.args.low_res_movie,
            keep_data=self.args.keep_movie_data
        )
        mmaker.run()
