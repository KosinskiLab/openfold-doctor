import logging
import os
import numpy as np

from openfold.np import protein
from openfold.utils.tensor_utils import tensor_tree_map


logger = logging.getLogger(__file__)
logger.setLevel(level=logging.DEBUG)

class Singleton(type):
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


class Doctor(metaclass=Singleton):
    def __init__(self):
        self.feature_dict = None
        self._processed_feature_dict = None
        self.output_name = None
        self.output_dir = None
        self.feature_processor = None
        self.config_preset = None
        self.multimer_ri_gap = None
        #self.subtract_plddt = None
        self.cif_output = None
        self.nseq = 0
        self.in_use = False

    @property
    def processed_feature_dict(self):
        return self._processed_feature_dict

    @processed_feature_dict.setter
    def processed_feature_dict(self, val):
        self._processed_feature_dict = tensor_tree_map(lambda x: np.array(x.cpu()), val)
        assert self._processed_feature_dict is not None


    def intermediate_output(self, out, from_evoformer=False):

        out = tensor_tree_map(lambda x: np.array(x.cpu()), out)

        unrelaxed_protein = self.prep_intermediate_output(out)

        unrelaxed_file_suffix = f"{'_evoformer' if from_evoformer else ''}.pdb"  # _unrelaxed.pdb
        if self.cif_output:
            unrelaxed_file_suffix = f"{'_evoformer' if from_evoformer else ''}.cif"  # _unrelaxed.cif
        unrelaxed_output_path = os.path.join(
            self.output_dir, f'{self.output_name}_{self.nseq:04d}{unrelaxed_file_suffix}'
        )

        with open(unrelaxed_output_path, 'w') as fp:
            if self.cif_output:
                fp.write(protein.to_modelcif(unrelaxed_protein))
            else:
                fp.write(protein.to_pdb(unrelaxed_protein))

        logger.info(f"Output written to {unrelaxed_output_path}...")

        self.nseq += 1

        # if model.save_outputs:
        #     output_dict_path = os.path.join(
        #         model.output_dir, f'{model.output_name}_output_dict.pkl'
        #     )
        #     with open(output_dict_path, "wb") as fp:
        #         pickle.dump(out, fp, protocol=pickle.HIGHEST_PROTOCOL)
        #
        #     logger.info(f"Model output written to {output_dict_path}...")


    def prep_intermediate_output(self, out):
        # Prep protein metadata
        template_domain_names = []
        template_chain_index = None
        if self.feature_processor.config.common.use_templates and "template_domain_names" in self.feature_dict:
            template_domain_names = [
                t.decode("utf-8") for t in self.feature_dict["template_domain_names"]
            ]

            # This works because templates are not shuffled during inference
            template_domain_names = template_domain_names[
                                    :self.feature_processor.config.predict.max_templates
                                    ]

            if "template_chain_index" in self.feature_dict:
                template_chain_index = self.feature_dict["template_chain_index"]
                template_chain_index = template_chain_index[
                                       :self.feature_processor.config.predict.max_templates
                                       ]

        no_recycling = self.feature_processor.config.common.max_recycling_iters
        remark = ', '.join([
            f"no_recycling={no_recycling}",
            f"max_templates={self.feature_processor.config.predict.max_templates}",
            f"config_preset={self.config_preset}",
        ])

        # For multi-chain FASTAs
        ri = self.feature_dict["residue_index"]
        chain_index = (ri - np.arange(ri.shape[0])) / self.multimer_ri_gap
        chain_index = chain_index.astype(np.int64)
        cur_chain = 0
        prev_chain_max = 0
        batch = self.processed_feature_dict
        for i, c in enumerate(chain_index):
            if c != cur_chain:
                cur_chain = c
                prev_chain_max = i + cur_chain * self.multimer_ri_gap

            batch["residue_index"][i] -= prev_chain_max

        unrelaxed_protein = protein.from_prediction(
            features=batch,
            result=out,
            remove_leading_feature_dimension=False,
            #b_factors=None,
            #remark=remark,
            #parents=template_domain_names,
            #parents_chain_index=template_chain_index,
        )

        return unrelaxed_protein

dr = Doctor()
