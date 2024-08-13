import logging
import os
import numpy

from openfold.np import protein

logger = logging.getLogger(__file__)
logger.setLevel(level=logging.DEBUG)

nseq = 0

def intermediate_output(model, out, batch):
    #from openfold.utils.script_utils import prep_output
    global nseq

    feature_dict = model.feature_dict
    feature_processor = model.feature_processor
    config_preset = model.config_preset
    multimer_ri_gap = model.multimer_ri_gap
    #subtract_plddt = model.subtract_plddt

    unrelaxed_protein = prep_intermediate_output(out, batch, feature_dict, feature_processor, config_preset, multimer_ri_gap)

    unrelaxed_file_suffix = "_unrelaxed.pdb"
    if model.cif_output:
        unrelaxed_file_suffix = "_unrelaxed.cif"
    unrelaxed_output_path = os.path.join(
        model.output_dir, f'{model.output_name}_{nseq:04d}{unrelaxed_file_suffix}'
    )

    with open(unrelaxed_output_path, 'w') as fp:
        if model.cif_output:
            fp.write(protein.to_modelcif(unrelaxed_protein))
        else:
            fp.write(protein.to_pdb(unrelaxed_protein))

    logger.info(f"Output written to {unrelaxed_output_path}...")

    nseq += 1

    # if model.save_outputs:
    #     output_dict_path = os.path.join(
    #         model.output_dir, f'{model.output_name}_output_dict.pkl'
    #     )
    #     with open(output_dict_path, "wb") as fp:
    #         pickle.dump(out, fp, protocol=pickle.HIGHEST_PROTOCOL)
    #
    #     logger.info(f"Model output written to {output_dict_path}...")


def prep_intermediate_output(out, batch, feature_dict, feature_processor, config_preset, multimer_ri_gap):
    # Prep protein metadata
    template_domain_names = []
    template_chain_index = None
    if feature_processor.config.common.use_templates and "template_domain_names" in feature_dict:
        template_domain_names = [
            t.decode("utf-8") for t in feature_dict["template_domain_names"]
        ]

        # This works because templates are not shuffled during inference
        template_domain_names = template_domain_names[
                                :feature_processor.config.predict.max_templates
                                ]

        if "template_chain_index" in feature_dict:
            template_chain_index = feature_dict["template_chain_index"]
            template_chain_index = template_chain_index[
                                   :feature_processor.config.predict.max_templates
                                   ]

    no_recycling = feature_processor.config.common.max_recycling_iters
    remark = ', '.join([
        f"no_recycling={no_recycling}",
        f"max_templates={feature_processor.config.predict.max_templates}",
        f"config_preset={config_preset}",
    ])

    # For multi-chain FASTAs
    ri = feature_dict["residue_index"]
    chain_index = (ri - numpy.arange(ri.shape[0])) / multimer_ri_gap
    chain_index = chain_index.astype(numpy.int64)
    cur_chain = 0
    prev_chain_max = 0
    for i, c in enumerate(chain_index):
        if c != cur_chain:
            cur_chain = c
            prev_chain_max = i + cur_chain * multimer_ri_gap

        batch["residue_index"][i] -= prev_chain_max

    unrelaxed_protein = protein.from_prediction(
        features=batch,
        result=out,
        b_factors=plddt_b_factors,
        remove_leading_feature_dimension=False,
        remark=remark,
        parents=template_domain_names,
        parents_chain_index=template_chain_index,
    )

    return unrelaxed_protein
