import logging
logger = logging.getLogger(__file__)
logger.setLevel(level=logging.DEBUG)

class AttnExporter:
    def __init__(self, model, args, output_dir):
        self.model = model
        self.args = args

        # set callbacks
        [block.msa_att_row.mha.register_forward_hook(self._attn_row_hook) for block in self.model.evoformer.blocks]
        [block.msa_att_col._msa_att.mha.register_forward_hook(self._attn_col_hook) for block in self.model.evoformer.blocks if not block.no_column_attention]

    def _attn_col_hook(self, attn_block, input, output):
        logger.debug("ciao col!")
    

    def _attn_row_hook(self, attn_block, input, output):
        logger.debug("ciao row!")
