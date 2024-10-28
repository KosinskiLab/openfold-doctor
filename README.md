# OpenfoldDoctor

**OpenfoldDoctor** is an extension of the [OpenFold](https://github.com/aqlaboratory/openfold) project, designed to enhance protein folding simulations with advanced inspection capabilities. With OpenFold Doctor, users can gain deeper insights into the folding process through various export functionalities.

## 🚀 Features

- **Export intermediate protein structures**: Capture and export all intermediate states of protein structures during the folding simulation.
- **Folding progression movie**: Generate a "movie" showing the trajectory of protein folding, providing a visual representation of the entire process.
- **MSA and pair representations**:
  - **MSA**: Visualize Multiple Sequence Alignment (MSA) data through heatmaps.
  - **Pair representation**: Export heatmaps of pairwise interactions at each iteration of the simulation.
- **Attention mechanism visualization**: Export heatmaps showing row-wise and column-wise attention.

## 📦 Installation

To get started with OpenFold Doctor, follow these steps:

### 1. **Clone this OpenFold repository**

Begin by cloning the OpenFold repository from GitHub:

```bash
git clone https://github.com/KosinskiLab/openfold-doctor.git
```

### 2. Checkout the `dr-dev branch`

Move into the cloned OpenFold Doctor directory, and switch to the `dr-dev` branch:

```bash
cd openfold
git checkout dr-dev
```

### 3. Proceed with OpenFold installation
Continue with the [original OpenFold installation process (`pl_upgrades` branch)](https://github.com/aqlaboratory/openfold/blob/pl_upgrades/README.md), however installing the conda environment as in the OpenFold Doctor `environment.yml` file.

Ensure you follow all the steps outlined in the OpenFold installation guide to set up the environment correctly (e.g. using CUDA 12, as per the `pl_upgrades` branch requirements) with the new dependencies introduced by OpenFold Doctor.

## 📈 Usage

This section covers how to utilize all the extensions added to the `pl_upgrades` branch of the official OpenFold repository through OpenfoldDoctor.

### 1. **Exporting intermediate protein structures**

**Description**: Capture and save all intermediate protein structures generated during the inference process.

**How to**:

When executing the folding simulation, use the `--intermediate_structures_export` flag to enable the export of intermediate structures in `.cif` format.

  ```bash
  python run_openfold.py [your usual openfold flags] --intermediate_structures_export
  ```

To generate a **movie of the protein structure evolution** during the inference, you can use the following additional flags:
- `--protein_movie`: enables the movie generation
- `--frame_duration_seconds N`: duration of each frame in seconds (optional; default: 1)
- `--low_res_movie`: exports low resolution png frames for a low resolution movie (optional; default: False)
- `--keep_movie_data`: preserve all intermediate files used in the movie generation process, i.e., png frames and pdb files (optional; default: False)

Example:
  ```bash
  python run_openfold.py [your usual openfold flags] --protein_movie --frame_durations_seconds 1.5 --low_res_movie --keep_movie_data
  ```

`--intermediate_structures_export` can be omitted when using `--protein_movie`. 

### 2. **Exporting MSA and pair representations**

**Description**: Capture the MSA and pair representations before and after they are processed by the evoformer stack. Export them as heatmaps, and optionally generate movies showing their evolution during the inference process.

**How to**:

When executing the folding simulation, use the `--representation_export` flag to enable the export of intermediate structures as png heatmaps. The files will be saved in two separate folders -- `msa` and `pair` -- under your main output folder.

  ```bash
  python run_openfold.py [your usual openfold flags] --representation_export
  ```

To generate a **movie of the MSA and pair representation evolution** during the inference, you can use the following additional flag:
- `--representations_movie`: enables the movie generation

The movies generated for the MSA and pair representation will be saved under the `msa` and `pair` folders, respectively, together with the heatmap png images.

| [<img src="https://github.com/user-attachments/assets/a3f8de36-963d-4434-9a6c-a32935dbfaa4" />](https://github.com/user-attachments/assets/a3f8de36-963d-4434-9a6c-a32935dbfaa4) |
|:--:| 
| 1D3Z (Ubiquitin) MSA (top) and pair (bottom) representation heatmaps before and after the first two recycles. |


### 3. **MSA coverage**

**Description**: Export a png chart of the MSA coverage.

**How to**: Use the `--plot_msa_coverage` flag to save the MSA coverage plot in the `alignment` folder under the main output folder.

  ```bash
  python run_openfold.py [your usual openfold flags] --plot_msa_coverage
  ```

| [<img src="https://github.com/user-attachments/assets/373d8961-5015-414a-aa71-8a883fa783c1" width="300" />](https://github.com/user-attachments/assets/373d8961-5015-414a-aa71-8a883fa783c1) [<img src="https://github.com/user-attachments/assets/cca0af4b-e8c7-4751-8404-f93d412f1b3c" width="300" />](https://github.com/user-attachments/assets/cca0af4b-e8c7-4751-8404-f93d412f1b3c) |
|:--:| 
| CASP14 target T1082 (left) and 1D3Z (Ubiquitin; right) MSA coverage plots |


### 4. **Visualization of attention**

**Description**: Export heatmaps of the row-wise and column-wise attention mechanisms.

**How to**: Use the `--attention_export` flag to save the plots in the `attn` folder under the main output folder.

  ```bash
  python run_openfold.py [your usual openfold flags] --attention_export
  ```

| [<img src="https://github.com/user-attachments/assets/d490e8d8-6c05-4f64-a1c0-bbf5eb6fc0e2" />](https://github.com/user-attachments/assets/d490e8d8-6c05-4f64-a1c0-bbf5eb6fc0e2) |
|:--:| 
| CASP14 target T1082 MSA column-wise attention heatmaps |

### 5. **MSA fasta**

**Description**: Export MSA fasta before and after filtering and embedding.

**How to**: Use the `--msa_fasta_export` flag to save the fasta files in the `alignment` folder under the main output folder.

  ```bash
  python run_openfold.py [your usual openfold flags] --msa_fasta_export
  ```


### 6. **A complete example using all features**

  ```bash
   python run_pretrained_openfold.py \
        $INPUT_FASTA_DIR \
        $TEMPLATE_MMCIF_DIR \
        --output_dir $OUTPUT_DIR \
        --config_preset model_1_ptm \
        --uniref90_database_path $BASE_DATA_DIR/uniref90/uniref90.fasta \
        --mgnify_database_path $BASE_DATA_DIR/mgnify/mgy_clusters_2018_12.fa \
        --pdb70_database_path $BASE_DATA_DIR/pdb70/pdb70 \
        --uniclust30_database_path $BASE_DATA_DIR/uniclust30/uniclust30_2018_08/uniclust30_2018_08 \
        --bfd_database_path $BASE_DATA_DIR/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
        --model_device "cuda:0" \
        --hhblits_binary_path $HHSUITE_BASE_DIR/hhblits \
        --hhsearch_binary_path $HHSUITE_BASE_DIR/hhsearch \
        --kalign_binary_path $KALIGN_PATH \
        --save_outputs \
        --cif_output \
        --skip_relaxation \
        --protein_movie \
        --frame_duration_seconds 0.5 \
        --low_res_movie \
        --keep_movie_data \
        --representation_movie \
        --attention_export \
        --plot_msa_coverage \
        --msa_fasta_export
  ```

## 🔗 Additional Resources
- **OpenFold repository**: https://github.com/aqlaboratory/openfold
- **OpenFold Doctor releases**: https://github.com/KosinskiLab/openfold/releases
- **Issues and support**: If you encounter any issues, feel free to open an issue on the [OpenfoldDoctor GitHub repository](https://github.com/lgiannantoni/openfold/issues).

## 📝 Contributing
Contributions are welcome! Please fork the repository and submit a pull request with your enhancements or bug fixes.
