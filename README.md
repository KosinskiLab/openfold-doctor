# OpenfoldDoctor

**OpenfoldDoctor** is an extension of the [OpenFold](https://github.com/aqlaboratory/openfold) project, designed to enhance protein folding simulations with advanced inspection capabilities. With OpenfoldDoctor, users can gain deeper insights into the folding process through various export functionalities.

## 🚀 Features

- **Export intermediate protein structures**: Capture and export all intermediate states of protein structures during the folding simulation.
- **Folding progression movie**: Generate a "movie" showing the trajectory of protein folding, providing a visual representation of the entire process.
- **MSA and pair representations**:
  - **MSA**: Visualize Multiple Sequence Alignment (MSA) data through heatmaps.
  - **Pair representation**: Export heatmaps of pairwise interactions at each iteration of the simulation.

## 📦 Installation

To get started with OpenfoldDoctor, follow these steps:

### 1. **Clone the original OpenFold repository**

Begin by cloning the original OpenFold repository from GitHub:

```bash
git clone https://github.com/aqlaboratory/openfold.git
```

### 2. Navigate to the OpenFold Directory
Move into the cloned OpenFold directory:

```bash
cd openfold
```

### 3. Download the OpenfoldDoctor patch
Download the `openfold_doctor.patch` file from the OpenfoldDoctor releases. Ensure you save the patch file to a known location on your local machine.

### 4. Checkout the `pl_upgrades branch`
Switch to the `pl_upgrades` branch and reset it to the specific commit `3bec3e9b2d1e8bdb83887899102eff7d42dc2ba9`:

```bash
git checkout pl_upgrades
git reset --hard 3bec3e9b2d1e8bdb83887899102eff7d42dc2ba9
```

### 5. Apply the OpenfoldDoctor patch
Apply the downloaded patch using to incorporate OpenfoldDoctor's enhancements:

```bash
git am --3way /path/to/openfold_doctor.patch
```

**Note**: Replace `/path/to/` with the actual path to your downloaded `openfold_doctor.patch` file.

### 6. Proceed with OpenFold installation
Continue with the [original OpenFold installation process](https://github.com/aqlaboratory/openfold/blob/pl_upgrades/README.md).

Ensure you follow all the steps outlined in the OpenFold installation guide to set up the environment correctly with the new dependencies introduced by OpenfoldDoctor.

## 📈 Usage

This section covers how to utilize all the extensions added to the `pl_upgrades` branch of the official OpenFold repository through OpenfoldDoctor.

### 1. **Exporting intermediate protein structures**

**Description**: Capture and save all intermediate protein structures generated during the inference process.

**How to**:

When executing the folding simulation, use the `--use_doctor` flag to enable the export of intermediate structures.

  ```bash
  python run_openfold.py [your usual openfold flags] --use_doctor
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

The `--protein_movie` flag automatically enables also the `--use_doctor flag`.

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

### 3. **A complete example using all features**

  ```bash
   62 python run_pretrained_openfold.py \
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
        --representation_movie
  ```

## 🔗 Additional Resources
- **OpenFold repository**: https://github.com/aqlaboratory/openfold
- **OpenfoldDoctor releases**: https://github.com/lgiannantoni/openfold/releases
- **Issues and support**: If you encounter any issues, feel free to open an issue on the [OpenfoldDoctor GitHub repository](https://github.com/lgiannantoni/openfold/issues).

## 📝 Contributing
Contributions are welcome! Please fork the repository and submit a pull request with your enhancements or bug fixes.
