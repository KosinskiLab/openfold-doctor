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
After installation, you can use OpenfoldDoctor's features to inspect and visualize the protein folding process. Refer to the documentation and examples provided in the repository for detailed usage instructions.

## 🔗 Additional Resources
- **OpenFold repository**: https://github.com/aqlaboratory/openfold
- **OpenfoldDoctor releases**: https://github.com/lgiannantoni/openfold/releases
- **Issues and support**: If you encounter any issues, feel free to open an issue on the [OpenfoldDoctor GitHub repository](https://github.com/lgiannantoni/openfold/issues).

## 📝 Contributing
Contributions are welcome! Please fork the repository and submit a pull request with your enhancements or bug fixes.
