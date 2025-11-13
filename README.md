# ST_SpatialFactor

Code and data for Projected Bayesian Spatial Factor Models


## Prerequisites

- Anaconda or Miniconda
- Julia 1.10 or later
- C compiler (gcc, clang, or MSVC)
- Git

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/yourusername/ST_SpatialFactor.git
cd ST_SpatialFactor
```

### 2. Create Conda Environment

```bash
# Create a new Conda environment
conda create -n spatial_factor python=3.11 -y
conda activate spatial_factor

# Install required Python packages
conda install -c conda-forge numpy pandas matplotlib scipy -y
conda install -c conda-forge scanpy anndata h5py -y
conda install -c bioconda scikit-learn -y

# Install Jupyter and Ipykernel
conda install -c conda-forge jupyter notebook ipykernel -y

# Optional: Install additional spatial transcriptomics packages
pip install squidpy 
pip install GraphST
pip install POT

# Install R with specific version
conda install -c conda-forge r-base=4.4.3 -y


# Install required R packages
conda install -c conda-forge \
    r-Matrix\
    r-ggplot2 \
    r-seurat \
    r-mclust \
    r-fields \
    -y

# For packages not available in conda-forge, install from R
R -e 'install.packages(c(
    "Matrix",
    "GPvecchia",
    "reticulate",
    "ggforce"
), repos="https://cloud.r-project.org/", dependencies = TRUE)'
```

Note: Some R packages might take a while to install due to compilations.

### 3. Install Julia (Skip if already installed)

> **Note:** This project uses Julia 1.10 or later. The latest stable version is recommended.

Choose one of the installation methods below:

#### Option A: Using Homebrew (macOS)
```bash
brew install julia
```

#### Option B: Direct Download (Recommended)
- Visit the [official Julia downloads page](https://julialang.org/downloads/)
- Download the appropriate version for your operating system
- Follow the platform-specific installation instructions

#### Option C: Using Juliaup (Version Manager)
```bash
# For macOS
brew install juliaup
juliaup add release
juliaup default release

# For other platforms, see https://github.com/JuliaLang/juliaup
```

#### Option D: Conda Installation (May not work on all systems)
```bash
conda install -c conda-forge julia -y
```

### 4. Install Julia Packages

```bash
# Create a new Julia project
julia --project -e 'using Pkg; Pkg.activate(".")'

# Install required packages
julia install_packages.jl

# Install all dependencies
julia --project -e 'using Pkg; Pkg.instantiate()'
```

The installation may take a few minutes depending on your internet connection.

### 5. Compile the C Library

```bash
# Compile the C code (serial mode)
julia utils/compile_c_lib.jl

# Optional: Compile with OpenMP for improved performance (if your machine supports OpenMP)
julia utils/compile_c_lib.jl --openmp
```

Note: Use the `--openmp` flag if your machine supports OpenMP and you want to enable parallel processing. 
This can significantly improve performance on multi-core systems by leveraging parallel computation 
for nearest neighbor calculations. If unsure, use the default serial compilation.

### 6. Verify Installation

```bash
# Test the installation
julia --project test_installation.jl
```

### 7. Register Julia kernel with Jupyter

```bash
julia --project -e 'using IJulia; IJulia.installkernel("Julia (ST_SpatialFactor)", "--project=.")'
```

## Project Structure

- `sim1/`: Simulation code and data for Simulation I
  - `projects/`: Jupyter notebooks for data preparation and analysis
    - (Step 1) `datapre.ipynb`: Code for generating data for Section 5
    - (Step 2) `BSLMC_diag_Sigma.ipynb`: Implementation and results for Gibbs+post in Section 5.1
    - (Step 3) `BSLMC_diag_Sigma_Proj.ipynb`: Implementation and results for ProjMC² in Section 5.1
    - (Step 4) `compare.ipynb`: Code for generating Table 3 in Section 5.1
- `sim/`: Simulation code and data for Simulation II and Sensitivity Analysis
  - `projects/`: Jupyter notebooks for data preparation and analysis
    - (Step 1) `datapre.ipynb`: Code for generating data for Simulation II
    - (Step 2) `BSLMC_diag_Sigma.ipynb`: Implementation and results for Gibbs+post for Simulation II
    - (Step 3) `BSLMC_diag_Sigma_Proj.ipynb`: Implementation and results for ProjMC² for Simulation II
    - (Step 4) `compare.ipynb`: Code for generating Table 8 in Supplement S.3
    - (Step 5) `Sensitive_analysis.ipynb`: Sensitivity analysis in Section 5.2
- `RDA/`: Real data analysis code and data
   - `projects/`: Code for data preparation and analysis
      - (Step 1) `readdata.R`: Code for reading data
      - (Step 2) `EDA.R`: Code for data preparation
      - (Step 3) `Kidney_nondisease_ProjMC.ipynb`: Code for ProjMC² for RDA
      - (Optional) `GraphST_test.ipynb`: Code for GraphST comparison, need to install GraphST using `pip install GraphST`
      - (Optional) `STAGATE_test.ipynb`: Code for STAGATE comparison
      - (Step 4) `cluster.R`: Code for summary
- `utils/`: Utility functions and C code
  - `julia-R-nn-ccall2/`: C code for nearest neighbor calculations
  - `util2.j`: Julia utility functions
  - `compile_c_lib.jl`: Script to compile the C library
  - `run_nn.jl`: Script to run the nearest neighbor calculations
- `install_packages.jl`: Script to install required Julia packages
- `test_installation.jl`: Script to test the installation

## Usage

### 1. Activate the Environment

```bash
conda activate spatial_factor
```

### 2. Running the Analysis

#### Simulation Studies
1. Launch Jupyter from your Conda environment:
   ```bash
   jupyter notebook
   ```
2. In the Jupyter interface:
   - Navigate to `sim/projects/datapre.ipynb`
   - Select the "Julia (ST_SpatialFactor)" kernel from the kernel menu
   - This ensures your notebook has access to both your Conda environment and your Julia project
   - Run the notebooks in sequence as numbered in the project structure

#### Real Data Analysis
1. Download the non-diseased kidney data from [10x Genomics](https://www.10xgenomics.com/datasets/human-kidney-preview-data-xenium-human-multi-tissue-and-cancer-panel-1-standard)

2. Run the code in `./projects` in the following order:
   - `readdata.R`
   - `EDA.R`
   - `Kidney_nondisease_ProjMC.ipynb`
   - `cluster.R`

3. For comparison with STAGATE:
   Since STAGATE_pyG has compatibility issues with Python 3.11, you'll need to create a separate environment. Follow these steps:

   ```bash
   # Create and activate a new environment with Python 3.9
   conda create -n stagate_env python=3.9 -y
   conda activate stagate_env

   # Install PyTorch and related packages
   pip install requests
   pip install torch torchvision torchaudio
   pip install torch-geometric

   # Install PyTorch Geometric dependencies
   pip install torch-sparse torch-scatter torch-cluster torch-spline-conv -f https://data.pyg.org/whl/torch-2.1.0+cpu.html

   # Install additional required packages
   conda install -c conda-forge scanpy pytorch -y
   conda install -c conda-forge jupyter notebook ipykernel -y
   ```

   After setting up the environment, follow the complete installation instructions at [STAGATE Documentation](https://stagate.readthedocs.io/en/latest/Installation_pyG.html) to ensure all dependencies are properly installed.

4. Optional Comparisons:
   - For GraphST comparison: Run `GraphST_test.ipynb`
   - For STAGATE comparison: Run `STAGATE_test.ipynb`

