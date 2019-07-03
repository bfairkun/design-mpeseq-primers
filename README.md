# design-mpeseq-primers
For designing RT primers with blast search, tm constraints, etc
See `dag.svg` and `config.yaml` to get idea of the basic workflow of primer design.

## Setup

To ensure all contributors are using the same computational environment, we use
[conda][] to manage software dependencies (made possible by the [bioconda][] and
[conda-forge][] projects). Please complete the following steps to replicate the
computing environment. Note that this is only guaranteed to work on a Linux-64
based architecture, but in theory should be able to work on macOS as well. All
commands shown below are intended to be run in a Bash shell from the root of the
project directory.

1. Install Git and register for an account on GitHub

1. Download and install Miniconda ([instructions](https://conda.io/miniconda.html))

1. Clone this repository (or your personal fork) using `git clone`
    ```
    git clone https://github.com/bfairkun/design-mpeseq-primers.git
    ```

1. Create the conda environment "my_PrimerDesign_env" using `environment.yaml`
    ```
    conda env create --file environment.yaml python=3.6
    ```

1. To use the conda environment, you must first activate it by running `source
activate my_PrimerDesign_env`. This will override your default settings for Python, and
various other software packages. When you are done working on this project, you
can either logout of the current session or deactivate the environment by
running `source deactivate`.

1. edit config and set parameters as needed
```
vim config.yaml
```
1. execute workflow.  Use `snakemake -n` to do a dry-run test. If using cluster, edit
`cluster-config.yaml` and execute the snakemake process and all snakemake-derived
jobs on cluster via the included sbatch script via `sbatch snakemake.sbatch`.
```
#Execute snakemake
snakemake
```

[bioconda]: https://bioconda.github.io
[conda]: https://conda.io/docs/
[conda-forge]: https://conda-forge.org/
[workflowr]: https://github.com/jdblischak/workflowr
