# GroupProject SequenceBioinformatics

This project aims to replicate the results of:

**Han, H., Wang, Z. & Zhu, S.**  
*Benchmarking metagenomic binning tools on real datasets across sequencing platforms and binning modes.*  
**Nat Commun 16, 2865 (2025).** [https://doi.org/10.1038/s41467-025-57957-6](https://doi.org/10.1038/s41467-025-57957-6)

A tool wrapper is made available by the authors at: [https://github.com/htaohan/databinning](https://github.com/htaohan/databinning)

---

## Installation

After cloning this repository (or after updating tools), run the following in a Bash shell:

source scripts/environment/install_and_activate.sh

This will install all dependencies in the conda environment and activates it.

If you install, update or deinstall tools run:

bash scripts/environment/save_environment.sh

---

## Downloading data

To download the .sra files run:

bash scripts/00_download.sh

This will fetch all .sra files (202, approx. 2 Tb) from the accesion codes in the data/subdir/*.txt files into the respective data/subdir folders. Rerunning this script only downloads files, if they are missing or incomplete.