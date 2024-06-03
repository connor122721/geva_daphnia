# Estimate the Age of Variants using Phased-SNP Data from *Daphnia pulex* genomes

This repository provides scripts to estimate the age of genetic variants using phased single nucleotide polymorphism (SNP) data. It leverages the GEVA (Genealogical Estimation of Variant Age) software developed by Peter K. Albers and colleagues. GEVA is a powerful tool that uses coalescent-based methods to estimate the age of variants from genotype data, useful for VCF formatted data.

## Installation

To get started, clone the Geva repository and make the executable:

```bash
git clone https://github.com/pkalbers/geva
cd geva
make .
```

# References
Albers, P. K., & McVean, G. (2020). Dating genomic variants and shared ancestry in population-scale sequencing data. *PLoS biology*, 18(1), e3000586. [https://doi.org/10.1371/journal.pbio.3000586]
