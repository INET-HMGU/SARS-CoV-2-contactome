# GWAS analysis of the covid interactome

In this project a protein-protein interaction map of all covid proteins with human host proteins was generated.

The goal of this project is to use GWAS data in combination PPI networks to understand the comorbidities of COVID-19.

For this we make use of the magma gene set enrichment workflow and apply it to the GTEx GWAS data of 114 studies.


## Install magma
```
mkdir -p packages/magma
cd packages/magma
wget https://ctg.cncr.nl/software/MAGMA/prog/magma_v1.07b_static.zip
unzip magma_v1.07b_static.zip 
git clone https://github.com/Kyoko-wtnb/FUMA_scRNA_data.git
cd ../..
```

Dowload the reference data
```
cd data/current
mkdir magma
cd magma
wget https://ctg.cncr.nl/software/MAGMA/aux_files/NCBI37.3.zip
unzip NCBI37.3.zip 
wget https://ctg.cncr.nl/software/MAGMA/aux_files/NCBI38.zip
unzip NCBI38.zip
wget https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip
unzip g1000_eur.zip
cd ../../../..
```

## Get the GWAS data
```
cd data/current
mkdir gtex_gwas_data
cd gtex_gwas_data
wget -O harmonized_imputed_gwas.tar https://zenodo.org/record/3629742/files/harmonized_imputed_gwas.tar?download=1
wget -O gwas_metadata.txt  https://zenodo.org/record/3629742/files/gwas_metadata.txt?download=1
tar -xf harmonized_imputed_gwas.tar
cd ../../../..
```

## Run the analysis

Knit the file 'R/magma.Rmd' to generate a report and individual results files in 'results/current/'.
