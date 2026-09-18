RSAT-NF is a nextflow pipeline to implementation of the workflow by Ksouri et al: https://academic.oup.com/plphys/article/185/3/1242/6094791

The pipeline uses available docker/singularity containers to allow for reproducible execution

# Use example

1.- Download the example input data from this repository, or use your own generated data.
2.- Generate a samplesheet to select the input data
3.- Generate a config file to configurate nextflow run. RSAT has some specific requirements when using singularity, which is already prepared in the the example config file. Use it as template to customize it to your needs
4.- Run!

## 1.- Downloading the data

This repository holds some data generated to be used as an example:

- Pregenerated modules
- Default values for analysis in config file

The reference sequence of the interested organisms, must be downloaded externally. This project is tested using _Prunus Prunica_ reference genome from ensembl: It can be downloaded from "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/prunus_persica/cdna/Prunus_persica.Prunus_persica_NCBIv2.cdna.all.fa.gz" . The sequence and gtf file The pipeline offers a way to automatically download and (re)use 

In config, results variable indicates where the folder to save the generated data.

## 2.- Generate a module samplesheet

When providing pregenerated modules, a samplesheet for them must be generated: an example is privided in the example folder folder as 

## 3.- Generating or updating the config file

Parameters that must be set to run the pipeline are;

## 4.- Run

```
nextflow run main.nf -c nextflow.config
```

## HPC compatibility

The default config comes with slurm profile, aiming to add compatibility with HPC clusters. By default, uses -qos=short to use the short queue of your HPC, which may or may not exist as a your option. Check and update the queue before running the pipeline
