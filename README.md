# ctg-sc-adt-rna-10x 
## Nextflow pipeline for preprocessing of 10x chromium cite-seq with rna+adt data with cellranger. 

- Only needs one samplesheet with description of each sample, and the feature.csv list with antibody tags  
- Supports mm10 and hg38 references, but can also be run with custom reference genome and annotation (must be added via nextflow.config). See custom genome below.
- Must be run project-wise (one project (rna+adt) per run)

## USAGE

1. For reproducibility, create a project / run directory in which you clone this repository, containing the nextflow.config, nf-pipeline, bin
2. Clone and build the Singularity container for this pipeline: https://github.com/perllb/ctg-sc-adt-rna-10x/tree/master/container/sc-adt-rna-10x.v6
3. Edit the nextflow.config file to fit your project and system. Set current directory as `basedir`. (from where you execute the pipeline).
4. Edit your samplesheet to match the example samplesheet. See section `SampleSheet` below
5. Run pipeline 
```
nohup nextflow run pipe-sc-adt-rna-10x.nf > log.pipe-sc-adt-rna-10x.txt &
```


## Pipeline steps:

Cellranger version: cellranger v6.0 

* `parse samplesheets`: Creates samplesheets (one for RNA, and one for ADT) for demux based on the input samplesheet. 
* `generate library csv`: Creates library.csv file based on input samplesheet. One .csv per matched RNA and ADT sample.
* `Demultiplexing` (cellranger mkfastq): Converts raw basecalls to fastq, and demultiplex samples based on index (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/6.0/using/mkfastq). Does this separately for RNA and ADT (since they often have different index types (dual/single)
* `FastQC`: FastQC calculates quality metrics on raw sequencing reads (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). MultiQC summarizes FastQC reports into one document (https://multiqc.info/).
* `Align` + `Counts` + `Feature Barcoding` (cellranger count): Aligns fastq files to reference genome, counts genes for each cell/barcode, and quantifies ADT features per barcode - Then performs secondary analysis such as clustering and generates the cloupe files (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis).
* `Aggregation` (cellranger aggr): Automatically creates the input csv pointing to molecule_info.h5 files for each sample to be aggregated and executes aggregation (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/aggregate). 
* `Cellranger count metrics` (bin/ctg-sc-adt-rna-count-metrics-concat.py): Collects main count metrics (#cells and #reads/cell etc.) from each sample and collect in table (**UPDATE**)
* `multiQC`: Compile fastQC and cellranger count metrics in multiqc report
* `md5sum`: md5sum of all generated files


## Output:
* ctg-PROJ_ID-output
    * `qc`: Quality control output. 
        * cellranger metrics: Main metrics summarising the count / cell output 
        * fastqc output (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
        * multiqc output: Summarizing FastQC output and demultiplexing (https://multiqc.info/)
    * `fastq`: Contains raw fastq files from cellranger mkfastq.
    * `count-cr`: Cellranger count output. Here you find gene/cell count matrices, secondary analysis output, and more. Please go to (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis) for more information on the output files.
    * `summaries`: 
        * web-summary files which provide an overview of essential metrics from the 10x run. 
        * cloupe files which can be used to explore the data interactively in the Loupe browser (https://support.10xgenomics.com/single-cell-gene-expression/software/visualization/latest/what-is-loupe-cell-browser)  
    * `aggregate`:
        * Output from cellranger aggregation. 
    * `ctg-md5.PROJ_ID.txt`: text file with md5sum recursively from output dir root    


## Samplesheet requirements:

Note: no header! only the rows shown below, starting with the column names.

 | Sample_ID | index | Sample_Project | Sample_Species | Sample_Lib | Sample_Pair | 
 | --- | --- | --- | --- | --- | --- | --- | --- |
 | Sr1 | SI-GA-D9 | proj_2021_012 | human | rna | 1 |
 | Sr2 | SI-GA-H9 | proj_2021_012 | human | rna | 2 |
 | Sadt1 | SI-GA-C9 | proj_2021_013 | mouse | adt | 1 |
 | Sadt2 | SI-GA-C9 | proj_2021_013 | mouse | adt | 2 |

- The nf-pipeline takes the following Columns from samplesheet to use in channels:

- Sample_ID : 'Sample_Name' will be ignored
- index : Must use index ID if dual. For single, index sequence works too.
- Sample_Project : Project ID
- Sample_Species : 'human'/'mouse'/'custom' - if custom, see below how to edit the config file
- Sample_Lib : 'rna'/'adt'
- Sample_Pair : To match the rna sample with the corresponding adt sample. e.g. in the example above, sample 'Sr1' is the rna library, that should be matched with 'Sadt1' which is the adt library of the sample

## Container
- `sc-adt-rna-10x`: For 10x sc-adt-rna-seq. Based on cellranger v6.
https://github.com/perllb/ctg-sc-adt-rna-10x/tree/master/container/sc-adt-rna-10x.v6

Build container:
NOTE: Environment.yml file has to be in current working directory
```
sudo -E singularity build sc-adt-rna-10x.v6.sif sc-adt-rna-10x.v6-build
```

Add path to .sif in nextflow.config

## Custom genome 

If custom genome (not hg38 or mm10) is used

1. Set "Sample_Species" column to 'custom' in samplesheet:

Example:
 | Sample_ID | index | Sample_Project | Sample_Species | Sample_Lib | Sample_Pair | 
 | --- | --- | --- | --- | --- | --- | --- | --- |
 | Sr1 | SI-GA-D9 | proj_2021_012 | **custom** | rna | 1 |
 | Sr2 | SI-GA-H9 | proj_2021_012 | **custom** | rna | 2 |
 | Sadt1 | SI-GA-C9 | proj_2021_013 | **custom** | adt | 1 |
 | Sadt2 | SI-GA-C9 | proj_2021_013 | **custom** | adt | 2 |
 
 2. In nextflow.config, set 
 `custom_genome=/PATH/TO/CUSTOMGENOME`
 
### Add custom genes (e.g. reporters) to cellranger annotation

Use the `ctg-cellranger-add2ref` script. 

https://github.com/perllb/ctg-cellranger-add2ref


