# ctg-sc-cite-seq-10x 
## Nextflow pipeline for processing of 10x chromium cite-seq with rna+(adt/hto) data with cellranger. 

- Analyze 10x cite-seq in one pipeline. 
- Demux -> cellranger -> QC
- Supports ADT/HTO and RNA libraries sequenced on same flowcell.
- Supports different indexing of RNA and ADT/HTO library (e.g. RNA dual and ADT/HTO single). See `Handle dual and single indexing in same sequencing run` for more info.


1. Clone and build the Singularity container for this pipeline: https://github.com/perllb/ctg-sc-cite-seq-10x/tree/master/container
2. Prepare the feature ref csv. See section `Feature reference` below
3. Edit your samplesheet to match the example samplesheet. See section `SampleSheet` below
4. Edit the nextflow.config file to fit your project and system. 
5. Run pipeline 
```
nohup nextflow run pipe-sc-cite-seq-10x.nf > log.pipe-sc-cite-seq-10x.txt &
```

## Input files

1. Samplesheet (see `SampleSheet` section below)
2. Feature reference csv (see `Feature Reference` section below)

### 1. Samplesheet (CTG_SampleSheet.sc-cite-seq-10x.csv):

 | Sample_ID | index | Sample_Project | Sample_Species | Sample_Lib | Sample_Pair | 
 | --- | --- | --- | --- | --- | --- | 
 | Sr1 | SI-GA-D9 | proj_2021_012 | human | rna | 1 |
 | Sr2 | SI-GA-H9 | proj_2021_012 | human | rna | 2 |
 | Sadt1 | SI-GA-C9 | proj_2021_013 | mouse | adt | 1 |
 | Sadt2 | SI-GA-C9 | proj_2021_013 | mouse | adt | 2 |

- The nf-pipeline takes the following Columns from samplesheet to use in channels:

- `Sample_ID` : ID of sample. Sample_ID can only contain a-z, A-Z and "_".  E.g space and hyphen ("-") are not allowed! If 'Sample_Name' is present, it will be ignored. 
- `index` : Must use index ID (10x ID) if dual index. For single index, the index sequence works too.
- `Sample_Project` : Project ID. E.g. 2021_033, 2021_192.
- `Sample_Species` : Only 'human'/'mouse'/'custom' are accepted. If species is not human or mouse, set 'custom'. This custom reference genome has to be specified in the nextflow config file. See below how to edit the config file.
- `Sample_Lib` : 'rna'/'adt'. Specify whether sample is RNA or ADT library. 
- `Sample_Pair` : To match the rna sample with the corresponding adt sample. e.g. in the example above, sample 'Sr1' is the rna library, that should be matched with 'Sadt1' which is the adt library of the sample


### 2. Feature reference (feature.ref.csv)
Csv that declares the molecule structure and unique Feature Barcode sequence of each feature present in your experiment 

See https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis for more info

Example (TotalSeq A): 
| id | name | read | pattern | sequence | feature_type | 
| --- | --- | --- | --- | --- | --- | 
| CD235a | CD235a | R2 | ^(BC) | AGAGTATGTATGGGA | Antibody Capture | 
| CD33 | CD33 | R2 | ^(BC) | TAACTCAGGGCCTAT | Antibody Capture | 
| CD71 | CD71 | R2 | ^(BC) | CCGTGTTCCTCATTA | Antibody Capture | 
| CD11b | CD11b | R2 | ^(BC) | GACAAGTGATCTGCA | Antibody Capture | 
| CD45RA | CD45RA | R2 | ^(BC) | TCAATCCTTCCGCTT | Antibody Capture | 
| CD34 | CD34 | R2 | ^(BC) | GCAGAAATCTCCCTT | Antibody Capture | 
| CD49d | CD49d | R2 | ^(BC) | CCATTCAACTTCCGG | Antibody Capture | 
| CD45 | CD45 | R2 | ^(BC) | TCCCTTGCGATTTAC | Antibody Capture | 


### CSV format templates

#### 1. Samplesheet : `CTG_SampleSheet.sc-cite-seq-10x.csv`
```
Sample_ID,index,Sample_Project,Sample_Species,Sample_Lib,Sample_Pair
Si1,SI-GA-D9,2021_012,human,rna,1
Si2,SI-GA-H9,2021_012,mouse,rna,2
Sample1,SI-GA-C9,2021_013,human,adt,1
Sample2,SI-GA-C9,2021_013,mouse,adt,2
``` 

#### 2. Feature reference : `feature.ref.csv`
```
id,name,read,pattern,sequence,feature_type
CD235a,CD235a,R2,^(BC),AGAGTATGTATGGGA,Antibody Capture
CD33,CD33,R2,^(BC),TAACTCAGGGCCTAT,Antibody Capture
CD71,CD71,R2,^(BC),CCGTGTTCCTCATTA,Antibody Capture
CD11b,CD11b,R2,^(BC),GACAAGTGATCTGCA,Antibody Capture
CD45RA,CD45RA,R2,^(BC),TCAATCCTTCCGCTT,Antibody Capture
CD34,CD34,R2,^(BC),GCAGAAATCTCCCTT,Antibody Capture
CD49d,CD49d,R2,^(BC),CCATTCAACTTCCGG,Antibody Capture
CD45,CD45,R2,^(BC),TCCCTTGCGATTTAC,Antibody Capture

``` 

## Pipeline steps:

Cellranger version: cellranger v6.0 

* `parse samplesheets`: Creates samplesheets (one for RNA, and one for ADT/HTO) for demux based on the input samplesheet. 
* `generate library csv`: Creates library.csv file based on input samplesheet. One .csv per matched RNA and ADT/HTO sample.
* `Demultiplexing` (cellranger mkfastq): Converts raw basecalls to fastq, and demultiplex samples based on index (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/6.0/using/mkfastq). Does this separately for RNA and ADT/HTO (since they often have different index types (dual/single)
* `FastQC`: FastQC calculates quality metrics on raw sequencing reads (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). MultiQC summarizes FastQC reports into one document (https://multiqc.info/).
* `Align` + `Counts` + `Feature Barcoding` (cellranger count): Aligns fastq files to reference genome, counts genes for each cell/barcode, and quantifies ADT/HTO features per barcode - Then performs secondary analysis such as clustering and generates the cloupe files (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis).
* `Aggregation` (cellranger aggr): Automatically creates the input csv pointing to molecule_info.h5 files for each sample to be aggregated and executes aggregation (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/aggregate). 
* `Cellranger count metrics` (bin/ctg-sc-cite-seq-count-metrics-concat.py): Collects main count metrics (#cells and #reads/cell etc.) from each sample and collect in table (**UPDATE**)
* `multiQC`: Compile fastQC and cellranger count metrics in multiqc report
* `md5sum`: md5sum of all generated files


## Handle dual and single indexing in same sequencing run

If your RNA and ADT/HTO libraries have different indexing it can be handled as following:

#### RNA dual - ADT/HTO single
In nextflow.config, set 
```
// bcl2fastq arguments
bcl2fastqarg_rna = "" 
bcl2fastqarg_adt = "--use-bases-mask=Y28n*,I6n*,N10,Y90n*" 
// Index type ('dual' or 'single')
index_rna = "dual"
index_adt = "single"	

```
	
## Container
- `sc-cite-seq-10x`: For 10x sc-cite-seq. Based on cellranger v6.
https://github.com/perllb/ctg-sc-cite-seq-10x/tree/master/container

Build container:
NOTE: Environment.yml file has to be in current working directory
```
sudo -E singularity build sc-cite-seq-10x.sif sc-cite-seq-10x-builder 
```

Add path to .sif in nextflow.config

## Output:
* ctg-PROJ_ID-output
    * `qc`: Quality control output. 
        * cellranger metrics: Main metrics summarising the count / cell output 
        * fastqc output (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
        * multiqc output: Summarizing FastQC output and demultiplexing (https://multiqc.info/)
    * `fastq`: Contains raw fastq files from cellranger mkfastq.
    * `count-cr`: Cellranger count output. Here you find gene/cell count matrices, feature quantification, secondary analysis output, and more. See (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis) for more information on the output files.
    * `summaries`: 
        * web-summary files which provide an overview of essential metrics from the 10x run. 
        * cloupe files which can be used to explore the data interactively in the Loupe browser (https://support.10xgenomics.com/single-cell-gene-expression/software/visualization/latest/what-is-loupe-cell-browser)  
    * `aggregate`:
        * Output from cellranger aggregation. 
    * `ctg-md5.PROJ_ID.txt`: text file with md5sum recursively from output dir root    



## Custom genome 

If custom genome (not hg38 or mm10) is used

1. Set "Sample_Species" column to 'custom' in samplesheet:

Example:
 | Sample_ID | index | Sample_Project | Sample_Species | Sample_Lib | Sample_Pair | 
 | --- | --- | --- | --- | --- | --- | 
 | Sr1 | SI-GA-D9 | proj_2021_012 | **custom** | rna | 1 |
 | Sr2 | SI-GA-H9 | proj_2021_012 | **custom** | rna | 2 |
 | Sadt1 | SI-GA-C9 | proj_2021_013 | **custom** | adt | 1 |
 | Sadt2 | SI-GA-C9 | proj_2021_013 | **custom** | adt | 2 |
 
 2. In nextflow.config, set 
 `custom_genome=/PATH/TO/CUSTOMGENOME`
 
### Add custom genes (e.g. reporters) to cellranger annotation

Use the `ctg-cellranger-add2ref` script. 

https://github.com/perllb/ctg-cellranger-add2ref

