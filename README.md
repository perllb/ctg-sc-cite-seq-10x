# ctg-sc-cite-seq-10x 
## Nextflow pipeline for processing of 10x cite-seq. RNA+(ADT/HTO) data with cellranger. 

- Run one project at a time only
- Demux -> cellranger -> QC -> delivery
- Supports ADT/HTO and RNA libraries sequenced on same flowcell.
- Supports different indexing of RNA and ADT/HTO library (e.g. RNA dual 10bp and ADT/HTO single 6bp). See `Handle dual and single indexing in same sequencing run` for more info.
- TotalSeq feature IDs used for ADT/HTO in experiment can be specified in samplesheet. Pipeline creates feature reference csv from IDs. See `2. Feature reference (feature.ref.csv)` and `1. Samplesheet` sections below.
- If antibody features are not in reference, a custom feature-ref csv can be created and specified with drivers -f tag.
- Currently, BioLegend TotalSeq ADT A,B and C as well as TotalSeq HTO A are in references, located at `/projects/fs1/shared/references/cellranger_totalseq`.


1. Clone and build the Singularity container for this pipeline: https://github.com/perllb/ctg-sc-cite-seq-10x/tree/master/container
2. If not using reference totalSeq antibodies, prepare the feature ref csv. See section `Feature reference` below
3. Edit your samplesheet to match the example samplesheet. See section `SampleSheet` below
4. Edit the nextflow.config file to fit your project and system. 
5. Run pipeline 
```
nohup nextflow run pipe-sc-cite-seq-10x.nf > log.pipe-sc-cite-seq-10x.txt &
```
## Driver
- Driver for pipeline located in `.../shared/ctg-pipelines/ctg-sc-cite-seq-10x/v2.0/sc-cite-seq-10x-driver`

## Input files

The following files must be in the runfolder to start pipeline successfully.

1. Samplesheet  (see `SampleSheet` section below)
2. (OPTIONAL): Feature reference csv (see `Feature Reference` section below) if not using reference features on lsens4 `/projects/fs1/shared/references/cellranger_totalseq`.

### 1. Samplesheet (CTG_SampleSheet.sc-cite-seq-10x.csv):

- Can have other names than CTG_SampleSheet.**sc-cite-seq-10x**.csv - then specify which sheet when starting driver: e.g. `sc-cite-seq-10x-driver -s CTG_SampleSheet.2022_102.csv`

#### Example sheet
```
[Header]
metaid,2021_067_citeseqTest
antibodies,"ADT_A0574,ADT_A0052,ADT_A0394,ADT_A0161,ADT_A0063,ADT_A0576,ADT_A0054,ADT_A0048"
email,per.a@med.lu.se
autodeliver,y
[Data]
Lane,Sample_ID,index,Sample_Species,Sample_Project,Sample_Lib,Sample_Pair
,EFS_21_022,SI-TT-D5,human,2021_067,rna,1
,con_21_023,SI-TT-E5,human,2021_067,rna,2
,rnaEFS_ADT,ACAGTG,human,2021_067,adt,1
,con_ADT,TGACCA,human,2021_067,adt,2
```
#### [Header] section
Two optional entries. They can both be skipped, but recommended to use both:
- `email` : Email to customer (or ctg staff) that should retrieve email with qc and deliver info upon completion of pipeline. Note: only lu emails works (e.g. @med.lu.se or @lth.se.
- `autodeliver` : set to `y` if email should be sent automatically upon completion. Otherwise, set to `n`.
- `metaid` : optional. set to create working folders with this name. otherwise, it will be created based on runfolder name/date.
- `antibodies` : optional. define antibodies used in experiment (across all samples). The sc-cite-seq-10x-driver will use these IDs to extract all info from cellranger_totalseq references, and create the "features" csv file needed in count analysis.
  - recommended if using antibodies that are defined in the totalSeq human cocktail csv files (/projects/fs1/shared/references/cellranger_totalseq). 
  	- Set IDs of antibodies that were used in experiment (and will be used to create the feature reference file for `count` analysis.)
  	- IMPORTANT: They MUST match the IDs of the totalSeq human cocktail references on lsens 
	- The list of antibodies should be comma-separated, and best if quoted (should also work without quote, but not tested.
  - Alternative is to
  	- create a feature.ref.csv file, and add it to runfolder. It must have the standard required cellranger format (see: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis). Note: The sc-cite-seq-10x-driver will first look for a `feature.ref.csv` file in runfolder. If not found, it will look for `antibodies,` field in header section. 
  	- The last alternative is to run the driver with -f flag, pointing to any file that serve as feature-ref csv. e.g. `sc-cite-seq-10x-driver -f /path/to/feature.ref.csv`

#### [Data] section

 | Sample_ID | index | Sample_Species | Sample_Project | Sample_Lib | Sample_Pair | 
 | --- | --- | --- | --- | --- | --- | 
 | Sr1 | SI-GA-D9 | human | 2022_022 | rna | 1 | 
 | Sr2 | SI-GA-H9 | mouse | 2022_022 | rna | 2 | 
 | Sadt1 | SI-GA-C9 | human | 2022_022 | adt | 1 | 
 | Sadt2 | SI-GA-C9 | mouse | 2022_022 | adt | 2 | 

- The nf-pipeline takes the following Columns from samplesheet to use in channels:

- `Sample_ID` : ID of sample. Sample_ID can only contain a-z, A-Z and "_".  E.g space and hyphen ("-") are not allowed! If 'Sample_Name' is present, it will be ignored. 
- `index` : Must use index ID (10x ID) if dual index. For single index, the index sequence works too.
- `Sample_Project` : Project ID. E.g. 2021_033, 2021_192.
- `Sample_Species` : Only 'human'/'mouse'/'custom' are accepted. If species is not human or mouse, set 'custom'. This custom reference genome has to be specified in the nextflow config file. See below how to edit the config file.
- `Sample_Lib` : 'rna'/'adt'. Specify whether sample is RNA or ADT library. Note - even if it is and HTO experiment, use `adt` for HTO.
- `Sample_Pair` : To match the rna sample with the corresponding adt sample. e.g. in the example above, sample 'Sr1' is the rna library, that should be matched with 'Sadt1' which is the adt library of the sample

### Samplesheet template

#### Samplesheet name : `CTG_SampleSheet.sc-cite-seq-10x.csv`
```
[Header]
email,per.a@med.lu.se
autodeliver,y
metaid,2021_067_citeseqTest
antibodies,"ADT_A0574,ADT_A0052,ADT_A0394,ADT_A0161,ADT_A0063,ADT_A0576,ADT_A0054,ADT_A0048"
[Data]
Lane,Sample_ID,index,Sample_Species,Sample_Project,Sample_Lib,Sample_Pair,email,autodeliver
,EFS_21_022,SI-TT-D5,human,2021_067,rna,1
,con_21_023,SI-TT-E5,human,2021_067,rna,2
,rnaEFS_ADT,ACAGTG,human,2021_067,adt
,con_ADT,TGACCA,human,2021_067,adt,2
```

### 2. Feature reference (feature.ref.csv)

#### NB: 
- Recommended is to use the `antibodies,` declarartion in the header: see examples above. 
- The references are found in `/projects/fs1/shared/references/cellranger_totalseq`. 
- The IDs must match "id" fields in these references
- The driver will use these ids to extract the information needed for the feature.ref.csv file 
- The driver will create the feature ref file and add it to nextflow.config to use in pipeline.

#### If not using antibodies, list in header
- Csv that declares the molecule structure and unique Feature Barcode sequence of each feature present in your experiment. 
- This should either be added in runfolder (must be named /path/to/runfolder/**feature.ref.csv**).
- Or it could be specified by driver with -f flag: sc-cite-seq-10x-driver -f /path/to/feature.ref.csv. 


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


### Feature reference template 
#### Filename : `feature.ref.csv`
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

RNA and ADT/HTO libraries must often have different indexing. It is handled by the pipeline by:
- It looks up the length of `adt` index, and setting the --use-bases-mask accordingly. If adt index is found to be 6 bases, it will set --use-bases-mask=Y28n*,I6n*,N10,Y90n* during mkfastq_adt. 
- By default, it will assume that RNA sample indices are dual, and ADT indices are single. It will thus set --filter-single-index during mkfastq_adt, and --filter-dual-index during mkfastq_rna. 
	
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

