ATAC-seq-analysis using ENCODE atac-seq-pipeline


@title: atac-seq-analysis using ENCODE atac-seq-pipeline from Kundaje lab

@date: 2018-01-22

@author: Xinlong

@update: 2019-01-10

## INSTALL atac-seq-pipeline

source: <https://github.com/kundajelab/atac-seq-pipeline/blob/master/docs/tutorial_local_conda.md>

**Note: make sure you have already installed *Java8* and *Conda* in your working environment.**

### 1. Git clone this pipeline

```bash
$ git clone https://github.com/ENCODE-DCC/atac-seq-pipeline
$ cd atac-seq-pipeline
```

### 2. Download `cromwell`

`cromwell` website: [cromwell](https://github.com/broadinstitute/cromwell)

```bash
$ wget https://github.com/broadinstitute/cromwell/releases/download/34/cromwell-34.jar
$ chmod +rx cromwell-34.jar
```

### 3. Install `Conda` dependencies

It will create two conda environments `CONDA_ENV`和`CONDA_ENV_PY3`.

```bash
$ bash uninstall_dependencies.sh  # to remove any old existing pipeline env
$ bash install_dependencies.sh
```

> ##### How to install atac-seq-pipeline conda envs under sepcific path:
>
> * Edit `.condarc` file under `$HOME`:
>
>   ```bash
>   channels:
>     - bioconda
>     - conda-forge
>     - defaults
>   
>   # add below code for indicating where conda new envs will be installed
>   envs_dirs:
>     - /mnt/nfs/data/Kian_Lab/Xinlong/install/miniconda3/envs/
> 
>   ```
>
> 
>
> * Modify `install_dependencies.sh` file:
>
>   用vim打开`install_dependencies.sh`，在下面的环境变量`CONDA_ENV`和`CONDA_ENV_PY3`后添加安装路径信息：
>
>   ```bash
>   CONDA_ENV=encode-atac-seq-pipeline
>   CONDA_ENV_PY3=encode-atac-seq-pipeline-python3
> 
>   # add below code for indicating dependency installation directory
>   # manually set installation directory
>   INSTALL_DIR=/your/install/directory/${CODNA_ENV}
>   INSTALL_DIR_PY3=/your/install/directory/${CONDA_ENV_PY3}
>   ```
>
> 

### 4. Download pre-built genome or build customized genome

#### 1) Download pre-built genome database for mm10

```bash
$ wget https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/@@download/mm10_no_alt_analysis_set_ENCODE.fasta.gz
$ tar xvf test_genome_database_hg38_atac.tar
```

Other genome download path is listed in the following section: INPUT JSON preparation --- Reference genome

#### 2) Build genome database

source: <https://github.com/kundajelab/atac-seq-pipeline/blob/master/docs/build_genome_database.md>

**After installing `Conda` and `Conda dependencies` from above steps**, choose `[GENOME]` from `hg19`, `hg38`, `mm9` and `mm10` and specify a destination directory. This will take several hours. We recommend not to run this installer on a login node of your cluster. It will take >8GB memory and >2h time.

```bash
$ bash conda/install_genome_data.sh [GENOME] [DESTINATION_DIR]
```

Find a **TSV file** on the destination directory and use it for `"atac.genome_tsv"` in  input JSON

#### 3) build genome database for your own genome:

Follow the manual from this [LINK](https://github.com/kundajelab/atac-seq-pipeline/blob/master/docs/build_genome_database.md).

##  Run the automated pipeline

### 1. Activate (or deactivate) working environment

```bash
#To activate this environment, use:
$ source activate bds_atac
#To deactivate an active environment, use:
$ conda deactivate
```

### 2. Run ATAC-seq pipeline

```bash
$ source activate encode-atac-seq-pipeline # IMPORTANT!
$ INPUT=examples.json # check the configuration of INPUT JSON file in the following section!
$ java -jar -Dconfig.file=backends/backend.conf cromwell-34.jar run atac.wdl -i ${INPUT}
```

> It will take about an hour. You will be able to find all outputs on `cromwell-executions/atac/[RANDOM_HASH_STRING]/`. See  output details in following sections.

## Input JSON preparation

source: <https://github.com/kundajelab/atac-seq-pipeline/blob/master/docs/input.md>

An input JSON file includes all input parameters and metadata for running pipelines:

1) Reference genome (hg38, mm10, hg19, ...) and genome specific parameters (indices, ...).
2) Input data file paths/URIs (FASTQs, BAMs, TAG-ALIGNs, ...).
3) Pipeline parameters.
4) Resource for instances/jobs.

### For DNANexus CLI users

dxWDL (DNANexus CLI for WDL) does not support definition of task level variables with a prefix `atac.` in an input JSON file. Therefore, `atac.[TASK_NAME].[VAR_NAME]` should be replaced with `[TASK_NAME].[VAR_NAME]`. Simply remove a prefix `atac.` for task level variables. BUT DO NOT REMOVE it for workflow level variables. For example, `atac.qc_report.name` is a task (task `qc_report` in a workflow `atac`) level variable so it should be replaced with `qc_report.name`. But `atac.genome_tsv` is a workflow (`atac`) level variable, so you need to keep it the same. This is the only difference between DNANexus CLI and other platforms.

### Reference genome

We currently support 4 genomes. You can also [build a genome database for your own genome](build_genome_database.md).

| genome | source | built from                                                   |
| ------ | ------ | ------------------------------------------------------------ |
| hg38   | ENCODE | [GRCh38_no_alt_analysis_set_GCA_000001405](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz) |
| mm10   | ENCODE | [mm10_no_alt_analysis_set_ENCODE](https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/@@download/mm10_no_alt_analysis_set_ENCODE.fasta.gz) |
| hg19   | UCSC   | [GRCh37/hg19](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/referenceSequences/male.hg19.fa.gz) |
| mm9    | UCSC   | [mm9, NCBI Build 37](<http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.2bit>) |

Choose one TSV file for `"atac.genome_tsv"` in your input JSON. `[GENOME]` should be `hg38`, `mm10`, `hg19` or `mm9`.

| platform              | path/URI                                                     |
| --------------------- | ------------------------------------------------------------ |
| Google Cloud Platform | `gs://encode-pipeline-genome-data/[GENOME]_google.tsv`       |
| DNANexus (CLI)        | `dx://project-BKpvFg00VBPV975PgJ6Q03v6:data/pipeline-genome-data/[GENOME]_dx.tsv` |
| DNANExus (Web)        | Choose `[GENOME]_dx.tsv` from [here](https://platform.dnanexus.com/projects/BKpvFg00VBPV975PgJ6Q03v6/data/pipeline-genome-data) |
| Stanford Sherlock     | `genome/scg/[GENOME]_scg.tsv`                                |
| Stanford SCG          | `genome/sherlock/[GENOME]_sherlock.tsv`                      |
| Local/SLURM/SGE       | You need to [build a genome database](build_genome_database.md). |

### Input data file

Choose any data type (FASTQ, BAM, nodup/filtered BAM, TAG-ALIGN and PEAK) you want and DO NOT define arrays for other types. For FASTQs and their corresponding adapter arrays, we provide two ways to define them since DNANexus web UI supports up to an 1-dim array. Choose between 3-dim `fastqs` or 1-dim `fastqs_rep[REP_ID]_R[READ_END_ID]` according to your preference. The pipeline supports up to 6 replicates.

- `"atac.fastqs"` : 3-dimensional array with FASTQ file path/URI.

  - 1st dimension: replicate ID: **最外层数（第一层数组）组是样本重复数**

  - 2nd dimension: merge ID (this dimension will be reduced after merging FASTQs)： **第二层数组是之后会merge到一起的数据**

  - 3rd dimension: endedness ID (0 for SE and 0,1 for PE): **最内层数组（第三层）是PE或者SE, 其中SE是一维数组，PE是二维数组**

    ```bash
    array[ [ [[MERGE1_R1],               # PE: R1          
              [MERGE1_R2]],              # PE: R2 # merge 1 for 2nd dimension             
              [[MERGE2_R1], [MERGE2_R2]]          # merge 2 for 2nd dimension
           ],                                           # rep1 for 1st demension
 
           [],                                          # rep2 for 1st demension
 
           []                                          # rep3 for 1st demension
    ]
 
    ```

 

- `"atac.fastqs_rep1_R1"` : Array of FASTQ file to be merged for rep1-R1.

- `"atac.fastqs_rep1_R2"` : Array of FASTQ file to be merged for rep1-R2. Do not define if your FASTQ is single ended.

- `"atac.fastqs_rep2_R1"` : Array of FASTQ file to be merged for rep2-R1. Do not define if you don't have replicate 2.

- `"atac.fastqs_rep2_R2"` : Array of FASTQ file to be merged for rep2-R2. Do not define if you don't have replicate 2.

- `"atac.fastqs_rep3_R1"` : Array of FASTQ file to be merged for rep3-R1. Do not define if you don't have replicate 3.

- `"atac.fastqs_rep3_R2"` : Array of FASTQ file to be merged for rep3-R2. Do not define if you don't have replicate 3.

- `"atac.fastqs_rep4_R1"` : Array of FASTQ file to be merged for rep4-R1. Do not define if you don't have replicate 4.

- `"atac.fastqs_rep4_R2"` : Array of FASTQ file to be merged for rep4-R2. Do not define if you don't have replicate 4.

- `"atac.bams"` : Array of raw (unfiltered) BAM file path/URI.

  - 1st dimension: replicate ID

- `"atac.nodup_bams"` : Array of filtered (deduped) BAM file path/URI.

  - 1st dimension: replicate ID

- `"atac.tas"` : Array of TAG-ALIGN file path/URI.

  - 1st dimension: replicate ID

- `"atac.peaks"` : Array of NARROWPEAK file path/URI.

  - 1st dimension: replicate ID

- `"atac.peaks_pr1"` : Array of NARROWPEAK file path/URI for 1st self pseudo replicate of replicate ID.

  - 1st dimension: replicate ID

- `"atac.peaks_pr2"` : Array of NARROWPEAK file path/URI for 2nd self pseudo replicate of replicate ID.

  - 1st dimension: replicate ID

- `"atac.peak_ppr1"` : NARROWPEAK file path/URI for pooled 1st pseudo replicates.

- `"atac.peak_ppr2"` : NARROWPEAK file path/URI for pooled 2nd pseudo replicates.

- `"atac.peak_pooled"` : NARROWPEAK file path/URI for pooled replicate.

If starting from peaks then always define `"atac.peaks"`. Define `"atac.peaks_pr1"`, `"atac.peaks_pr2"`, `"atac.peak_pooled"`, `"atac.peak_ppr1"` and `"atac.peak_ppr2"` according to the following rules:

```
if num_rep>1:
    if true_rep_only: peak_pooled,
    else: peaks_pr1[], peaks_pr2[], peak_pooled, peak_ppr1, peak_ppr2
else:
    if true_rep_only: "not the case!"
    else: peaks_pr1[], peaks_pr2[]
```

### Pipeline parameters

1. General

   Pipeline type (ATAC-Seq or DNase-Seq) : The only difference between two types is TN5 shifting for TAG-ALIGN outputs.

   - `"atac.pipeline_type` : `atac` for ATAC-Seq. `dnase` for DNase-Seq.

   Input data endedness.

   - `"atac.paired_end"` : Set it as `true` if input data are paired end, otherwise `false`.

   Other optional settings.

   - `"atac.align_only"` : (optional) Disable all downstream analysis (peak calling, ...) after mapping.
   - `"atac.multimapping"` : (optional) Multimapping reads.
   - `"atac.true_rep_only"` : (optional) Set it as `true` to disable all analyses (including IDR, naive-overlap and reproducibility QC) related to pseudo replicates. This flag suppresses `"atac.enable_idr"`.
   - `"atac.disable_xcor` : (optional) Disable cross-correlation analysis.
   - `"atac.qc_report.name"` : (optional) Name of sample.
   - `"atac.qc_report.desc"` : (optional) Description for sample.

2. Adapter trimmer settings

   Structure/dimension of `"atac.adapters` must match with that of `"atac.fastqs"`. If no adapters are given then do not define `"atac.adapters"` in `input.json`. If some adapters are known then define them in `"atac.adapters"` and leave other entries empty (`""`) while keeping the same structure/dimension as in `"atac.fastqs"`. All undefined/non-empty adapters will be trimmed without auto detection.

   - `"atac.trim_adapter.auto_detect_adapter"` : (optional) Set it as `true` to automatically detect/trim adapters for empty entries in `"atac.adapters"`. There will be no auto detection for non-empty entries it. If `"atac.adapters"` is not defined then all adapters will be detected/trimmed for all fastqs.
   - `"atac.trim_adapter.min_trim_len"` : (optional) Minimum trim length for `cutadapt -m`.
   - `"atac.trim_adapter.err_rate"` : (optional) Maximum allowed adapter error rate for `cutadapt -e`.

3. Bowtie2 settings (remove a prefix `atac.` for DNANexus CLI).

   - `"atac.bowtie2.score_min"` : (optional) Min. acceptable alignment score function w.r.t read length.

4. Filter/dedup (post-alignment) settings (remove a prefix `atac.` for DNANexus CLI).

   - `"atac.filter.dup_marker"` : (optional) Dup marker. Choose between `picard` (default) and `sambamba`.
   - `"atac.filter.mapq_thresh"` : (optional) Threshold for low MAPQ reads removal.
   - `"atac.filter.no_dup_removal"` : (optional) No dup reads removal when filtering BAM.

5. BAM-2-TAGALIGN settings (remove a prefix `atac.` for DNANexus CLI).

   Pipeline filters out chrM reads by default.

   - `"atac.bam2ta.regex_grep_v_ta"` : (optional) Perl-style regular expression pattern to remove matching reads from TAGALIGN (default: `chrM`).
   - `"atac.bam2ta.subsample"` : (optional) Number of reads to subsample TAGALIGN. Subsampled TAGALIGN will be used for all downstream analysis (MACS2, IDR, naive-overlap).

6. Cross correlation analysis settings (remove a prefix `atac.` for DNANexus CLI).

   - `"atac.xcor.subsample"` : (optional) Number of reads to subsample TAGALIGN.

7. MACS2 settings

   **DO NOT DEFINE MACS2 PARAMETERS IN `"atac.macs2"` SCOPE**. All MACS2 parameters must be defined in `"atac"` scope.

   - `"atac.cap_num_peak"` : (optional) Cap number of raw peaks called from MACS2.
   - `"atac.pval_thresh"` : (optional) P-value threshold.
   - `"atac.smooth_win"` : (optional) Size of smoothing window.

8. IDR settings

   **DO NOT DEFINE IDR PARAMETERS IN `"atac.idr"` SCOPE**. All IDR parameters must be defined in `"atac"` scope.

   - `"atac.enable_idr"` : (optional) Set it as `true` to enable IDR on raw peaks.
   - `"atac.idr_thresh"` : (optional) IDR threshold.

9. ATAQC (annotation based analysis) settings

   - `"atac.disable_ataqc"` : (optional) Set it as `true` to disable ATAQC.

### Resource

**RESOURCES DEFINED IN AN INPUT JSON ARE PER TASK**. For example, if you have FASTQs for 2 replicates (2 tasks) and set `cpu` for `bowtie2` task as 4 then total number of cpu cores to map FASTQs is 2 x 4 = 8.

CPU (`cpu`), memory (`mem_mb`) settings are used for submitting jobs to cluster engines (SGE and SLURM) and Cloud platforms (Google Cloud Platform, AWS, ...). VM instance type on cloud platforms will be automatically chosen according to each task's `cpu` and `mem_mb`. Number of cores for tasks without `cpu` parameter is fixed at 1.

- `"atac.trim_adapter.cpu"` : (optional) Number of cores for `trim_adapter` (default: 2).
- `"atac.bowtie2.cpu"` : (optional) Number of cores for `bowtie2` (default: 4).
- `"atac.filter.cpu"` : (optional) Number of cores for `filter` (default: 2).
- `"atac.bam2ta.cpu"` : (optional) Number of cores for `bam2ta` (default: 2).
- `"atac.xcor.cpu"` : (optional) Number of cores for `xcor` (default: 2).
- `"atac.trim_adapter.mem_mb"` : (optional) Max. memory limit in MB for `trim_adapter` (default: 10000).
- `"atac.bowtie2.mem_mb"` : (optional) Max. memory limit in MB for `bowtie2` (default: 20000).
- `"atac.filter.mem_mb"` : (optional) Max. memory limit in MB for `filter` (default: 20000).
- `"atac.bam2ta.mem_mb"` : (optional) Max. memory limit in MB for `bam2ta` (default: 10000).
- `"atac.spr.mem_mb"` : (optional) Max. memory limit in MB for `spr` (default: 12000).
- `"atac.xcor.mem_mb"` : (optional) Max. memory limit in MB for `xcor` (default: 10000).
- `"atac.macs2_mem_mb"` : (optional) Max. memory limit in MB for `macs2` (default: 16000).

Disks (`disks`) is used for Cloud platforms (Google Cloud Platforms, AWS, ...).

- `"atac.trim_adapter.disks"` : (optional) Disks for `trim_adapter` (default: "local-disk 100 HDD").
- `"atac.bowtie2.disks"` : (optional) Disks for `bowtie2` (default: "local-disk 100 HDD").
- `"atac.filter.disks"` : (optional) Disks for `filter` (default: "local-disk 100 HDD").
- `"atac.bam2ta.disks"` : (optional) Disks for `bam2ta` (default: "local-disk 100 HDD").
- `"atac.xcor.disks"` : (optional) Disks for `xcor` (default: "local-disk 100 HDD").
- `"atac.macs2_disks"` : (optional) Disks for `macs2` (default: "local-disk 100 HDD").

Walltime (`time`) settings (for SGE and SLURM only).

- `"atac.trim_adapter.time_hr"` : (optional) Walltime for `trim_adapter` (default: 24).
- `"atac.bowtie2.time_hr"` : (optional) Walltime for `bowtie2` (default: 48).
- `"atac.filter.time_hr"` : (optional) Walltime for `filter` (default: 24).
- `"atac.bam2ta.time_hr"` : (optional) Walltime for `bam2ta` (default: 6).
- `"atac.xcor.time_hr"` : (optional) Walltime for `xcor` (default: 6).
- `"atac.macs2_time_hr"` : (optional) Walltime for `macs2` (default: 24).

 

## Output specification

All output filenames keep prefixes from corresponding input filenames. For example. If you have started from `REP1.fastq.gz` and `REP2.fastq.gz` then corresponding alignment log for each replicate has a filename of `REP1.flagstat.qc` and `REP2.flagstat.qc`, respectively.

Final HTML report (`qc.html`) and QC json (`qc.json`) files do not have any prefix.

1. `DNANexus`: If you choose to use `dxWDL` and run pipelines on DNANexus platform, then output will be stored on the specified output directory without any subdirectories.
2. `Cromwell`: Otherwise `Cromwell` will store outputs for each task under `cromwell-executions/[WORKFLOW_ID]/call-[TASK_NAME]/shard-[IDX]`. For all tasks except `idr` and `overlap`, `[IDX]` means a zero-based index for each replicate but for tasks `idr` and `overlap` it stands for a zero-based index for all possible pair of replicates. For example, you have 3 replicates and all possible combination of two replicates are `[(rep1,rep2), (rep1,rep3), (rep2,rep3)]`. Therefore, `call-idr/shard-2` should be an output directory for the pair of replicate 2 and 3.

For more details, refer to the file table section in an HTML report generated by the pipeline. Files marked as (E) are outputs to be uploaded during ENCODE accession.

| task            | filename                   | description                                                  |
| --------------- | -------------------------- | ------------------------------------------------------------ |
| merge_fastq     | merge_fastqs_R?_*.fastq.gz | Merged FASTQ                                                 |
| trim_fastq      | *.trim_*bp.fastq.gz        | Trimmed FASTQ (R1 only)                                      |
| bwa             | * .bam                     | Raw BAM                                                      |
| bwa             | * .bai                     | BAI for Raw BAM                                              |
| bwa             | * .flagstat.qc             | Samtools flagstat log for raw BAM                            |
| filter          | * .nodup.bam               | Filtered/deduped BAM                                         |
| filter          | * .nodup.flagstat.qc       | Samtools flagstat log for filtered/deduped BAM               |
| filter          | * .dup.qc                  | Picard/sambamba markdup log                                  |
| filter          | * .pbc.qc                  | PBC QC log                                                   |
| bam2ta          | * .tagAlign.gz             | TAG-ALIGN generated from filtered BAM                        |
| bam2ta          | * .N.tagAlign.gz           | Subsampled (N reads) TAG-ALIGN generated from filtered BAM   |
| bam2ta          | * .tn5.tagAlign.gz         | TN5-shifted TAG-ALIGN                                        |
| spr             | * .pr1.tagAlign.gz         | 1st pseudo-replicated TAG-ALIGN                              |
| spr             | * .pr2.tagAlign.gz         | 2nd pseudo-replicated TAG-ALIGN                              |
| pool_ta         | * .tagAlign.gz             | Pooled TAG-ALIGN from all replciates                         |
| fingerprint     | * .jsd.qc                  | DeepTools fingerprint log                                    |
| fingerprint     | * .png                     | DeepTools fingerprint plot                                   |
| choose_ctl      | ctl_for_rep*.tagAlign.gz   | Chosen control for each IP replicate                         |
| xcor            | * .cc.plot.pdf             | Cross-correlation plot PDF                                   |
| xcor            | * .cc.plot.png             | Cross-correlation plot PNG                                   |
| xcor            | * .cc.qc                   | Cross-correlation analysis score log                         |
| xcor            | * .cc.fraglen.txt          | Estimated fragment length                                    |
| macs2           | * .narrowPeak.gz           | NARROWPEAK                                                   |
| macs2           | * .bfilt.narrowPeak.gz     | Blacklist-filtered NARROWPEAK                                |
| macs2           | * .pval.signal.bigwig      | p-val signal BIGWIG                                          |
| macs2           | * .fc.signal.bigwig        | fold enrichment signal BIGWIG                                |
| macs2           | * .frip.qc                 | Fraction of read (TAG-ALIGN) in peaks (NARROWPEAK)           |
| spp             | * .regionPeak.gz           | SPP NARROWPEAK(REGIONPEAK)                                   |
| spp             | * .bfilt.regionPeak.gz     | Blacklist-filtered REGIONPEAK                                |
| spp             | * .frip.qc                 | Fraction of read (TAG-ALIGN) in peaks (REGIONPEAK)           |
| idr             | * .*Peak.gz                | IDR NARROWPEAK                                               |
| idr             | * .bfilt.*Peak.gz          | Blacklist-filtered IDR NARROWPEAK                            |
| idr             | * .txt.png                 | IDR plot PNG                                                 |
| idr             | * .txt.gz                  | Unthresholded IDR output                                     |
| idr             | * .log                     | IDR STDOUT log                                               |
| idr             | * .frip.qc                 | Fraction of read (TAG-ALIGN) in peaks (IDR NARROWPEAK)       |
| overlap         | * .*Peak.gz                | Overlapping NARROWPEAK                                       |
| overlap         | * .bfilt.*Peak.gz          | Blacklist-filtered overlapping NARROWPEAK                    |
| overlap         | * .frip.qc                 | Fraction of read (TAG-ALIGN) in peaks (overlapping NARROWPEAK) |
| reproducibility | * .reproducibility.qc      | Reproducibililty QC log                                      |
| reproducibility | optimal_peak.gz            | Optimal final peak file                                      |
| reproducibility | conservative_peak.gz       | Conservative final peak file                                 |
| qc_report       | qc.html                    | Final HTML QC report                                         |
| qc_report       | qc.json                    | Final QC JSON                                                |