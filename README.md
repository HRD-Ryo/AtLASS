# AtLASS
AtLASS (Attention-based bidirectional LSTM for Annotating Splice Sites) is fully automated learning and predicting tool for Splice Sites. <br>
AtLASS requires Genome file in fasta format and RNA-seq files in fastq format or BAM file.

## Hardware requirements
- 48 GB of RAM
- 16 GB of GPU memory

## Installation
You need to install Singularity on your system.
```
$ singularity build --fakeroot Atlass.sif singularity.def
```

## Usage

You can get usage information by the following command. <br>
```
$ singularity run --nv Atlass.sif -h
```
Use of alias is recommended, for example, by adding the following to your ~/.bashrc.
```bash
alias Atlass='singularity run --nv Atlass.sif -h'
```
### Select modes
AtLASS has three modes depending on combination of input files.

#### 1. Leaning and Prediction mode
This mode is the standard. <br>
RNA mapping is performed on the input genome by HISAT2, and models are trained using the annotated splice sites. <br>
Finally, the splice sites of the input genome are predicted. <br>
Command examples are as follows.
```
$ Atlass -g <genome> -1 <rna_1.fastq> -2 <rna_2.fastq> -c <cpu> -gpu
$ Atlass -g <genome> -b <bam> -gpu
```

#### 2. Prediction only mode
This mode uses pre-trained models to predict the input genome. <br>
Command example is as follows.
```
$ Atlass -g <genome> -m1 <model1> -m2 <model2> -gpu
```

#### 3. Learning on the other species and Prediction mode
This mode trains models on the genome of other species and predicts the input genome. <br>
Evidence for splice sites in other species genome can be provided by GenBank file, RNA-seq files, or BAM file. <br>
Command examples are as follows.
```
$ Atlass -g <genome> -o-g <other_genome> -o-gbf <other_GenBank> -gpu
$ Atlass -g <genome> -o-g <other_genome> -o-1 <other_rna_1.fastq> -o-2 <other_rna_2.fastq> -c <cpu> -gpu
$ Atlass -g <genome> -o-g <other_genome> -o-b <other_bam> -gpu
```

## Output
By default, `Atlass_out_${year+month+day}/` directory is created under the current directory, and all output files are created in it. <br>
The main output files to be created are as follows.
- __Atlass.tsv__ : Final prediction results in TSV format. <br>Each column shows:
  1. Scaffold name
  2. Position number of splice site
  3. Intron start (`forward`) or Intron end (`revcom`)
  4. Probability
  5. The results of validation with RNA-seq mapping and one of four categories (mode 1 only): `True Positive`, `False Negative`, `False Positive`, and `Atlass Prediction`. `AtLASS Prediction` should be noted because it is a predicted result in genomic regions not covered by RNA-seq. 

  |  Scaffold  |  Position  |  Start/End  |  Probability  |  Category  |
  | :----: | :----: | :----: | :----: | :----: |
  |  NODE_1  |  124  |  forward  |  0.999  |  True Positive  |
  |  NODE_1  |  211  |  revcom  |  0.998  |  True Positive  |
  |  NODE_1  |  1288  |  revcom  |  0.922  |  False Positive  |
  |  NODE_1  |  3286  |  revcom  |  0.998  |  AtLASS Prediction  |
  |  NODE_1  |  22752  |  forward  |  N/A  |  False Negative  |

- __intron.tsv__ and __exon.tsv__ : TSV files of intron and exon regions
- __model1/2_dir__ : The trained model state files are saved. These can be used in the mode 2 run.
- __filtered.bam__ : The BAM file after mapping and filtering, which can be used with the --BAM flag in the next run.
- __IntronTrain1/2.log__ : Logfile of each model training.
- __Std Error__ : The end time of each step is output. Example is follows.
```
  [2022-07-29 04:15:03]  AtLASS-0.0.1 mode:1 start.
  [2022-07-29 04:15:03]  Outdir is made at "/home/harada/AtLASS_test/Atlass_out_20220729".
  [2022-07-29 05:35:45]  IntronExtract Finish.
  [2022-07-29 05:35:46]  DataMake1 Finish.
  [2022-07-30 16:29:41]  IntronTrain Finish.
  [2022-07-31 22:30:44]  GenomePred Finish.
  [2022-07-31 23:13:58]  DataMake2 Finish.
  [2022-08-01 10:13:17]  IntronTrain Finish.
  [2022-08-02 16:15:28]  GenomePred Finish.
  [2022-08-02 22:38:21]  AtLASS-0.0.1 mode:1 Finish.
```

## Flags
```
usage: Atlass [-h] [-o OUTDIR] [-g GENOME] [-b BAM] [-1 FASTQ_1] [-2 FASTQ_2] [-c CPU] [-gpu] [-f] [-m1 MODEL1]
              [-m2 MODEL2] [-o-g OTHER_GENOME] [-o-gbf OTHER_GENBANK] [-o-b OTHER_BAM] [-o-1 OTHER_FASTQ_1]
              [-o-2 OTHER_FASTQ_2] [-v]

optional arguments:
  -h, --help            show this help message and exit
  -o OUTDIR, --outdir OUTDIR
                        Output directory 
                        Default is ./output_${year+month+day}
  -g GENOME, --genome GENOME
                        Genome fasta file to be predicted
  -b BAM, --bam BAM, --BAM BAM, --sam BAM, --SAM BAM
                        BAM or SAM file for learning 
                        HISAT2 is expected mapper
  -1 FASTQ_1, --fastq_1 FASTQ_1
                        Fastq or Fastq.gz file 1 for learning
  -2 FASTQ_2, --fastq_2 FASTQ_2
                        Fastq or Fastq.gz file 2 for learning
  -c CPU, --cpu CPU, --CPU CPU
                        Number of cpu in mapping phase 
                        Default is 1
  -gpu, --gpu, --GPU    Run with GPU device in learning and predicting phase 
                        Default is False, recommended is True 
                        CPU mode will use ALL CPUs
  -f, --fast, --FAST    Enable fast model 
                        Fast mode uses smaller models and may result in less accuracy
  -m1 MODEL1, --model1 MODEL1
                        Pre-trained model_1 file
  -m2 MODEL2, --model2 MODEL2
                        Pre-trained model_2 file
  -o-g OTHER_GENOME, --other_genome OTHER_GENOME
                        Genome fasta file of other species for learning
  -o-gbf OTHER_GENBANK, --other_GenBank OTHER_GENBANK
                        GenBank file of other species for learning
  -o-b OTHER_BAM, --other_bam OTHER_BAM
                        BAM or SAM file of other species for learning
  -o-1 OTHER_FASTQ_1, --other_fastq_1 OTHER_FASTQ_1
                        Fastq or Fastq.gz file 1 of other species for learning
  -o-2 OTHER_FASTQ_2, --other_fastq_2 OTHER_FASTQ_2
                        Fastq or Fastq.gz file 2 of other species for learning
  -v, --version         Show version and exit
```

## License
See the [LICENSE](LICENSE) file for details.

## Publication and Citation
Comming soon
