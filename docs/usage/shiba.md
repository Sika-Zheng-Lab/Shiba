# Shiba usage

``` bash
usage: shiba.py [-h] [-p PROCESS] [-s START_STEP] [-v] config

Shiba v0.5.5 - Pipeline for identification of differential RNA splicing

Step 1: bam2gtf.py
    - Assembles transcript structures based on mapped reads using StringTie2.
Step 2: gtf2event.py
    - Converts GTF files to event format.
Step 3: bam2junc.py
    - Extracts junction reads from BAM files.
Step 4: psi.py
    - Calculates PSI values and perform differential analysis
Step 5: expression.py
    - Analyzes gene expression.
Step 6: pca.py
    - Performs PCA.
Step 7: plots.py
    - Generates plots from results.

positional arguments:
  config                Config file in yaml format

options:
  -h, --help            show this help message and exit
  -p PROCESS, --process PROCESS
                        Number of processors to use (default: 1)
  -s START_STEP, --start-step START_STEP
                        Start the pipeline from the specified step (default: 0, run all steps)
  -v, --verbose         Verbose mode
```

Check the [Manual](../manual/diff_splicing_bulk.md/#1-prepare-inputs) to learn how to prepare the `config`.

The `Shiba` command will run the following steps sequentially:

1. `bam2gtf.py`
2. `gtf2event.py`
3. `bam2junc.py`
4. `psi.py`
5. `expression.py`
6. `pca.py`
7. `plots.py`

---

## Step1: `bam2gtf.py`

``` bash
usage: bam2gtf.py [-h] -i INPUT -r REFERENCE -o OUTPUT [-p PROCESSORS] [-v]

Pipeline for transcript assembly using StringTie2

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Experiment table
  -r REFERENCE, --reference REFERENCE
                        Reference GTF file
  -o OUTPUT, --output OUTPUT
                        Assembled GTF file
  -p PROCESSORS, --processors PROCESSORS
                        Number of processors to use (default: 1)
  -v, --verbose         Verbose output
```

## Step2: `gtf2event.py`

``` bash
usage: gtf2event.py [-h] -i GTF [-r REFERENCE_GTF] -o OUTPUT [-p NUM_PROCESS] [-v]

Extract alternative splicing events from GTF file

optional arguments:
  -h, --help            show this help message and exit
  -i GTF, --gtf GTF     Input GTF file
  -r REFERENCE_GTF, --reference-gtf REFERENCE_GTF
                        Reference GTF file
  -o OUTPUT, --output OUTPUT
                        Output directory
  -p NUM_PROCESS, --num-process NUM_PROCESS
                        Number of processors to use
  -v, --verbose         Verbose output
```

Example of `EVENT_SE.txt`:

``` bash
event_id  pos_id  exon  intron_a  intron_b  intron_c  strand  gene_id  gene_name  label
SE_1  SE@GL456354.1@84521-85111@83560-85765  GL456354.1:84521-85111  GL456354.1:83560-84521  GL456354.1:85111-85765  GL456354.1:83560-85765  -  ENSMUSG00000094337  Gm3286  annotated
SE_2  SE@chr10@100080857-100080940@100080130-100087347  chr10:100080857-100080940  chr10:100080130-100080857  chr10:100080940-100087347  chr10:100080130-100087347  +  ENSMUSG00000019966  Kitl  annotated
SE_3  SE@chr10@100485051-100485125@100478022-100487162  chr10:100485051-100485125  chr10:100478022-100485051  chr10:100485125-100487162  chr10:100478022-100487162  -  ENSMUSG00000036676  Tmtc3  annotated
SE_4  SE@chr10@100485051-100485185@100478022-100487162  chr10:100485051-100485185  chr10:100478022-100485051  chr10:100485185-100487162  chr10:100478022-100487162  -  ENSMUSG00000036676  Tmtc3  annotated
SE_5  SE@chr10@100495641-100495661@100494954-100495823  chr10:100495641-100495661  chr10:100494954-100495641  chr10:100495661-100495823  chr10:100494954-100495823  +  ENSMUSG00000019971  Cep290  annotated
SE_6  SE@chr10@100578315-100578431@100577358-100583914  chr10:100578315-100578431  chr10:100577358-100578315  chr10:100578431-100583914  chr10:100577358-100583914  -  ENSMUSG00000046567  4930430F08Rik  annotated
SE_7  SE@chr10@100582263-100582322@100578431-100583914  chr10:100582263-100582322  chr10:100578431-100582263  chr10:100582322-100583914  chr10:100578431-100583914  -  ENSMUSG00000046567  4930430F08Rik  annotated
SE_8  SE@chr10@100594537-100594656@100592429-100595035  chr10:100594537-100594656  chr10:100592429-100594537  chr10:100594656-100595035  chr10:100592429-100595035  +  ENSMUSG00000056912  1700017N19Rik  annotated
SE_9  SE@chr10@100610596-100610715@100609254-100612429  chr10:100610596-100610715  chr10:100609254-100610596  chr10:100610715-100612429  chr10:100609254-100612429  +  ENSMUSG00000056912  1700017N19Rik  annotated
...
```

## Step3: `bam2junc.py`

``` bash
usage: bam2junc.py [-h] -i INPUT -r RI_EVENT -o OUTPUT [-p PROCESSORS] [-a ANCHOR] [-m MIN_INTRON] [-M MAX_INTRON] [-s STRAND] [-v]

Pipeline for processing junction read counts.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Experiment table
  -r RI_EVENT, --ri_event RI_EVENT
                        Intron retention event file
  -o OUTPUT, --output OUTPUT
                        Output junction read counts file
  -p PROCESSORS, --processors PROCESSORS
                        Number of processors to use (default: 1)
  -a ANCHOR, --anchor ANCHOR
                        Minimum anchor length (default: 8)
  -m MIN_INTRON, --min_intron MIN_INTRON
                        Minimum intron size (default: 70)
  -M MAX_INTRON, --max_intron MAX_INTRON
                        Maximum intron size (default: 500000)
  -s STRAND, --strand STRAND
                        Strand specificity (default: XS)
  -v, --verbose         Verbose output
```

## Step4: `psi.py`

``` bash
usage: psi.py [-h] [-p NUM_PROCESS] [-g GROUP] [-f FDR] [-d PSI] [-r REFERENCE] [-a ALTERNATIVE] [-m MINIMUM_READS] [-i] [-t] [--onlypsi] [--onlypsi-group] [--excel] [-v] junctions event output

PSI calculation for alternative splicing events

positional arguments:
  junctions             A bed file of Junction read counts generated by bam2junc.sh
  event                 Directory that contains text files of alternative splicing events generated by gtf2event.py
  output                Directory for output files

optional arguments:
  -h, --help            show this help message and exit
  -p NUM_PROCESS, --num-process NUM_PROCESS
                        Number of processors to use (default: 1)
  -g GROUP, --group GROUP
                        Group information for detecting differential events (default: None)
  -f FDR, --fdr FDR     FDR for detecting differential events (default: 0.05)
  -d PSI, --psi PSI     Threshold of delta PSI for detecting differential events (default: 0.1)
  -r REFERENCE, --reference REFERENCE
                        Reference group for detecting differential events (default: None)
  -a ALTERNATIVE, --alternative ALTERNATIVE
                        Alternative group for detecting differential events (default: None)
  -m MINIMUM_READS, --minimum-reads MINIMUM_READS
                        Minumum value of total reads for each junction for detecting differential events (default: 10)
  -i, --individual-psi  Print PSI for individual samples to output files (default: False)
  -t, --ttest           Perform Welch's t-test between reference and alternative group (default: False)
  --onlypsi             Just calculate PSI for each sample, not perform statistical tests (default: False)
  --onlypsi-group       Just calculate PSI for each group, not perform statistical tests (Overrides --onlypsi when used together) (default: False)
  --excel               Make result files in excel format (default: False)
  -v, --verbose         Verbose output (default: False)
```

## Step5: `expression.py`

``` bash
usage: expression.py [-h] -i INPUT -g REFERENCE -o OUTPUT [-r REFGROUP] [-a ALTGROUP] [-p PROCESSORS] [-v]

RNA expression analysis using featureCounts and DESeq2.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Experiment table
  -g REFERENCE, --reference REFERENCE
                        Reference GTF file
  -o OUTPUT, --output OUTPUT
                        Output directory
  -r REFGROUP, --refgroup REFGROUP
                        Reference group for differential expression analysis
  -a ALTGROUP, --altgroup ALTGROUP
                        Alternative group for differential expression analysis
  -p PROCESSORS, --processors PROCESSORS
                        Number of processors to use (default: 1)
  -v, --verbose         Increase output verbosity
```

## Step6: `pca.py`

``` bash
usage: pca.py [-h] [--input-tpm INPUT_TPM] [--input-psi INPUT_PSI] [-g GENES] [-o OUTPUT] [-v]

Principal Component Analysis for matrix of gene expression and splicing

optional arguments:
  -h, --help            show this help message and exit
  --input-tpm INPUT_TPM
                        Input TPM file (default: None)
  --input-psi INPUT_PSI
                        Input PSI file (default: None)
  -g GENES, --genes GENES
                        Number of highly-variable genes or splicing events to calculate PCs (default: 3000)
  -o OUTPUT, --output OUTPUT
                        Output directory (default: None)
  -v, --verbose         Verbose output (default: False)
```

## Step7: `plots.py`

``` bash
usage: plots.py [-h] [-i INPUT] [-e EXPERIMENT_TABLE] [-s SHIBA_COMMAND] [-o OUTPUT] [-v]

Make plots for alternative splicing events

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Directory that contains result files (default: None)
  -e EXPERIMENT_TABLE, --experiment-table EXPERIMENT_TABLE
                        Experiment table file (default: None)
  -s SHIBA_COMMAND, --shiba-command SHIBA_COMMAND
                        Shiba command (default: None)
  -o OUTPUT, --output OUTPUT
                        Directory for output files (default: None)
  -v, --verbose         Verbose output (default: False)
```
