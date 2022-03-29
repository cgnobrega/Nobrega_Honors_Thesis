# Running RNA-seq pipeline using singularity image

***Note** .bashrc was not automatically running, so we had to do `source .bashrc` so be sure to check if bashrc is running*

## Locating the images
```
ls /compbio/software/sif
```

In the above directory is the `.sif` files for all the RNA seq processing software. The extension `.sif` is for singularity container. In order to get the most recent version of the software, there are symbolic links for each software to grab the most up to date version. 

## Getting help on running an image
```
singularity run /compbio/software/sif/fastqc_v_latest.sif fastqc -h
```

This command will run the singularity image for fastqc with the help text. 

## Script baiscs
`fastqc_singularity` will give you the script help:
```
---------------------------------
---------- SCRIPT HELP ----------

COMMAND: fastqc_singularity.sh

USAGE
-----
SHOW HELP:
./fastqc_singularity.sh -h, --help

RUN COMMAND SCRIPT:
./fastqc_singularity.sh -p, --param_file <path_to_param_file>/<your_param_file>.txt

EXPORT PARAM FILE TEMPLATE:
./fastqc_singularity.sh -x, --export_param_template [optional_template_file_name]
```

To run script you need a param file which is created with this command:
```
fastqc_singularity -x
```
Which will output a file called: `fastqc_parameter_template.txt`. However, you can add your own descriptive filename to have easier to read output. 

**You need a param file for each sample,** so there will be 18 param files for each software for this project. Give them meaningful names and store them in a directory to keep organized.

### Param file

Will look slightly different for each program (see end of this file). Use `nano` to edit the param file and put in appropriate path names and file names. 

Use the following convention for the read files: (replace SLxxxxxx with the sample #)

```
left_read_file=jcoffman_007/SLxxxxxx_R1.fastq.gz
right_read_file=jcoffman_007/SLxxxxxx_R2.fastq.gz
```
Use a meaningful output directory sturcture for `output_directory=`. Do not modify anything else in the file unless necessary and do not put spaces or extra `/` characters. 

***Note for `trimgalore`: be sure that container image either doesn't have `_`, check path*** 

### Run the script

Now that you have your param file, you can run the script using: 
```
fastqc_singularity -p fastqc_parameter_template.txt
```

## What is the script doing? 

To get an idea of what each script acutally does take a look at the code: 
```
cat /compbio/software/compbio_command_scripts/fastqc/fastqc_singularity.sh
```

## Order of go

1. `fastqc` goes first. 18 param files
2. `trimgalore` goes next. Will run `cutadapt` first and `fastqc` again. 18 param files
3. `kallisto_index` goes only **once**.
*:question: we ran it during demo meeting, that should mean I'm all set?*
5. `kallisto_quant` 18 param files, see below for param file notes 

# Param files contd. 

## FastQC:

```
left_read_file=jcoffman_007/SL139772_R1.fastq.gz
right_read_file=jcoffman_007/SL139772_R2.fastq.gz        #This parameter is only required if samples are paired-end.
output_directory=fastqc_out
root_working_directory=/compbio

#------------------------------------------------------------------- #
# DO NOT MODIFY PARAMETERS BELOW UNLESS YOU KNOW WHAT YOU ARE DOING! #
#------------------------------------------------------------------- #

fastqc_options=""
fastqc_flags=""

container_image=fastqc_v_latest.sif
```
**Notes**:
* `output_directory=` should be sample name, and maybe be in parent directory to hold all fastqc output


## Trim Galore: 
```
left_read_file=jcoffman_007/SL139772_R1.fastq.gz
right_read_file=jcoffman_007/SL139772_R2.fastq.gz       #This parameter is only required if samples are "paired-end".
trim_galore_output_directory=trimgalore_out
root_working_directory=/compbio

#--------------------------------------------------------------------#
# DO NOT MODIFY PARAMETERS BELOW UNLESS YOU KNOW WHAT YOU ARE DOING! #
#--------------------------------------------------------------------#

trim_galore_options="--fastqc_args \"--noextract\" --length_1 35 --length_2 35 --stringency 1 --length 20 --quality 20"
trim_galore_flags="--paired --illumina --retain_unpaired"

container_image=trimgalore_v_latest.sif
```
**Notes**:
* `trim_galore_output_directory` should be sample name, and maybe be in parent directory to hold all trim galore output
* `container_image` should not have and underscore (double check path)

## Kallisto Index
```
input_file=ensembl_91_zebrafish/danio_rerio-joined.fa
output_directory=kallisto_index #This directory must exist. Kallisto will not create it.
index=ensembl_91_zebrafish      #Filename for the kallisto index to be constructed
root_working_directory=/compbio

#------------------------------------------------------------------- #
# DO NOT MODIFY PARAMETERS BELOW UNLESS YOU KNOW WHAT YOU ARE DOING! #
#------------------------------------------------------------------- #

kallisto_options="--kmer-size=25"
kallisto_flags=""

container_image=kallisto_v_latest.sif
```
**Notes**: 
* only run **once**
* `input_file` uses joined
* `output_directory` should match begining `kallisto_index` directory for kallisto quant (see below)

## Kallisto Quant
```
left_read_file=jcoffman_007/SL139772_R1.fastq.gz
right_read_file=jcoffman_007/SL139772_R2.fastq.gz       #Required if reads are paired-end
output_directory=SL139772
kallisto_index=kallisto_index/ensembl_91_zebrafish              #Filename for the kallisto index constructed with kallisto index
threads=8
root_working_directory=/compbio

#------------------------------------------------------------------- #
# DO NOT MODIFY PARAMETERS BELOW UNLESS YOU KNOW WHAT YOU ARE DOING! #
#------------------------------------------------------------------- #

kallisto_options=""
kallisto_flags=""

container_image=kallisto_v_latest.sif
```
**Notes**:
* add parent directory to output directory 
* note `kallisto_index`

## Tximport
  * `tximport` needs transcript to gene map so it can make both a transcript 
    and gene count matrix
  * there already exists one, but it has only transcript id without 
    version number which will mismatch with `kallisto` output and make 
    `tximport` fail.
  * so before running `tximport` we need to remake the tx2gene file based 
    the gtf file. The steps for how this was done are below:
    
### tx2gene file
  * `.gtf` file has 9 columns. The 3rd column describes the entry type, and the 
    9th column has all the information (attributes) about the entry in a semi-colon
    delimited list. We are only interested in rows where the 3rd column says "trasncript"
  * the final tx2gene file will be a text file with two columns tab delimited where the 
    first column is the `gene_id` and the second column is the `transcript_id` with the 
    `transcript_version`.
  * the first step is to take all the rows from the `.gtf` file that are transcripts. This 
    was done with the following code:
    ```
     grep -v "#" ensembl_91_zebrafish/Danio_rerio.GRCz10.91.gtf | cut -f 3,9 | awk '{ if($1 == "transcript") print $0;}' > transcripts_Danio_rerio.GRCz10.91.txt
    ```
  * this puts all the transcript entries into a file called `transcripts_Danio_rerio.GRCz10.91.txt`
    with only rows 3 and 9 (the rows we are interested in)
  * next, to create the final tx2gene file, the following bash script was used:
    ```
    #!/bin/bash
    . /home/ubuntu/.profile
    touch v_Danio_rerio.GRCz10.91.tx2gene.txt
    
    while IFS= read -r line
    do
    ATR="$(echo $line | cut -f2)"
    GENE_ID="$(echo $ATR | cut -d ';' -f1 | grep -o '".*"' | sed 's/"//g')"
    TX_ID="$(echo $ATR | cut -d ';' -f3 | grep -o '".*"' | sed 's/"//g')"
    TX_V="$(echo $ATR | cut -d ';' -f4 | grep -o '".*"' | sed 's/"//g')"
    TX_F="$TX_ID.$TX_V"
    echo -e $GENE_ID'\t'$TX_F >> v_Danio_rerio.GRCz10.91.tx2gene.txt
    done < transcripts_Danio_rerio.GRCz10.91.txt
    ```
  * this generated the file `v_Danio_rerio.GRCz10.91.tx2gene.txt` which is the 
    final tx2gene file. Here is the first 10 lines:
    ```
    ENSDARG00000104632      ENSDART00000166186.2
    ENSDARG00000100660      ENSDART00000166174.2
    ENSDARG00000098417      ENSDART00000157825.1
    ENSDARG00000100422      ENSDART00000172566.1
    ENSDARG00000100422      ENSDART00000169187.1
    ENSDARG00000100422      ENSDART00000169667.2
    ENSDARG00000102128      ENSDART00000167290.2
    ENSDARG00000102128      ENSDART00000169805.1
    ENSDARG00000102128      ENSDART00000158376.2
    ENSDARG00000102128      ENSDART00000179581.1
    ```

# Param file scripts 

### FastQC

```
#!/bin/bash
. /home/ubuntu/.profile
PRGM='_fastqc_param'
for FILE in ./jcoffman_007/SL*R1*;
do
LEFT=${FILE:2:23}1.fastq.gz;
RIGHT=${FILE:2:23}2.fastq.gz;
PARAM_NAME=${FILE:15:8}$PRGM.txt;
OUT=${FILE:15:8};

/compbio/software/compbio_command_scripts/fastqc/fastqc_singularity.sh -x ./fastqc_param_files/$PARAM_NAME;

sed -i -e "s|left_read_file=|left_read_file=$LEFT|" ./fastqc_param_files/$PARAM_NAME;
sed -i -e "s|right_read_file=|right_read_file=$RIGHT|" ./fastqc_param_files/$PARAM_NAME;

mkdir fastqc_out/$OUT;
sed -i -e "s|output_directory=|output_directory=fastqc_out/$OUT|" ./fastqc_param_files/$PARAM_NAME;

# echo "left read file: $LEFT "
# echo "right read file: $RIGHT "
# echo "param file name: $PARAM_NAME "
# echo "out directory: $OUT "

/compbio/software/compbio_command_scripts/fastqc/fastqc_singularity.sh -p ./fastqc_param_files/$PARAM_NAME
done
```

### Trim Galore

```
#!/bin/bash
PRGM='_trimgalore_param'
for FILE in ./jcoffman_007/SL*R1*;
do
LEFT=${FILE:2:23}1.fastq.gz;
RIGHT=${FILE:2:23}2.fastq.gz;
PARAM_NAME=${FILE:15:8}$PRGM.txt;
OUT=${FILE:15:8};

/compbio/software/compbio_command_scripts/trim_galore/trim_galore_singularity.sh -x ./trimgalore_param_files/$PARAM_NAME;

sed -i -e "s|left_read_file=|left_read_file=$LEFT|" ./trimgalore_param_files/$PARAM_NAME;
sed -i -e "s|right_read_file=|right_read_file=$RIGHT|" ./trimgalore_param_files/$PARAM_NAME;

mkdir trimgalore_out/$OUT;
sed -i -e "s|trim_galore_output_directory=|trim_galore_output_directory=trimgalore_out/$OUT|" ./trimgalore_param_files/$PARAM_NAME;

sed -i -e "s|trim_galore_v_latest.sif|trimgalore_v_latest.sif|" ./trimgalore_param_files/$PARAM_NAME

# echo "left read file: $LEFT "
# echo "right read file: $RIGHT "
# echo "param file name: $PARAM_NAME "
# echo "out directory: $OUT "

/compbio/software/compbio_command_scripts/trim_galore/trim_galore_singularity.sh -p ./trimgalore_param_files/$PARAM_NAME
done
```

### Kallisto Quant

```
#!/bin/bash
PRGM='_kallisto_quant_param'
for FILE in ./jcoffman_007/SL*R1*;
do
LEFT=${FILE:2:23}1.fastq.gz;
RIGHT=${FILE:2:23}2.fastq.gz;
PARAM_NAME=${FILE:15:8}$PRGM.txt;
OUT=${FILE:15:8};

/compbio/software/compbio_command_scripts/kallisto/kallisto_quant_singularity.sh -x ./kallisto_quant_param_files/$PARAM_NAME;

sed -i -e "s|left_read_file=|left_read_file=$LEFT|" ./kallisto_quant_param_files/$PARAM_NAME;
sed -i -e "s|right_read_file=|right_read_file=$RIGHT|" ./kallisto_quant_param_files/$PARAM_NAME;

mkdir kallisto_quant_out/$OUT;
sed -i -e "s|output_directory=|output_directory=kallisto_quant_out/$OUT|" ./kallisto_quant_param_files/$PARAM_NAME;

sed -i -e "s|kallisto_index=|kallisto_index=kallisto_index/ensembl_91_zebrafish|" ./kallisto_quant_param_files/$PARAM_NAME;

sed -i -e "s|threads=|threads=8|" ./kallisto_quant_param_files/$PARAM_NAME;

# echo "left read file: $LEFT "
# echo "right read file: $RIGHT "
# echo "param file name: $PARAM_NAME "
# echo "out directory: $OUT "

/compbio/software/compbio_command_scripts/kallisto/kallisto_quant_singularity.sh -p ./kallisto_quant_param_files/$PARAM_NAME
done
```
















