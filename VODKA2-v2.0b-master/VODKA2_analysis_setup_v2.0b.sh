#!/bin/bash

set -e

usage=$'\n\t\033[1m*** DVG ANALYSIS SETUP ***\033[0m\n
Setup VODKA and BLAST anaysis.\n
Usage:\n\t\033[1m'$(basename "$0")$'\033[0m [-h] \033[1m-f files.txt -d bt2_index -p project\033[0m [-s SUFFIX]
where:
\t \033[1m-f <files.txt>\t\tTXT file containing the list of fastq files\033[0m (1 per sample, with the relative path)
\t \033[1m-d <bt2_index>\tName of VODKA DB index for Bowtie2\033[0m (make sure already exists in VODKA installation folder)
\t \033[1m-p <project>\t\tProject name to label files\033[0m
(option) -i <folder>\t\tVODKA2 installation folder (default: ./VODKA2-v2.0b-master)
(option) -s <suffix>\t\tsuffix of FASTQ files (anything after samplename, default: _merged_nohuman.fq.gz)
(option) -l <on|off>\t\tLSF on or off (default: on)
\t -h  \t\t\tshow this help text and exit\n'

# default values
suffix="_nonhumanG.fq.gz"
lsf="on"
vodkainstall="$(pwd)/VODKA2-v2.0b-master"

options=':hf:d:p:i:s:l:'
while getopts $options option; do
  case "$option" in
    h) echo "$usage"; exit;;
    f) samplesfile=$OPTARG;;
    d) vodkabt2=$(basename "$OPTARG");;
    i) vodkainstall=$OPTARG;;
    p) project=$OPTARG;;
    s) suffix=$OPTARG;;
    l) lsf=$OPTARG;;
    :) echo "$usage"; echo -e "\033[5m\033[4mERROR\033[0m\033[5m:\033[0m missing argument for $OPTARG\n"; exit 1;;
   \?) echo "$usage"; echo -e "\033[5m\033[4mERROR\033[0m\033[5m:\033[0m illegal option -$OPTARG\n"; exit 1;;
  esac
done

# mandatory arguments
if [ ! "$samplesfile" ] || [ ! "$vodkabt2" ] || [ ! "$project" ] ; then
  echo "$usage"
  echo -e "\033[5m\033[4mERROR\033[0m\033[5m:\033[0m options -f, -d, and -p must be provided\n"
  exit 1
fi

if [ $lsf == "on" ] || [ $lsf == "off" ]; then
  echo "LSF is ${lsf} for this run!"
else
  echo -e "\033[5m\033[4mERROR\033[0m\033[5m:\033[0m options -l is ambiguous (please specify 'on' for enabling job submission with LSF or 'off')\n"
  exit 1
fi

echo
echo -e "\033[1m\t* Setup DVG analysis *\033[0m"
echo

# project folder
vodka_folder="${project}_vodka_output"
# path to bowtie2 index
vodka_bt2_idx="${vodkainstall}/bt2_index"
# path to VODKA scripts
vodka_scripts="${vodkainstall}/scripts"

# suffix of VODKA result filenames
vodka_suffix="${vodkabt2}_RESULTS.txt"

# extract DVGtype, genome length and read length from VODKA index name
gl=$(echo ${vodkabt2} | rev | cut -f 2 -d "." | rev)
rl=$(echo ${vodkabt2} | rev | cut -f 1 -d "." | rev)
ty=$(echo ${vodkabt2} | rev | cut -f 3 -d "." | rev)

vodka_blast_script=VODKA2_blast_v2.0b.sh
if [ ${ty} == "DEL" ]; then
  analysis_script=delDVG_analysis_${project}.sh
  FLvirus=$(basename ${vodkabt2} .DEL.${gl}.${rl})
  DVGtype=DEL
else
  analysis_script=cbDVG_analysis_${project}.sh
  FLvirus=$(basename ${vodkabt2} .${gl}.${rl})
  DVGtype=CB
fi

# create results folder on Acive
mkdir -p ${vodka_folder}/vodka_output_${DVGtype}

## SETUP BLAST ANALYSIS FOLDER ###
echo "mkdir ${vodka_folder}/${project}_${DVGtype}_dvg; cd ${vodka_folder}/${project}_${DVGtype}_dvg" > vodka_blast_${project}_${DVGtype}.sh
echo "cp ${vodka_bt2_idx}/${FLvirus}.fasta ." >> vodka_blast_${project}_${DVGtype}.sh
cmd="makeblastdb -in ${FLvirus}.fasta -out ${FLvirus} -dbtype nucl"
if [ $lsf == "on" ]; then
  cmd="LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /emna/rsv -q general -a 'docker(ncbi/blast:latest)' -K ${cmd}"
fi
echo ${cmd} >> vodka_blast_${project}_${DVGtype}.sh
echo "cp ${vodka_scripts}/vodka_blast/${vodka_blast_script} ." >> vodka_blast_${project}_${DVGtype}.sh
echo "cp ${vodka_scripts}/vodka_results_to_fasta/vodka_results_to_fasta_v3.sh ." >> vodka_blast_${project}_${DVGtype}.sh
echo "chmod 777 *sh" >> vodka_blast_${project}_${DVGtype}.sh

# print the list of samples
echo -e "\033[4mList of samples:\033[0m"
echo
while read i; do
        #### sample filename parsing method may need to be adapted to data ########################
	# parsing method 1: delete suffix from samplename
	samplename=$(basename "${i}" ${suffix})
	# parsing method 2: extract info before first "_" and R1/R2
	#samplename=$(basename "${i}" ${suffix} | cut -d "_" -f 1)_$(basename ${i} | sed -e "s/.*\(R[12]\).*/\1/")
        # parsing method 3: extract info before first "_"
        #samplename=$(basename "${i}" ${suffix} | cut -d "_" -f 1)
        ###########################################################################################
        echo ${samplename}

	## VODKA ANALYSIS ##

	# create sample analysis folders
	alignment="${vodka_folder}/vodka_output_${DVGtype}/${samplename}_vodka_output/alignment"
	results="${vodka_folder}/vodka_output_${DVGtype}/${samplename}_vodka_output/results"
	mkdir -p ${alignment} ${results}

	# VODKA STEP 1: Alignment against VODKA DVG database
	cmd1="bowtie2 -p 30 -x ${vodka_bt2_idx}/${vodkabt2} --local -U \"${i}\" -S ${alignment}/${samplename}.${vodkabt2}.sam --no-sq --no-unal --mp 0,0 &"
	if [ $lsf == "on" ]; then
	  cmd1="bsub -o ${project}_step1_aligndvg.out -e ${project}_step1_aligndvg.err -J ${samplename}_step1 -g /emna/rsv -n 12 -R 'rusage[mem=150GB]' -q general -a 'docker(biocontainers/bowtie2:v2.4.1_cv1)' -K ${cmd1}"
	fi
	echo ${cmd1} >> vodka_step1_${project}_${DVGtype}.sh

	# VODKA STEP 2: Extract reads matching DVG ref
	cmd2="perl ${vodka_scripts}/search1.pl ${alignment}/${samplename}.${vodkabt2}.sam ${alignment}/${samplename}.${vodkabt2}.search_output.txt ${gl} ${rl} &"
	if [ $lsf == "on" ]; then
	  cmd2="bsub -g /emna/rsv -q general -a 'docker(perl)' -o ${project}_step2_search.out -e ${project}_step2_search.err -K ${cmd2}"
	fi
	echo ${cmd2} >> vodka_step2_${project}_${DVGtype}.sh

	# VODKA STEP 3: generate DVG candidates fasta file (for step4) and prepare results files (for step5)
	cmd3="perl ${vodka_scripts}/organize.pl ${vodka_folder}/vodka_output_${DVGtype}/${samplename}_vodka_output ${vodkabt2} ${samplename} &"
	if [ $lsf == "on" ]; then
	  cmd3="bsub -g /emna/rsv -q general -a 'docker(perl)' -o ${project}_step3_organize.out -e ${project}_step3_organize.err -K ${cmd3}"
	fi
	echo ${cmd3} >> vodka_step3_${project}_${DVGtype}.sh

	## BLAST VERIFICATION ##
	echo "cp ${results}/${vodka_suffix} ${samplename}_${vodka_suffix}" >> vodka_blast_${project}_${DVGtype}.sh
        cmd4="./${vodka_blast_script}  -r ${FLvirus} -v ${samplename}_${vodka_suffix} -s _${vodka_suffix} -t ${DVGtype} &"
	if [ $lsf == "on" ]; then
	  cmd4="LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /emna/rsv -q general -a 'docker(ncbi/blast:latest)' -o ${samplename}_blast.out -e ${samplename}_blast.err -K ${cmd4}"
	fi
	echo ${cmd4} >> vodka_blast_${project}_${DVGtype}.sh

done < ${samplesfile}

echo "wait" >> vodka_step1_${project}_${DVGtype}.sh
echo "wait" >> vodka_step2_${project}_${DVGtype}.sh
echo "wait" >> vodka_step3_${project}_${DVGtype}.sh
echo "wait" >> vodka_blast_${project}_${DVGtype}.sh
echo "rm ${FLvirus}.*" >> vodka_blast_${project}_${DVGtype}.sh
echo "cd .." >> vodka_blast_${project}_${DVGtype}.sh

echo "cd ${vodka_folder}" > vodka_report_${project}_${DVGtype}.sh
echo "echo -e \"sample\" > samplenames.tmp" >> vodka_report_${project}_${DVGtype}.sh
echo "ls ${vodka_folder}/${project}_${DVGtype}_dvg/*_blast.out | cut -d \"/\" -f 11 | sed \"s/_blast.out//\" >> samplenames.tmp" >> vodka_report_${project}_${DVGtype}.sh
echo "echo -e \"nb_junctions\tnb_dvg_reads\" > dvg_reads.tmp" >> vodka_report_${project}_${DVGtype}.sh
echo "grep Found ${vodka_folder}/${project}_${DVGtype}_dvg/*out | cut -d \" \" -f 4,14 --output-delimiter=\"	\" >> dvg_reads.tmp" >> vodka_report_${project}_${DVGtype}.sh
echo "paste samplenames.tmp dvg_reads.tmp > ${project}_${DVGtype}_dvg_reads.txt" >> vodka_report_${project}_${DVGtype}.sh
echo "rm samplenames.tmp dvg_reads.tmp" >> vodka_report_${project}_${DVGtype}.sh

if [ $lsf == "on" ]; then
  echo "export LSF_DOCKER_PRESERVE_ENVIRONMENT=false" >> ${analysis_script}
fi

cat vodka_step*_${project}_${DVGtype}.sh vodka_blast_${project}_${DVGtype}.sh vodka_report_${project}_${DVGtype}.sh >> ${analysis_script}
chmod 777 ${analysis_script}

echo
echo "Your script is ready!"
ls ${analysis_script}
echo

