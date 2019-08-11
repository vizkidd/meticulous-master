#!/bin/sh

#module load python
#module load fastqc

show_spinner()
{
  local -r pid="${1}"
  local -r delay='0.75'
  local spinstr='\|/-'
  local temp
  while ps a | awk '{print $1}' | grep -q "${pid}"; do
    temp="${spinstr#?}"
    printf " [%c]  " "${spinstr}"
    spinstr=${temp}${spinstr%"${temp}"}
    sleep "${delay}"
    printf "\b\b\b\b\b\b"
  done
  printf "    \b\b\b\b"
}

##FUNCTIONS:
paralellize_fastqc() ##Runs fastqc on each sample file
{
nohup fastqc -o "./" $1 &> "./"$1.out& ##fastqc & write output to directory of the bash script
proc_ids+=($!) ##Store process IDs
}

##STANDALONE SCRIPT/SUBSCRIPT
##ENTRYPOINT
##Check paramter count ($#)
if [ $# -eq 0 ]
 then
	"Please give path to fastq files followed by file extension(without .) for analysis"
	exit 1
fi

proc_ids=() ##Used to store process IDs
#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

##Set extension to look for
##Note: This is a stupid overhaul
if [ -z $2 ]
then
   ext="*"
else
   ext=$2
fi

cd "$1" ##cd to given dir

##if cd unsuccessful then dir doesn't exist
if [ $? != 0 ]
then
  echo "Check file path"
  exit 1
fi

##Print info
echo "Dir: $(pwd)"
echo "Ext: $ext"
echo "Files: $(ls)"

##Choose files by extension & run fastqc
for f in *.$ext; do
  echo "File Chosen-> $f"
  paralellize_fastqc $f ##Called for each file
done

##Wait for processes to complete
for id in ${proc_ids[@]}; do
  ##show_spinner $id
  wait $id
done

echo $?
multiqc "./" ##Run multiqc on directory

