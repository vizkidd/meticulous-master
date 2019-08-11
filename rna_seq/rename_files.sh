#!/bin/bash
sample_file=$1
folder_path=$2
line_count=$(wc -l < $1)
echo $line_count
echo $(ls $folder_path)
for SAMPLE in $(< $sample_file); 
do
pattern=$(echo "$SAMPLE" | awk -F "," {'print $1'})
sample_name=$(echo "$SAMPLE" | awk -F "," {'print $2'})
file=$(ls $folder_path | grep -i "$pattern")
ext=$(echo ${file#"$pattern"})
rename -f 's/oldname/newname/' $folder_path$file $folderpath/$sample_name$ext
#mv $folder_path$file $folderpath/$sample_name$ext
done
