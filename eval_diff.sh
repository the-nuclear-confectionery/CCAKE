#!/bin/bash

# clone BSQ and vUSPhydro2, then: cd ./BSQ

BSQ_dir="../BSQ/Hydro-master/"
EBE_dir="../v-USPhydro2/"
OUT_dir="../BSQ/textfiles/"
changes_file="../BSQ/textfiles/changes.txt"

# diff between whole directories
diff -q $BSQ_dir $EBE_dir >> full_diff_vUSPhydro2_BSQ-Hydro-master.txt

# list files that differ
diff -q $BSQ_dir $EBE_dir | grep "differ" >> changed_files_vUSPhydro2_BSQ-Hydro-master.txt

# save the differences into files
while read p; 
do echo "$p" ; 
diff $BSQ_dir$p $EBE_dir$p >> $OUT_dir"diff_"$p.txt
done < $changes_file
