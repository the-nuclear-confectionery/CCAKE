#!/bin/bash

# clone BSQ and EBE-vUSPhydro, then: cd ./BSQ

BSQ_dir="../BSQ/Hydro-master/"
EBE_dir="../EBE-vUSPhydro/EBE-Node/v-USPhydro/"
OUT_dir="../BSQ/textfiles/"
changes_file="../BSQ/textfiles/changes.txt"

# diff between whole directories
diff -q $BSQ_dir $EBE_dir >> full_diff_EBE-vUSPhydro_BSQ-Hydro-master.txt

# list files that differ
diff -q $BSQ_dir $EBE_dir | grep "differ" >> changed_files_EBE-vUSPhydro_BSQ-Hydro-master.txt

# save the differences into files
while read p; 
do echo "$p" ; 
diff $BSQ_dir$p $EBE_dir$p >> $OUT_dir"diff_"$p.txt
done < $changes_file
