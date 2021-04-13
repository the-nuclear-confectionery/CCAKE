#!/bin/bash

# 1) clone BSQ to get this file
# 2) clone v-USPhydro2  (folder to compare)
# 3) cd ./BSQ/
# 4) modify if needed next paths

BSQ_dir="../BSQ/Hydro-master/"
EBE_dir="../v-USPhydro2/"
OUT_dir="../BSQ/textfiles/"

changes_file="../BSQ/textfiles/changes.txt"

# 5) save this file and type ./evaluate_diff.sh
# 6) text files with diff inside OUT_dir

# diff between whole directories
diff -q $BSQ_dir $EBE_dir >> $OUT_dir"full_diff_vUSPhydro2_BSQ-Hydro-master.txt"

# list files that differ
diff -q $BSQ_dir $EBE_dir | grep "differ" >> $OUT_dir"changed_files_vUSPhydro2_BSQ-Hydro-master.txt"

# save the differences into files
while read p; 
do echo "$p" ; 
diff -b $BSQ_dir$p $EBE_dir$p >> $OUT_dir"diff_"$p.txt
done < $changes_file
