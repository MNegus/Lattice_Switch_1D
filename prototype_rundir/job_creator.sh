#!/bin/bash

OVERALL_TIME_MIN=0.001;
TIME_INC=0.001;

TEMP_MIN=0.1;
TEMP_MAX=1.0;
TEMP_INC=0.1;

WORK_DIRECTORY="/storage/molsim/phuhzr/Lattice_Switch_1D"


for POTENTIAL_NAME in "KT" "QUARTIC" "DIFF_WIDTH"
do
	DATA_DIRECTORY="/home/theory/phuhzr/Documents/URSS/data_storage/"$POTENTIAL_NAME;
	for JOB_NO in $(seq 1 10)
	do
		TIME_MIN=$(python -c "print $OVERALL_TIME_MIN + ($JOB_NO - 1) * 10 * $TIME_INC");
		TIME_MAX=$(python -c "print $TIME_MIN + 9 * $TIME_INC");

		jobfilename="runjob_"$POTENTIAL_NAME"_"$JOB_NO".pbs";
		echo "#!/bin/bash">$jobfilename;
		echo "#PBS -l nodes=1:ppn=1,mem=256mb,walltime=08:00:00">>$jobfilename;
		echo "#PBS -V">>$jobfilename;
		echo " ">>$jobfilename;
		echo "WORK_DIRECTORY="$WORK_DIRECTORY";">>$jobfilename;
		echo "DATA_DIRECTORY="$DATA_DIRECTORY";">>$jobfilename;
		echo "cd $WORK_DIRECTORY;">>$jobfilename;
		echo " ">>$jobfilename;
		echo "POTENTIAL_NAME="$POTENTIAL_NAME";">>$jobfilename;
		echo "JOB_NO="$JOB_NO";">>$jobfilename;
		echo "TIME_INC="$TIME_INC";">>$jobfilename;
		echo "TEMP_INC="$TEMP_INC";">>$jobfilename;
		echo " ">>$jobfilename;
		echo "TIME_MIN="$TIME_MIN";">>$jobfilename;
		echo "TIME_MAX="$TIME_MAX";">>$jobfilename;
		cat bottom_half.txt>>$jobfilename;
	done
done
