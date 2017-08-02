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
	for JOB_NO in $(seq 1 25)
	do
		TIME_MIN=$(python -c "print $OVERALL_TIME_MIN + ($JOB_NO - 1) * 4 * $TIME_INC");
		TIME_MAX=$(python -c "print $TIME_MIN + 3 * $TIME_INC");

		subjobfilename="subjob_"$POTENTIAL_NAME"_"$JOB_NO".sh";
		echo "#!/bin/bash">$subjobfilename;
		echo " ">>$subjobfilename
		echo "WORK_DIRECTORY="$WORK_DIRECTORY";">>$subjobfilename;
		echo "DATA_DIRECTORY="$DATA_DIRECTORY";">>$subjobfilename;
		echo "cd $WORK_DIRECTORY;">>$subjobfilename;
		echo " ">>$subjobfilename;
		echo "POTENTIAL_NAME="$POTENTIAL_NAME";">>$subjobfilename;
		echo "JOB_NO="$JOB_NO";">>$subjobfilename;
		echo "TIME_INC="$TIME_INC";">>$subjobfilename;
		echo "TEMP_INC="$TEMP_INC";">>$subjobfilename;
		echo " ">>$subjobfilename;
		echo "TIME_MIN="$TIME_MIN";">>$subjobfilename;
		echo "TIME_MAX="$TIME_MAX";">>$subjobfilename;
		echo "TEMP_MIN="$TEMP_MIN";">>$subjobfilename;
		echo "TEMP_MAX="$TEMP_MAX";">>$subjobfilename;
		cat bottom_half.txt>>$subjobfilename;
	done
done
