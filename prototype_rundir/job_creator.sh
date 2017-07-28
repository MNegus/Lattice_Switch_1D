#!/bin/bash

OVERALL_TIME_MIN=0.001;
TIME_INC=0.001;

TEMP_MIN=0.1;
TEMP_MAX=1.0;
TEMP_INC=0.1;



for POTENTIAL_NAME in "KT" "QUARTIC" "DIFF_WIDTH"
do
	for JOB_NO in $(seq 1 10)
	do
		TIME_MIN=$(python -c "print $OVERALL_TIME_MIN + ($JOB_NO - 1) * 10 * $TIME_INC");
		TIME_MAX=$(python -c "print $TIME_MIN + 9 * $TIME_INC");

		jobfilename="runjob_"$POTENTIAL_NAME"_"$JOB_NO".sh";
		echo "#!/bin/bash">$jobfilename;
		echo "POTENTIAL_NAME="$POTENTIAL_NAME>>$jobfilename;
		echo "JOB_NO="$JOB_NO>>$jobfilename;
		echo "TIME_INC="$TIME_INC>>$jobfilename;
		echo "TEMP_INC="$TEMP_INC>>$jobfilename;
		echo " ">>$jobfilename;
		echo "TIME_MIN="$TIME_MIN>>$jobfilename;
		echo "TIME_MAX="$TIME_MAX>>$jobfilename;
		cat bottom_half.txt>>$jobfilename;
	done
done