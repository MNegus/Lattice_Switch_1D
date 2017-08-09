#!/bin/bash

TEMP_MIN=0.1
TEMP_INC=0.1
TEMP_MAX=1.5

TIMESTEP_MIN=0.001
TIMESTEP_INC=0.010
TIMESTEP_MAX=0.101

DYNAMICS_TYPE="BAOAB_LIMIT"
TOT_STEPS=10000000
START_WELL=0
SWITCH_REGULARITY=10
X_MIN=-4
X_MAX=4
NOBINS=1000
MASS=1



for POTENTIAL_NAME in "KT" "QUARTIC" "DIFF_WIDTH"
do
	INPUT_FILENUM=1;

	for TEMP in $(seq $TEMP_MIN $TEMP_INC $TEMP_MAX)
	do
		for TIMESTEP in $(seq $TIMESTEP_MIN $TIMESTEP_INC $TIMESTEP_MAX)
		do
			input_filename="input_"$POTENTIAL_NAME"_"$INPUT_FILENUM".txt";
			echo $DYNAMICS_TYPE>$input_filename;
			echo $POTENTIAL_NAME>>$input_filename;
			echo $TOT_STEPS>>$input_filename;
			echo $START_WELL>>$input_filename;
			echo $SWITCH_REGULARITY>>$input_filename;
			echo $X_MIN>>$input_filename;
			echo $X_MAX>>$input_filename;
			echo $NOBINS>>$input_filename;
			echo $TEMP>>$input_filename;
			echo $MASS>>$input_filename;

			if [ "$DYNAMICS_TYPE" == "BAOAB_LIMIT" ]; then
				echo $TIMESTEP>>$input_filename;
			fi

			INPUT_FILENUM=$((INPUT_FILENUM+1));
			echo $INPUT_FILENUM;
		done
	done
done
