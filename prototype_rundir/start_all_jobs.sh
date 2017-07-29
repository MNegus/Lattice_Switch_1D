#!/bin/bash
for f in /storage/molsim/phuhzr/Lattice_Switch_1D/prototype_rundir*.pbs; do
	qsub $f
done