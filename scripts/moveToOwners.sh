#!/bin/bash
for i in $(squeue -u $USER -h -t PD -o %i)
do
scontrol update jobid=$i partition=owners
done
