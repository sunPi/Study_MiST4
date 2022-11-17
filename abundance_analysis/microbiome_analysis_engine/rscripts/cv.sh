#!/bin/bash

REPS=$1
DO_SPLIT=$2
# FILE=$2


# echo $REPS
echo Will run boosting for $REPS iterations!


for R in $(seq 1 $REPS)
do
    echo Cross-validation iteration $R...
    # mkdir -pv model_$R 
    Rscript roc_analysis.R $R $DO_SPLIT
done

echo "Set a name to move the results into a seperate folder (or press enter to skip): "
read OUT

mkdir $OUT
mv ./model* ./$OUT/