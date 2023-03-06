#!/bin/bash
# Filters Copycat (https://github.com/MariaNattestad/copycat) output
# to prepare it for visualization in Circa v1.2.3 (https://omgenomics.com/circa/)

DATA_DIR="../data"
OUT_DIR="out"

convert() {
    input="$DATA_DIR/copycat.$1.coverage.10kb.csv"
    output="$OUT_DIR/$1.coverage.circa.csv"
    threshold=$2

    echo "chromosome,start,end,unsegmented_coverage" > $output
    while IFS="," read -r rec_column1 rec_column2 rec_column3 rec_column4
    do
        mod=$(($rec_column2%200000))
        if [ "$mod" -eq "0" ]; then
            if (( $(echo "$rec_column4 < $threshold" |bc -l) )); then
                echo "$rec_column1,$rec_column2,$rec_column3,$rec_column4" >> $output
            else
                echo "$rec_column1,$rec_column2,$rec_column3,$threshold" >> $output
            fi
        fi
    done < <(tail -n +2 $input)
}

convert "pcrfree" 70
convert "pb" 60
convert "ont" 60
