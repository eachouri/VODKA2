#!/usr/bin/env bash

vodka_input=$1
dir=$(dirname ${vodka_input})
fastaB=$(basename ${vodka_input} .txt)_B.fa
fastaR=$(basename ${vodka_input} .txt)_R.fa

tail -n +2 ${vodka_input} | cut -f 1,2,10 | tr -d "<>" | sed 's/[ATCGN]*::://' | sed 's/\t/-/' | sed 's/^/>/g' | tr "\t" "\n" > ${dir}/${fastaB}
tail -n +2 ${vodka_input} | cut -f 1,2,10 | tr -d "<>" | sed 's/:::.*//' | sed 's/\t/-/' | sed 's/^/>/g' | tr "\t" "\n" > ${dir}/${fastaR}

exit 0
