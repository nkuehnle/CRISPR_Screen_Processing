#!/bin/bash
PROJ_NAME="${PWD##*/}" # Set default name to the name of the current directory
LIBRARY=Brunello # Set default library to Brunello

cd "$(dirname "$0")" || exit

while getopts n:l flag # Names of flags
do
    case "${flag}" in
        n) PROJ_NAME=${OPTARG};; # PROJ_NAME Flag -n
        l) LIBRARY=${OPTARG};; # LIBRARY Flag -l
    esac
done

DRUGZ_CONTROLS_LIST="$(echo Libraries/Controls/"$LIBRARY"_drugz_control_list.txt)"

while IFS=$'\t' read -r Output_Name Treatment Control
do
    python drugz.py -i "MAGeCK/$PROJ_NAME".count.txt -o DrugZ/"$Output_Name" -c "$Control" -x "$Treatment" -r "$DRUGZ_CONTROLS_LIST"
done < DrugZ_Tests_Table.txt
