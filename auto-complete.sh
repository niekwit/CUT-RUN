#!/bin/bash

#This script enables autocompletion for the CRISPR library options (-l) in crispr-pipeline.sh
#Add the following line to your .bashrc file: source /path/to/CRISPR-tools/auto-complete.sh

#finds CRISPR libary names
SCRIPT_DIR=$(find $HOME -type d -name "CUT-RUN")
genome_list=$(cat "$SCRIPT_DIR/config.yml" | shyaml get-values genome | grep -v ':' | tr "\n" " ")

#enables autocompletion of `-g` flag
function libs()
{
case $3 in
    -g) COMPREPLY=($(compgen -W "$genome_list" "${COMP_WORDS[$COMP_CWORD]}"));;
  esac
}

complete -F libs cutrun.sh
