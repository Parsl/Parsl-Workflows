#!/bin/bash


create_fast_q () {
    echo "Creating a mock fast_q file"
    shuf -i 1-1000000 -o mock.fastq
    echo "Done : mock.fastq"
}


create_fast_q
export PYTHONPATH=$PWD/..:$PYTHONPATH
alias python='python3.5'
alias python3='python3.5'
