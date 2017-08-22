# Parsl Aligner

Parsl Aligner skeleton for initial handoff.

## Requirements

* Make sure you have Python3.5.
* Install parsl package:

      pip3 install parsl


## How to run:

* Run the setup script to generate a mock fastq file
* Start ipcluster with command:

    ipcluster start -n 2

* Run the script:

    python3 par_align.py
