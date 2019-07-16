# Open Targets Tractability Pipeline

## Installation
Change to the directory containing this file

`pip install .`

Install cxoracle

Set the following environment variables:

`CHEMBL_DB=oracle://address:to@local.chembl` 

`CHEMBL_VERSION=25`

## Getting Started

Run the pipeline with the following command:

`run-ot-pipeline genes.csv`

Where `genes.csv` is a file with one Ensembl Gene ID per line with no headers

