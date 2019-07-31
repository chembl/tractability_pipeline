# Open Targets Tractability Pipeline

## Introduction

This pipeline has been developed to produce tractability data for list of input Ensembl Gene IDs. This implementation
is based on the public version of the GSK tractability pipeline, published [here](https://pubs.rsc.org/en/content/articlelanding/2018/md/c7md00633k#!divAbstract)

The pipeline produces a TSV file with one target per row. Small molecule tractability buckets are denoted with "Bucket_X",
antibody buckets with "Bucket_X_ab" and PROTAC buckets with "Bucket_X_PROTAC".

In addition to PROTAC tractability buckets, there is an additional "PROTAC_location_Bucket", which allows you to assess
whether a target's location is suitable for the PROTAC approach.

1. High confidence good location
2. Med confidence good location
3. High confidence grey location
4. Med condifence grey location
5. Unknown location
6. Med confidence bad location
7. High confidence bad location



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

