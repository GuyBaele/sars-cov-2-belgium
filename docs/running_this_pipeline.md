# Running the pipeline

## Data downloads

### Nextstrain fasta and metadata files

1. navigate to gisaid.org and login
2. select EpiCoV
3. select Downloads
4. scroll to the bottom where it says "Genomic Epidemiology"
  * select "FASTA", agree to the conditions and download
    * I download to "data/gisaid/", but this is just a temp. location
  * repeat for "metadata"
5. go to wherever the files are downloaded and for each file:
 `$ gunzip <fname>`
 * this will take a while for the fasta, usually around 5 minutes on my computer
6. after the files are extracted, rename them to `data/sequences.fasta` and `data/metadata.tsv`, respectively

### Variant lists

### Nextstrain background sequence list
