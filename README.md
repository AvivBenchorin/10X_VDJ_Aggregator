# 10X_VDJ_Aggregator
This directory contains Python scripts used to process output files from the 10X Genomics [`cellranger vdj`](https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/using/vdj) pipeline.

# Overview
 This tool is intended to help process single-cell immune profiling data outputed from 10X Genomics' [`cellranger vdj`](https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/using/vdj) pipeline. This directory contains:
 1)  `vdj_aggr.py`, which is used to aggregate the FASTA and annotation CSV files across multiple `vdj` runs. This script takes as input a configuration file with paths to the desired input FASTA and annotation files and reference CSV files with the transcripts to keep, and outputs an aggregated FASTA file and/or aggregated annotation CSV file. Users can specify to only aggregated the FASTA files or annotation CSV files if desired.
 2) `AIRR_processor.py`, which is used to add metadata to a [AIRR](https://docs.airr-community.org/en/stable/index.html) formatted file. This script takes as input an AIRR TSV file and a CSV file with metadata labels and the metadata to add specified cells, and outputs a modified AIRR file.


# `vdj_aggr.py` Usage

### Input/output FASTA and annotaiton CSV file formats
This script processes FASTA and annotation CSV files outputed by `cellranger vdj` as input. For example, you could use the `filtered_contig.fasta` and `filtered_contig_annotations.csv` outputs from `vdj` for the filtered contigs, or the `all_contig.fasta` and `all_contig_annotations.csv` outputs for the non-filterd contigs.

This script can outputs one FASTA file and/or one annotation CSV file, which contain the aggregated contents of the provided files.
### Running the program
Usage of `vdj_aggr.py` is as follows:
```
python3 vdj_aggr.py [-h] [-a ANNOTATION_OUTPUT] [-f FASTA_OUTPUT] [-v] config

positional arguments:
  config                CSV with the input paths and metadata label

optional arguments:
  -h, --help            show this help message and exit
  -a ANNOTATION_OUTPUT, --annotation_output ANNOTATION_OUTPUT
                        path to write annotation CSV file to
  -f FASTA_OUTPUT, --fasta_output FASTA_OUTPUT
                        path to write contigs fasta file to
  -v, --verbose         increase logging verbosity
```

The only required argument is the configuration file specifying the input FASTA and annotation CSV files to process.
### Configuration file format
The configuration file should be formatted as a CSV file. The first line should specify which type of files of aggregate and the labels of any added metadata to the annotation file, with the following structure: 

    <compare_annotation>,<compare_fasta>,[metadata_label_1],[metadata_label_2]...
In which `<compare_annotation>` and `<compare_fasta>` are booleans (TRUE or FALSE), specifying whether to compare annotation CSV files or FASTA files respectively. Any number of comma-separated metadata labels can be added, and the metadata columns will be added to the annotation file header in the order given.


Subsequent lines should specify the paths to the files from each `vdj` run, with the following structure assuming both FASTA and annotation CSV files will be processed:

    <run_id>,<reference_path>,<annotation_input>,<fasta_input>

The `<run_id>` is an id used to uniquely identify the barcodes of cells from each run, and would replace the `-1` or similar `-<number>` indicator at the end each barcode from that run. To not modify the original barcode, just leave this argument as an empty string (e.g. `,<reference_path>,...`).

The `<reference_path>` is a path to a reference CSV file, which contains all of the transcripts of cells to keep from this run, as well as the values of any added metadata columns for specified cells. The structure of this file is specified in the section below.

The `<annotation_input>` and `<fasta_input>` are the paths to the input annotation CSV and FASTA files respectively for this `vdj` run. If it is specified to only compare FASTA files or annotation CSV files, use the following structure:

    <run_id>,<reference_path>,<annotation_input OR fasta_input>.

Here are a few sample configuration files:
1) Aggregating both annotation CSV and FASTA files for two `vdj` runs, not adding metadata labels or unique `<run_id>` values:
```
TRUE,TRUE
,/path/to/sample1_reference.csv,/path/to/sample1_contigs_annotation.csv,/path/to/sample1_contigs.fasta
,/path/to/sample2_reference.csv,/path/to/sample2_contigs_annotation.csv,/path/to/sample2_contigs.fasta
```

2) Aggregating only annotation CSV files for three `vdj` runs, adding metadata labels and unique `<run_id>` values:
```
TRUE,FALSE,metadata1,metadata2
run1,/path/to/sample1_reference.csv,/path/to/sample1_contigs_annotation.csv
run2,/path/to/sample2_reference.csv,/path/to/sample2_contigs_annotation.csv
run3,/path/to/sample3_reference.csv,/path/to/sample3_contigs_annotation.csv
```

### Reference file format
The reference file should be formatted as a CSV file that contains all of the transcripts of cells to keep from its corresponding `vdj` run, as well as the values of added metadata columns to assign to each cell. The structure of the file should be as follows:
```
<transcript>,[metadata_value],[metadata_value]
<transcript>,[metadata_value],[metadata_value]
...
```
The `<transcript>` value is the transcript of the cell, while the  `[metadata_value]` are values for any metadata columns added to the annotation CSV file in the configuration file to assign to each the cell.

A sample reference file if no metadata columns were added would be:
```
AAAAAAAAAAAAAAAA
AAAAAAAAAAAAAATT
AAAAAAAAAAAGCGCG
AAAAAAAAAAGTCATC
...
```

A sample reference file if two metadata columns were added would be:
```
AAAAAAAAAAAAAAAA,10,6
AAAAAAAAAAAAAATT,9,4
AAAAAAAAAAAGCGCG,15,9
AAAAAAAAAAGTCATC,2,10
...
```

All cells whose transcripts are not included in the reference file will not be included in the aggregated output file, and each cell must be assigned a value for each added metadata column.

# `AIRR_processor.py` Usage

### Input/output AIRR file formats
This script processes [AIRR-formatted](https://docs.airr-community.org/en/stable/index.html) TSV files. An example AIRR file this can take in as input would be the `filtered_contig_heavy_germ-pass.tsv` file outputed by [Change-O](https://immcantation.readthedocs.io/en/stable/tutorials/10x_tutorial.html#assign-v-d-and-j-genes-and-define-clonal-groups) in the Immcantation pipeline.

### Running the program
Usage of `AIRR_processor.py` is as follows:
```
AIRR_processor.py [-h] [-o AIRR_OUTPUT] [-v] AIRR_input metadata_list

positional arguments:
  AIRR_input            path to input AIRR TSV file
  metadata_list         CSV with metadata labels, and the metadata to add to each transcript

optional arguments:
  -h, --help            show this help message and exit
  -o AIRR_OUTPUT, --AIRR_output AIRR_OUTPUT
                        path to write the AIRR TSV file to
  -v, --verbose         increase logging verbosity
```

The only required arguments are `AIRR_input` and `metadata_list`.
### Metadata file format
The metadata file should be formatted as a CSV file. The first line should specify the labels of metadata columns to add to the AIRR file. Subsequent lines should contains the transcripts of cells to assign metadata values to, and their associated metadata values. The metadata file will have the following structure:
```
<metadata1_label>,[metadata2_label],...
<transcript1>,<metadata1_value>,[metadata2_value],...
<transcript2>,<metadata1_value>,[metadata2_value],...
```

Any number of comma-separated metadata labels can be added, and the metadata columns will be added to the output AIRR file header in the order given.

A sample metadata file would be:
```
num_metadata,bool_metadata
AAAAAAAAAAAAAAAA,9,TRUE
AAAAAAAAAAAAAATT,5,FALSE
AAAAAAAAAAAGCGCG,10,TRUE
AAAAAAAAAAGTCATC,2, TRUE
...
```
# License
This program is licensed under the Apache License Version 2.0. A copy of the Apache 2.0 license can be found [here](https://github.com/AvivBenchorin/10X_VDJ_Aggregator/blob/main/LICENSE).

