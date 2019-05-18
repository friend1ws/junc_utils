# junc_utils

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Build Status](https://travis-ci.org/friend1ws/junc_utils.svg?branch=devel)](https://travis-ci.org/friend1ws/junc_utils)

Utility functions for analyzing splicing junction (generated by STAR, .SJ.out.tab files)

## Dependency

### Python
`pysam`, [`annot_utils>=0.3.0`](https://github.com/friend1ws/annot_utils) packages.

### Software
[htslib](http://www.htslib.org), [bedtools](http://bedtools.readthedocs.io/en/latest/)

## Install

```
pip install junc_utils
```

Alternatively, you can install from the source code
```
wget https://github.com/friend1ws/junc_utils/archive/v0.5.0.tar.gz
tar zxvf v0.5.0.tar.gz
cd junc_utils-0.5.0
python setup.py build install
```

This package has been tested on Python 2.7, 3.5, 3.6.

## Commands

### fitler

filter out splicing junctions outside specified conditions
```
junc_utils filter [-h] [--read_num_thres READ_NUM_THRES]
                         [--overhang_thres OVERHANG_THRES] [--keep_annotated]
                         [--pooled_control_file POOLED_CONTROL_FILE]
                         sample.SJ.out.tab output.txt
```

### annotate

annotate splicing junctions

```
junc_utils annotate [-h] [--grc] [--genome_id {hg19,hg38,mm10}]
                           [--junction_margin JUNCTION_MARGIN]
                           [--exon_margin EXON_MARGIN]
                           sample.SJ.out.tab output.txt
```

### merge_control

merge, compress and index the splicing junction list

```
junc_utils merge_control [-h] [--read_num_thres READ_NUM_THRES]
                                [--overhang_thres OVERHANG_THRES]
                                [--keep_annotated]
                                [--sample_num_thres SAMPLE_NUM_THRES]
                                junc_list.txt output_path
```

### associate

associate junctions with mutations or SVs

```
junc_utils associate [-h] [--grc] [--genome_id {hg19,hg38,mm10}]
                            [--donor_size donor_size]
                            [--acceptor_size acceptor_size]
                            [--reference reference.fa] [--debug]
                            [--mutation_format {vcf,anno}] [--only_dist]
                            [--only_dist_search_margin only_dist_search_margin]
                            [--sv] [--branchpoint]
                            [--branchpoint_size branchpoint_size]
                            annotated_junction.SJ.out.tab mutation.vcf.gz
                            output_file
```

