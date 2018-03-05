#! /usr/bin/env python

from run import *
import argparse
from version import __version__

def create_parser():

    parser = argparse.ArgumentParser(prog = "junc_utils")

    parser.add_argument("--version", action = "version", version = "%(prog)s " + __version__)

    subparsers = parser.add_subparsers()

    ##########
    # filter

    filter = subparsers.add_parser("filter",
                                   help = "filter out splicing junctions outside specified conditions")

    filter.add_argument("junc_file", metavar = "sample.SJ.out.tab", default = None, type = str,
                        help = "the path to star junction file")

    filter.add_argument("output_path", metavar = "output.txt", default = None, type = str,
                        help = "the path to the output file")

    filter.add_argument("--read_num_thres", type = int, default = 3,
                        help = "remove splicing junctions whose supporting numbers are below this value (default: %(default)s)")
                               
    filter.add_argument("--overhang_thres", type = int, default = 10,
                        help = "remove splicing junctions whose overhang sizes are this value (default: %(default)s)")
                               
    filter.add_argument("--keep_annotated", default = False, action = 'store_true',
                        help = "do not remove annotated splicing junctions")

    filter.add_argument("--pooled_control_file", default = None, type = str,
                        help = "the path to control data created by merge_control (default: %(default)s)")

    filter.set_defaults(func = filter_main)

    ##########
    # annotate

    annotate = subparsers.add_parser("annotate",
                                     help = "annotate junctions")

    annotate.add_argument("junc_file", metavar = "sample.SJ.out.tab", default = None, type = str,
                          help = "the path to star splicing junction file")

    annotate.add_argument("output_path", metavar = "output.txt", default = None, type = str,
                          help = "the path to the output file")

    # annotate.add_argument("annotation_dir", metavar = "annotation_dir", default = None, type = str,
    #                       help = "the path to the database directory")

    annotate.add_argument("--grc", default = False, action = 'store_true',
                          help = "convert chromosome names to Genome Reference Consortium nomenclature (default: %(default)s)")

    annotate.add_argument("--genome_id", choices = ["hg19", "hg38", "mm10"], default = "hg19",
                          help = "the genome id used for selecting UCSC-GRC chromosome name corresponding files (default: %(default)s)")

    annotate.add_argument("--junction_margin", type = int, default = 3,
                          help = "the allowed margin size for the difference between splicing breakpoints and exon-intron junctions")

    annotate.add_argument("--exon_margin", type = int, default = 30,
                          help = "the allowed margin size for the difference between splicing breakpoints and exonic region")

    annotate.set_defaults(func = annotate_main)

    ##########
    # associate

    associate = subparsers.add_parser("associate",
                                     help = "associate junctions with mutations or SVs")

    associate.add_argument("annotated_junction_file", metavar = "annotated_junction.SJ.out.tab", default = None, type = str,
                           help = "the path to the annotated splicing junction file")

    associate.add_argument("mutation_file", metavar = "mutation.vcf.gz", default = None, type = str,
                           help = "the path to the mutation file")

    associate.add_argument("output_file", metavar = "output_file", default = None, type = str, 
                           help = "the path to the output")

    # associate.add_argument("annotation_dir", metavar = "annotation_dir", type = str,
    #                        help = "the path to database directory")

    associate.add_argument("--grc", default = False, action = 'store_true',
                           help = "convert chromosome names to Genome Reference Consortium nomenclature (default: %(default)s)")

    associate.add_argument("--genome_id", choices = ["hg19", "hg38", "mm10"], default = "hg19",
                           help = "the genome id used for selecting UCSC-GRC chromosome name corresponding files (default: %(default)s)")

    associate.add_argument("--donor_size", metavar = "donor_size", default = "3,6", type = str,
                           help = "splicing donor site size (exonic region size, intronic region size) (default: %(default)s)")

    associate.add_argument("--acceptor_size", metavar = "acceptor_size", default = "6,1", type = str,
                           help = "splicing donor site size (intronic region size, exonic region size) (default: %(default)s)")

    associate.add_argument("--reference", metavar = "reference.fa", type = str,
                           help = "the path to the reference genomoe sequence (necessary for anno format)")

    associate.add_argument("--debug", default = False, action = 'store_true', help = "keep intermediate files")

    associate.add_argument("--mutation_format", choices=["vcf", "anno"], default = "vcf",
                           help = "Deprecated. This is not used any more. The format of mutation file vcf or annovar (tsv) format")

    # associate.add_argument("--is_edit_dist", default = False, action = 'store_true', 
    #                        help = "for SNV, only extract those changing edit distance to splicing motifs in consistent manners")

    # associate.add_argument("--skip_creation_indel", default = False, action = 'store_true',
    #                        help = "do not consider indels around splicing junctions (putative splicing motif creatin indels)")

    associate.add_argument('--only_dist', action='store_true',
                           help = "associate all mutations within the junctions plus specified margin (this is mainly for evaluation)")

    associate.add_argument("--only_dist_search_margin", metavar = "only_dist_search_margin", default = 100, type = int,
                           help = "margin size used when only_dist option is active (default: %(default)s)")

    associate.add_argument('--sv', action='store_true',
                           help = "analysis structural variation file")

    associate.add_argument('--branchpoint', action='store_true',
                           help = "include branch points to the association analysis")

    associate.add_argument("--branchpoint_size", metavar = "branchpoint_size", default = "3,2", type = str,
                           help = "splicing branchpoint size (intronic-side region size, exonic-side region size) (default: %(default)s)")

    associate.set_defaults(func = associate_main)

    ##########
    # merge control
    merge_control = subparsers.add_parser("merge_control",
                                          help = "merge, compress and index the splicing junction list (.SJ.out.tab file generated by STAR)")

    merge_control.add_argument("junc_list", metavar = "junc_list.txt", default = None, type = str,
                               help = "junction path list")

    merge_control.add_argument("output_path", default = None, type = str,
                               help = "the path of the output file")


    merge_control.add_argument("--read_num_thres", type = int, default = 3,
                               help = "remove splicing junctions whose supporting numbers are below this value (default: %(default)s)")

    merge_control.add_argument("--overhang_thres", type = int, default = 10,
                               help = "remove splicing junctions whose overhang sizes are this value (default: %(default)s)")

    merge_control.add_argument("--keep_annotated", default = False, action = 'store_true', 
                               help = "do not remove annotated splicing junctions")

    merge_control.add_argument("--sample_num_thres", type = int, default = 2,
                               help = "register splicing junctions at least shared by specified number of samples (default: %(default)s)")

    merge_control.set_defaults(func = merge_control_main)

    """
    # temporary remove this function, but want to recover in near future (2018/02/02)
    ##########
    # exonization_pair
    exonization_pair = subparsers.add_parser("exonization_pair",
                                          help = "get candidates for the other junction of an exonization supporting junction pair")

    exonization_pair.add_argument("half_exonizaiton_junction", metavar = "half_exonizaiton_junction", default = None, type = str,
                                  help = "a half of exonization junction pairs")

    exonization_pair.add_argument("junc_file", metavar = "sample.SJ.out.tab", default = None, type = str,
                                  help = "the path to star splicing junction file")

    exonization_pair.add_argument("output_file", metavar = "output_file", default = None, type = str,
                                  help = "the path to the output")

    exonization_pair.add_argument("--control_file", default = None, type = str,
                                  help = "the path to control data created by junc_utils merge_control (default: %(default)s)")

    exonization_pair.add_argument("--grc", default = False, action = 'store_true',
                                  help = "convert chromosome names to Genome Reference Consortium nomenclature (default: %(default)s)")

    exonization_pair.add_argument("--genome_id", choices = ["hg19", "hg38", "mm10"], default = "hg19",
                                  help = "the genome id used for selecting UCSC-GRC chromosome name corresponding files (default: %(default)s)")

    exonization_pair.add_argument("--read_num_thres", type = int, default = 2,
                                  help = "remove splicing junctions whose supporting numbers are below this value (default: %(default)s)")

    exonization_pair.add_argument("--min_new_intron_size", type = int, default = 30,
                                  help = "minimum new intron size")

    exonization_pair.add_argument("--max_new_exon_size", type = int, default = 500,
                                  help = "maximum new exon size")

    exonization_pair.add_argument("--boundary_margin", type = int, default = 3,
                                  help = "margin size when comparing splicing junctions with exon-intron boundaries")

    exonization_pair.set_defaults(func = exonization_pair_main)
    ##########
    """

    return parser

