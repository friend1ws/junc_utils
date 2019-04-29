#! /usr/bin/env python

import unittest
import os, tempfile, shutil, filecmp
import junc_utils
from check_download import *

class TestAssociate(unittest.TestCase):

    def setUp(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        check_download("https://storage.googleapis.com/friend1ws_package_data/common/GRCh37.fa", \
                       cur_dir + "/resource/reference_genome/GRCh37.fa")

        self.parser = junc_utils.parser.create_parser()


    def tearDown(self):
        pass 


    def test1(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        annotated_sj_file = cur_dir + "/data/annotate/CCLE-HCC1143-RNA-08.SJ.out.filter.annotate.tab"
        mutation_file = cur_dir + "/data/mutation/CCLE-HCC1143-DNA-08.genomon_mutation.result.txt" 
        output_file = tmp_dir + "/CCLE-HCC1143-RNA-08.SJ.out.filter.annotate.associate.tab"
        ref_genome = cur_dir + "/resource/reference_genome/GRCh37.fa"
        answer_file = cur_dir + "/data/associate/CCLE-HCC1143-RNA-08.SJ.out.filter.annotate.associate.tab"

        args = self.parser.parse_args(["associate", annotated_sj_file, mutation_file, output_file, \
                                       "--reference", ref_genome])
        args.func(args)

        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))

        shutil.rmtree(tmp_dir)


    def test2(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        annotated_sj_file = cur_dir + "/data/annotate/CCLE-HCC1143-RNA-08.SJ.out.filter.annotate.tab"
        mutation_file = cur_dir + "/data/mutation/CCLE-HCC1143-DNA-08.genomon_mutation.result.vcf"
        output_file = tmp_dir + "/CCLE-HCC1143-RNA-08.SJ.out.filter.annotate.associate.tab"
        answer_file = cur_dir + "/data/associate/CCLE-HCC1143-RNA-08.SJ.out.filter.annotate.associate.tab"

        args = self.parser.parse_args(["associate", annotated_sj_file, mutation_file, output_file])
        args.func(args)

        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))

        shutil.rmtree(tmp_dir)


    def test3(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        annotated_sj_file = cur_dir + "/data/annotate/CCLE-GSU-RNA-08.SJ.out.filter.annotate.tab"
        mutation_file = cur_dir + "/data/sv/CCLE-GSU-DNA-08.genomonSV.result.filt.txt"
        output_file = tmp_dir + "/CCLE-GSU-DNA-08.SJ.out.filter.annotate.associate.tab"
        # ref_genome = cur_dir + "/resource/reference_genome/GRCh37.fa"
        answer_file = cur_dir + "/data/associate/CCLE-GSU-RNA-08.SJ.out.filter.annotate.associate.tab"

        args = self.parser.parse_args(["associate", annotated_sj_file, mutation_file, output_file, \
                                       "--sv"])
        args.func(args)

        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))

        shutil.rmtree(tmp_dir)


if __name__ == "__main__":
    unittest.main()


