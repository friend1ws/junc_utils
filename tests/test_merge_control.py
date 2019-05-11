#! /usr/bin/env python

from __future__ import print_function
import unittest
import os, glob, tempfile, shutil, filecmp
import junc_utils
from check_download import *

class TestMergeControl(unittest.TestCase):

    def setUp(self):
        self.parser = junc_utils.parser.create_parser()


    def tearDown(self):
        pass 


    def test1(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        all_sj_file = glob.glob(cur_dir + "/data/*.out.tab")
        with open(tmp_dir + "/CCLE.SJ_out_tab.list.txt", 'w') as hout:
            for sj_file in sorted(all_sj_file):
                if os.path.basename(sj_file) in ["CCLE-HCC1143-RNA-08.SJ.out.tab", "CCLE-HCC1954-RNA-08.SJ.out.tab", "CCLE-MCF7-RNA-08.SJ.out.tab"]:
                    print(sj_file, file = hout)


        input_list_file = tmp_dir + "/CCLE.SJ_out_tab.list.txt"
        output_file = tmp_dir + "/merge_control.bed.gz"
        answer_file = cur_dir + "/data/merge_control/merge_control.bed.gz"
 
        args = self.parser.parse_args(["merge_control", input_list_file, output_file])
        args.func(args)

        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))

        shutil.rmtree(tmp_dir)

if __name__ == "__main__":
    unittest.main()


