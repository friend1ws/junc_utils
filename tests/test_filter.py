#! /usr/bin/env python

import unittest
import os, tempfile, shutil, filecmp
import junc_utils

class TestFilter(unittest.TestCase):

    def setUp(self):
        self.parser = junc_utils.parser.create_parser()


    def tearDown(self):
        pass 


    def test1(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        input_file = cur_dir + "/data/CCLE-HCC1143-RNA-08.SJ.out.tab"
        output_file = tmp_dir + "/CCLE-HCC1143-RNA-08.SJ.out.filter.tab"
        answer_file = cur_dir + "/data/filter/CCLE-HCC1143-RNA-08.SJ.out.filter.tab"
 
        args = self.parser.parse_args(["filter", input_file, output_file])
        args.func(args)

        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))

        shutil.rmtree(tmp_dir)


    def test2(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        input_file = cur_dir + "/data/CCLE-GSU-RNA-08.SJ.out.tab"
        output_file = tmp_dir + "/CCLE-GSU-RNA-08.SJ.out.filter.tab"
        answer_file = cur_dir + "/data/filter/CCLE-GSU-RNA-08.SJ.out.filter.tab"

        args = self.parser.parse_args(["filter", input_file, output_file])
        args.func(args)

        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))

        shutil.rmtree(tmp_dir)



if __name__ == "__main__":
    unittest.main()


