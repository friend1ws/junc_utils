#! /usr/bin/env python

import unittest
import os, tempfile, shutil, filecmp
import junc_utils

class TestAnnotate(unittest.TestCase):

    def setUp(self):
        self.parser = junc_utils.parser.create_parser()


    def tearDown(self):
        pass 


    def test1(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        input_file = cur_dir + "/data/filter/CCLE-HCC1143-RNA-08.SJ.out.filter.tab"
        output_file = tmp_dir + "/CCLE-HCC1143-RNA-08.SJ.out.filter.annotate.tab"
        answer_file = cur_dir + "/data/annotate/CCLE-HCC1143-RNA-08.SJ.out.filter.annotate.tab"
 
        args = self.parser.parse_args(["annotate", input_file, output_file, "--grc"])
        args.func(args)

        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))

        shutil.rmtree(tmp_dir)


    def test2(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        input_file = cur_dir + "/data/filter/CCLE-GSU-RNA-08.SJ.out.filter.tab"
        output_file = tmp_dir + "/CCLE-GSU-RNA-08.SJ.out.filter.annotate.tab"
        answer_file = cur_dir + "/data/annotate/CCLE-GSU-RNA-08.SJ.out.filter.annotate.tab"

        args = self.parser.parse_args(["annotate", input_file, output_file, "--grc"])
        args.func(args)

        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))

        shutil.rmtree(tmp_dir)


if __name__ == "__main__":
    unittest.main()


