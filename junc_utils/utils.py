#! /usr/bin/env python

from __future__ import print_function
import sys, gzip
import pysam

def convert_anno2vcf(input_file, output_file, reference, header = False):

    """
        convert annovar format to vcf format
        the difficulty is: for the annovar format, the reference base for insertion or deletion is removed,
        and we have to recover them...
    """

    hin = open(input_file, 'r')
    hout = open(output_file, 'w')

    """
    # get current date
    today = datetime.date.today()

    # print meta-information lines
    print >> hout, "##fileformat=#VFCv4.1"
    print >> hout, "##fileDate=" + today.strftime("%Y%m%d")
    # print "##source=ebcallv2.0"
    # print "##INFO=<ID=TD,Number=1,Type=Integer,Description=\"Tumor Depth\">"
    # print "##INFO=<ID=TV,Number=1,Type=Integer,Description=\"Tumor Variant Read Num\">"
    # print "##INFO=<ID=ND,Number=1,Type=Integer,Description=\"Normal Depth\">"
    # print "##INFO=<ID=NV,Number=1,Type=Integer,Description=\"Normal Variant Read Num\">"
    print >> hout, "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic Variation\">"
    """

    if header == True:
        print("##fileformat=VCFv4.1", file = hout)
        print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", file = hout)

    for line in hin:
        F = line.rstrip('\n').split('\t')
        if F[0].startswith('#'): continue
        if F[0] == "Chr" and F[1] == "Start" and F[2] == "End" and F[3] == "Ref" and F[4] == "Alt": continue

        pos, ref, alt = F[1], F[3], F[4]
        
        # insertion
        if F[3] == "-":
            # get the sequence for the reference base
            seq = ""    
            for item in pysam.faidx(reference, F[0] + ":" + str(F[1]) + "-" + str(F[1])):
                seq = seq + item.rstrip('\n')
            seq = seq.replace('>', '')
            seq = seq.replace(F[0] + ":" + str(F[1]) + "-" + str(F[1]), '')
            ref, alt = seq, seq + F[4]


        # deletion
        if F[4] == "-":
            # get the sequence for the reference base
            seq = ""    
            for item in pysam.faidx(reference, F[0] + ":" + str(int(F[1]) - 1) + "-" + str(int(F[1]) - 1)):
                seq = seq + item.rstrip('\n')
            seq = seq.replace('>', '')
            seq = seq.replace(F[0] + ":" + str(int(F[1]) - 1) + "-" + str(int(F[1]) - 1), '')
            pos, ref, alt = str(int(F[1]) - 1), seq + F[3], seq


        # QUAL = int(float(F[15]) * 10)
        QUAL = 60
        # INFO = "TD=" + F[9] + ";TV=" + F[10] + ";ND=" + F[13] + ";NV=" + F[14] + ";SOMATIC"
        INFO = "SOMATIC"

        print(F[0] + "\t" + pos + "\t.\t" + ref + "\t" + alt \
                + "\t" + str(QUAL) + "\t" + "PASS" + "\t" + INFO, file = hout)


    hin.close()
    hout.close()


def remove_vcf_header(input_file, output_file):

    hout = open(output_file, 'w')
    with open(input_file, 'r') as hin:
        for line in hin:
            line = line.rstrip('\n')
            if not line.startswith('#'): print(line, file = hout)

    hout.close()



def convert_genosv2bed(input_file, output_file):

    hout = open(output_file, 'w')
    num = 1
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[0].startswith('#'): continue
            if F[0] == "Chr_1" and F[1] == "Pos_1": continue
            chr1, pos1, dir1, chr2, pos2, dir2, inseq, sv_type = F[0], F[1], F[2], F[3], F[4], F[5], F[6], F[7]

            if chr1 != chr2: continue
            print('\t'.join([chr1, str(int(pos1) - 1), pos2, chr1, pos1, dir1, chr2, pos2, dir2, inseq, sv_type]), file = hout) 
            # start1, end1 = str(int(F[1]) - 1), F[1]
            # start2, end2 = str(int(F[4]) - 1), F[4]
            # dir1, dir2 = F[2], F[5]
            # name = "SV_" + str(num)
            # inseq = F[6] 

            # print >> hout, '\t'.join([chr1, start1, end1, chr2, start2, end2, name, inseq, dir1, dir2])

    hout.close()


def proc_star_junction(input_file, output_file, control_file, read_num_thres, overhang_thres, remove_annotated, convert_map_splice2):
    
    is_control = True if control_file is not None else False
    # if is_control: control_db = pysam.TabixFile(control_file)
    
    control_db = {}
    if is_control:
        with gzip.open(control_file, 'r') as hin:
            for line in hin:
                F = line.rstrip('\n').split('\t')
                key = F[0] + '\t' + F[1] + '\t' + F[2]
                control_db[key] = 1


    if read_num_thres is None: read_num_thres = 0
    if overhang_thres is None: overhang_thres = 0
    if remove_annotated is None: remove_annotated = False
    
    hout = open(output_file, 'w')
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = F[0] + '\t' + F[1] + '\t' + F[2]
            if remove_annotated == True and F[5] != "0": continue
            if int(F[6]) < read_num_thres: continue
            if int(F[8]) < overhang_thres: continue


            if key in control_db: continue
            """
            ##########
            # remove control files
            if is_control:
                tabixErrorFlag = 0
                try:
                    records = control_db.fetch(F[0], int(F[1]) - 5, int(F[1]) + 5)
                except Exception as inst:
                    # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
                    # tabixErrorMsg = str(inst.args)
                    tabixErrorFlag = 1

                control_flag = 0;
                if tabixErrorFlag == 0:
                    for record_line in records:
                        record = record_line.split('\t')
                        if F[0] == record[0] and F[1] == record[1] and F[2] == record[2]:
                            control_flag = 1

                if control_flag == 1: continue
            ##########
            """

            if convert_map_splice2:
                # convert to map-splice2 coordinate
                F[1] = str(int(F[1]) - 1)
                F[2] = str(int(F[2]) + 1)

            print('\t'.join(F), file = hout)

    hout.close()

