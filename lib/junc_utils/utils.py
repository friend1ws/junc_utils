#! /usr/bin/env python

import pysam, tabix


def proc_star_junction(input_file, output_file, control_file, read_num_thres, overhang_thres, remove_annotated, convert_map_splice2):
    
    is_control = True if control_file is not None else False
    if is_control: control_db = pysam.TabixFile(control_file)

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

            if F[1].startswith("2959542"):
                pass

            ##########
            # remove control files
            if is_control:
                tabixErrorFlag = 0
                try:
                    records = control_db.fetch(F[0], int(F[1]) - 5, int(F[1]) + 5)
                except Exception as inst:
                    print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
                    tabixErrorMsg = str(inst.args)
                    tabixErrorFlag = 1

                control_flag = 0;
                if tabixErrorFlag == 0:
                    for record_line in records:
                        record = record_line.split('\t')
                        if F[0] == record[0] and F[1] == record[1] and F[2] == record[2]:
                            control_flag = 1

                if control_flag == 1: continue
            ##########

            if convert_map_splice2:
                # convert to map-splice2 coordinate
                F[1] = str(int(F[1]) - 1)
                F[2] = str(int(F[2]) + 1)

            print >> hout, '\t'.join(F)

    hout.close()

