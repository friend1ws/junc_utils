#! /usr/bin/env python

from __future__ import print_function
import sys, re, subprocess
import annot_utils.gene, annot_utils.exon
import pysam

def check_splicing_junction_for_exonization(exonizaiton_junction, output_file, genome_id, is_grc, min_new_intron_size, boundary_margin):

    annot_utils.gene.make_gene_info(output_file + ".tmp.refGene.bed.gz", "refseq", genome_id, is_grc, True)
    annot_utils.exon.make_exon_info(output_file + ".tmp.refExon.bed.gz", "refseq", genome_id, is_grc, True)

    gene_tb = pysam.TabixFile(output_file + ".tmp.refGene.bed.gz")
    exon_tb = pysam.TabixFile(output_file + ".tmp.refExon.bed.gz")

    pos_match = re.match(r'([\w\d]+)\:(\d+)\-(\d+)', exonizaiton_junction)
    if pos_match is None:
        print("exonization_junction is not the right format", file = sys.stderr)
        sys.exit(1)

    schr, sstart, send = pos_match.group(1), int(pos_match.group(2)), int(pos_match.group(3))


    ##########
    is_same_gene, is_boundary, is_intron = False, False, False
    new_exon_info = []
    # check gene annotation for the side 1  
    tabixErrorFlag = 0
    try:
        records = gene_tb.fetch(schr, sstart, send)
    except Exception as inst:
        # print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
        # print >> sys.stderr, '\t'.join(F)
        tabixErrorFlag = 1

    if tabixErrorFlag == 0:
        for record_line in records:
            record = record_line.split('\t')
            if send > int(record[1]) and send < int(record[2]):
                is_same_gene = True
                tgene_chr, tgene_start, tgene_end, tgene, tstrand = record[0], int(record[1]), int(record[2]), record[3], record[5]

                # check exon information
                tabixErrorFlag2 = 0
                try:
                    records2 = exon_tb.fetch(tgene_chr, tgene_start, tgene_end)
                except Exception as inst:
                    tabixErrorFlag2 = 1

                tstarts, tends = [], []
                if tabixErrorFlag2 == 0:
                    for record_line2 in records2:
                        record2 = record_line2.split('\t')
                        if record2[3] == tgene and record2[5] == tstrand:
                            tstarts.append(int(record2[1]) + 1) # the last exonic positions
                            tends.append(int(record2[2])) # the first exonic positions

                    tstarts.sort()
                    tends.sort()

                    for i in range(len(tends) - 1):
                        if abs((sstart - 1) - tends[i]) <= boundary_margin and tstarts[i + 1] - (send + 1) - 1 >= min_new_intron_size:
                            if tstrand == '+':
                                new_exon_info.append((schr, tends[i], tstarts[i + 1], send + 1, '+', "acceptor"))
                            else:
                                new_exon_info.append((schr, tends[i], tstarts[i + 1], send + 1, '-', "donor"))

                        if abs((send + 1) - tstarts[i + 1]) <= boundary_margin and (sstart - 1) - tends[i] - 1 >= min_new_intron_size:
                            if tstrand == '+':
                                new_exon_info.append((schr, tstarts[i + 1], tends[i], sstart - 1, '+', "donor"))
                            else:
                                new_exon_info.append((schr, tstarts[i + 1], tends[i], sstart - 1, '-', "acceptor"))


    new_exon_info = list(set(new_exon_info))

    subprocess.call(["rm","-rf", output_file + ".tmp.refGene.bed.gz"])
    subprocess.call(["rm","-rf", output_file + ".tmp.refGene.bed.gz.tbi"])
    subprocess.call(["rm","-rf", output_file + ".tmp.refExon.bed.gz"])
    subprocess.call(["rm","-rf", output_file + ".tmp.refExon.bed.gz.tbi"])

    return new_exon_info



def check_opposite_junction(junc_file, exonization_info, output_file, control_file, read_num_thres, max_new_exon_size, boundary_margin):

    if control_file is not None:
        control_db = pysam.TabixFile(control_file)

    hout = open(output_file, 'w')

    with open(junc_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[0] != exonization_info[0][0]: continue

            new_exon_flag = 0
            for i in range(len(exonization_info)):
                if (exonization_info[i][4] == '+' and exonization_info[i][5] == "acceptor") or \
                   (exonization_info[i][4] == '-' and exonization_info[i][5] == "donor"):
                    if int(F[1]) - 1 >= exonization_info[i][3] and int(F[1]) - 1 <= exonization_info[i][3] + max_new_exon_size and \
                       abs(int(F[2]) + 1 - exonization_info[i][2]) <= boundary_margin and int(F[6]) >= read_num_thres:
                        new_exon_flag = 1
                else:
                    if int(F[2]) +  1 >= exonization_info[i][3] - max_new_exon_size and int(F[2]) + 1 <= exonization_info[i][3] and \
                        abs(int(F[1]) - 1 - exonization_info[i][2]) <= boundary_margin and int(F[6]) >= read_num_thres:
                        new_exon_flag = 1

            if new_exon_flag == 1:

                control_flag = 0
                if control_file is not None:
                    tabixErrorFlag = 0
                    try:
                        records = control_db.fetch(F[0], int(F[1]) - 5, int(F[1]) + 5)
                    except Exception as inst:
                        # prin >> sys.stderr, "%s: %s" % (type(inst), inst.args)
                        # tabixErrorMsg = str(inst.args)
                        tabixErrorFlag = 1

                    if tabixErrorFlag == 0:
                        for record_line in records:
                            record = record_line.split('\t')
                            if F[0] == record[0] and F[1] == record[1] and F[2] == record[2]:
                                control_flag = 1

                if control_flag == 0:
                    print('\t'.join(F), file = hout)

     
    hout.close()




