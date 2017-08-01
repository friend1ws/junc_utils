#! /usr/bin/env python

import sys, re, subprocess
import annot_utils.gene, annot_utils.exon
import pysam

def check_splicing_junction_for_exonization(exonizaiton_junction, output_file, genome_id, is_grc):

    annot_utils.gene.make_gene_info(output_file + ".tmp.refGene.bed.gz", "refseq", genome_id, is_grc, True)
    annot_utils.exon.make_exon_info(output_file + ".tmp.refExon.bed.gz", "refseq", genome_id, is_grc, True)

    gene_tb = pysam.TabixFile(output_file + ".tmp.refGene.bed.gz")
    exon_tb = pysam.TabixFile(output_file + ".tmp.refExon.bed.gz")

    pos_match = re.match(r'([\w\d]+)\:(\d+)\-(\d+)', exonizaiton_junction)
    if pos_match is None:
        print >> sys.stderr, "exonization_junction is not the right format"
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
                            tstarts.append(int(record2[1]))
                            tends.append(int(record2[2]))

                    tstarts.sort()
                    tends.sort()

                    for i in range(len(tends) - 1):
                        if abs(sstart - tends[i]) < 5 and send < tstarts[i + 1] - 100:
                            if tstrand == '+':
                                new_exon_info.append((schr, tends[i], tstarts[i + 1], send, '+', "acceptor"))
                            else:
                                new_exon_info.append((schr, tends[i], tstarts[i + 1], send, '-', "donor"))

                        if abs(send - tstarts[i + 1]) < 5 and sstart > tends[i] + 100:
                            if tstrand == '+':
                                new_exon_info.append((schr, tstarts[i + 1], tends[i], sstart, '+', "donor"))
                            else:
                                new_exon_info.append((schr, tstarts[i + 1], tends[i], sstart, '-', "acceptor"))


    new_exon_info = list(set(new_exon_info))

    subprocess.call(["rm","-rf", output_file + ".tmp.refGene.bed.gz"])
    subprocess.call(["rm","-rf", output_file + ".tmp.refGene.bed.gz.tbi"])
    subprocess.call(["rm","-rf", output_file + ".tmp.refExon.bed.gz"])
    subprocess.call(["rm","-rf", output_file + ".tmp.refExon.bed.gz.tbi"])

    return new_exon_info



def check_opposite_junction(junc_file, exonization_info, output_file):

    with open(junc_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[0] != exonization_info[0][0]: continue

            for i in range(len(exonization_info)):
                if (exonization_info[i][4] == '+' and exonization_info[i][5] == "acceptor") or \
                   (exonization_info[i][4] == '-' and exonization_info[i][5] == "donor"):
                    if int(F[1]) > exonization_info[i][3] - 10 and int(F[1]) < exonization_info[i][3] + 500 and \
                       abs(int(F[2]) - exonization_info[i][2]) < 5 and int(F[6]) >= 2:
                        print '\t'.join(F)
                else:
                    if int(F[2]) > exonization_info[i][3] - 500 and int(F[2]) < exonization_info[i][3] + 10 and \
                        abs(int(F[1]) - exonization_info[i][2]) < 5 and int(F[6]) >= 2:
                        print '\t'.join(F)
 



