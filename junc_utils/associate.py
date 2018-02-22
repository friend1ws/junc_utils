#! /usr/bin/env python

import sys, subprocess, pkg_resources
import pysam
import annot_utils.exon


def get_snv_junction2(input_file, output_file, mutation_file, donor_size, acceptor_size, 
                      genome_id, is_grc, is_branchpoint, branch_size):

    """
        a script for detecting candidate somatic substitutions causing splicing changes
        
        1. exon-skip:
           exon-intron junction region within spliced region
        
        2. splice-site-slip, pseudo-exon-inclusion
           in addition to exon-intron junction region, 
           we check for the region around non exon-intron junction break points 
        Maybe, we should check the exon-intron junction region within, 
        e.g., +-10000bp from the spliced region to search for unknown phenomena.
    """

    donor_size_exon, donor_size_intron = [int(x) for x in donor_size.split(',')]
    acceptor_size_intron, acceptor_size_exon = [int(x) for x in acceptor_size.split(',')]
    branch_size_intron, branch_size_exon = [int(x) for x in branch_size.split(',')]

    searchMargin1 = 100 
    searchMargin2 = 10

    splicingDonnorMotif = ["AG", "GTRAGT"]
    splicingAcceptorMotif = ["YYYYNCAG", "G"]

    annot_utils.exon.make_exon_info(output_file + ".tmp.refExon.bed.gz", "refseq", genome_id, is_grc, True)

    hout = open(output_file, 'w')

    mutation_tb = pysam.TabixFile(mutation_file)
    exon_tb = pysam.TabixFile(output_file + ".tmp.refExon.bed.gz")
    if is_branchpoint:
        if genome_id != "hg19":
            print >> sys.stderr, "branchpoint can be used only for hg19"
            sys.exit(1)
        if is_grc:
            branchpoint_bed = pkg_resources.resource_filename("annot_utils", "data/hg19/branchpoint_signal.grc.bed.gz")
        else:
            branchpoint_bed = pkg_resources.resource_filename("annot_utils", "data/hg19/branchpoint_signal.bed.gz")
        branch_tb = pysam.TabixFile(branchpoint_bed)

    header2ind = {}
    with open(input_file, 'r') as hin:
        
        # read header
        header = hin.readline().rstrip('\n').split('\t')
        for i in range(len(header)):
            header2ind[header[i]] = i

        # print header
        print >> hout, '\t'.join(header) + '\t' + "Mutation_Key" + '\t' + "Motif_Pos" + '\t' + "Mutation_Type" + '\t' + "Is_Canonical"

        for line in hin:
            F = line.rstrip('\n').split('\t') 

            sj_start = int(F[header2ind["SJ_2"]]) - 1
            sj_end = int(F[header2ind["SJ_3"]]) + 1

            if F[header2ind["Splicing_Class"]] not in ["exon-skip", "alternative-3'-splice-site", "alternative-5'-splice-site",
                                                       "intronic-alternative-3'-splice-site", "intronic-alternative-5'-splice-site"]: continue

            firstSearchRegion = [F[header2ind["SJ_1"]], sj_start, sj_end]
            splicingMotifRegions = []
            targetGene  =[]

            # we need to detect the non exon-intron junction break points
            # current procedure may be not perfect and be subject to change..

            gene1 = F[header2ind["Gene_1"]].split(';')
            gene2 = F[header2ind["Gene_2"]].split(';')
            junction1 = F[header2ind["Is_Boundary_1"]].split(';')   
            junction2 = F[header2ind["Is_Boundary_2"]].split(';')
            offset1 = F[header2ind["Offset_1"]].split(';')
            offset2 = F[header2ind["Offset_2"]].split(';')

            # just consider genes sharing the exon-intron junction with the breakpoints of splicings
            for i in range(0, len(gene1)):
                if junction1[i] != "*": targetGene.append(gene1[i])
            for i in range(0, len(gene2)):
                if junction2[i] != "*": targetGene.append(gene2[i])
            targetGene = list(set(targetGene))


            if F[header2ind["Splicing_Class"]] in ["alternative-3'-splice-site", "alternative-5'-splice-site",
                                                   "intronic-alternative-3'-splice-site", "intronic-alternative-5'-splice-site"]:

                # for non exon-intron junction breakpoints
                for i in range(0, len(gene1)):

                    if junction1[i] == "*" and junction2[i] == "s": # splicing donor motif, plus direction
                        firstSearchRegion[1] = sj_start - searchMargin1
                        motif_start = sj_start - donor_size_exon + 1 - int(offset2[i])
                        motif_end = sj_start + donor_size_intron - int(offset2[i])
                        splicingMotifRegions.append((F[header2ind["SJ_1"]], motif_start, motif_end, "donor", "+", 0))
                    if junction1[i] == "*" and junction2[i] == "e": # splicing acceptor motif, minus direction
                        firstSearchRegion[1] = sj_start - searchMargin1 
                        motif_start = sj_start - acceptor_size_exon + 1 - int(offset2[i])
                        motif_end = sj_start + acceptor_size_intron - int(offset2[i])
                        splicingMotifRegions.append((F[header2ind["SJ_1"]], motif_start, motif_end, "acceptor", "-", 0))
                    if junction1[i] == "s" and junction2[i] == "*": # splicing donor motif, minus direction
                        firstSearchRegion[2] = sj_end + searchMargin1
                        motif_start = sj_end - donor_size_intron - int(offset1[i])
                        motif_end = sj_end + donor_size_exon - 1 - int(offset1[i])
                        splicingMotifRegions.append((F[header2ind["SJ_1"]], motif_start, motif_end, "donor", "-", 0))
                    if junction1[i] == "e" and junction2[i] == "*": # # splicing acceptor motif, plus direction
                        firstSearchRegion[2] = sj_end + searchMargin1
                        motif_start = sj_end - acceptor_size_intron - int(offset1[i])
                        motif_end = sj_end + acceptor_size_exon - 1 - int(offset1[i])
                        splicingMotifRegions.append((F[header2ind["SJ_1"]], motif_start, motif_end, "acceptor", "+", 0))


            ##########
            # rough check for the mutation between the spliced region
            tabixErrorFlag1 = 0
            try:
                mutation_lines = mutation_tb.fetch(firstSearchRegion[0], firstSearchRegion[1], firstSearchRegion[2])
            except Exception as inst:
                # print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
                # print >> sys.stderr, '\t'.join(F)
                tabixErrorFlag1 = 1

            # if there are some mutaions
            if tabixErrorFlag1 == 0 and mutation_lines is not None:

                # chr_name = grch2ucsc[F[header2ind["SJ_1"]]] if F[header2ind["SJ_1"]] in grch2ucsc else F[header2ind["SJ_1"]] 
                chr_name = F[header2ind["SJ_1"]] 
                # check the exons within the spliced regions
                tabixErrorFlag2 = 0
                try:
                    exon_lines = exon_tb.fetch(chr_name, firstSearchRegion[1], firstSearchRegion[2])
                except Exception as inst:
                    # print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
                    # print >> sys.stderr, '\t'.join(F)
                    tabixErrorFlag2 = 1

                # first, add the exon-intron junction for detailed check region list
                if tabixErrorFlag2 == 0:
                    for exon_line in exon_lines:
                        exon = exon_line.split('\t')
                        if exon[3] not in targetGene: continue
                        if exon[5] == "+":
                            # splicing acceptor motif, plus direction
                            splicingMotifRegions.append((exon[0], int(exon[1]) - acceptor_size_intron + 1, int(exon[1]) + acceptor_size_exon, "acceptor", "+", 1))
                            # splicing donor motif, plus direction
                            splicingMotifRegions.append((exon[0], int(exon[2]) - donor_size_exon + 1, int(exon[2]) + donor_size_intron, "donor", "+", 1))
                        if exon[5] == "-":
                            # splicing donor motif, minus direction 
                            splicingMotifRegions.append((exon[0], int(exon[1]) - donor_size_intron + 1, int(exon[1]) + donor_size_exon, "donor", "-", 1))
                            # splicing acceptor motif, minus direction
                            splicingMotifRegions.append((exon[0], int(exon[2]) - acceptor_size_exon + 1, int(exon[2]) + acceptor_size_intron, "acceptor", "-", 1))


                # check branchpoint within the spliced regions
                if is_branchpoint:
                    tabixErrorFlag3 = 0
                    try:
                        branche_lines = branch_tb.fetch(chr_name, firstSearchRegion[1], firstSearchRegion[2])
                    except Exception as inst:
                        print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
                        print >> sys.stderr, '\t'.join(F)
                        tabixErrorFlag3 = 3

                    if tabixErrorFlag3 == 0:
                        for branch_line in branche_lines:
                            branch = branch_line.split('\t')
                            if branch[5] == "+":
                                splicingMotifRegions.append((branch[0], int(branch[2]) - branch_size_intron, int(branch[2]) + branch_size_exon - 1, "branchpoint", "+", 1))
                            else:
                                splicingMotifRegions.append((branch[0], int(branch[2]) - branch_size_exon + 1, int(branch[2]) + branch_size_intron, "branchpoint", "-", 1))


                splicingMotifRegions = list(set(splicingMotifRegions))

                # compare each mutation with exon-intron junction regions and non-exon-intorn junction breakpoints.
                for mutation_line in mutation_lines:
                    mutation = mutation_line.split('\t')
                    RegMut = []
                    for reg in splicingMotifRegions:

                        # insertion or deletion (just consider the disruption of splicing motifs now)
                        # if (len(mutation[3]) > 1 or len(mutation[4]) > 1) and reg[5] == 1:
                        if (len(mutation[3]) > 1 or len(mutation[4]) > 1):

                            indel_start = int(mutation[1]) + 1
                            indel_end = int(mutation[1]) + len(mutation[3]) - 1 if len(mutation[3]) > 1 else indel_start
                            if indel_start <= reg[2] and reg[1] <= indel_end:
                           
                                if reg[3] in ["acceptor", "donor"]: 
                                    if reg[3] == "acceptor" and reg[4] == "+": canonical_start_pos = reg[2] - acceptor_size_exon - 1
                                    if reg[3] == "acceptor" and reg[4] == "-": canonical_start_pos = reg[1] + acceptor_size_exon 
                                    if reg[3] == "donor" and reg[4] == "+": canonical_start_pos = reg[1] + donor_size_exon 
                                    if reg[3] == "donor" and reg[4] == "-": canonical_start_pos = reg[2] - donor_size_exon - 1

                                    is_canonical = "non-canonical"
                                    if len(mutation[3]) > 1: # deletion
                                        if indel_start <= canonical_start_pos + 1 and canonical_start_pos <= indel_end:
                                            is_canonical = "canonical"
                                    else: # insertion
                                        if indel_start == canonical_start_pos + 1:
                                            is_canonical = "canonical"

                                # branch point
                                else:
                                    if reg[4] == "+": canonical_start_pos = reg[2] - branch_size_exon + 1
                                    if reg[4] == "-": canonical_start_pos = reg[1] + branch_size_exon - 1
 
                                    is_canonical = "non-canonical"
                                    if len(mutation[3]) > 1: # deletion
                                        if indel_start <= canonical_start_pos and canonical_start_pos <= indel_end:
                                            is_canonical = "canonical"
                                    else: # insertion
                                        pass





                                if reg[5] == 0: RegMut.append([reg, "splicing " + reg[3] + " creation", is_canonical])
                                if reg[5] == 1: RegMut.append([reg, "splicing " + reg[3] + " disruption", is_canonical])
 
                        # base substitution
                        if len(mutation[3]) == 1 and len(mutation[4]) == 1 and reg[1] <= int(mutation[1]) <= reg[2]:
                            motifSeq = ""
                            if reg[3] == "acceptor": motifSeq = splicingAcceptorMotif[0] + splicingAcceptorMotif[1]                  
                            if reg[3] == "donor": motifSeq = splicingDonnorMotif[0] + splicingDonnorMotif[1]

                            is_canonical = "non-canonical"
                            if reg[3] in ["acceptor", "donor"]:
                                if reg[3] == "acceptor" and reg[4] == "+": canonical_start_pos = reg[2] - acceptor_size_exon - 1
                                if reg[3] == "acceptor" and reg[4] == "-": canonical_start_pos = reg[1] + acceptor_size_exon 
                                if reg[3] == "donor" and reg[4] == "+": canonical_start_pos = reg[1] + donor_size_exon 
                                if reg[3] == "donor" and reg[4] == "-": canonical_start_pos = reg[2] - donor_size_exon - 1

                                is_canonical = "canonical" if canonical_start_pos <= int(mutation[1]) <= canonical_start_pos + 1 else "non-canonical"

                            # branchpoint 
                            else:
                                if reg[4] == "+": canonical_start_pos = reg[2] - branch_size_exon + 1 
                                if reg[4] == "-": canonical_start_pos = reg[1] + branch_size_exon - 1 

                                is_canonical = "canonical" if canonical_start_pos == int(mutation[1]) else "non-canonical"


                            if reg[5] == 0: RegMut.append([reg, "splicing " + reg[3] + " creation", is_canonical])
                            if reg[5] == 1: RegMut.append([reg, "splicing " + reg[3] + " disruption", is_canonical])

                    for item in RegMut:
                        mut_print_str = ','.join([mutation[i] for i in [0, 1, 3, 4]])
                        print >> hout, '\t'.join(F) + '\t' + mut_print_str + '\t' + F[header2ind["SJ_1"]] + ':' + str(item[0][1]) + '-' + str(item[0][2]) + ',' + item[0][4] + '\t' + item[1] + '\t' + item[2]


    hout.close()

    subprocess.call(["rm", "-rf", output_file + ".tmp.refExon.bed.gz"])
    subprocess.call(["rm", "-rf", output_file + ".tmp.refExon.bed.gz.tbi"])



def get_snv_junction_only_dist(input_file, output_file, mutation_file, annotation_dir, searchMargin):

    ref_exon_bed = annotation_dir + "/refExon.bed.gz"
    grch2ucsc_file = annotation_dir + "/grch2ucsc.txt"

    # relationship between CRCh and UCSC chromosome names
    grch2ucsc = {}
    with open(grch2ucsc_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            grch2ucsc[F[0]] = F[1]

    hout = open(output_file, 'w')
    mutation_tb = tabix.open(mutation_file)

    header2ind = {}
    with open(input_file, 'r') as hin:
        
        # read header
        header = hin.readline().rstrip('\n').split('\t')
        for i in range(len(header)):
            header2ind[header[i]] = i

        # print header
        print >> hout, '\t'.join(header) + '\t' + "Mutation_Key" + '\t' + "Motif_Pos" + '\t' + "Mutation_Type" + '\t' + "Is_Canonical"

        for line in hin:
            F = line.rstrip('\n').split('\t') 

            sj_start = int(F[header2ind["SJ_2"]]) - 1
            sj_end = int(F[header2ind["SJ_3"]]) + 1

            if F[header2ind["Splicing_Class"]] not in ["exon-skip", "alternative-3'-splice-site", "alternative-5'-splice-site",
                                                       "intronic-alternative-3'-splice-site", "intronic-alternative-5'-splice-site"]: continue

            firstSearchRegion = [F[header2ind["SJ_1"]], sj_start, sj_end]

            ##########
            # rough check for the mutation between the spliced region
            tabixErrorFlag1 = 0
            try:
                mutations = mutation_tb.query(firstSearchRegion[0], firstSearchRegion[1] - searchMargin, firstSearchRegion[2] + searchMargin)
            except Exception as inst:
                # print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
                # print >> sys.stderr, '\t'.join(F)
                tabixErrorFlag1 = 1

            # if there are some mutaions
            if tabixErrorFlag1 == 0 and mutations is not None:

                for mutation in mutations:

                    mut_print_str = ','.join([mutation[i] for i in [0, 1, 3, 4]])
                    print >> hout, '\t'.join(F) + '\t' + mut_print_str + '\t' + "---,---" + '\t' + "---" + '\t' + "---"


    hout.close()

def get_sv_junction(input_file, output_file, sv_file, genome_id, is_grc):
# def get_sv_junction(input_file, output_file, mutation_file, annotation_dir):

    from collections import namedtuple
    SV = namedtuple('SV', ('chr1', 'pos1', 'dir1', 'chr2', 'pos2', 'dir2', 'inseq', 'sv_type'))

    """
        a script for detecting candidate structural variations (currently just deletions causing splicing changes)
    """

    sv_comp_margin = 20
    exon_comp_margin = 10


    hout = open(output_file, 'w')

    sv_tb = pysam.TabixFile(sv_file)
    # annot_utils.exon.make_exon_info(output_file + ".tmp.refExon.bed.gz", "refseq", genome_id, is_grc, True)
    # exon_tb = pysam.TabixFile(output_file + ".tmp.refExon.bed.gz")

    header2ind = {}
    with open(input_file, 'r') as hin:

        # read header
        header = hin.readline().rstrip('\n').split('\t')
        for i in range(len(header)):
            header2ind[header[i]] = i

        # print header
        print >> hout, '\t'.join(header) + '\t' + "SV_Key" + '\t' + "SV_Type" + '\t' + "Dist_To_Junc1" + '\t' + "Dist_To_Junc2"


        for line in hin:
            F = line.rstrip('\n').split('\t')

            sj_chr = F[header2ind["SJ_1"]]
            sj_start = int(F[header2ind["SJ_2"]]) - 1
            sj_end = int(F[header2ind["SJ_3"]]) + 1

            gene1 = F[header2ind["Gene_1"]].split(';')
            gene2 = F[header2ind["Gene_2"]].split(';')
            junction1 = F[header2ind["Is_Boundary_1"]].split(';')
            junction2 = F[header2ind["Is_Boundary_2"]].split(';')


            """
            # just consider exon skipping genes
            for i in range(0, len(gene1)):
                if junction1[i] != "*" and junction2[i] != "*": 
                    targetGene.append(gene1[i])
                    targetGene.append(gene2[i])
            targetGene = list(set(targetGene))
            """

            associated_sv_list = []
            ##########
            # rough check for the mutation between the spliced region
            tabixErrorFlag1 = 0
            try:
                sv_list = sv_tb.fetch(sj_chr, sj_start - sv_comp_margin, sj_end + sv_comp_margin)
            except ValueError as inst:
                # print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
                tabixErrorFlag1 = 1

            # if there are some mutaions
            if tabixErrorFlag1 == 0 and sv_list is not None:
        
                for sv_cur in sv_list:

                    sv_curs = sv_cur.split('\t')
                    SV_cur = SV(chr1 = sv_curs[3], pos1 = int(sv_curs[4]), dir1 = sv_curs[5], chr2 = sv_curs[6], pos2 = int(sv_curs[7]), dir2 = sv_curs[8], inseq = sv_curs[9], sv_type = sv_curs[10])
                    if SV_cur.sv_type != "deletion": continue # For now, we just consider deletion
            

                    # the SV should be deletion confied within spliced junction
                    if  SV_cur.chr1 == sj_chr and SV_cur.chr2 == sj_chr and \
                      sj_start - sv_comp_margin <= SV_cur.pos1 and SV_cur.pos2 <= sj_end + sv_comp_margin:

                        """ 
                        # the splicing junction should be shared by SV breakpoint or exon-intron junction
                        junc_flag1 = 0
                        for i in range(0, len(gene1)):
                            if junction1[i] != "*": junc_flag1 = 1
                        
                        junc_flag2 = 0
                        for i in range(0, len(gene2)): 
                            if junction2[i] != "*": junc_flag2 = 1 
                          
                        if junc_flag1 == 1 or abs(sj_start - int(mutation[2])) <= sv_comp_margin and \
                          junc_flag2 == 1 or abs(sj_end - int(mutation[5])) <= sv_comp_margin:
                            mutation_sv.append('\t'.join(mutation))
                        """

                        # if F[header2ind["Splicing_Class"]] in ["exon-skip", "splice-site-slip", "pseudo-exon-inclusion"]: 
                        if F[header2ind["Splicing_Class"]] in ["exon-skip", "alternative-3'-splice-site", "alternative-5'-splice-site",
                                                               "intronic-alternative-3'-splice-site", "intronic-alternative-5'-splice-site"]:

                            associated_sv_list.append(SV_cur)
                        elif abs(sj_start - SV_cur.pos1) <= sv_comp_margin and abs(sj_end - SV_cur.pos2) <= sv_comp_margin:
                            associated_sv_list.append(SV_cur)

                            

            for SV_cur in associated_sv_list:

                junc_to_dist1 = SV_cur.pos1 - sj_start
                junc_to_dist2 = sj_end - SV_cur.pos2

                sv_key = ','.join([SV_cur.chr1, str(SV_cur.pos1), SV_cur.dir1, SV_cur.chr2, str(SV_cur.pos2), SV_cur.dir2, SV_cur.inseq])

                print >> hout, '\t'.join(F) + '\t' + sv_key + '\t' + SV_cur.sv_type + '\t' + str(junc_to_dist1) + '\t' + str(junc_to_dist2)

    hout.close()

