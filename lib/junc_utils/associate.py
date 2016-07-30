#! /usr/bin/env python

import sys
import tabix
import numpy


def get_snv_junction(input_file, output_file, mutation_file, annotation_dir):

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

    searchMargin1 = 30
    searchMargin2 = 10

    splicingDonnorMotif = ["AG", "GTRAGT"]
    splicingAcceptorMotif = ["YYYYNCAG", "G"]

    nuc2vec = {'A': [1, 0, 0, 0], 'C': [0, 1, 0, 0], 'G': [0, 0, 1, 0], 'T': [0, 0, 0, 1], \
               'W': [1, 0, 0, 1], 'S': [0, 1, 1, 0], 'M': [1, 1, 0, 0], 'K': [0, 0, 1, 1], \
               'R': [1, 0, 1, 0], 'Y': [0, 1, 0, 1], 'B': [0, 1, 1, 1], 'D': [1, 0, 1, 1], \
               'H': [1, 1, 0, 1], 'V': [1, 1, 1, 0], 'N': [1, 1, 1, 1]}


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
    exon_tb = tabix.open(ref_exon_bed)

    header2ind = {}
    with open(input_file, 'r') as hin:
        
        # read header
        header = hin.readline().rstrip('\n').split('\t')
        for i in range(len(header)):
            header2ind[header[i]] = i

        # print header
        print >> hout, '\t'.join(header) + '\t' + "Mutation_Info" + '\t' + "Splicing_Motif_Pos" + '\t' + "Splicing_Mutation_Type" + '\t' + "Is_Cannonical"

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


            # just consider genes sharing the exon-intron junction with the breakpoints of splicings
            for i in range(0, len(gene1)):
                if junction1[i] != "*": targetGene.append(gene1[i])
            for i in range(0, len(gene2)):
                if junction2[i] != "*": targetGene.append(gene2[i])
            targetGene = list(set(targetGene))


            if F[header2ind["Splicing_Class"]] in ["alternative-3'-splice-site", "alternative-5'-splice-site",
                                                   "intronic-alternative-3'-splice-site", "intronic-alternative-5'-splice-site"]:

                # for non exon-intron junction breakpoints
                if "*" in junction1 and "s" in junction2: # splicing donor motif, plus direction
                    firstSearchRegion[1] = firstSearchRegion[1] - searchMargin1
                    splicingMotifRegions.append((F[header2ind["SJ_1"]], sj_start - len(splicingDonnorMotif[0]) + 1, sj_start + len(splicingDonnorMotif[1]), "donor", "+", 0))
                if "*" in junction1 and "e" in junction2: # splicing acceptor motif, minus direction
                    firstSearchRegion[1] = firstSearchRegion[1] - searchMargin1
                    splicingMotifRegions.append((F[header2ind["SJ_1"]], sj_start - len(splicingAcceptorMotif[1]) + 1, sj_start + len(splicingAcceptorMotif[0]), "acceptor", "-", 0))
                if "s" in junction1 and "*" in junction2: # splicing donor motif, minus direction
                    firstSearchRegion[2] = firstSearchRegion[2] + searchMargin1
                    splicingMotifRegions.append((F[header2ind["SJ_1"]], sj_end - len(splicingDonnorMotif[1]), sj_end + len(splicingDonnorMotif[0]) - 1, "donor", "-", 0))
                if "e" in junction1 and "*" in junction2: # # splicing acceptor motif, plus direction
                    firstSearchRegion[2] = firstSearchRegion[2] + searchMargin1
                    splicingMotifRegions.append((F[header2ind["SJ_1"]], sj_end - len(splicingAcceptorMotif[0]), sj_end + len(splicingAcceptorMotif[1]) - 1, "acceptor", "+", 0))


            ##########
            # rough check for the mutation between the spliced region
            tabixErrorFlag1 = 0
            try:
                mutations = mutation_tb.query(firstSearchRegion[0], firstSearchRegion[1], firstSearchRegion[2])
            except Exception as inst:
                # print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
                # print >> sys.stderr, '\t'.join(F)
                tabixErrorFlag1 = 1

            # if there are some mutaions
            if tabixErrorFlag1 == 0 and mutations is not None:

                chr_name = grch2ucsc[F[header2ind["SJ_1"]]] if F[header2ind["SJ_1"]] in grch2ucsc else F[header2ind["SJ_1"]] 
                # check the exons within the spliced regions
                tabixErrorFlag2 = 0
                try:
                    exons = exon_tb.query(chr_name, firstSearchRegion[1], firstSearchRegion[2])
                except Exception as inst:
                    # print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
                    # print >> sys.stderr, '\t'.join(F)
                    tabixErrorFlag2 = 1

                # first, add the exon-intron junction for detailed check region list
                if tabixErrorFlag2 == 0:
                    for exon in exons:
                        if exon[3] not in targetGene: continue
                        if exon[5] == "+":
                            # splicing acceptor motif, plus direction
                            splicingMotifRegions.append((exon[0], int(exon[1]) - len(splicingAcceptorMotif[0]) + 1, int(exon[1]) + len(splicingAcceptorMotif[1]), "acceptor", "+", 1))
                            # splicing donor motif, plus direction
                            splicingMotifRegions.append((exon[0], int(exon[2]) - len(splicingDonnorMotif[0]) + 1, int(exon[2]) + len(splicingDonnorMotif[1]), "donor", "+", 1))
                        if exon[5] == "-":
                            # splicing donor motif, minus direction 
                            splicingMotifRegions.append((exon[0], int(exon[1]) - len(splicingDonnorMotif[1]) + 1, int(exon[1]) + len(splicingDonnorMotif[0]), "donor", "-", 1))
                            # splicing acceptor motif, minus direction
                            splicingMotifRegions.append((exon[0], int(exon[2]) - len(splicingAcceptorMotif[1]) + 1, int(exon[2]) + len(splicingAcceptorMotif[0]), "acceptor", "-", 1))


                splicingMotifRegions = list(set(splicingMotifRegions))

                # compare each mutation with exon-intron junction regions and non-exon-intorn junction breakpoints.
                for mutation in mutations:
                    RegMut = []
                    for reg in splicingMotifRegions:

                        # insertion or deletion (just consider the disruption of splicing motifs now)
                        if (len(mutation[3]) > 1 or len(mutation[4]) > 1) and reg[5] == 1:
                            indel_start = int(mutation[1]) + 1
                            indel_end = int(mutation[1]) + len(mutation[3]) - 1 if len(mutation[3]) > 1 else indel_start
                            if indel_start <= reg[2] and reg[1] <= indel_end:
                            
                                is_canonical = "non-canonical" 
                                if reg[3] == "acceptor" and reg[4] == "+" and indel_start <= reg[2] - 1 and reg[1] + 6 <= indel_end: is_canonical = "canonical"
                                if reg[3] == "acceptor" and reg[4] == "-" and indel_start <= reg[2] - 6 and reg[1] + 1 <= indel_end: is_canonical = "canonical" 
                                if reg[3] == "donor" and reg[4] == "+" and indel_start <= reg[2] - 4 and reg[1] + 2 <= indel_end: is_canonical = "canonical" 
                                if reg[3] == "donor" and reg[4] == "-" and indel_start <= reg[2] - 2 and reg[1] + 4 <= indel_end: is_canonical = "canonical" 


                                RegMut.append([reg, "splicing " + reg[3] + " disruption", is_canonical])
                        
                        # base substitution
                        if len(mutation[3]) == 1 and len(mutation[4]) == 1 and reg[1] <= int(mutation[1]) <= reg[2]:
                            motifSeq = ""
                            if reg[3] == "acceptor": motifSeq = splicingAcceptorMotif[0] + splicingAcceptorMotif[1]                  
                            if reg[3] == "donor": motifSeq = splicingDonnorMotif[0] + splicingDonnorMotif[1]

                            if reg[4] == "-":
                                complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', \
                                              'W': 'S', 'S': 'W', 'M': 'K', 'K': 'M', \
                                              'R': 'Y', 'Y': 'R', 'B': 'V', 'D': 'H', \
                                              'H': 'D', 'V': 'B', 'N': 'N'}
                                motifSeq = "".join(complement.get(base) for base in reversed(motifSeq)) 


                            vecAtMut = nuc2vec[motifSeq[int(mutation[1]) - int(reg[1])]]
                            editDistDiff = numpy.dot(vecAtMut, nuc2vec[mutation[4]]) - numpy.dot(vecAtMut, nuc2vec[mutation[3]]) 

                            is_canonical = "non-canonical"
                            if reg[3] == "acceptor" and reg[4] == "+" and reg[1] + 6 <= int(mutation[1]) <= reg[2] - 1: is_canonical = "canonical"
                            if reg[3] == "acceptor" and reg[4] == "-" and reg[1] + 1 <= int(mutation[1]) <= reg[2] - 6: is_canonical = "canonical"
                            if reg[3] == "donor" and reg[4] == "+" and reg[1] + 2 <= int(mutation[1]) <= reg[2] - 4: is_canonical = "canonical"
                            if reg[3] == "donor" and reg[4] == "-" and reg[1] + 4 <= int(mutation[1]) <= reg[2] - 2: is_canonical = "canonical"

                            if editDistDiff > 0 and reg[5] == 0: RegMut.append([reg, "splicing " + reg[3] + " creation", is_canonical])
                            if editDistDiff < 0 and reg[5] == 1: RegMut.append([reg, "splicing " + reg[3] + " disruption", is_canonical])

                    for item in RegMut:
                        print >> hout, '\t'.join(F) + '\t' + ','.join(mutation) + '\t' + F[header2ind["SJ_1"]] + ':' + str(item[0][1]) + '-' + str(item[0][2]) + ',' + item[0][4] + '\t' + item[1] + '\t' + item[2]


    hout.close()



def get_sv_junction(input_file, output_file, mutation_file, annotation_dir):

    """
        a script for detecting candidate structural variations (currently just deletions causing splicing changes
        
    """

    sv_comp_margin = 20
    exon_comp_margin = 10


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
    exon_tb = tabix.open(ref_exon_bed)

    header2ind = {}
    with open(input_file, 'r') as hin:

        # read header
        header = hin.readline().rstrip('\n').split('\t')
        for i in range(len(header)):
            header2ind[header[i]] = i

        # print header
        print >> hout, '\t'.join(header) + '\t' + "SV_Info"


        for line in hin:
            F = line.rstrip('\n').split('\t')


            # if F[3] not in ["within-gene", "exon-skip", "exon-exon-junction", "spliced-chimera", "unspliced-chimera"]: continue

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

            mutation_sv = []
            ##########
            # rough check for the mutation between the spliced region
            tabixErrorFlag1 = 0
            try:
                mutations = mutation_tb.query(F[header2ind["SJ_1"]], sj_start - sv_comp_margin, sj_end + sv_comp_margin)
            except Exception as inst:
                # print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
                # print >> sys.stderr, '\t'.join(F)
                tabixErrorFlag1 = 1

            # if there are some mutaions
            if tabixErrorFlag1 == 0 and mutations is not None:
        
                for mutation in mutations:
                    # the SV should be deletion and the SV should be confied within spliced junction
                    if mutation[8] == '+' and mutation[9] == '-' and mutation[0] == F[0] and mutation[3] == F[0] and \
                      sj_start - sv_comp_margin <= int(mutation[2]) and int(mutation[5]) <= sj_end + sv_comp_margin:

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

                        if F[header2ind["Splicing_Class"]] in ["exon-skip", "splice-site-slip", "pseudo-exon-inclusion"]: 
                            mutation_sv.append('\t'.join(mutation))
                        elif abs(sj_start - int(mutation[2])) <= sv_comp_margin and abs(sj_end - int(mutation[5])) <= sv_comp_margin:
                            mutation_sv.append('\t'.join(mutation))


            for mutation in mutation_sv:
                muts = mutation.split('\t')
                print >> hout, '\t'.join(F) + '\t' + muts[0] + ',' + muts[2] + ',' + muts[8] + ',' + \
                                                     muts[3] + ',' + muts[5] + ',' + muts[9] + ',' + muts[7]

    hout.close()

