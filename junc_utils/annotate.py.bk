#! /usr/bin/env python

import pysam 
import sys, subprocess
import annot_utils.gene, annot_utils.exon, annot_utils.coding

def spliced_coding_size(gene1, gene2, sj_chr, sj_start, sj_end, ref_coding_tb, exon_margin):

    tabixErrorFlag = 0
    try:
        records = ref_coding_tb.fetch(sj_chr, sj_start - exon_margin - 1, sj_end + exon_margin)
    except Exception as inst:
        tabixErrorFlag = 1

    coding_size = 0
    if tabixErrorFlag == 0:
        for record_line in records:
            record = record_line.split('\t')
            if record[3] in [gene1, gene2] and record[4] == "coding":

                # sj_start overlaps with the current coding region
                if int(record[1]) <= sj_start + exon_margin and sj_start - exon_margin <= int(record[2]):
                    # sj_end also overlaps with the current coding region
                    if int(record[1]) <= sj_end + exon_margin and sj_end - exon_margin <= int(record[2]):
                        coding_size = sj_end - sj_start - 1
                    else:
                        coding_size = coding_size + int(record[2]) - sj_start
                # sj_end overlaps with the current coding region
                elif int(record[1]) <= sj_end + exon_margin and sj_end - exon_margin <= int(record[2]):
                    coding_size = coding_size + sj_end - 1 - int(record[1])
                else:
                    coding_size = coding_size + int(record[2]) - int(record[1])

    return coding_size

 
# def annot_junction(input_file, output_file, annotation_dir, junction_margin, exon_margin):
def annot_junction(input_file, output_file, junction_margin, exon_margin, genome_id = "hg19", is_grc = True):
    """
        The purpose of this script is to classify splicing changes
        mainly by comparing the two breakpoints with the exon-intorn junction of genes
        within the database.
        Also, we generate the sequence arrond the breakpoints, which will be helpful
        for checking the authenticity of the splicing and evaluating the relationships
        with the somatic mutations.

        here is the classification categories:
        1. known (The splicing pattern is included in the database)
        (the start and end breakpoints are next exon-intron junctions of the same gene) 
        2. exon skipping 
        (the start and end breakpoints are exon-intron junctions of the same gene,
         but not the next ones)
        3. splice-site slip
        (one of the two breakpoints is an exon-intron junction and the other is within the 30bp exon of the same gene)
        4. pseudo-exon inclusion
        (one of the two break points is an exon-intron junction and the other is located in the same gene, but more than 30bp from exons of the gene)
        5. other
        (neighter of the two breakpoins are exon-intron junction, but located in the same gene)
        6. chimeric (spliced)
        7. chimeric (un-spliced)


        The algorithm for the annotation is as follows
        1. for both breakpoints, list up the exon-intron junctions matching to the breakpoints
        2. for both breakpoints, list up the exons within 30bp from the breakpoints
        3. for both breakpoints, list up the genes matching to the breakpoints
        4. summarize the above results and induce the annotation from them
        5. get the sequence arround the breakpoints.

    """

    annot_utils.gene.make_gene_info(output_file + ".tmp.refGene.bed.gz", "refseq", genome_id, is_grc, True)
    annot_utils.exon.make_exon_info(output_file + ".tmp.refExon.bed.gz", "refseq", genome_id, is_grc, True)
    annot_utils.coding.make_coding_info(output_file + ".tmp.refCoding.bed.gz", "refseq", genome_id, is_grc, True)
    
    hout = open(output_file, 'w') 
    gene_tb = pysam.TabixFile(output_file + ".tmp.refGene.bed.gz")
    exon_tb = pysam.TabixFile(output_file + ".tmp.refExon.bed.gz")
    coding_tb = pysam.TabixFile(output_file + ".tmp.refCoding.bed.gz") 

    # read header
    with open(input_file, 'r') as hin:
        line = hin.readline()
        F = line.rstrip('\n').split('\t')
        print >> hout, '\t'.join(["SJ_" + str(i) for i in range(1, len(F) + 1)]) + '\t' + \
                       "Splicing_Class" + '\t' + "Is_Inframe" + '\t' + "Gene_1" + '\t' + "Exon_Num_1" + '\t' + \
                       "Is_Boundary_1" + '\t' + "Offset_1" + '\t' + "Gene_2" + '\t' + "Exon_Num_2" + '\t' + \
                       "Is_Boundary_2" + '\t' + "Offset_2"

    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            chr_name = F[0]
            sj_start = int(F[1]) - 1
            sj_end = int(F[2]) + 1
            ##########
            # check gene annotation for the side 1  
            tabixErrorFlag = 0
            try:
                records = gene_tb.fetch(chr_name, sj_start - 1, sj_start)
            except Exception as inst:
                # print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
                # print >> sys.stderr, '\t'.join(F)
                tabixErrorFlag = 1

            gene1 = [];
            if tabixErrorFlag == 0:
                for record_line in records:
                    record = record_line.split('\t')
                    gene1.append(record[3])

            gene1 = list(set(gene1))
            ##########

            ##########
            # check gene annotation for the side 2  
            tabixErrorFlag = 0
            try:
                records = gene_tb.fetch(chr_name, sj_end - 1, sj_end)
            except Exception as inst:
                # print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
                # print >> sys.stderr, '\t'.join(F)
                tabixErrorFlag = 1
                
            gene2 = [];
            if tabixErrorFlag == 0:
                for record_line in records:
                    record = record_line.split('\t')
                    gene2.append(record[3])
                    
            gene2 = list(set(gene2))
            ##########

            ##########
            # check exon and junction annotation for the side 1  
            tabixErrorFlag = 0
            try:
                records = exon_tb.fetch(chr_name, sj_start - exon_margin - 1, sj_start + exon_margin)
            except Exception as inst:
                # print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
                # print >> sys.stderr, '\t'.join(F)
                tabixErrorFlag = 1
            
            exon1 = {}
            junction1 = {}
            offset1 = {}
            if tabixErrorFlag == 0:
                for record_line in records:
                    record = record_line.split('\t')
                    exon1[record[3]] = int(record[4])
                    if abs(sj_start - int(record[2])) < junction_margin:
                        offset1[record[3]] = sj_start - int(record[2])
                        if record[5] == "+": junction1[record[3]] = "e"
                        if record[5] == "-": junction1[record[3]] = "s"
            ##########

            ##########
            # check exon and junction annotation for the side 2
            tabixErrorFlag = 0
            try:
                records = exon_tb.fetch(chr_name, sj_end - exon_margin - 1, sj_end + exon_margin)
            except Exception as inst:
                # print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
                # print >> sys.stderr, '\t'.join(F)
                tabixErrorFlag = 1

            exon2 = {}
            junction2 = {}
            offset2 = {}
            if tabixErrorFlag == 0:
                for record_line in records:
                    record = record_line.split('\t')
                    exon2[record[3]] = int(record[4])
                    if abs(sj_end - 1 - int(record[1])) < junction_margin:
                        offset2[record[3]] = sj_end - 1 - int(record[1])
                        if record[5] == "+": junction2[record[3]] = "s"
                        if record[5] == "-": junction2[record[3]] = "e"
            ##########


            spliceClass = ""
            in_frame = "---"
            checkGenes = list(set(gene1 + gene2))
            ##########
            # check for know junction
            passGene = []
            for gene in checkGenes:
                if gene in gene1 and gene in gene2 and gene in junction1 and gene in junction2:
                    if junction1[gene] == "e" and junction2[gene] == "s" and exon2[gene] - exon1[gene] == 1: passGene.append(gene)
                    if junction2[gene] == "e" and junction1[gene] == "s" and exon1[gene] - exon2[gene] == 1: passGene.append(gene)

            if len(passGene) > 0: spliceClass = "known"

            ##########
            # check for exon skip
            if spliceClass == "":
                passGene = []
                inframe_gene = []
                for gene in checkGenes:
                    if gene in gene1 and gene in gene2 and gene in junction1 and gene in junction2:
                        if (junction1[gene] == "e" and junction2[gene] == "s" and exon2[gene] - exon1[gene] > 1) or \
                           (junction2[gene] == "e" and junction1[gene] == "s" and exon1[gene] - exon2[gene] > 1): 
                            passGene.append(gene)

                            sc_size = spliced_coding_size(gene, None, chr_name, sj_start, sj_end, coding_tb, exon_margin)
                            if sc_size != 0 and sc_size % 3 == 0:
                                inframe_gene.append(gene)

                if len(passGene) > 0: spliceClass = "exon-skip"
                if len(inframe_gene) > 0: in_frame = "in-frame"

            ##########
            # check for alternative-3'-splice-site
            if spliceClass == "":
                passGene = []
                inframe_gene = []
                for gene in checkGenes:
                    if gene in gene1 and gene in gene2:
                        if (gene in junction1 and junction1[gene] == "e" and gene not in junction2 and gene in exon2) or \
                           (gene in junction2 and junction2[gene] == "e" and gene not in junction1 and gene in exon1):
                            passGene.append(gene)

                            sc_size = spliced_coding_size(gene, None, chr_name, sj_start, sj_end, coding_tb, exon_margin)
                            if sc_size != 0 and sc_size % 3 == 0:
                                inframe_gene.append(gene)

                if len(passGene) > 0: spliceClass = "alternative-3'-splice-site"
                if len(inframe_gene) > 0: in_frame = "in-frame"


            # check for alternative-5'-splice-site
            if spliceClass == "":
                passGene = []
                inframe_gene = []
                for gene in checkGenes:
                    if gene in gene1 and gene in gene2:
                        if (gene in junction1 and junction1[gene] == "s" and gene not in junction2 and gene in exon2) or \
                           (gene in junction2 and junction2[gene] == "s" and gene not in junction1 and gene in exon1):
                            passGene.append(gene)

                            sc_size = spliced_coding_size(gene, None, chr_name, sj_start, sj_end, coding_tb, exon_margin)
                            if sc_size != 0 and sc_size % 3 == 0:
                                inframe_gene.append(gene)
            
                if len(passGene) > 0: spliceClass = "alternative-5'-splice-site"
                if len(inframe_gene) > 0: in_frame = "in-frame"

            ##########
            # check for intronic-alternative-3'-splice-site
            if spliceClass == "":
                passGene = []
                inframe_gene = []
                for gene in checkGenes:
                    if gene in gene1 and gene in gene2:
                        if (gene in junction1 and junction1[gene] == "e" and gene not in junction2 and gene not in exon2) or \
                           (gene in junction2 and junction2[gene] == "e" and gene not in junction1 and gene not in exon1):
                            passGene.append(gene)
            
                if len(passGene) > 0: spliceClass = "intronic-alternative-3'-splice-site"

            # check for intronic-alternative-5'-splice-site
            if spliceClass == "":
                passGene = []
                inframe_gene = []
                for gene in checkGenes:
                    if gene in gene1 and gene in gene2:
                        if (gene in junction1 and junction1[gene] == "s" and gene not in junction2 and gene not in exon2) or \
                           (gene in junction2 and junction2[gene] == "s" and gene not in junction1 and gene not in exon1):
                            passGene.append(gene)
                
                if len(passGene) > 0: spliceClass = "intronic-alternative-5'-splice-site"

            ##########
            # within-exon
            if spliceClass == "":
                passGene = []
                inframe_gene = []
                for gene in checkGenes:
                    if gene in gene1 and gene in gene2 and gene in exon1 and gene in exon2:
                        if exon1[gene] == exon2[gene]:
                            passGene.append(gene)

                            sc_size = spliced_coding_size(gene, None, chr_name, sj_start, sj_end, coding_tb, exon_margin)
                            if sc_size != 0 and sc_size % 3 == 0:
                                inframe_gene.append(gene)

                if len(passGene) > 0: spliceClass = "within-exon"
                if len(inframe_gene) > 0: in_frame = "in-frame"

            ##########
            # check for exon-exon-junction
            if spliceClass == "":
                passGene = []
                inframe_gene = []
                for gene in checkGenes:
                    if gene in gene1 and gene in gene2:
                        if gene in exon1 and gene in exon2:
                            passGene.append(gene)

                            sc_size = spliced_coding_size(gene, None, chr_name, sj_start, sj_end, coding_tb, exon_margin)
                            if sc_size != 0 and sc_size % 3 == 0:
                                inframe_gene.append(gene)
     
                if len(passGene) > 0: spliceClass = "exon-exon-junction"
                if len(inframe_gene) > 0: in_frame = "in-frame"

            ##########
            # check for within-gene 
            if spliceClass == "":
                passGene = []
                for gene in checkGenes:
                    if gene in gene1 and gene in gene2 and gene: passGene.append(gene)   
            
                if len(passGene) > 0: spliceClass = "within-gene"


            ##########
            # check for spliced-chimera 
            if spliceClass == "":
                passGene = []
                for g1 in gene1:
                    for g2 in gene2:
                        if (g1 in junction1 and junction1[g1] == "s" and g2 in junction2 and junction2[g2] == "e") or \
                           (g1 in junction1 and junction1[g1] == "e" and g2 in junction2 and junction2[g2] == "s"): 
                            passGene.append(g1 + ',' + g2)

                            """
                            sc_size = spliced_coding_size(gene, None, chr_name, sj_start, sj_end, coding_tb, exon_margin)
                            if sc_size != 0 and sc_size % 3 == 0:
                                inframe_gene.append(gene)
                            """

                if len(passGene) > 0: spliceClass = "spliced-chimera"
                if len(inframe_gene) > 0: in_frame = "in-frame"

            ##########
            # check for unspliced-chimera 
            if spliceClass == "":
                passGene = []
                for g1 in gene1:
                    for g2 in gene2:
                        passGene.append(g1 + ',' + g2)

                if len(passGene) > 0: spliceClass = "unspliced-chimera"


            if spliceClass == "": spliceClass = "other"
            

            # summarize the exon and junction information for display
            geneInfo1 = []
            exonInfo1 = []
            junctionInfo1 = []
            offsetInfo1 = []
            if len(gene1) > 0:
                for g1 in sorted(gene1):
                    if g1 not in passGene: continue 
                    geneInfo1.append(g1)
                    if g1 in exon1: 
                        exonInfo1.append(str(exon1[g1]))
                    else:
                        exonInfo1.append("*")

                    if g1 in junction1:
                        junctionInfo1.append(junction1[g1])
                    else:
                        junctionInfo1.append("*")

                    if g1 in offset1:
                        offsetInfo1.append(str(offset1[g1]))
                    else:
                        offsetInfo1.append("*")

            if len(geneInfo1) == 0: 
                geneInfo1.append("---")
                exonInfo1.append("---")
                junctionInfo1.append("---")
                offsetInfo1.append("---")

            geneInfo2 = []
            exonInfo2 = []
            junctionInfo2 = []
            offsetInfo2 = []
            if len(gene2) > 0:
                for g2 in sorted(gene2):
                    if g2 not in passGene: continue
                    geneInfo2.append(g2)
                    if g2 in exon2:
                        exonInfo2.append(str(exon2[g2]))
                    else:
                        exonInfo2.append("*")
                    
                    if g2 in junction2:
                        junctionInfo2.append(junction2[g2])
                    else:
                        junctionInfo2.append("*")

                    if g2 in offset2:
                        offsetInfo2.append(str(offset2[g2]))
                    else:
                        offsetInfo2.append("*")

            if len(geneInfo2) == 0:
                geneInfo2.append("---")
                exonInfo2.append("---")
                junctionInfo2.append("---")
                offsetInfo2.append("---")

         

            print >> hout, '\t'.join(F) + '\t' + spliceClass + '\t' + in_frame + '\t' + '\t'.join([';'.join(geneInfo1), ';'.join(exonInfo1), ';'.join(junctionInfo1), ';'.join(offsetInfo1), ';'.join(geneInfo2), ';'.join(exonInfo2), ';'.join(junctionInfo2), ';'.join(offsetInfo2)])
     

    hout.close()


    subprocess.call(["rm", "-rf", output_file + ".tmp.refGene.bed.gz"])
    subprocess.call(["rm", "-rf", output_file + ".tmp.refGene.bed.gz.tbi"]) 
    subprocess.call(["rm", "-rf", output_file + ".tmp.refExon.bed.gz"]) 
    subprocess.call(["rm", "-rf", output_file + ".tmp.refExon.bed.gz.tbi"])
    subprocess.call(["rm", "-rf", output_file + ".tmp.refCoding.bed.gz"]) 
    subprocess.call(["rm", "-rf", output_file + ".tmp.refCoding.bed.gz.tbi"])


