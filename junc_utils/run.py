#! /usr/bin/env python

import os, subprocess
import utils

def filter_main(args):

    utils.proc_star_junction(args.junc_file, args.output_path, args.pooled_control_file,
                             args.read_num_thres, args.overhang_thres, not args.keep_annotated, False)


def annotate_main(args):
   
    import annotate 
    annotate.annot_junction(args.junc_file, args.output_path, args.junction_margin, args.exon_margin, args.genome_id, args.grc)


def associate_main(args):

    import associate

    mutation_file = args.mutation_file
    output_file = args.output_file
    # is_anno = True if args.mutation_format == "anno" else False
    is_anno = False if mutation_file.endswith(".vcf") or mutation_file.endswith(".vcf.gz") else True 
    # control_file = args.ctrl
    is_sv = True if args.sv else False
    is_debug = True if args.debug else False 
    reference_genome = args.reference

    if not is_sv and is_anno and reference_genome is None:
        print >> sys.stderr, "When the mutation file format is annovar format, reference genome should be specified"
        sys.exit(1)

    output_dir = os.path.dirname(output_file)
    if output_dir != "" and not os.path.exists(output_dir):
       os.makedirs(output_dir)

    ##########
    # processing mutation file

    if not is_sv: 

        # convert mutation data to vcf (if --anno is on)
        if is_anno == True:
            utils.convert_anno2vcf(mutation_file, output_file + ".mutran_tmp.unsorted.vcf", reference_genome)
        else:
            utils.remove_vcf_header(mutation_file, output_file + ".mutran_tmp.unsorted.vcf")

        hout = open(output_file + ".mutran_tmp.vcf", 'w')
        s_ret = subprocess.call(["sort", "-k1,1", "-k2,2n", output_file + ".mutran_tmp.unsorted.vcf"], stdout = hout)
        hout.close()

        if s_ret != 0:
            print >> sys.stderr, "Error in sorting vcf file"
            sys.exit(1)

        s_ret = subprocess.call(["bgzip", "-f", output_file + ".mutran_tmp.vcf"])
        if s_ret != 0:
            print >> sys.stderr, "Error in bgzip compression"
            sys.exit(1)

        s_ret = subprocess.call(["tabix", "-p", "vcf", output_file + ".mutran_tmp.vcf.gz"])
        if s_ret != 0:
            print >> sys.stderr, "Error in tabix indexing"
            sys.exit(1)

    else:
        
        utils.convert_genosv2bed(mutation_file, output_file + ".mutran_tmp.unsorted.bedpe")

        hout = open(output_file + ".mutran_tmp.bedpe", 'w') 
        s_ret = subprocess.call(["sort", "-k1,1", "-k2,2n", "-k3,3n", "-k4,4", "-k5,5n", "-k6,6n", output_file + ".mutran_tmp.unsorted.bedpe"], stdout = hout)
        hout.close()
        
        if s_ret != 0:
            print >> sys.stderr, "Error in sorting bedpe file"
            sys.exit(1)
    
        s_ret = subprocess.call(["bgzip", "-f", output_file + ".mutran_tmp.bedpe"])
        if s_ret != 0:
            print >> sys.stderr, "Error in bgzip compression"
            sys.exit(1)
    
        
        s_ret = subprocess.call(["tabix", "-p", "bed", output_file + ".mutran_tmp.bedpe.gz"])
        if s_ret != 0:
            print >> sys.stderr, "Error in tabix indexing"
            sys.exit(1)

    ##########

    ##########
    # associate mutation and junction
    if not is_sv:
        """
        associate.get_snv_junction(args.annotated_junction_file,
                                   output_file,
                                   output_file + ".mutran_tmp.vcf.gz",
                                   args.annotation_dir, True, args.skip_creation_indel)
                                   # args.annotation_dir, args.is_edit_dist, args.skip_creation_indel)
        """

        if not args.only_dist:
            associate.get_snv_junction2(args.annotated_junction_file,
                                        output_file,
                                        output_file + ".mutran_tmp.vcf.gz",
                                        args.donor_size, args.acceptor_size,
                                        args.genome_id, args.grc,
                                        args.branchpoint, args.branchpoint_size)
        else:
            associate.get_snv_junction_only_dist(args.annotated_junction_file, 
                                                 output_file, 
                                                 output_file + ".mutran_tmp.vcf.gz",
                                                 args.annotation_dir, args.only_dist_search_margin)

    else:
        associate.get_sv_junction(args.annotated_junction_file,
                                  output_file,
                                  output_file + ".mutran_tmp.bedpe.gz",               
                                  args.genome_id, args.grc)
    ##########

    if is_debug != True:
        if not is_sv:
            subprocess.call(["rm", "-rf", output_file + ".mutran_tmp.unsorted.vcf"])
            subprocess.call(["rm", "-rf", output_file + ".mutran_tmp.vcf.gz"])
            subprocess.call(["rm", "-rf", output_file + ".mutran_tmp.vcf.gz.tbi"])
        else:
            subprocess.call(["rm", "-rf", output_file + ".mutran_tmp.unsorted.bedpe"])
            subprocess.call(["rm", "-rf", output_file + ".mutran_tmp.bedpe.gz"])
            subprocess.call(["rm", "-rf", output_file + ".mutran_tmp.bedpe.gz.tbi"])

        subprocess.call(["rm", "-rf", output_file + ".mutran_tmp.junction.txt"])
        subprocess.call(["rm", "-rf", output_file + ".mutran_tmp.junction.annot.txt"])



def merge_control_main(args):

    input_file_list = args.junc_list
    output_file = args.output_path
    read_num_thres = args.read_num_thres
    overhang_thres = args.overhang_thres
    remove_annotated = False if args.keep_annotated else True
    sample_num_thres = args.sample_num_thres

    # make directory for output if necessary
    if os.path.dirname(output_file) != "" and not os.path.exists(os.path.dirname(output_file)):
        os.makedirs(os.path.dirname(output_file))

    hin = open(input_file_list, 'r')
    hout = open(output_file + ".unsorted", 'w')

    with open(input_file_list, 'r') as hin:
        for line in hin:

            junction_file = line.rstrip('\n')
            with open(junction_file, 'r') as hin2:
                for line2 in hin2:

                    F = line2.rstrip('\n').split('\t')
                    if remove_annotated == True and F[5] != "0": continue
                    if int(F[6]) < read_num_thres: continue 
                    if int(F[8]) < overhang_thres: continue 
         
                    # convert to map-splice2 coordinate
                    # F[1] = str(int(F[1]) - 1)
                    # F[2] = str(int(F[2]) + 1)
         
                    print >> hout, '\t'.join(F)
                

    hout = open(output_file + ".sorted", 'w')
    s_ret = subprocess.call(["sort", "-k1,1", "-k2,2n", "-k3,3n", output_file + ".unsorted"], stdout = hout)
    hout.close()

    if s_ret != 0:
        print >> sys.stderr, "Error in sorting merged junction file"
        sys.exit(1)


    hout = open(output_file + ".merged", 'w')
    with open(output_file + ".sorted", 'r') as hin:
        temp_key = ""
        temp_read_num = []
        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = F[0] + '\t' + F[1] + '\t' + F[2]
            read_num = int(F[6])
            if key != temp_key:
                if temp_key != "":
                    if len(temp_read_num) >= sample_num_thres:
                        print >> hout, temp_key + '\t' + ','.join([str(x) for x in sorted(temp_read_num)])
                temp_key = key
                temp_read_num = []
            else:
                temp_read_num.append(read_num)

        if temp_key != "":
            if len(temp_read_num) >= sample_num_thres:
                print >> hout, temp_key + '\t' + ','.join([str(x) for x in sorted(temp_read_num)])



    hout = open(output_file, 'w')
    s_ret = subprocess.call(["bgzip", "-f", "-c", output_file + ".merged"], stdout = hout)
    hout.close()

    if s_ret != 0:
        print >> sys.stderr, "Error in compression merged junction file"
        sys.exit(1)


    s_ret = subprocess.call(["tabix", "-p", "bed", output_file])
    if s_ret != 0:
        print >> sys.stderr, "Error in indexing merged junction file"
        sys.exit(1)

    subprocess.call(["rm", "-f", output_file + ".unsorted"])
    subprocess.call(["rm", "-f", output_file + ".sorted"])
    subprocess.call(["rm", "-f", output_file + ".merged"])



def exonization_pair_main(args):

    # investigate the input splicing junction
    # 1. whether the both ends are located within the same gene
    # 2. whether one of the ends is at the authentic splicing junction
    # 3. whether one of the ends is located on intronic region
    # 4. the opposite authentic splicing junction to the one of the ends of input splicing junction

    # exonization_info = [chromosome, authentic exon-intron boundary, opposite exon-intron boundary, new splice site, strand, donor or acceptor]

    import exonization_pair

    exonization_info = exonization_pair.check_splicing_junction_for_exonization(args.half_exonizaiton_junction, args.output_file, 
                         args.genome_id, args.grc, args.min_new_intron_size, args.boundary_margin)
    
    # print exonization_info

    if len(exonization_info) == 0: 
        open(args.output_file, 'a').close()
    else:
        exonization_pair.check_opposite_junction(args.junc_file, exonization_info, args.output_file, args.control_file,
                                                 args.read_num_thres, args.max_new_exon_size, args.boundary_margin)

 

