from Bio import SeqIO
import re
import argparse

# ============= main part ==========>>>>>>>>>>>>>>>>>>>>>>
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Get all the related insertions and deletions information based on the CRISPResso2 allele table result file")
    parser.add_argument("-ifa", "--input_fasta", help="The input reference fasta file", required=True)
    parser.add_argument("-iat", "--input_allele_table", help="The input CRISPResso2 allele table file", required=True)
    parser.add_argument("-event_indel", "--output_event_indel", help="Output each reference index exist event position plus total indel information", required=True)
    parser.add_argument("-idx_indel", "--output_idx_indel", help="Output the only the idx indel info, not plus each idx", required=True)
    parser.add_argument("-ins", "--output_ins_info", help="Output each insertion event position, length and count", required=True)
    parser.add_argument("-del", "--output_del_info", help="Output each deletion event position, length and count", required=True)

    ARGS = parser.parse_args()

    ### ==== identify params and load files =====###
    in_fasta = ARGS.input_fasta
    in_allele = ARGS.input_allele_table
    out_event_indel = ARGS.output_event_indel
    out_idx_indel = ARGS.output_idx_indel
    out_ins = ARGS.output_ins_info
    out_del = ARGS.output_del_info


    ### =======================================>>>>>>
    ### ============ Main analysis part ========= ###
    ### =======================================>>>>>>
    ref_dict = {}
    for seq_record in SeqIO.parse(in_fasta, "fasta"):
        ref_dict[seq_record.id] = str(seq_record.seq)


    ref_idx_indel = {}
    ins_info = {}
    del_info = {}

    with open(in_allele) as f:
        rows = f.readlines()[1:]

        for i in ref_dict.keys():
            ref_seq = ref_dict[i]

        def get_ins_del_info_separate(align_read_seq,align_ref_seq):
            ### ====== insertion information ====== ###
            l = 0
            ins_span = [(m.start()+1,m.end()) for m in re.finditer(r'((-){1,})', align_ref_seq)]

            for i in range(len(ins_span)):
                if i == 0:
                    idx = ins_span[0][0]
                    ins_seq = str(align_read_seq[ins_span[0][0]-1:ins_span[0][1]])
                    ins_key = str(idx) + "_" +str((ins_span[0][1]-ins_span[0][0]+1))+"bp"+"_"+ins_seq
                    if ins_key in ins_info.keys():
                        ins_info[ins_key] = int(ins_info[ins_key])+int(align_dict[align_ref+"_"+row_idx])
                    else:
                        ins_info[ins_key] = int(align_dict[align_ref+"_"+row_idx])
                    l += ins_span[0][1]-ins_span[0][0]+1

                if i > 0:
                    idx = ins_span[i][0]-l-1
                    ins_seq = str(align_read_seq[ins_span[i][0] - 1:ins_span[i][1]])
                    ins_key = str(idx) + "_" +str((ins_span[i][1]-ins_span[i][0]+1))+"bp"+"_"+ins_seq
                    if ins_key in ins_info.keys():
                        ins_info[ins_key] = int(ins_info[ins_key]) + int(align_dict[align_ref+"_"+row_idx])
                    else:
                        ins_info[ins_key] = int(align_dict[align_ref+"_"+row_idx])
                    l += ins_span[i][1]-ins_span[i][0]+1

            ### ====== deletion information ====== ###
            del_span = [(m.start()+1,m.end()) for m in re.finditer(r'((-){1,})', align_read_seq)]
            for i in range(len(del_span)):
                del_key = str(del_span[i][0])+"-"+str(del_span[i][1]) + "_" + str((del_span[i][1]-del_span[i][0]+1)) + "bp"
                if del_key in del_info.keys():
                    del_info[del_key] = int(del_info[del_key]) + int(align_dict[align_read+"_"+row_idx])
                else:
                    del_info[del_key] = int(align_dict[align_read+"_"+row_idx])


        def get_ins_seq_only(align_read_seq, align_ref_seq):
            ins_span = [(m.start() + 1, m.end()) for m in re.finditer(r'((-){1,})', align_ref_seq)]

            final_ins_seq = []
            for i in range(len(ins_span)):
                if i == 0:
                    ins_seq = str(align_read_seq[ins_span[0][0] - 1:ins_span[0][1]])
                    final_ins_seq.append(ins_seq)

                if i > 0:
                    ins_seq = str(align_read_seq[ins_span[i][0] - 1:ins_span[i][1]])
                    final_ins_seq.append(ins_seq)
            return final_ins_seq


        def each_idx_indel_info(align_read_seq,align_ref_seq):
            ### ========= deal insertion ========= ###
            ins_c = 0
            del_c = 0

            ins_idx_lst = list()
            l = 0

            ins_span = [(m.start()+1,m.end()) for m in re.finditer(r'((-){1,})', align_ref_seq)]
            for i in range(len(ins_span)):
                if i == 0:
                    idx = ins_span[0][0]
                    ins_seq = str(align_read_seq[ins_span[0][0] - 1:ins_span[0][1]])
                    ins_idx_lst.append(idx)
                    l += ins_span[0][1]-ins_span[0][0]+1

                if i > 0:
                    idx = ins_span[i][0]-l-1
                    ins_seq = str(align_read_seq[ins_span[i][0] - 1:ins_span[i][1]])
                    ins_idx_lst.append(idx)
                    l += ins_span[i][1]-ins_span[i][0]+1

            if ref_idx in ins_idx_lst:
                ins_c += int(align_dict[align_ref+"_"+row_idx])
    #             ins_seq.append(align_ref+"_"+align_dict[align_ref+"_"+row_idx])
            else:
                ins_c = ins_c

            ### ========= deal deletion ========= ###

            del_pos = [m.start()+1 for m in re.finditer(r'((-))', align_read_seq)]
            if ref_idx in del_pos:
                del_c += int(align_dict[align_read+"_"+row_idx])
            else:
                del_c = del_c

            return ins_c,del_c


        #### ======================================>>>>>
        ### Get insertion and deletion separate info
        #### ======================================>>>>>
        align_dict = {}
    #     test_ins_dict = {}
        for line in rows:
            data = line.strip().split("\t")
            align_read = data[0]
            align_ref = data[1]
            count = data[7]
            row_idx = str(data[-1])

            align_dict[align_read+"_"+row_idx]=count
            align_dict[align_ref+"_"+row_idx]=count

            ### align_read exit "-"
            if "-" in align_read and "-" not in align_ref:
                if align_ref[0] != "-":
                    get_ins_del_info_separate(align_read,align_ref)

                if align_ref[0] == "-" and align_ref[-1] == "-":
                    contain_rm_span = [(m.start()+1,m.end()) for m in re.finditer(r'((-){1,})', align_ref)]
                    new_align_read = align_read[contain_rm_span[0][1]:contain_rm_span[-1][0]-1]
                    new_align_ref = align_ref[contain_rm_span[0][1]:contain_rm_span[-1][0]-1]

                    get_ins_del_info_separate(new_align_read,new_align_ref)

            ### align_ref exit "-"
            elif "-" in align_ref and "-" not in align_read:
                if align_ref[0] != "-":
                    get_ins_del_info_separate(align_read, align_ref)

                if align_ref[0] == "-" and align_ref[-1] == "-":
                    contain_rm_span = [(m.start() + 1, m.end()) for m in re.finditer(r'((-){1,})', align_ref)]
                    new_align_read = align_read[contain_rm_span[0][1]:contain_rm_span[-1][0] - 1]
                    new_align_ref = align_ref[contain_rm_span[0][1]:contain_rm_span[-1][0] - 1]

                    get_ins_del_info_separate(new_align_read, new_align_ref)

            ### align_read and align_ref exit "-" simultaneously
            elif "-" in align_read and align_ref:
                if align_ref[0] != "-":
                    ### cal ins: rm "-" in align_read
                    rm_gap_align_read = align_read.replace("-","N")
                    get_ins_del_info_separate(rm_gap_align_read, align_ref)
                    get_ins_seq = get_ins_seq_only(rm_gap_align_read,align_ref)
                    # print(get_ins_seq)

                    rm_ins_align_read = None
                    rm_ins_align_read = align_read.replace(get_ins_seq[0], "")

                    for i in range(1, len(get_ins_seq)):
                        rm_ins_align_read = rm_ins_align_read.replace(get_ins_seq[i], "")

                    ### cal del: rm ins_seq in align_read and rm "-" in align_ref
                    rm_gap_align_ref = align_ref.replace("-","")
                    get_ins_del_info_separate(rm_ins_align_read,rm_gap_align_ref)

                if align_ref[0] == "-" and align_ref[-1] == "-":
                    contain_rm_span = [(m.start() + 1, m.end()) for m in re.finditer(r'((-){1,})', align_ref)]
                    new_align_read = align_read[contain_rm_span[0][1]:contain_rm_span[-1][0] - 1]
                    new_align_ref = align_ref[contain_rm_span[0][1]:contain_rm_span[-1][0] - 1]

                    ### cal ins: rm "-" in align_read
                    rm_gap_align_read = new_align_read.replace("-", "N")
                    get_ins_del_info_separate(rm_gap_align_read, new_align_ref)
                    get_ins_seq = get_ins_seq_only(rm_gap_align_read, new_align_ref)
                    # print(get_ins_seq)

                    if not get_ins_seq:
                        rm_ins_align_read = align_read
                    else:
                        rm_ins_align_read = None
                        rm_ins_align_read = align_read.replace(get_ins_seq[0], "")

                        for i in range(1, len(get_ins_seq)):
                            rm_ins_align_read = rm_ins_align_read.replace(get_ins_seq[i], "")

                    ### cal del: rm ins_seq in align_read and rm "-" in align_ref
                    rm_gap_align_ref = new_align_ref.replace("-", "")
                    get_ins_del_info_separate(rm_ins_align_read, rm_gap_align_ref)


        #### ======================================>>>>>
        ### Get each reference index position indel info
        #### ======================================>>>>>

        for ref_idx in range(1,len(ref_seq)+2):
            f_ins_c = 0
            f_del_c = 0

            ins_seq = list()
            for line in rows:
                data = line.strip().split("\t")
                align_read = data[0]
                align_ref = data[1]
                count = data[7] 
                row_idx = str(data[-1])

                ### align_read exit "-"
                if "-" in align_read and "-" not in align_ref:
                    if align_ref[0] != "-":
                        f_ins_c += each_idx_indel_info(align_read, align_ref)[0]
                        f_del_c += each_idx_indel_info(align_read, align_ref)[1]

                    ### ========= if align_ref start with "-", indicate that reads sequencing carry random N ========= ###
                    if align_ref[0] == "-" and align_ref[-1] == "-":
                        contain_rm_span = [(m.start() + 1, m.end()) for m in re.finditer(r'((-){1,})', align_ref)]
                        new_align_read = align_read[contain_rm_span[0][1]:contain_rm_span[-1][0] - 1]
                        new_align_ref = align_ref[contain_rm_span[0][1]:contain_rm_span[-1][0] - 1]

                        f_ins_c += each_idx_indel_info(new_align_read, new_align_ref)[0]
                        f_del_c += each_idx_indel_info(new_align_read, new_align_ref)[1]

                ### align_ref exit "-"
                elif "-" in align_ref and "-" not in align_read:
                    if align_ref[0] != "-":
                        f_ins_c += each_idx_indel_info(align_read, align_ref)[0]
                        f_del_c += each_idx_indel_info(align_read, align_ref)[1]

                    ### ========= if align_ref start with "-", indicate that reads sequencing carry random N ========= ###
                    if align_ref[0] == "-" and align_ref[-1] == "-":
                        contain_rm_span = [(m.start() + 1, m.end()) for m in re.finditer(r'((-){1,})', align_ref)]
                        new_align_read = align_read[contain_rm_span[0][1]:contain_rm_span[-1][0] - 1]
                        new_align_ref = align_ref[contain_rm_span[0][1]:contain_rm_span[-1][0] - 1]

                        f_ins_c += each_idx_indel_info(new_align_read, new_align_ref)[0]
                        f_del_c += each_idx_indel_info(new_align_read, new_align_ref)[1]

                ### align_read and align_ref exit "-" simultaneously
                elif "-" in align_read and align_ref:
                    if align_ref[0] != "-":
                        ### cal ins: rm "-" in align_read
                        rm_gap_align_read = align_read.replace("-", "N")
                        get_ins_seq = get_ins_seq_only(rm_gap_align_read, align_ref)
                        f_ins_c += each_idx_indel_info(rm_gap_align_read, align_ref)[0]

                        ### cal del: rm ins_seq in align_read and rm "-" in align_ref
                        rm_gap_align_ref = align_ref.replace("-", "")

                        rm_ins_align_read = None
                        rm_ins_align_read = align_read.replace(get_ins_seq[0], "")

                        for i in range(1, len(get_ins_seq)):
                            rm_ins_align_read = rm_ins_align_read.replace(get_ins_seq[i], "")

                        f_del_c += each_idx_indel_info(rm_ins_align_read, rm_gap_align_ref)[1]

                    if align_ref[0] == "-" and align_ref[-1] == "-":
                        contain_rm_span = [(m.start() + 1, m.end()) for m in re.finditer(r'((-){1,})', align_ref)]
                        new_align_read = align_read[contain_rm_span[0][1]:contain_rm_span[-1][0] - 1]
                        new_align_ref = align_ref[contain_rm_span[0][1]:contain_rm_span[-1][0] - 1]

                        ### cal ins: rm "-" in align_read
                        rm_gap_align_read = new_align_read.replace("-", "N")
                        get_ins_seq = get_ins_seq_only(rm_gap_align_read, new_align_ref)
                        f_ins_c += each_idx_indel_info(rm_gap_align_read, new_align_ref)[0]

                        ### cal del: rm ins_seq in align_read and rm "-" in align_ref
                        rm_gap_align_ref = new_align_ref.replace("-", "")

                        if not get_ins_seq:
                            rm_ins_align_read = align_read
                        else:
                            rm_ins_align_read = None
                            rm_ins_align_read = align_read.replace(get_ins_seq[0], "")

                            for i in range(1, len(get_ins_seq)):
                                rm_ins_align_read = rm_ins_align_read.replace(get_ins_seq[i], "")

                        f_del_c += each_idx_indel_info(rm_ins_align_read, rm_gap_align_ref)[1]

            ref_idx_indel[str(ref_idx)+"_ins-count"]=f_ins_c
            ref_idx_indel[str(int(ref_idx)+1)+"_del-count"]=f_del_c


    ### ====== output each ref index position insertion and deletion info ====== ###
    out_indel_info = open(out_event_indel,"w")
    out_indel_info.write("Idx"+"\t"+"Event_type"+"\t"+"c_Reads"+"\n")

    print(len(ref_idx_indel))
    for i in ref_idx_indel.keys():
        if int(i.split("_")[0]) <= len(ref_idx_indel)/2:
            if int(i.split("_")[0]) == 1:
                out_indel_info.write(i.split("_")[0]+"\t"+i.split("_")[1]+"\t"+str(ref_idx_indel[i])+"\n")
            elif int(i.split("_")[0]) > 1:
                out_indel_info.write(str(int(i.split("_")[0])-1)+"\t"+i.split("_")[1]+"\t"+str(ref_idx_indel[i])+"\n")
    out_indel_info.close()

    ### ====== output the insertion position and sequences length ====== ###
    out_ins_info = open(out_ins,"w")
    out_ins_info.write("Ins_Ref_idx"+"\t"+"Length_ins(bp)"+"\t"+"c_Reads"+"\t"+"Ins_seq"+"\n")
    for i in ins_info.keys():
        ins_ref_idx = i.split("_")[0]
        ins_l = i.split("_")[1].split("bp")[0]
        c_reads = str(ins_info[i])
        seqs = str(i.split("_")[2])
        if int(ins_ref_idx) == 1:
            out_ins_info.write(ins_ref_idx+"\t"+ins_l+"\t"+c_reads+"\t"+seqs+"\n")
        elif int(ins_ref_idx) > 1:
            out_ins_info.write(str(int(ins_ref_idx)-1)+"\t"+ins_l+"\t"+c_reads+"\t"+seqs+"\n")
    out_ins_info.close()

    ### ====== output the deletion region and sequences length ====== ###
    out_del_info = open(out_del,"w")
    out_del_info.write("Del_Ref_Start"+"\t"+"Del_Ref_End"+"\t"+"Length_del(bp)"+"\t"+"c_Reads"+"\n")
    for i in del_info.keys():
        del_start = i.split("_")[0].split("-")[0]
        del_end = i.split("_")[0].split("-")[1]
        del_l = i.split("_")[1].split("bp")[0]
        c_reads = str(del_info[i])
        out_del_info.write(del_start+"\t"+del_end+"\t"+del_l+"\t"+c_reads+"\n")
    out_del_info.close()

    ### ====== output the only idx indel results ====== ###
    idx_ins_count = {}
    idx_del_count = {}

    with open(out_ins, "r") as f:
        rows = f.readlines()[1:]
        for line in rows:
            line = line.strip().split("\t")

            idx = line[0]
            cread = int(line[2])

            if idx in idx_ins_count.keys():
                idx_ins_count[idx] = idx_ins_count[idx] + cread
            elif idx not in idx_ins_count.keys():
                idx_ins_count[idx] = cread

    with open(out_del, "r") as f:
        rows = f.readlines()[1:]
        for line in rows:
            line = line.strip().split("\t")

            idx = line[0]
            cread = int(line[3])

            if idx in idx_del_count.keys():
                idx_del_count[idx] = idx_del_count[idx] + cread
            elif idx not in idx_del_count.keys():
                idx_del_count[idx] = cread

    final_count = open(out_idx_indel, "w")
    final_count.write("idx" + "\t" + "ins_count" + "\t" + "del_count" + "\n")

    for i in ref_dict.keys():
        ref_seq = ref_dict[i]
    ref_len = len(ref_seq)+1

    for i in range(1, ref_len):
        if str(i) in idx_ins_count.keys():
            if str(i) in idx_del_count.keys():
                final_count.write(str(i) + "\t" + str(idx_ins_count[str(i)]) + "\t" + str(idx_del_count[str(i)]) + "\n")
            elif str(i) not in idx_del_count.keys():
                final_count.write(str(i) + "\t" + str(idx_ins_count[str(i)]) + "\t" + str(0) + "\n")
        elif str(i) not in idx_ins_count.keys():
            if str(i) in idx_del_count.keys():
                final_count.write(str(i) + "\t" + str(0) + "\t" + str(idx_del_count[str(i)]) + "\n")
            elif str(i) not in idx_del_count.keys():
                final_count.write(str(i) + "\t" + str(0) + "\t" + str(0) + "\n")

    final_count.close()
