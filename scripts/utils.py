'''
WorkFlow for gutphage annotation
Author:
    MAO Sihong, Xianglab
Date:
    2022/11/10
'''
import os
import os.path as osp

import numpy as np
import pandas as pd

'''
use viga for anntation of phage
eg:
cd /mnt/g/viga_dev/
python3 /mnt/g/viga_dev/VIGA_temp.py --input /mnt/c/Users/Maosihong/Desktop/BG16_2/BG16_2.fasta 
                                    --diamonddb databases/RefSeq_Viral_DIAMOND/refseq_viral_proteins 
                                    --blastdb databases/RefSeq_Viral_BLAST/refseq_viral_proteins 
                                    --hmmerdb databases/pvogs_rvdb/pvogs_vogs_RVDB.hmm 
                                    --rfamdb databases/rfam/Rfam.cm --modifiers modifiers_c11.txt
'''

def make_viga_command(workdir, fastafile, vigapath="/mnt/g/viga_dev",modifier="c11") -> str:
    chdir = "cd {}".format(vigapath)
    vigacommand = "python3.8 VIGA_temp.py --input {} --diamonddb {} --blastdb {} --hmmerdb {} --rfamdb {} --modifiers {}".format(osp.join(workdir, fastafile),
                                                                                                                               "databases/RefSeq_Viral_DIAMOND/refseq_viral_proteins",
                                                                                                                               "databases/RefSeq_Viral_BLAST/refseq_viral_proteins",
                                                                                                                               "databases/pvogs_rvdb/pvogs_vogs_RVDB.hmm",
                                                                                                                               "databases/rfam/Rfam.cm",
                                                                                                                               "modifiers_{}.txt".format(modifier))
    print("VIGA command like:\n==========================")
    print(vigacommand)
    return chdir + "\n" + vigacommand
#print(make_viga_command("/mnt/c/Users/Maosihong/Desktop/BG16_2", "BG16_2.fasta"))

'''
use HMMsearch for anntation of phage
eg:
cd /mnt/c/Users/Maosihong/Desktop/gut_phage_all_fasta/7W
hmmsearch -o 7W.pfam.out --tblout 7W.pfam.tbl 
        --domtblout 7W.pfam.dom.tbl --pfamtblout 7W.pfam.ptbl 
        --noali --notextw /mnt/g/Pfam/Pfam-A.hmm 7W_annotated.genes.faa 
'''

def make_HMM_PFAM_command(workdir, faafile, pfampath = "/mnt/g/Pfam/Pfam-A.hmm") -> str:
    chdir = "cd {}".format(workdir)
    prefix = workdir.split("/")[-1]
    pfamcommand = "hmmsearch -o {} --tblout {} --domtblout {} --pfamtblout {} --noali --notextw {} {}".format("{}.pfam.out".format(prefix),
                                                                                                              "{}.pfam.tbl".format(prefix),
                                                                                                              "{}.pfam.dom.tbl".format(prefix),
                                                                                                              "{}.pfam.ptbl".format(prefix),
                                                                                                              pfampath,
                                                                                                              faafile)
    print("PFAM HMMsearch command like:\n==========================")
    print(pfamcommand)
    return chdir + "\n" + pfamcommand

'''
use HHsuite for anntation of phage
eg:
ffindex_from_fasta -s your_multifasta.ff{data,index} your_multifasta.faa
hhsearch_omp -i /mnt/c/Users/Maosihong/Desktop/BG16_2/c11/BG16_2_annotated.genes -d /mnt/g/phrogs_hhsuite_db/phrogs -o /mnt/c/Users/Maosihong/Desktop/BG16_2/BG16_2.phrogs.test -blasttab /mnt/c/Users/Maosihong/Desktop/BG16_2/BG16_2.phrogs.tsv.test
'''

def make_HMM_PHROG_command(workdir, faafile, phrogspath="/mnt/g/phrogs_hhsuite_db/phrogs") -> str:
    chdir = "cd {}".format(workdir)
    ffindexfile = ".".join(faafile.split(".")[0:2])
    prefix = workdir.split("/")[-1]
    ffindexcommand = "ffindex_from_fasta -s %s.ff{data,index} %s"%(ffindexfile, faafile)
    hhsearchcommand = "hhsearch_omp -i {} -d {} -o {} -blasttab {}".format(ffindexfile, phrogspath, prefix + ".phrogs", prefix + ".phrogs.tsv")
    print("PHROG HMMsearch command like:\n==========================")
    total_command = ffindexcommand + "\n" + hhsearchcommand
    print(total_command)
    return chdir + "\n" + total_command

def parse_pfam_result(pfamresult):
    resultfile = open(pfamresult, "r")
    # header = "target_name        target_accession  query_name           query_accession    E-value_full  score_full  bias_full   E-value_dom  score_dom  bias_dom   exp reg clu  ov env dom rep inc description of target"
    content = [x for x in resultfile.readlines() if not x.startswith("#")]
    def make_parse(row):
        row_split = row.strip().split()
        # print(row_split)
        parse_dict = {"LOC": row_split[0], "pfam_anno": [row_split[2]], "pfam_ID": [row_split[3]], "e-value":[float(row_split[4])]}
        return parse_dict
    total_dict = {}
    for row in content:
        parse_dict = make_parse(row)
        total_dict[parse_dict["LOC"]] = total_dict.setdefault(parse_dict["LOC"], {})
        total_dict[parse_dict["LOC"]]["pfam_anno"] = total_dict[parse_dict["LOC"]].setdefault("pfam_anno", []) + parse_dict["pfam_anno"]
        total_dict[parse_dict["LOC"]]["pfam_ID"] = total_dict[parse_dict["LOC"]].setdefault("pfam_ID", []) + parse_dict["pfam_ID"]
        total_dict[parse_dict["LOC"]]["e-value"] = total_dict[parse_dict["LOC"]].setdefault("e-value", []) + parse_dict["e-value"]
    return total_dict

def parse_phrog_result(phrogresult,phroganno):
    resultfile = open(phrogresult, "r")
    content = [x for x in resultfile.readlines() if len(x.split("\t")) == 12]
    def make_parse(row):
        if not row.startswith("LOC"):
            row = row[1:]
        row_split = row.strip().split("\t")
        parse_dict = {"LOC": row_split[0], "phrog_ID": row_split[1], "phrog_identity": float(row_split[2]), "e-value": float(row_split[-2])}
        return parse_dict
    total_dict = {}
    for row in content:
        parse_dict = make_parse(row)
        total_dict[parse_dict["LOC"]] = total_dict.setdefault(parse_dict["LOC"], {"phrog_ID": "", "phrog_identity": 0, "e-value": np.Inf})
        # if parse_dict["LOC"]=="LOC_1_24":
        #     print(parse_dict)
        #     print(total_dict[parse_dict["LOC"]])
        if parse_dict["e-value"] < total_dict[parse_dict["LOC"]]["e-value"]:
            total_dict[parse_dict["LOC"]]["phrog_ID"] = parse_dict["phrog_ID"]
            total_dict[parse_dict["LOC"]]["phrog_identity"] = parse_dict["phrog_identity"]
            total_dict[parse_dict["LOC"]]["e-value"] = parse_dict["e-value"]
    phroganno_df = pd.read_csv(phroganno, sep = "\t", header=0)
    # print(total_dict)
    # print(phroganno_df)
    phroganno_df["phrog"] = phroganno_df["phrog"].apply(lambda x: "phrog_" + str(x))
    # print(phroganno_df)
    for loc in total_dict.keys():
        # print(loc)
        lookup_anno = phroganno_df[phroganno_df["phrog"] == total_dict[loc]["phrog_ID"]]["annot"]
        # print(lookup_anno)
        # print(total_dict[loc])
        total_dict[loc]["phrog_anno"] = lookup_anno.values[0]
    return total_dict

