'''
WorkFlow for gutphage annotation
Author:
    MAO Sihong, Yang Yaoyu; Xianglab, Tsinghua University
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
    vigacommand = "python3 {}/VIGA_temp.py --input {} --diamonddb {} --blastdb {} --hmmerdb {} --rfamdb {} --modifiers {}".format(workdir,osp.join(workdir, fastafile),
                                                                                                                               os.path.join(vigapath,"DIAMOND_Refseq/refseq_viral_proteins"),
                                                                                                                               os.path.join(vigapath,"BLAST_Refseq/refseq_viral_proteins"),
                                                                                                                               os.path.join(vigapath,"HMMER_VOGs/pvogs_vogs_RVDB.hmm"),
                                                                                                                               os.path.join(vigapath,"Rfam/Rfam.cm"),
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

'''
use ESMFOLD and FOLDSEEK for anntation of phage
'''

def make_ESMFOLD_command(workdir, faafile, esmfoldpath="/mnt/g/esmfold/")-> str:
    os.system("mkdir "+esmfoldpath)
    chdir = "cd {}".format(workdir)
    prefix = workdir.split("/")[-1]
    esmfoldcommand="python esmfold_inference.py -i {}  -o {}".format(faafile,esmfoldpath)
    print("ESMFOLD command like:\n==========================")
    print(esmfoldcommand)
    return chdir + "\n" + esmfoldcommand

def make_FOLDSEEK_command(workdir,esmfoldpath, foldseekpath="/mnt/g/foldseek/pdbdb/pdb")-> str:
    chdir = "cd {}".format(workdir)
    prefix = workdir.split("/")[-1]
    foldseekcommand="foldseek easy-search {} {} {} tmp --format-output 'query,target,theader,pident,alntmscore,qtmscore,ttmscore,evalue,prob,rmsd'".format(esmfoldpath,foldseekpath,prefix+ ".foldseek.tsv")
    print("FOLDSEEK command like:\n==========================")
    print(foldseekcommand)
    return chdir + "\n" + foldseekcommand


###########################Result processing###############################
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

def parse_foldseek_result(foldseekresult):
    resultfile=open(foldseekresult,"r")
    content = [x for x in resultfile.readlines()]
    def make_parse(row):
        row_split = row.strip().split("\t")
        parse_dict = {"LOC": row_split[0], "foldseek_anno": [row_split[2]], "foldseek_ID": [row_split[1]],"foldseek_qtm": [float(row_split[5])],"foldseek_ttm": [float(row_split[6])],"foldseek_average_tm": [ (float(row_split[5]) + float(row_split[6])) / 2 ],"foldseek_evalue":[float(row_split[7])],"foldseek_rmsd":[float(row_split[9])]}
        return parse_dict
    total_dict = {}
    for row in content:
        if (float(row.strip().split("\t")[5])+ float(row.strip().split("\t")[6]))<1 or float(row.strip().split("\t")[7])>=0.01:
            continue
        else:
            parse_dict = make_parse(row)
            total_dict[parse_dict["LOC"]] = total_dict.setdefault(parse_dict["LOC"], {})
            total_dict[parse_dict["LOC"]]["foldseek_anno"] = total_dict[parse_dict["LOC"]].setdefault("foldseek_anno", []) + parse_dict["foldseek_anno"]
            total_dict[parse_dict["LOC"]]["foldseek_ID"] = total_dict[parse_dict["LOC"]].setdefault("foldseek_ID", []) + parse_dict["foldseek_ID"]
            total_dict[parse_dict["LOC"]]["foldseek_qtm"] = total_dict[parse_dict["LOC"]].setdefault("foldseek_qtm", []) + parse_dict["foldseek_qtm"]
            total_dict[parse_dict["LOC"]]["foldseek_ttm"] = total_dict[parse_dict["LOC"]].setdefault("foldseek_ttm", []) + parse_dict["foldseek_ttm"]
            total_dict[parse_dict["LOC"]]["foldseek_average_tm"] = total_dict[parse_dict["LOC"]].setdefault("foldseek_average_tm", []) + parse_dict["foldseek_average_tm"]
            total_dict[parse_dict["LOC"]]["foldseek_evalue"] = total_dict[parse_dict["LOC"]].setdefault("foldseek_evalue", []) + parse_dict["foldseek_evalue"]
            total_dict[parse_dict["LOC"]]["foldseek_rmsd"] = total_dict[parse_dict["LOC"]].setdefault("foldseek_rmsd", []) + parse_dict["foldseek_rmsd"]
    return total_dict
        
def make_annotation_file(anno_result,outable):
    d0=pd.read_table("../resources/vog.annotations.tsv",header=0)
    dic0= dict(zip(d0["GroupName"],d0["ConsensusFunctionalDescription"]))
    df=pd.read_table(anno_result,header=0)
    df['BLAST_anno'] = df.apply(lambda row: f"[Homology BLAST #Rank1#]{row['Description']}" if row['Source'] != 'NO_HIT' else '', axis=1)
    df['VOG'] = df['HMMer'].map(dic0).fillna('')  # 用'Not Found'替换NaN
    df["VOG_anno"]=df.apply(lambda row: f"[hmm profile #Rank2#]{row['HMMer']} {row['VOG']}" if row['HMMer'] != 'NO' else '', axis=1)
    df["PHROG_anno"]=df.apply(lambda row: f"[hmm profile #Rank2#]{row['phrog_ID']} {row['phrog_anno']}" if row['phrog_ID'] != '' else '', axis=1)
    df["foldseekID_1"]=df['foldseek_ID'].apply(lambda x: x.split(',') if pd.notna(x) and x != '' else [])
    df["foldseekanno_1"]=df['foldseek_anno'].apply(lambda x: x.split(',') if pd.notna(x) and x != '' else [])
    df['FOLDSEEK_anno'] = df.apply(lambda row: ["[Structure alignment #Rank3#]"+a +" "+ b for a, b in zip(row['foldseekID_1'], row['foldseekanno_1'])], axis=1)
    df['ANNO'] = df.apply(lambda row: [x for x in [row['BLAST_anno'],row['VOG_anno'],row['PHROG_anno']] + row['FOLDSEEK_anno'] if x != '' and x != []], axis=1)
    df['GPnote_anno'] = df['ANNO'].apply(lambda x: '; '.join([f"Annot{i+1}={item}" for i, item in enumerate(x)]) if len(x) > 0 else '')
    dfnew=df[["Protein ID","GPnote_anno"]]
    dfnew.to_csv(outable,sep="\t",index=None)
    return outable

def make_gbk_file(outable,gbk):
    gpdic=pd.read_table(outable,header=0)
    gpdic=gpdic.fillna("Hypothetical protein")
    gpdic2=gpdic.set_index('Protein ID')['GPnote_anno'].to_dict()
    f11=open(gbk)
    name=gbk.split("_annotated")[0]
    fo=open(gbk.replace("gbk","_gpnote.gbk"),"w")
    for line in f11:
        if "feature	gene" in line:
            newline=line.replace("LOC_1_",name+"gp")
            newline=newline.replace("LOC_1",name)
            fo.write(newline)
        elif "feature	CDS" in line:
            id=line.split("protein_id=")[1].split(";")[0]
            fullid=name+"_"+id
            nline=re.sub(r'locus_tag=[^;\n]*',"locus_tag="+id.replace("LOC_1_",name+"gp"),line)
            nline=re.sub(r'note=[^;\n]*','note='+str(gpdic2[fullid]),nline)
            if len(str(gpdic2[fullid]).split(";"))>1:
                newline=re.sub(r'product=[^;\n]*','product='+str(gpdic2[fullid].split(";")[0]),nline)
            else:
                newline=re.sub(r'product=[^;\n]*','product='+str(gpdic2[fullid]),nline)
            newline=newline.replace("LOC_1",name)
            fo.write(newline)
        else:
            newline=line.replace("LOC_1",name)
            fo.write(newline)
    fo.close()
    f11.close()
    return

def make_gff_file(outable,gff):
    gpdic=pd.read_table(outable,header=0)
    gpdic=gpdic.fillna("Hypothetical protein")
    gpdic2=gpdic.set_index('Protein ID')['GPnote_anno'].to_dict()
    f11=open(gff)
    name=gff.split("_annotated")[0]
    fo=open(gff.replace("gff","_gpnote.gff"),"w")
    for line in f11:
        if "feature	gene" in line:
            newline=line.replace("LOC_1_",name+"gp")
            newline=newline.replace("LOC_1",name)
            fo.write(newline)
        elif "feature	CDS" in line:
            id=line.split("protein_id=")[1].split(";")[0]
            fullid=name+"_"+id
            nline=re.sub(r'locus_tag=[^;\n]*',"locus_tag="+id.replace("LOC_1_",name+"gp"),line)
            nline=re.sub(r'note=[^;\n]*','note='+str(gpdic2[fullid]),nline)
            if len(str(gpdic2[fullid]).split(";"))>1:
                newline=re.sub(r'product=[^;\n]*','product='+str(gpdic2[fullid].split(";")[0]),nline)
            else:
                newline=re.sub(r'product=[^;\n]*','product='+str(gpdic2[fullid]),nline)
            newline=newline.replace("LOC_1",name)
            fo.write(newline)
        else:
            newline=line.replace("LOC_1",name)
            fo.write(newline)
    fo.close()
    f11.close()
    return 
