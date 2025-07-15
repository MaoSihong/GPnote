'''
WorkFlow for gutphage annotation
Author:
    MAO Sihong, Yang Yaoyu; Xianglab, Tsinghua University
Date:
    2022/11/10
'''
import numpy as np
import pandas as pd

'''
three dependant software/tools are needed for running the script for full anno functions of this workflow
    they are :
                hhsuite for HMM FUNCTION INFERENCE with PHROGs database
                hmmsearch for HMM FUNCTION INFERENCE with PFAM database
                VIGA (cite from https://github.com/EGTortuero/viga.git) for workflow with DIAMOND/BLAST+ anno \
                    HMM FUNCTION INFERENCE with PVOGs and VOGs databases
'''
import utils
import argparse
import os
import os.path as osp
import time
import shutil
logo=r''' 
    ___
   /   \ 
  /__ __\   GGGG   PPPP   N   N  OOO TTTTTTT EEEEE
  |\/ \/|  G       P   P  NN  N O   O   T    E
  | ___ |  G  GG   PPPP   N N N O   O   T    EEEE
  \     /  G   G   P      N  NN O   O   T    E
   \___/   GGGG    P      N   N  OOO    T    EEEEE
  __ | __
 /  / \  \ 
/  /   \  \
'''

if __name__ == "__main__":

    print(logo)
    print("Xiangye Lab, Tsinghua university")
    '''
    assign the parameters for run the pipeline
    '''
    parser = argparse.ArgumentParser(description='WorkFlow for GPnote (Gut phage annotation)')
    parser.add_argument('--workdir', type=str, default="/mnt/c/Users/Maosihong/Desktop/Xianglab/gut_phage/BL19_2",
                        help='workpath that contains the raw fasta file of the phage sequence')
    parser.add_argument('--fastafile', type=str, default="BL19_2.fasta",
                        help='fasta file name')
    parser.add_argument('--vigapath', type=str, default="/mnt/g/viga_dev",
                        help='vigapath on your computer or cluster')
    parser.add_argument('--pfampath', type=str, default="/mnt/g/Pfam/Pfam-A.hmm",
                        help='pfampath on your computer or cluster')
    parser.add_argument('--phrogspath', type=str, default="/mnt/g/phrogs_hhsuite_db/phrogs",
                        help='phrogspath on your computer or cluster')
    parser.add_argument('--esmfoldpath', type=str, default="/mnt/g/esmfold/",
                        help='esmfoldpath on your computer or cluster')
    parser.add_argument('--foldseekpath', type=str, default="/mnt/g/foldseek/pdbdb/pdb",
                        help='foldseekpath on your computer or cluster')
    
    args = parser.parse_args()
    '''
    generate the command and make sh script for the pipeline
    AKA stage 1
    '''
    
    prefix = args.fastafile.split(".")[0]
    vigacommand = utils.make_viga_command(args.workdir, args.fastafile, vigapath=args.vigapath)
    pfamcommand = utils.make_HMM_PFAM_command(args.workdir, prefix + "_annotated.genes.faa", pfampath=args.pfampath)
    phrogscommand = utils.make_HMM_PHROG_command(args.workdir, prefix + "_annotated.genes.faa", phrogspath=args.phrogspath)
    esmfoldcommand = utils.make_ESMFOLD_command(args.workdir, prefix + "_annotated.genes.faa", esmfoldpath=args.esmfoldpath)
    foldseekcommand= utils.make_FOLDSEEK_command(args.workdir,esmfoldpath=args.esmfoldpath, foldseekpath=args.foldseekpath)
    shfilepath= osp.join(args.workdir, "{}.workflow.sh".format(prefix))
    shfile = open(shfilepath, "w")
    shfile.write(vigacommand + "\n")
    shfile.write(pfamcommand + "\n")
    shfile.write(phrogscommand + "\n")
    shfile.write(esmfoldcommand + "\n")
    shfile.write(foldseekcommand + "\n")
    shfile.close()
    shutil.copy2("../resources/vog.annotations.tsv",osp.join(args.workdir,"vog.annotations.tsv"))
    shutil.copy2("../resources/modifiers_c11.txt",osp.join(args.workdir,"modifiers_c11.txt"))

    print("WORKFLOW STARTS!")
    st = time.time()
    os.system("bash {}".format(shfilepath))
    print("WORKFLOW STAGE1 ENDS WITHIN {:.2f}".format(time.time()-st))

    '''
    extract and prune the result and generate the read-friendly result 
    AKA stage 2
    the useful files include:
        {prefix}_annotated.csv from viga
        {prefix}.pfam.tbl from pfam HMM
        {prefix}.phrogs.tsv.ffdata from phrog HMM
        {prefix}.foldseek.tsv for foldseek
    '''
    os.chdir(args.workdir)
    try:
        os.rename(osp.join(args.workdir, prefix + "_annotated.csv"), osp.join(args.workdir, prefix + "_annotated.tsv"))
    except:
        print("no need to rename")
    temp_df = pd.read_csv(prefix + "_annotated.tsv", sep="\t")
    pfam_dict = utils.parse_pfam_result(prefix + ".pfam.tbl")
    phrog_dict = utils.parse_phrog_result(prefix + ".phrogs.tsv.ffdata", osp.join("/".join(args.phrogspath.split("/")[:-1]), "phrog_annot_v4.tsv"))
    foldseek_dict = utils.parse_foldseek_result(prefix + ".foldseek.tsv")
    print(temp_df)
    print(type(temp_df))
    # print(pfam_dict)
    # print(phrog_dict)
    pfamid = []
    pfamanno = []
    pfamev = []
    phrogid = []
    phroganno = []
    phrogidentity = []
    phrogev = []
    foldseekid = []
    foldseekanno=[]
    foldseekqtm = []
    foldseekttm = []
    foldseekatm=[]
    foldseekev = []
    foldseekrmsd = []
    for loc in temp_df["Protein ID"]:
        pfam_query = pfam_dict.setdefault(loc, {"pfam_ID": [], "pfam_anno": [], "e-value": []})
        phrog_query = phrog_dict.setdefault(loc, {"phrog_ID": "", "phrog_anno": "", "phrog_identity": np.nan, "e-value": np.nan})
        foldseek_query=foldseek_dict.setdefault(loc, {"foldseek_ID": [],"foldseek_anno": [], "foldseek_qtm": [], "foldseek_ttm": [], "foldseek_average_tm": [], "foldseek_evalue": [], "foldseek_rmsd": []})
        pfamid.append(",".join(pfam_query["pfam_ID"]))
        pfamanno.append(",".join(pfam_query["pfam_anno"]))
        pfamev.append(",".join([str(x) for x in pfam_query["e-value"]]))
        phrogid.append(phrog_query["phrog_ID"])
        phroganno.append(phrog_query["phrog_anno"])
        phrogidentity.append(str(phrog_query["phrog_identity"]))
        phrogev.append(str(phrog_query["e-value"]))
        foldseekid.append(",".join(foldseek_query["foldseek_ID"]))
        foldseekanno.append(",".join(foldseek_query["foldseek_anno"]))
        foldseekqtm.append(",".join([str(x) for x in foldseek_query["foldseek_qtm"]]))
        foldseekttm.append(",".join([str(x) for x in foldseek_query["foldseek_ttm"]]))
        foldseekatm.append(",".join([str(x) for x in foldseek_query["foldseek_average_tm"]]))
        foldseekev.append(",".join([str(x) for x in foldseek_query["foldseek_evalue"]]))
        foldseekrmsd.append(",".join([str(x) for x in foldseek_query["foldseek_rmsd"]]))
    temp_df.insert(len(temp_df.columns), "pfam_ID", pfamid)
    temp_df.insert(len(temp_df.columns), "pfam_anno", pfamanno)
    temp_df.insert(len(temp_df.columns), "pfam_e-value", pfamev)
    temp_df.insert(len(temp_df.columns), "phrog_ID", phrogid)
    temp_df.insert(len(temp_df.columns), "phrog_anno", phroganno)
    temp_df.insert(len(temp_df.columns), "phrog_identity", phrogidentity)
    temp_df.insert(len(temp_df.columns), "phrog_e-value", phrogev)
    temp_df.insert(len(temp_df.columns), "foldseek_ID", foldseekid)
    temp_df.insert(len(temp_df.columns), "foldseek_anno", foldseekanno)
    temp_df.insert(len(temp_df.columns), "foldseek_qtm", foldseekqtm)
    temp_df.insert(len(temp_df.columns), "foldseek_ttm", foldseekttm)
    temp_df.insert(len(temp_df.columns), "foldseek_average_tm", foldseekatm)
    temp_df.insert(len(temp_df.columns), "foldseek_evalue", foldseekev)
    temp_df.insert(len(temp_df.columns), "foldseek_rmsd", foldseekrmsd)
    print(temp_df)
    temp_df.to_csv(prefix + "_anno_result.tsv", sep="\t", index=False)
    #move tempfiles in viga src directory to avoid multi treatment
    tempfiles = [x for x in os.listdir(args.vigapath) if (x.endswith(".fasta") or x.endswith(".fna") or x.endswith(".ptt") or x.endswith(".csv") or x.startswith("logfile"))]
    print(tempfiles)
    for file in tempfiles:
        os.remove(osp.join(args.vigapath, file))
        print("{} removed from {}".format(file, args.vigapath))
    print([x for x in os.listdir(args.vigapath) if (x.endswith(".fasta") or x.endswith(".fna") or x.endswith(".ptt") or x.endswith(".csv") or x.startswith("logfile"))] == [])
    print("WORKFLOW ENDS WITHIN {:.2f} secs".format(time.time() - st))
    make_annotation_file(prefix + "_anno_result.tsv",prefix + "_gpnote.txt")
    make_gbk_file(prefix + "_gpnote.txt",prefix+"_annotated.gbk")
    make_gff_file(prefix + "_gpnote.txt",prefix+"_annotated.gff")
