'''
WorkFlow for gutphage annotation
Author:
    MAO Sihong, Xianglab, THU
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

if __name__ == "__main__":
    '''
    assign the parameters for run the pipeline
    '''
    parser = argparse.ArgumentParser(description='WorkFlow for Gut Phage Annotation')
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
    args = parser.parse_args()
    '''
    generate the command and make sh script for the pipeline
    AKA stage 1
    '''
    prefix = args.fastafile.split(".")[0]
    vigacommand = utils.make_viga_command(args.workdir, args.fastafile, vigapath=args.vigapath)
    pfamcommand = utils.make_HMM_PFAM_command(args.workdir, prefix + "_annotated.genes.faa", pfampath=args.pfampath)
    phrogscommand = utils.make_HMM_PHROG_command(args.workdir, prefix + "_annotated.genes.faa", phrogspath=args.phrogspath)
    shfilepath= osp.join(args.workdir, "{}.workflow.sh".format(prefix))
    shfile = open(shfilepath, "w")
    shfile.write(vigacommand + "\n")
    shfile.write(pfamcommand + "\n")
    shfile.write(phrogscommand + "\n")
    shfile.close()
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
    '''
    os.chdir(args.workdir)
    try:
        os.rename(osp.join(args.workdir, prefix + "_annotated.csv"), osp.join(args.workdir, prefix + "_annotated.tsv"))
    except:
        print("no need to rename")
    temp_df = pd.read_csv(prefix + "_annotated.tsv", sep="\t")
    pfam_dict = utils.parse_pfam_result(prefix + ".pfam.tbl")
    phrog_dict = utils.parse_phrog_result(prefix + ".phrogs.tsv.ffdata", osp.join("/".join(args.phrogspath.split("/")[:-1]), "phrog_annot_v4.tsv"))
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
    for loc in temp_df["Protein ID"]:
        pfam_query = pfam_dict.setdefault(loc, {"pfam_ID": [], "pfam_anno": [], "e-value": []})
        phrog_query = phrog_dict.setdefault(loc, {"phrog_ID": "", "phrog_anno": "", "phrog_identity": np.nan, "e-value": np.nan})
        pfamid.append(",".join(pfam_query["pfam_ID"]))
        pfamanno.append(",".join(pfam_query["pfam_anno"]))
        pfamev.append(",".join([str(x) for x in pfam_query["e-value"]]))
        phrogid.append(phrog_query["phrog_ID"])
        phroganno.append(phrog_query["phrog_anno"])
        phrogidentity.append(str(phrog_query["phrog_identity"]))
        phrogev.append(str(phrog_query["e-value"]))
    temp_df.insert(len(temp_df.columns), "pfam_ID", pfamid)
    temp_df.insert(len(temp_df.columns), "pfam_anno", pfamanno)
    temp_df.insert(len(temp_df.columns), "pfam_e-value", pfamev)
    temp_df.insert(len(temp_df.columns), "phrog_ID", phrogid)
    temp_df.insert(len(temp_df.columns), "phrog_anno", phroganno)
    temp_df.insert(len(temp_df.columns), "phrog_identity", phrogidentity)
    temp_df.insert(len(temp_df.columns), "phrog_e-value", phrogev)
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