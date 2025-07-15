'''
WorkFlow for gutphage annotation
Author:
    MAO Sihong, Xianglab
Date:
    2022/11/10
'''
import multiprocessing.pool as Pool
import os
import argparse
import os.path as osp

def single_process(prefix, batchdir):
    print("NOW WORK ON {} annotation".format(prefix))
    try:
        workdir = osp.join(batchdir, prefix)
        os.mkdir(workdir)
        os.rename(osp.join(batchdir, prefix + ".fasta"), osp.join(workdir, prefix + ".fasta"))
    except:
        print("no need to move {}".format(prefix + ".fasta"))
    if not prefix + ".fasta" in os.listdir(workdir):
        return False
    os.system(
        "python3 ./gpnote.py --workdir {} --fastafile {}".format(
            workdir, prefix + ".fasta"))
    return True

parser = argparse.ArgumentParser(description='WorkFlow for Gut Phage Annotation batch in a directory')
parser.add_argument('--batchdir', type=str, default="/mnt/c/Users/Maosihong/Desktop/Xianglab/gut_phage/",
                    help='batchdir that contains the all raw fasta file of the phage sequence')
parser.add_argument('--multicore', type=int, default=8,
                    help='cores used for batch running')
args = parser.parse_args()

fastafiles = [x for x in os.listdir(args.batchdir) if x.endswith(".fasta")]
fastadirs = [x for x in os.listdir(args.batchdir) if not x.endswith(".fasta")]
allprefix = []
for fasta in fastafiles:
    allprefix.append(fasta.split(".")[0])
for fasta in fastadirs:
    allprefix.append(fasta.split(".")[0])
print("BATCH PLAN:")
print(allprefix)
process_status = []
for prefix in allprefix:
    status = single_process(prefix, args.batchdir)
    process_status.append(status)
print("=========== batch log ===========")
for i in range(len(allprefix)):
    print("{} : {}".format(allprefix[i], ("True" if process_status[i] else "False")))
print("=========== batch log ===========")
