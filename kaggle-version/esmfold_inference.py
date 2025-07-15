from transformers import AutoTokenizer, EsmForProteinFolding
from transformers.models.esm.openfold_utils.protein import to_pdb, Protein as OFProtein
from transformers.models.esm.openfold_utils.feats import atom14_to_atom37
import torch
import argparse
from Bio import SeqIO
import os

def esmfold(test_protein):
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
    model = model.cuda()
    model.esm = model.esm.half()
    torch.backends.cuda.matmul.allow_tf32 = True
    model.trunk.set_chunk_size(64)
    tokenized_input = tokenizer([test_protein], return_tensors="pt", add_special_tokens=False)['input_ids']
    tokenized_input = tokenized_input.cuda()
    with torch.no_grad():
        outputs = model(tokenized_input)
    return outputs

def convert_outputs_to_pdb(outputs,outdir,pdbname):
    final_atom_positions = atom14_to_atom37(outputs["positions"][-1], outputs)
    outputs = {k: v.to("cpu").numpy() for k, v in outputs.items()}
    final_atom_positions = final_atom_positions.cpu().numpy()
    final_atom_mask = outputs["atom37_atom_exists"]
    pdbs = []
    for i in range(outputs["aatype"].shape[0]):
        aa = outputs["aatype"][i]
        pred_pos = final_atom_positions[i]
        mask = final_atom_mask[i]
        resid = outputs["residue_index"][i] + 1
        pred = OFProtein(
            aatype=aa,
            atom_positions=pred_pos,
            atom_mask=mask,
            residue_index=resid,
            b_factors=outputs["plddt"][i],
            chain_index=outputs["chain_index"][i] if "chain_index" in outputs else None,
        )
        pdbs.append(to_pdb(pred))
        full_pdb = "\n".join(pdbs)
        with open(os.path.join(outdir,pdbname), "w") as pdb_file:
            pdb_file.write(full_pdb)
    return pdbs

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='ESMFOLD protein structure prediction')
    parser.add_argument('-i', type=str, default="BL19_2.fasta",
                            help='input fasta name')
    parser.add_argument('-o', type=str, default="BL19_2.pdb",
                            help='outputput pdb name')
    args = parser.parse_args()
    for record in SeqIO.parse(args.i, "fasta"):
        seq_id = record.id
        sequence = str(record.seq)
        outputs = esmfold(sequence)
        pdbs= convert_outputs_to_pdb(outputs,args.o,seq_id+".pdb")

