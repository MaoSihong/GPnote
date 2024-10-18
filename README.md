# GPnote
  Phage annotation pipeline taking advantage of protein structure comparison
# General
![GPnote Pipeline Overview](./images/GPnote_workflow.png)

  We developed a workflow called Gut Phage genome annotation note (GPnote) for phage genome annotation. GPnote is an automated de novo genome annotation pipeline that integrates both the sequence-based and structure-based predictions of gene functions. 

  we employed **VIGA** (Gonzalez-Tortuero et al., 2018) as our annotation start point(which provides the BLAST Hit of RefSeqViral and VOGs/pVOGs Hit by profile HMM search) and further incorprated the HMM profile comparison Hits in **PHROGs** via hhsuite.Also we did hmmsearch against PFam database for more function indications. Besides all these sequence-based strategies, we achieved advanced structure-based annotation by **ESM-Fold** (Lin Zeming et al., 2023) structure prediction and **FoldSeek** (Michel van Kempen et al., 2023) 3Di structure search against the **AFUniProt50** database.
