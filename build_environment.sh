#!/bin/bash
conda create -n GPnote
conda activate GPnote
conda install -y -c bioconda -c conda-forge bcbio-gff numpy pandas scipy lastz Aragorn infernal piler-cr Prodigal diamond blast hmmer hhsuite foldseek
pip install biopython==1.79
pip install "fair-esm[esmfold]"
pip install 'dllogger @ git+https://github.com/NVIDIA/dllogger.git'
pip install 'openfold @ git+https://github.com/aqlaboratory/openfold.git'