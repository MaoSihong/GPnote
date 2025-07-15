#!/bin/bash
#download GPnote sequence database
wget  https://cloud.tsinghua.edu.cn/seafhttp/files/a25123e2-ff0c-45b9-96aa-3f66dbdbd9f0/GPnote_seuqence_db.tar.gz 
tar zxvf GPnote_seuqence_db.tar.gz
cd GPnote_db/
mkdir foldseek_db
#download GPnote structure database
conda activate GPnote
foldseek databases Alphafold/UniProt50 foldseek_db/afdb tmp

