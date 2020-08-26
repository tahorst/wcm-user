#! /usr/bin/env bash

project_root=$(git rev-parse --show-toplevel)
regulon_db_dir="${project_root}/nca/data/regulon-db"
mkdir -p $regulon_db_dir
cd $regulon_db_dir

wget http://regulondb.ccg.unam.mx/menu/download/datasets/files/network_tf_gene.txt -O tf_genes.tsv
wget http://regulondb.ccg.unam.mx/menu/download/datasets/files/network_sigma_gene.txt -O sigma_genes.tsv
