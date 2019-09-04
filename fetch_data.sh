#!/usr/bin/env bash


#myfile=$(date +"%F_all_Fungi")
myfile=$(date +"%F_endophytic_fungi")


### download necessary scripts
# NCBI analysis
if test -f ./scripts/search_ncbi_gi_by_term.py
then
  echo "./scripts/search_ncbi_gi_by_term.py found"
else
  git clone https://github.com/jgmv/ncbi_data_analysis.git scripts
fi

# bash scripts
if test -f ./scripts/removeTaxonTag.sh
then
  echo "./scripts/removeTaxonTag.sh found"
else
  git clone https://github.com/jgmv/bash_seq_analysis.git scripts
fi

# system-wide access to scripts
. scripts/*.sh
export PATH="$PATH:scripts"


### fetch data
# retrieve GIs for all fungal ITS from GenBank
#python scripts/search_ncbi_gi_by_term.py -o $myfile.gi 'Fungi[Organism] and (internal transcribed spacer OR ITS OR ITS1 OR ITS2) not chromosome not Neocallimastigales'

# retrieve GIs for all fungal endophytes ITS from GenBank
python scripts/search_ncbi_gi_by_term.py -o $myfile.gi 'Fungi[Organism] and (internal transcribed spacer OR ITS OR ITS1 OR ITS2) not chromosome and 100:2000[Sequence Length] and endoph*[ALL]'

# get metadata from GIs (takes time, better use tool in https://www.ncbi.nlm.nih.gov/sites/batchentrez)
python fetch_gb.py -o $myfile.gb $myfile.gi

# convert GenBank data to table
python get_metadata_from_gb.py -o $myfile.csv $myfile.gb

# get fasta sequences from gb file
python get_fasta_from_gb.py -o $myfile.fasta $myfile.gb

mkdir -p files
mv $myfile* files


### identify ITS sequences
# download last version of UNITE ITS reference dataset, mothur release
# check https://unite.ut.ee/repository.php
if test -f ./data/UNITE_sh_dynamic.tax
then
  echo "UNITE database found"
else
  mkdir -p data
  wget -P data https://files.plutof.ut.ee/public/orig/56/25/5625BDC830DC246F5B8C7004220089E032CC33EEF515C76CD0D92F25BDFA9F78.zip
  unzip data/*.zip
  rm data/*.zip
  mv data/UNITE*_sh_dynamic.fasta data/UNITE_sh_dynamic.fasta
  mv data/UNITE*_sh_dynamic.tax data/UNITE_sh_dynamic.tax
  rm data/UNITEv*
fi

# identify sequences with mothur
mothur "#classify.seqs(fasta=files/20190726_endophytic_Fungi.fasta,\
  template=data/UNITEv6_sh_dynamic_s.fasta,\
  taxonomy=data/UNITEv6_sh_dynamic_s.tax, cutoff=60, probs=T)" 
removeTaxonTag files/20190726_endophytic_Fungi.UNITEv6_sh_dynamic_s.wang.taxonomy \
  files/taxonomy_boot.csv
cp files/taxonomy_boot.csv files/taxonomy.csv
sed -i 's/([^()]*)//g' files/taxonomy.csv


### end
