Collection of HMMs for the protein transporters from http://www.tcdb.org/

Software Requirements

* hmmer-3.1b2
* mafft v7.310
* cdhit v4.7
* R 3.4.0 (with ggplot2)

Steps to generate the HMMs.

1. Download the sequences from TCDB, in the folder fasta/all, create it if needed.
  We used two sets of sequences a) http://www.tcdb.org/public/tcdb and the fasta sequences associated to the IDs in b) http://www.tcdb.org/public/refseq.tsv, all these sequences are stored in file fasta/all/tcdb_plus_refseq.fasta

2. Prepare a file that maps the sequence identitifer to the TCDB #TC, up to the four first places. This file should look like this (mapping/seq2family.txt):
 1001796365	4.F.1.1
 1002048002	2.A.4.2
 1002048004	2.A.4.5
 1092917156	1.C.4.5
 159164663	1.A.59.1
 222437377	1.M.1.3
 24418985	1.G.3.1
 254071257	1.A.50.4
 327179101	1.A.92.1
 357380318	1.A.58.1

You can use this recipe to generate the file (From the folder where tcdb.fasta and refseq.tsv are located) i.e., fasta/all:

  grep ">" tcdb.fasta |sed 's/>//'|sed 's/  / /g'|cut -f 1,2 -d' '|sed 's/ /	/' |sed 's/	\(\([0-9]*\.[A-Z]*\.[0-9]*\.[0-9]*\)\.[0-9]*\)/	\2/' > seq2family_tcdb.temp.txt
  cut -f 1,2 refseq.tsv |sed 's/	\(\([0-9]*\.[A-Z]*\.[0-9]*\.[0-9]*\)\.[0-9]*\)/	\2/' > seq2family_refseq.temp.txt
  mkdir ../../mapping/
  cat seq2family_refseq.temp.txt seq2family_tcdb.temp.txt |sort -u > ../../mapping/seq2family.txt
  rm seq2family_tcdb.temp.txt seq2family_refseq.temp.txt

4. Generate one fasta file per group, with all the sequence in that group. Run the following in the base dir.
  ./scripts/getSequencesPerGroup.pl -f fasta/all/tcdb_plus_refseq.fasta -i mapping/seq2family.txt
  mv tcdb_group_* fasta/seeds/

5. Run the script createHMMFromDataSplit.pl for each of the fasta files in (fasta/seeds/). 
  for seedFile in $(ls -1 *.fasta); do seedName=${seedFile/tcdb_group_};seedName=${seedName/\.fasta};../../scripts/createHMMFromDataSplit.pl --truepositives $seedFile --all ../all/tcdb_plus_refseq.fasta --name $seedName; done

