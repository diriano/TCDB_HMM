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

3. Generate one fasta file per group, with all the sequence in that group. Run the following in the base dir.
  ./scripts/getSequencesPerGroup.pl -f fasta/all/tcdb_plus_refseq.fasta -i mapping/seq2family.txt
  mv tcdb_group_* fasta/seeds/

4. Run the script createHMMFromDataSplit.pl for each of the fasta files in (fasta/seeds/). 
  for seedFile in $(ls -1 *.fasta); do seedName=${seedFile/tcdb_group_};seedName=${seedName/\.fasta};../../scripts/createHMMFromDataSplit.pl --truepositives $seedFile --all ../all/tcdb_plus_refseq.fasta --name $seedName; done

 This will take a few hours, go for a walk on the park!

5. The run of createHMMFromDataSplit.pl will create a *.selectGA.tbl file for each family, concatenate them in a single file:

  cat *selectGA.tbl > all.selGA.tbl #This is ran from within the fasta/seeds folder

 Then we need to check which families have a well defined GA, that is where the minimal score for TP is higher than the maximal score for TN. For these families we will keep the generated HMM and for the remaining the HMM will be discarded. The families with well defined GA, will have the string 'OK' in the column StatusGA in the file all.selGA.tbl
 
 Counting the number of families with not defined GA:

  cat *select*|grep -v Model|grep NotOK|wc -l #Results in 902 families. You results migth vary.

 And the number of families with defined GA:

  cat *select*|grep -v Model|grep -v NotOK|wc -l #Results in 1636 families

 Adding the GA line to the hmm file for the families with well defined GA:

  ../../scripts/addGA2HMMs.pl -g all.selGA.tbl #From within the fasta/seeds folder

 That will generate the *.mod.hmm files, which can be concatenated for later use.

  mkdir ../../HMMs
  cat *mod.hmm > ../../HMMs/TCDB_HMMs.hmm #From within the fasta/seeds folder

6. We try further with the families with not defined GA. We will cluster at 90% identity and for each cluster we run the whole process again.

  ../scripts/clusterSequences.pl -g all.selGA.tbl -i 0.9 #From within the fasta/seeds folder
  cd seeds_cdhit_0.9
  for seedFile in $(ls -1 *.fasta); do seedName=${seedFile/tcdb_group_};seedName=${seedName/\.fasta};../../../scripts/createHMMFromDataSplit.pl --truepositives $seedFile --all ../../all/tcdb_plus_refseq.fasta --name $seedName; done
  cat *selectGA.tbl > all.selGA.tbl
  ../../scripts/addGA2HMMs.pl -g all.selGA.tbl
  cat *mod.hmm >> ../../../HMMs/TCDB_HMMs.hmm
  
7. We try again with the families with not defined GA. We will cluster at 95% identity and for each cluster we run the whole process again.

  ../../scripts/clusterSequences.pl -g all.selGA.tbl -i 0.95 #From within the fasta/seeds/seeds_cdhit_0.9 folder
  cd seeds_cdhit_0.95
  for seedFile in $(ls -1 *.fasta); do seedName=${seedFile/tcdb_group_};seedName=${seedName/\.fasta};../../../../scripts/createHMMFromDataSplit.pl --truepositives $seedFile --all ../../../all/tcdb_plus_refseq.fasta --name $seedName; done
  cat *selectGA.tbl > all.selGA.tbl
  ../../scripts/addGA2HMMs.pl -g all.selGA.tbl
  cat *mod.hmm >> ../../../HMMs/TCDB_HMMs.hmm

8. One more try, clustering at 99%

  ../../scripts/clusterSequences.pl -g all.selGA.tbl -i 0.95 #From within the fasta/seeds/seeds_cdhit_0.9/seeds_cdhit_0.95/ folder
