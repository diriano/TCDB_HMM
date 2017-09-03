#!/usr/bin/perl

use strict;
use warnings;
use List::Util 'shuffle';
use Bio::DB::Fasta;
use Bio::SeqIO;
use Getopt::Long;
use File::Basename;

#Global vars
my $version="0.1";
my $license='';
my $help='';
my $truePositives='';
my $allSeqs='';
my $nameHmm='';
my %id_class;
my $debug=0;

#Get input from user
GetOptions(
    'license|l'        => \$license,
    'help|h|?'         => \$help,
    'truepositives|t=s'=> \$truePositives,
    'all|a=s'          => \$allSeqs,
    'name|n=s'         => \$nameHmm,
    'debug|d=i'        => \$debug
);

#Check input fRom user
if ($help){
    &usage();
    exit(1);
}
if ($license){
    &license();
    exit(1);
}
if (!$nameHmm){
    print STDERR "FATAL:  You must provide a string that will be used as the NAME of the new HMM.\n";
    &usage();
    exit(1);
}
if (!-s $truePositives){
    print STDERR "FATAL:  You must provide a multifasta file with the sequences of true positives.\n";
    &usage();
    exit(1);
}
if(!-s $allSeqs){
    print STDERR "FATAL:  You must provide a multifasta file with sequences of true positives and true negatives.\n";
    &usage();
    exit(1);
}

my ($allSeqsFilename,$allSeqsPath)=fileparse($allSeqs);

my $db  = Bio::DB::Fasta->new($truePositives);
my @ids = $db->get_all_primary_ids;

my $n = int(scalar(@ids)* 0.5);
my @ids_shuffled=shuffle(@ids);

my $outfileTraining=$truePositives.".trainingseqs.fasta";
my $outfileTesting =$truePositives.".testingseqs.fasta";
my $trainingSeqs_out = Bio::SeqIO->new(
                              -file   => ">$outfileTraining",
                              -format => 'fasta',
                                 );
#my $testingSeqs_out = Bio::SeqIO->new(
#                              -file   => ">$outfileTesting",
#                              -format => 'fasta',
#                                 );

#If we have at least 4 sequences in the original file,
# make a data split (original goal of this script) and 
# create an HMM from the training set (half of the input seqs)
#
my $trainingSeqs_aln='';

if(@ids_shuffled>=2){
 my $n=@ids_shuffled==3 ? 2 : $n;
 for (my $i=0; $i<@ids_shuffled;$i++){
  my $seqObj  = $db->get_Seq_by_id($ids_shuffled[$i]);
  if($i < $n){ 
   $id_class{$ids_shuffled[$i]}='TruePositiveTraining';
   $trainingSeqs_out->write_seq($seqObj);
  }
  else{
   $id_class{$ids_shuffled[$i]}='TruePositiveTest';
 #  $testingSeqs_out->write_seq($seqObj);
  }
 }
 if(@ids_shuffled>=3){
  ## Only aling sequences when there are at least 2 TruePositiveTraining sequences
  ###Align training sequence
  print STDOUT "Aligning sequences with MAFFT ".localtime()."\n";
  
  $trainingSeqs_aln   =$truePositives.".trainingseqs_".$$.".aln";
  my $trainingSeqs_alnLog=$truePositives.".trainingseqs_".$$.".aln.log";
  my $mafftCmd="mafft --auto $outfileTraining > $trainingSeqs_aln 2>$trainingSeqs_alnLog";
  system("$mafftCmd");
 }
 else{
  $trainingSeqs_aln =$outfileTraining;
 }
}
else{
 $id_class{$ids_shuffled[0]}='TruePositiveTraining';
 $trainingSeqs_aln = $truePositives;
} 
###Count number of sequences in training set
open trainingSeqs, $outfileTraining;
my @headersTraining=grep /^>/, <trainingSeqs>;
close trainingSeqs;

###Build HMM
print STDOUT "Creating HMM based on $trainingSeqs_aln".localtime()."\n";
my $hmmFile=$nameHmm."_".$$.".hmm";
my $hmmFileLog=$nameHmm."_".$$.".hmm.log";
my $hmmbuildCmd="hmmbuild -o $hmmFileLog -n $nameHmm $hmmFile $trainingSeqs_aln";
print $hmmbuildCmd."\n" if $debug >4;
system("$hmmbuildCmd");

###Run hmmsearch with a low threshold, keep the domtblout output
print STDOUT "Running hmmsearch on full set of proteins ".localtime()."\n";
my $hmmsearchOut=$allSeqsFilename."_vs_".$nameHmm."_$$.hmmsearch.domtblout";
my $hmmsearchCmd="hmmsearch -o /dev/null --domtblout $hmmsearchOut -E 100000 --domE 100000  --incE 100000 --incdomE 100000 -Z 100000 --domZ 100000 $hmmFile $allSeqs";
print $hmmsearchCmd."\n" if $debug >4;
system("$hmmsearchCmd");

###Parsing hmmsearch domtblout output
my %resHmmsearch;
print STDOUT "Parsing hmmsearch output ".localtime()."\n";
my $hmmsearchOutShort=$allSeqsFilename."_vs_".$nameHmm."_$$.hmmsearch.tbl";
open OUT, ">$hmmsearchOutShort";
open DOMTBL, $hmmsearchOut;
while(<DOMTBL>){
 next if /^#/;
 my $class='TrueNegative';
 #Pick the score and Evalue for the best domain hit.
 my ($target,undef,undef,undef,undef,undef,$evaluef,$scoref,undef,$dn,undef,undef,$evalue,$score,undef)=split(/ +/);
 if($id_class{$target}){$class=$id_class{$target}}
 if($dn == 1 ){#This get the data only for the best domain, first in results list
  $resHmmsearch{$target}="$target\t$score\t$evalue\t$class";
 # print OUT "$target\t$score\t$evalue\t$class\n";
 }
}
foreach my $tp(@ids){
 if($resHmmsearch{$tp}){
  print OUT $resHmmsearch{$tp}."\n";
 }
 else{
  print OUT "$tp\t-9999999\t9999999\t".$id_class{$tp}."\n";
 }
}

foreach my $pre_tn(keys %resHmmsearch){
 my @s=split(/\t/,$resHmmsearch{$pre_tn});
 if($s[-1] eq 'TrueNegative'){
  print OUT $resHmmsearch{$pre_tn}."\n";
 }
}

close OUT;
close DOMTBL;
unlink $hmmsearchOut unless $debug > 10;

###Generating R graph
print STDOUT "Generating R graph ".localtime()."\n";
my $rTitle=$nameHmm."_".$$;

my ($filename,$filepath) = fileparse($hmmsearchOutShort);
my $rscript=generateRscript($hmmsearchOutShort,$rTitle,$filepath,scalar(@ids),scalar(@headersTraining));
my $RscriptOut=$rTitle.".R";
open OUT, ">$RscriptOut";
print OUT $rscript."\n";
close OUT;

system("Rscript $RscriptOut");

unlink $rscript;


sub generateRscript{
 my $tbl       =shift;
 my $title     =shift;
 my $path      =shift;
 my $nSeqsTotal=shift;
 my $nSeqsTrain=shift;
 my $selectGA=$title.".selectGA.tbl";
 my $pdfFileScores =$title.".Scores.R.pdf";
 my $pdfFileEvalues =$title.".Evalues.R.pdf";
 my $rscript = <<"RSCRIPT";
library(ggplot2)
setwd("$path")
hmmsearch<-read.table("$tbl")
colnames(hmmsearch)<-c('SeqID','Score','Evalue','Class')
checkEval<-mean(hmmsearch[ which(hmmsearch\$Class == 'TruePositiveTraining' | hmmsearch\$Class == 'TruePositiveTest'),'Evalue'])
title="$title"
if(checkEval == 0){
 title = paste("$title", "\\nEvalue has pseudocounts (+1e-323), actual value is 0 for TP", sep=" ")
 hmmsearch\$Evalue<-hmmsearch\$Evalue+1e-323
}
if(length(hmmsearch[ which(hmmsearch\$Class == 'TrueNegative'),'Score']) == 0 ){
 #Dealing with the case where there are no True Negatives, the maxTN will be set as the mininum score achieved by a true positive
 maxTN<-min(hmmsearch[,'Score'])
}else{
 maxTN<-max(hmmsearch[ which(hmmsearch\$Class == 'TrueNegative'),'Score'])
}
minTP<-min(hmmsearch[ which(hmmsearch\$Class == 'TruePositiveTraining' | hmmsearch\$Class == 'TruePositiveTest'),'Score'])
putativeGA=minTP-((minTP-maxTN)/2)
if(maxTN == -9999999){
 statusGA='NotOK'
}else if(minTP == -9999999){
 statusGA='NotOK'
}else if(minTP >= maxTN){
 statusGA='OK'
}else{
 statusGA='NotOK'
}
res<-cbind("$title",maxTN,minTP,putativeGA,statusGA,"$nSeqsTotal","$nSeqsTrain")
colnames(res)<-c('Model','maxTN','minTP','pGA','statusGA','NumberSeqsTotal','NumberSeqsTraining') 
write.table(res,file="$selectGA",row.names=F,sep="\t",quote=F)
pdf(file="$pdfFileScores",width=12)
ggplot(hmmsearch,aes(x=SeqID,y=Score,colour=Class))+
  geom_point() +
  xlab("Protein ID") +
  ylab("Score") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_hline(yintercept=maxTN,colour='blue',linetype=2) +
  geom_hline(yintercept=minTP,colour='blue',linetype=2) +
  geom_hline(yintercept=putativeGA,colour='red') +
  ggtitle(title)
dev.off()
pdf(file="$pdfFileEvalues",width=12)
ggplot(hmmsearch,aes(x=SeqID,y=Evalue,colour=Class))+
  geom_point() +
  xlab("Protein ID") +
  ylab("Evalue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_log10() +
  ggtitle(title)
dev.off()

RSCRIPT
}
sub usage{
    print STDERR "$0 version $version, Copyright (C) 2016-2017 Diego Mauricio Riaño Pachón\n";
    print STDERR "$0 comes with ABSOLUTELY NO WARRANTY; for details type `$0 -l'.\n";
    print STDERR "This is free software, and you are welcome to redistribute it under certain conditions;\n";
    print STDERR "type `$0 -l' for details.\n";
    print STDERR <<EOF;
NAME
    $0   Creates a new HMM based on a list of true positives (only uses half of the sequences) and tested it against a set of proteins that includes both true postives and true negatives.

USAGE
    $0 --truepositives positives.fasta --all all.fasta --name Test

OPTIONS
    --truepositives -t    Multifasta file with true positives.         REQUIRED
    --all           -a    Multifasta file with true positives and 
                          true negatives.                              REQUIRED
    --name          -n    Name for the new HMM.                        REQUIRED
    --help          -h    This help.
    --license       -l    License.

EOF
}

sub license{
    print STDERR <<EOF;

Copyright (C) 2016-2017 Diego Mauricio Riaño Pachón
https://diriano.github.io/
e-mail: diriano\@gmail.com

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
EOF
exit;
}

