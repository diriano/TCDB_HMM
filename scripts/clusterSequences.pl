#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Path qw(make_path remove_tree);
use Bio::DB::Fasta;
use Bio::SeqIO;

#Global vars
my $version="0.1";
my $license='';
my $help='';
my $selectGA='';
my $ident='';
my %seq2cluster;
my $debug=0;

#Get input from user
GetOptions(
    'license|l'        => \$license,
    'help|h|?'         => \$help,
    'selectGA|g=s'     => \$selectGA,
    'ident|i=f'        => \$ident,
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
if (!-s $selectGA){
    print STDERR "FATAL:  You must provide the file all.selectGA.txt.\n";
    &usage();
    exit(1);
}
if(!$ident){
    print STDERR "FATAL:  You must provide a float that will be used for clusterin in CD-HIT. Must be between 0 and 1.\n";
    &usage();
    exit(1);
}
if($ident > 1){
    print STDERR "FATAL:  You must provide a float that will be used for clusterin in CD-HIT. Must be between 0 and max 1.\n";
    &usage();
    exit(1);
}

#Check if output dir exist. Delete if yes, otherwise create it
my $outDir='seeds_cdhit_'.$ident;
if(-d $outDir){
 #Remove previous results and then re-create the directory
 remove_tree($outDir);
 make_path("$outDir/tmp");
}
else{
 make_path("$outDir/tmp");
# mkdir $outDir;
}

#Select families, and their sequences for clustering, these with undefined GA in the previous iteration
open selectGA, $selectGA;
while(<selectGA>){
 chomp;
 next if /^Model/;
 my ($model,$maxTN,$minTP,$pGA,$statusGA,$NumberSeqsTotal,$NumberSeqsTraining)=split(/\t/);
 my ($family,undef)=split(/_/,$model);
 next if $statusGA eq 'OK';
 my $familySeqs='tcdb_group_'.$family.'.fasta';
 my $familyClusteredFasta=$outDir.'/tmp/tcdb_group_'.$family.'_cdhit'.$ident.'.fasta';
 my $familyClusterClstr=$outDir.'/tmp/tcdb_group_'.$family.'_cdhit'.$ident.'.fasta.clstr';
 my $cdhitLog=$outDir.'/tmp/tcdb_group_'.$family.'_cdhit'.$ident.'.cdhit.log';
 my $cdhitCmd="cd-hit -M 0 -c $ident -d 0 -i $familySeqs -o $familyClusteredFasta > $cdhitLog";
 system("$cdhitCmd");
 print "$familyClusterClstr\n" if $debug >15;
 open CLUSTR,$familyClusterClstr;
 my $clusterId=0;
 while(<CLUSTR>){
  chomp;
  if(/^>/){
   $clusterId++
  }
  else{
   my @fields=split(/,/);
   my @fields2=split(/ /,$fields[1]);
   my $seqId=$fields2[1];
   $seqId=~s/^>//;
   $seqId=~s/\.+$//;
   $seq2cluster{$family}{$clusterId}{$seqId}=1;
  }
 }
 print $_."\t$familySeqs\n" if $debug > 10;
}
close selectGA;

#generate new seed files
foreach my $family(keys %seq2cluster){
 my $familySeqs='tcdb_group_'.$family.'.fasta';
 my $dbinx=indexdb($familySeqs);
 foreach my $clusterId(keys %{$seq2cluster{$family}}){
  my $cluster=sprintf("%04d", $clusterId);
  my $outfile=$outDir.'/tcdb_group_'.$family.'_cdhit'.$ident.'_cluster'.$cluster.'.fasta'; 
  my $outfile_obj = Bio::SeqIO->new('-file'   => ">$outfile",
                                    '-format' => 'Fasta');
  foreach my $seqId(keys %{$seq2cluster{$family}{$clusterId}}){
   my $seq_obj=$dbinx->get_Seq_by_id($seqId);
   $outfile_obj->write_seq($seq_obj);
  }
 }
}

sub indexdb{
 my $file=shift;
# my $reinx=1;
 my $reinx=0;
 my $dbinx = Bio::DB::Fasta->new("$file", -reindex=>$reinx);
}

sub usage{
    print STDERR "$0 version $version, Copyright (C) 2016-2017 Diego Mauricio Riaño Pachón\n";
    print STDERR "$0 comes with ABSOLUTELY NO WARRANTY; for details type `$0 -l'.\n";
    print STDERR "This is free software, and you are welcome to redistribute it under certain conditions;\n";
    print STDERR "type `$0 -l' for details.\n";
    print STDERR <<EOF;
NAME
    $0   Take the file all.selectGA.txt and select the families with undefined GA, then cluster them using cd-hit, so that each cd-hit group would become a new seed for createHMMFromDataSplit.pl

USAGE
    $0 --selectGA all.selectGA.txt --ident 0.90

OPTIONS
    --selectGA      -g    all.selectGA.txt file.                      REQUIRED
    --ident         -i    Identity for clustering. Between 0 and 1.   REQUIRED
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

