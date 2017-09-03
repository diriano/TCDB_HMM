#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;

#Global vars
my $version="0.1";
my $license='';
my $help='';
my $GAfile='';
my $debug=0;

#Get input from user
GetOptions(
    'license|l'        => \$license,
    'help|h|?'         => \$help,
    'ga_file|g=s'      => \$GAfile,
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
if(!-s $GAfile){
    print STDERR "FATAL:  You must provide a file with potential GAs.\n";
    &usage();
    exit(1);
}

open(GAFILE, $GAfile);

while(<GAFILE>) {
 chomp;
 next if (/^Model/);
 my ($Model,$maxTN,$minTP,$pGA,$statusGA)=split(/\t/);
 my ($Family,undef)=split(/_/,$Model);
 my $hmmFile = $Model.".hmm";
 #Checking whether GA is well defined, and if so update HMM;
 if($statusGA eq 'OK'){
  addGA2HMM($hmmFile,$pGA);
 }
 elsif($statusGA eq 'NotOK'){
  warn "GA is not well defined in family $Family, $Model, file: $hmmFile\n"
 }
 else {
  warn "Something wrong with $Model\t$Family\t$maxTN\t$minTP\t$pGA\n";
 }
}
close(GAFILE);

sub addGA2HMM{
 my $oldFile=shift;
 my $pGA    =shift;
 open(OLDHMM,"<",$oldFile);
 my $newFile = $oldFile;
 $newFile =~ s/\.hmm/\.mod\.hmm/g;
 open(NEWHMM,">",$newFile);
 while(<OLDHMM>){
  chomp;
  if(/^CKSUM/){
   print NEWHMM "$_\nGA    $pGA $pGA\n";
  } else {
   print NEWHMM "$_\n";
  }
 }
 close(NEWHMM);
 close(OLDHMM);
}

sub usage{
    print STDERR "$0 version $version, Copyright (C) 2017 Diego Mauricio Riaño Pachón, Renato Augusto Correa dos Santos\n";
    print STDERR "$0 comes with ABSOLUTELY NO WARRANTY; for details type `$0 -l'.\n";
    print STDERR "This is free software, and you are welcome to redistribute it under certain conditions;\n";
    print STDERR "type `$0 -l' for details.\n";
    print STDERR <<EOF;
NAME
    $0   Modifies an HMM to add the GA line. It will generate a list of the HMMs where the GA is not defined.

USAGE
    $0 --ga_file all.selectedGA.tbl

OPTIONS
    --ga_file       -g    File with selected GAs.         REQUIRED
    --help          -h    This help.
    --license       -l    License.

EOF
}

sub license{
    print STDERR <<EOF;

Copyright (C) 2017 Diego Mauricio Riaño Pachón, Renato Augusto Correa dos Santos
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
