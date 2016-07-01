#!/usr/bin/perl
 use strict;
 use warnings;
 use Getopt::Long;
# this is for the new data sent on June 29, 2016



##############  read options  #########################

my ($ref, $file, $out);

# default values
$ref="";
$file="";
$out="";

GetOptions 
(
  'r=s'=>\$ref,
  'f=s'=>\$file,
  'o=s'=>\$out,
  );


if ($ref eq "" || $file eq "") {die "Missing one or both file name(s)!\n 
    USAGE: perl step1_find_position.pl -r <GTF file> -f <BED file> -o <name of output file>\n
    NOTE: output file name is not required. So far the reference only support the gtf file from flybase. i.e. Drosophila malanogaster\n";} # check parameters


############# step1, get the reference (transcript.gtf) #######################

 print "loading reference\n  Warning: make sure your RAM is 24GB or larger. Else your computer may freeze\n";

my $gene; my $on; my $st_int; my %utr3=(); my %utr5=(); my %intron=(); my %cds=(); my %ncrna=();

my %utr3_trans; my %utr5_trans; my %intron_trans; my %cds_trans; my %ncrna_trans; my %trans;

open (hand1, $ref) or die $!;

while (<hand1>)        {
  chomp;
  my @a = split /\t/; 
  my @b = split /"/, $a[8];
  #print "the gene is $a[0]\t$a[6]\t$a[8]\n";

  if ($a[2] eq 'gene') {
     $gene = $b[1];
     $on = 0;
     next  }
    
 
  if ($a[2] eq 'ncRNA')  {
     for (($a[3]-1)..$a[4])  {
        $ncrna{$a[0]}{$_} = $gene; }
     $ncrna_trans{$b[5]} =  $a[4] - $a[3] + 1;
     $trans{$b[5]} = $gene; 
     next  }
  
  if ($a[2] eq 'mRNA')   {
     $intron_trans{$b[5]} = $a[4]-$a[3]+1;
     $trans{$b[5]} = $gene;
     next    }

  if ($on ==1 and $a[2] eq 'exon')  {
     for ($st_int..($a[3]-2))   {  #gtf offset is 1, the ref is bed, so position will be -2
        $intron{$a[0]}{$_} .= "$b[5],";   }
             }  
   

  if ($a[2] eq "5UTR")  {
     for (($a[3]-1)..($a[4]-1))  {
        $utr5{$a[0]}{$_} .= "$b[5],";   }
     $intron_trans{$b[5]} -= $a[4] - $a[3] + 1;
     $utr5_trans{$b[5]} += $a[4] - $a[3] + 1; }  
   

  elsif ($a[2] eq "3UTR")  {
     for (($a[3]-1)..($a[4]-1))  {
        $utr3{$a[0]}{$_} .= "$b[5],";   }
     $intron_trans{$b[5]} -= $a[4] - $a[3] + 1;
     $utr3_trans{$b[5]} += $a[4] - $a[3] + 1;}  
   

  elsif ($a[2] eq "CDS")  {
     for (($a[3]-1)..($a[4]-1))  {
        $cds{$a[0]}{$_} .= "$b[5],";   }
     $intron_trans{$b[5]} -= $a[4] - $a[3] + 1;
     $cds_trans{$b[5]} += $a[4] - $a[3] + 1;}  
    

  elsif ($a[2] eq "exon")  {
    $st_int = $a[4]; $on =1; }    # do not -1 because offset is 1

                          }
    
close hand1;

open (hand2, ">gene_stat.txt");
open (hand3, ">ncRNA_stat.txt");

print hand2 "transcript_id\tgene\tUTR5\tCDS\tUTR3\tIntron\n";

foreach my $id (keys %trans)   {

  if (exists $ncrna_trans{$id})  {
     print hand3 "$id\t$trans{$id}\t$ncrna_trans{$id}\n";  }
  elsif (exists $utr3_trans{$id} and exists $utr5_trans{$id} and exists $cds_trans{$id} and exists $intron_trans{$id}) {
     print hand2 "$id\t$trans{$id}\t$utr5_trans{$id}\t$cds_trans{$id}\t$utr3_trans{$id}\t$intron_trans{$id}\n";  }   }

close hand2; 

print "exporting results\n";


  open (hand2, $file) or die $!;
  open (hand3, ">$out\_out.txt");

  print hand3 "Chrom\tstart\tstrand\treading\tUTR3\tUTR5\tCDS\tintron\tncRNA\n";

  while (<hand2>)    {
    chomp;
    my ($chr, $st, $end, $str) =split /\t/;
    $chr =~ s/^chr//;
    my $utr3 = search_ref($chr, $st, \%utr3);
    my $utr5 = search_ref($chr, $st, \%utr5);
    my $cds = search_ref($chr, $st, \%cds);
    my $intron = search_ref($chr, $st, \%intron);
    my $ncrna = search_ref($chr, $st, \%ncrna);
  
    print hand3 "$chr\t$st\t$str\t1\t$utr3\t$utr5\t$cds\t$intron\t$ncrna\n"; }
 
close hand3;

########### step2 summarize the output  #####################  

system("perl summary.pl ".$out) == 0 or die "summary.pl failed\n";
    


sub search_ref        {
  my ($ch, $start, $ref) = @_;
  my $o = 'NA';

  
  if (exists $$ref{$ch}{$start})  {
     $o = $$ref{$ch}{$start}; }
  return $o;  
    } 
    

