#!/usr/bin/perl
 use strict;
 use warnings;

my ($sam)=@ARGV;

my $gene; my %trans; my %utr3=(); my %utr5=(); my %intron=(); my %cds=(); my %ncrna = ();

open (hand1, "$sam\_out.txt") or die $!;

while (<hand1>)     {
  next if /^Chrom/;
  chomp;
  #print "$_\n";
  my @a = split /\t/;
  
  # for 3 utr
  if ($a[4] ne "NA")   {
   my @b = split /,/, $a[4];
   foreach my $i (@b)   {
     $utr3{$i} += $a[3]; 
     $trans{$i} ++ ;     } }

  # for 5 utr
  if ($a[5] ne "NA")   {
   my @b = split /,/, $a[5];
   foreach my $i (@b)   {
     $utr5{$i} += $a[3];  
     $trans{$i} ++ ;    } }
   

   # for cds
  if ($a[6] ne "NA")   {
   my @b = split /,/, $a[6];
   foreach my $i (@b)   {
     $cds{$i} += $a[3]; 
     $trans{$i} ++ ;     } }
   
 # for intron
  if ($a[7] ne "NA")   {
   my @b = split /,/, $a[7];
   foreach my $i (@b)   {
     $intron{$i} += $a[3];  
     $trans{$i} ++ ;    } }

 # for ncRNA
  if ($a[8] ne "NA")   {
   my @b = split /,/, $a[8];
   foreach my $i (@b)   {
     $ncrna{$i} += $a[3]; 
     $trans{$i} ++ ;     } }
     }
   
close hand1;

open (hand2, ">$sam\_gene.txt");

print hand2 "Gene\tUTR5\tCDS\tUTR3\tIntron\tncRNA\n";

foreach my $ge (sort keys %trans)   {
  if (exists $utr5{$ge})  {
    print hand2 "$ge\t$utr5{$ge}";}
  else {
    print hand2 "$ge\t0"; }

  if (exists $cds{$ge})  {
    print hand2 "\t$cds{$ge}";}
  else {
    print hand2 "\t0"; }
 
  if (exists $utr3{$ge})  {
    print hand2 "\t$utr3{$ge}";}
  else {
    print hand2 "\t0"; } 

  if (exists $intron{$ge})  {
    print hand2 "\t$intron{$ge}";}
  else {
    print hand2 "\t0"; }       

  if (exists $ncrna{$ge})  {
    print hand2 "\t$ncrna{$ge}\n";}
  else {
    print hand2 "\t0\n"; }     }   

close hand2;   



