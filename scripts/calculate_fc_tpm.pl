#!/usr/bin/env perl 

## simple script that takes a matrix of TPM values for 17 in vitro conditions
## and calculates fold change of any condition to other 16
## outputting genes if FC >= 4 (log2FC >=2 )

use strict; 
use warnings; 
use Data::Dumper; 

my $tsv = shift @ARGV; 
my $tag = $tsv; 
$tag =~ s/\.tsv//g; 

open TSV,"<",$tsv or die "$!"; 

my $header = <TSV>; 
chomp $header; 
my @h = split /\t/,$header; 
shift @h; 
my $H = {};

while (<TSV>) { 
  chomp; 
  my @tpms = split /\t/; 
  my $id = shift @tpms;
  for (my $i = 0; $i < scalar @tpms; $i++) { 
    my $ave = ave_but_one($i,@tpms); 
    my $fc = ($ave == 0) ? 1 : $tpms[$i]/$ave;
    $H->{$i}->{$id} = $fc if ($fc >= 4); 
  } 
}  

print STDERR Dumper $H; 
 
for (my $i = 0; $i < scalar @h; $i++) { 
  my $name = $tag.".".$h[$i].".list"; 
  open NAME,">",$name or die "$!"; 
  foreach my $id (keys %{$H->{$i}}) { 
    print NAME "$id\n"; 
  } 
  close NAME;
} 
 
sub ave_but_one { 
  my ($index, @tpms) = @_;
  my $sum = 0; 
  for (my $i = 0; $i < scalar @tpms; $i++) { 
    $sum += $tpms[$i] if ($i != $index); 
  } 
  my $ave = $sum/(scalar @tpms - 1);
  return $ave; 
} 

