#!/usr/bin/env perl
############################################################################
# SCRIPT NAME:	bed2graph
# DESCRIPTION:	parse bed file, generate cypher to create a graph (Neo4j)
#
# DATE WRITTEN: 2024-05-03
# WRITTEN BY:   Martin Maiers
#
############################################################################
use strict;    
use warnings;  

my $file = "cat ../KIR2/assemble_results.bed|";
open FILE, $file or die "$!: $file";

my %I; #ids
my %N; #notes
my %E; #edges

# run for selected ids
$I{"cA01~tA01_GU182338.1/1155286"}++;
$I{"cB01~tB01_GU182339.1/1239295"}++;
my $cypher = "CREATE ";
my @clauses;

my $prev="";
my $curid="";
while(<FILE>) {
  chomp;
  my ($id,$start, $stop, $blockname) = split /	/;	# hard tab
  next unless defined $id;
  next if $id=~/^#/;

  # comment this out to run for all ids

  #next unless defined $I{$id};

  my ($blocknum) = (split /:/, $blockname)[0];
  $N{"b$blocknum"}++;
  if ($id ne $curid) {
    $curid = $id;
  }  else {
    if ($prev) {
      $E{"b$blocknum"}{$prev}++;
    }
  }
  $prev = "b$blocknum";
}
foreach my $n (keys %N) {
  push @clauses,"($n:Block {num: \"$n\"})";
}

foreach my $from (keys %E) {
  foreach my $to (keys %{$E{$from}}) {
    push @clauses,"($from)-[:LINK]-> ($to)";
  }
} 
print $cypher, join (',', @clauses), "\n\n";;

exit 0;

