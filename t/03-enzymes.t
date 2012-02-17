#! /usr/bin/perl -T

use Test::More tests => 3;

BEGIN
{
  use_ok( 'Bio::GeneDesign::RestrictionEnzymes' ) || print "Bail out!\n";
  diag( "Testing Bio::GeneDesign::RestrictionEnzymes $Bio::GeneDesign::RestrictionEnzymes::VERSION, Perl $], $^X" );
}

use Bio::GeneDesign::Basic qw(configure_GeneDesign);
use Bio::GeneDesign::RestrictionEnzymes qw(:all);
my $tGD = configure_GeneDesign("t/");
my $orf = "ATGGACAGATCTTGGAAGCAGAAGCTGAACCGCGACACCGTGAAGCTGACCGAGGTGATGACCTGGA";
   $orf .= "GAAGACCCGCCGCTAAATGGTTTTATACTTTAATTAATGCTAATTATTTGCCACCATGCCCACCCG";
   $orf .= "ACCACCAAGATCACCGGCAGCAACAACTACCTGAGCCTGATCAGCCTGAACATCAACGGCCTGAAC";
   $orf .= "AGCCCCATCAAGCGGCACCGCCTGACCGACTGGCTGCACAAGCAGGACCCCACCTTCTGTTGCCTC";
   $orf .= "CAGGAGACCCACCTGCGCGAGAAGGACCGGCACTACCTGCGGGTGAAGGGCTGGAAGACCATCTTT";
   $orf .= "CAGGCCAACGGCCTGAAGAAGCAGGCTGGCGTGGCCATCCTGATCAGCGACAAGATCGACTTCCAG";
   $orf .= "CCCAAGGTGATCAAGAAGGACAAGGAGGGCCACTTCATCCTGATCAAGGGCAAGATCCTGCAGGAG";
   $orf .= "GAGCTGAGCATTCTGAACATCTACGCCCCCAACGCCCGCGCCGCCACCTTCATCAAGGACACCCTC";
   $orf .= "GTGAAGCTGAAGGCCCACATCGCTCCCCACACCATCATCGTCGGCGACCTGAACACCCCCCTGAGC";
   $orf .= "AGTGA";

my $rpositions = {
  TspRI => {'592' => 'AGCAGTG'}, PspOMI => {}, AlwI => {'450' => 'GATCC'},
  BsrI => {'227' => 'ACTGG'}, TaqI => {'386' => 'TCGA'}, BsrBI => {},
  AcuI => {'535' => 'CTGAAG', '343' => 'CTGAAG'}, AccI => {},
  AarI => {'274' => 'CACCTGC'}, MseI => {'100' => 'TTAA', '96' => 'TTAA'},
  AciI => {'503' => 'CCGC', '210' => 'GCGG', '303' => 'GCGG', '498' => 'CCGC', '73' => 'CCGC', '216' => 'CCGC', '76' => 'CCGC', '29' => 'CCGC'},
  AhdI => {}, BsmI => {'470' => 'GCATTC'}, MscI => {'362' => 'TGGCCA'},
  BbvCI => {}, RsaI => {}, AatII => {}, BtrI => {}, BmgBI => {}, PspXI => {},
  PspGI => {'61' => 'CCTGG', '264' => 'CCAGG'}, MlyI => {}, BmtI => {},
  BtsCI => {'366' => 'CATCC', '432' => 'CATCC'}, BsrDI => {},
  NlaIII => {'121' => 'CATG'}, BseYI => {}, BtsI => {'593' => 'GCAGTG'},
  Bpu10I => {'588' => 'CCTGAGC', '162' => 'CCTGAGC'}
};

#TESTING define_sites()
my $tRES = define_sites($tGD->{enzyme_list});

#TESTING define_site_status()
my $tSITE_STATUS = define_site_status($orf, $tRES);
my $rSITE_STATUS = {
  PspOMI => 0, TspRI => 1, BsrI => 1, TaqI => 1, AcuI => 2, AccI => 0,
  AarI => 1, AciI => 8, AhdI => 0, BsmI => 1, BbvCI => 0, MscI => 1, RsaI => 0,
  AatII => 0, BtrI => 0, BmgBI => 0, PspXI => 0, MlyI => 0, PspGI => 2,
  BsrBI => 0, BmtI => 0, AlwI => 1, BtsCI => 2, BsrDI => 0, BseYI => 0,
  NlaIII => 1, Bpu10I => 2, MseI => 2, BtsI => 1
};
is_deeply($tSITE_STATUS, $rSITE_STATUS, "define site status");

#TESTING siteseeker()
subtest "site seeker" => sub
{
  plan tests => scalar(keys %$tSITE_STATUS);
  foreach my $enz (keys %$tSITE_STATUS)
  {
    $tpositions = siteseeker($orf, $tRES->{$enz});
    is_deeply($tpositions, $rpositions->{$enz}, "find site $enz");
  }
};

