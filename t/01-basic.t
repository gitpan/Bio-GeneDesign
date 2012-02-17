#! /usr/bin/perl -T

use Test::More tests => 7;

BEGIN
{
  use_ok( 'Bio::GeneDesign::Basic' ) || print "Bail out!\n";
  diag( "Testing Bio::GeneDesign::Basic $Bio::GeneDesign::Basic::VERSION, Perl $], $^X" );
}

use Bio::GeneDesign::Basic qw(:all);

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
my $busted = "ATGGAYMGNWSNTGGAARCARAARYTNAAYMGNGAYACNGTNAARYTNACNGARGTNATGACNT";
   $busted .= "GGMGNMGNCCNGCNGCNAARTGGTTYTAYACNYTNATHAAYGCNAAYTAYYTNCCNCCNTGYC";
   $busted .= "CNCCNGAYCAYCARGAYCAYMGNCARCARCARYTNCCNGARCCNGAYCARCCNGARCAYCARM";
   $busted .= "GNCCNGARCARCCNCAYCARGCNGCNCCNCCNGAYMGNYTNGCNGCNCARGCNGGNCCNCAYY";
   $busted .= "TNYTNYTNCCNCCNGGNGAYCCNCCNGCNMGNGARGGNCCNGCNYTNCCNGCNGGNGARGGNY";
   $busted .= "TNGARGAYCAYYTNWSNGGNCARMGNCCNGARGARGCNGGNTGGMGNGGNCAYCCNGAYCARM";
   $busted .= "GNCARGAYMGNYTNCCNGCNCARGGNGAYCARGARGGNCARGGNGGNCCNYTNCAYCCNGAYC";
   $busted .= "ARGGNCARGAYCCNGCNGGNGGNGCNGARCAYWSNGARCAYYTNMGNCCNCARMGNCCNMGNM";
   $busted .= "GNCAYYTNCAYCARGGNCAYCCNMGNGARGCNGARGGNCCNCAYMGNWSNCCNCAYCAYCAYM";
   $busted .= "GNMGNMGNCCNGARCAYCCNCCNGARCARTRR";
my $shortorf = "ATGGACAGATCTTGGAAGCAGAAGCTGAACCGC";
my $shortbusted = "ABCDGHKMNRSTVWXY";
my $shortpep = "MDRSWKQKLNRDTVKLTEVMTWR*";
my $rlshortorf =  "ATGGACAGAT";
my $shortmessy = "\"ATGGAC>fastaTCttg ecZAGGWXYNNRB";
my $shortmessypep = "MDRSWKxbjQKLNRDTVBZKLTEVMTW!& R*";
my $fasta = ">holy crap, this sequence is awesome and starts with ATGG\nATGGAC";
   $fasta .= "AGATCTTGGAAGCAGAAGCTGAACCGC\n";

#TESTING configure_GeneDesign()
my $tGD = configure_GeneDesign("t/");
my $rGD = {this_server => "127.0.0.1", codon_dir => "t/codon_tables",
  enzyme_list => "t/enzymes/test", VERSION => "3.05",
  ORGANISMS => {Saccharomyces_cerevisiae => 't/codon_tables/Saccharomyces_cerevisiae.rscu'}, 
  CODONTABLES => {Standard => "t/codon_tables/Standard.ct"}};
is_deeply($tGD, $rGD, "configuring GeneDesign");

#TESTING count()
subtest "count" => sub
{
  plan tests => 2;

  my $torfcount = count($orf);
  my $rorfcount = {A => 159, T => 95, C => 197, G => 149, R => 0, Y => 0,
    W => 0, S => 0, M => 0, K => 0, B => 0, D => 0, H => 0, V => 0, U => 0,
    N => 0, length => 600, "?" => 0, d => 600, n => 0, GCp => 58, ATp => 42};
  is_deeply($torfcount, $rorfcount, "count non ambiguous bases");

  my $tbustcount = count($busted);
  my $rbustcount = {A => 91, T => 31, C => 125, G => 113, R => 43, Y => 54,
    W => 4, S => 4, M => 21, K => 0, B => 0, D => 0, H => 1, V => 0, U => 0,
    N => 113, length => 600, "?" => 0, d => 360, n => 240, GCp => 60,
    ATp => 40};
  is_deeply($tbustcount, $rbustcount, "count ambiguous bases");
};

#TESTING melt()
subtest "melt" => sub
{
  plan tests => 4;

  my $tmelt1 = melt($rlshortorf, 1);
  my $rmelt1 = 28;
  is ($tmelt1, $rmelt1, "simple Tm");

  my $tmelt2 = melt($shortorf, 2);
  my $rmelt2 = 60.5695687386446;
  is ($tmelt2, $rmelt2, "baldwin Tm");

  my $tmelt3 = melt($shortorf, 3);
  my $rmelt3 = 62.8422960113718;
  is ($tmelt3, $rmelt3, "primer3 Tm");

  my $tmelt4 = melt($shortorf, 4);
  my $rmelt4 = 67.3152820409049;
  is ($tmelt4, $rmelt4, "nearest neighbor thermodynamics Tm");
};

#TESTING complement()
subtest "complementing" => sub
{
  plan tests => 3;

  my $tfro = complement($orf, 1);
  my $rfro = "TCACTGCTCAGGGGGGTGTTCAGGTCGCCGACGATGATGGTGTGGGGAGCGATGTGGGCCTTCA";
     $rfro .= "GCTTCACGAGGGTGTCCTTGATGAAGGTGGCGGCGCGGGCGTTGGGGGCGTAGATGTTCAGAA";
     $rfro .= "TGCTCAGCTCCTCCTGCAGGATCTTGCCCTTGATCAGGATGAAGTGGCCCTCCTTGTCCTTCT";
     $rfro .= "TGATCACCTTGGGCTGGAAGTCGATCTTGTCGCTGATCAGGATGGCCACGCCAGCCTGCTTCT";
     $rfro .= "TCAGGCCGTTGGCCTGAAAGATGGTCTTCCAGCCCTTCACCCGCAGGTAGTGCCGGTCCTTCT";
     $rfro .= "CGCGCAGGTGGGTCTCCTGGAGGCAACAGAAGGTGGGGTCCTGCTTGTGCAGCCAGTCGGTCA";
     $rfro .= "GGCGGTGCCGCTTGATGGGGCTGTTCAGGCCGTTGATGTTCAGGCTGATCAGGCTCAGGTAGT";
     $rfro .= "TGTTGCTGCCGGTGATCTTGGTGGTCGGGTGGGCATGGTGGCAAATAATTAGCATTAATTAAA";
     $rfro .= "GTATAAAACCATTTAGCGGCGGGTCTTCTCCAGGTCATCACCTCGGTCAGCTTCACGGTGTCG";
     $rfro .= "CGGTTCAGCTTCTGCTTCCAAGATCTGTCCAT";
  is ($tfro, $rfro, "reverse complement");

  my $tfro2 = complement($orf);
  my $rfro2 = "TACCTGTCTAGAACCTTCGTCTTCGACTTGGCGCTGTGGCACTTCGACTGGCTCCACTACTGG";
     $rfro2 .= "ACCTCTTCTGGGCGGCGATTTACCAAAATATGAAATTAATTACGATTAATAAACGGTGGTAC";
     $rfro2 .= "GGGTGGGCTGGTGGTTCTAGTGGCCGTCGTTGTTGATGGACTCGGACTAGTCGGACTTGTAG";
     $rfro2 .= "TTGCCGGACTTGTCGGGGTAGTTCGCCGTGGCGGACTGGCTGACCGACGTGTTCGTCCTGGG";
     $rfro2 .= "GTGGAAGACAACGGAGGTCCTCTGGGTGGACGCGCTCTTCCTGGCCGTGATGGACGCCCACT";
     $rfro2 .= "TCCCGACCTTCTGGTAGAAAGTCCGGTTGCCGGACTTCTTCGTCCGACCGCACCGGTAGGAC";
     $rfro2 .= "TAGTCGCTGTTCTAGCTGAAGGTCGGGTTCCACTAGTTCTTCCTGTTCCTCCCGGTGAAGTA";
     $rfro2 .= "GGACTAGTTCCCGTTCTAGGACGTCCTCCTCGACTCGTAAGACTTGTAGATGCGGGGGTTGC";
     $rfro2 .= "GGGCGCGGCGGTGGAAGTAGTTCCTGTGGGAGCACTTCGACTTCCGGGTGTAGCGAGGGGTG";
     $rfro2 .= "TGGTAGTAGCAGCCGCTGGACTTGTGGGGGGACTCGTCACT";
  is ($tfro2, $rfro2, "complement");

  my $tdets = complement($busted, 1);
  my $rdets = "YYAYTGYTCNGGNGGRTGYTCNGGNCKNCKNCKRTGRTGRTGNGGNSWNCKRTGNGGNCCYTC";
     $rdets .= "NGCYTCNCKNGGRTGNCCYTGRTGNARRTGNCKNCKNGGNCKYTGNGGNCKNARRTGYTCNS";
     $rdets .= "WRTGYTCNGCNCCNCCNGCNGGRTCYTGNCCYTGRTCNGGRTGNARNGGNCCNCCYTGNCCY";
     $rdets .= "TCYTGRTCNCCYTGNGCNGGNARNCKRTCYTGNCKYTGRTCNGGRTGNCCNCKCCANCCNGC";
     $rdets .= "YTCYTCNGGNCKYTGNCCNSWNARRTGRTCYTCNARNCCYTCNCCNGCNGGNARNGCNGGNC";
     $rdets .= "CYTCNCKNGCNGGNGGRTCNCCNGGNGGNARNARNARRTGNGGNCCNGCYTGNGCNGCNARN";
     $rdets .= "CKRTCNGGNGGNGCNGCYTGRTGNGGYTGYTCNGGNCKYTGRTGYTCNGGYTGRTCNGGYTC";
     $rdets .= "NGGNARYTGYTGYTGNCKRTGRTCYTGRTGRTCNGGNGGRCANGGNGGNARRTARTTNGCRT";
     $rdets .= "TDATNARNGTRTARAACCAYTTNGCNGCNGGNCKNCKCCANGTCATNACYTCNGTNARYTTN";
     $rdets .= "ACNGTRTCNCKRTTNARYTTYTGYTTCCANSWNCKRTCCAT";
  is ($tdets, $rdets, "busted reverse complement");
};

#TESTING regres()
subtest "regular expression generation" => sub
{
  plan tests => 2;

  my $tbreg = regres($shortbusted, 1);
  my $rbreg = "A[BCGKSTY]C[ADGKRTW]G[ACHMTWY][GKT][ACM][ABCDGHKMNRSTVWY][AGR]";
     $rbreg .= "[CGS]T[ACGMRSV][ATW][X][CTY]";
  is ($tbreg, $rbreg, "nucleotide regexg");

  my $tpreg = regres($shortpep, 2);
  my $bpreg = "MDRSWKQKLNRDTVKLTEVMTWR[\*]";
  is ($tpreg, $bpreg, "protein regexg");
};

#TESTING cleanup()
subtest "cleanup" => sub
{
  plan tests => 4;

	my $tcleanst = cleanup($shortmessy, 0);
	my $rcleanst = "ATGGACATATCTTGCAGG";
	is ($tcleanst, $rcleanst, "strict nucleotide clean");
	
	my $tclean = cleanup($shortmessy, 1);
	my $rclean = "ATGGACASTATCTTGCAGGWYNNRB";
	is ($tclean, $rclean, "nucleotide clean");
	
	my $tcleanpep = cleanup($shortmessypep, 2);
  is ($tcleanpep, $shortpep, "peptide clean");

	$tcleanfasta = cleanup($fasta, 1);
  is ($tcleanfasta, $shortorf, "FASTA clean");
};

