#! /usr/bin/perl -T

use Test::More tests => 22;

BEGIN
{
  use_ok( 'Bio::GeneDesign::Codons' ) || print "Bail out!\n";
  diag( "Testing Bio::GeneDesign::Codons $Bio::GeneDesign::Codons::VERSION, "
      . "Perl $], $^X" );
}

use Bio::GeneDesign::Basic qw(configure_GeneDesign);
use Bio::GeneDesign::Codons qw(:all);

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
my $shortamb = "ABGCDT";
my $shortorf = "ATGGACAGATCTTGGAAGCAGAAGCTGAACCGC";
my $notorf = "ATGGGATAGTGAAATATG";
my $nostop = "ATGGGTAACCAATTA";

# TESTING define_codons();
my $tcodons = define_codons();
my $rcodons = [qw(AAA AAT AAC AAG ATA ATT ATC ATG ACA ACT ACC ACG AGA AGT
  AGC AGG TAA TAT TAC TAG TTA TTT TTC TTG TCA TCT TCC TCG TGA TGT TGC TGG CAA
  CAT CAC CAG CTA CTT CTC CTG CCA CCT CCC CCG CGA CGT CGC CGG GAA GAT GAC GAG
  GTA GTT GTC GTG GCA GCT GCC GCG GGA GGT GGC GGG)];
is_deeply( $tcodons, $rcodons, "define codons" );

# TESTING define_codon_table()
my $tCT = define_codon_table("Standard", $tGD);
$rCT = {TTT => "F", TTC => "F", TTA => "L", TTG => "L", CTT => "L", CTC => "L",
  CTA => "L", CTG => "L", ATT => "I", ATC => "I", ATA => "I", ATG => "M",
  GTT => "V", GTC => "V", GTA => "V", GTG => "V", TCT => "S", TCC => "S",
  TCA => "S", TCG => "S", CCT => "P", CCC => "P", CCA => "P", CCG => "P",
  ACT => "T", ACC => "T", ACA => "T", ACG => "T", GCT => "A", GCC => "A",
  GCA => "A", GCG => "A", TAT => "Y", TAC => "Y", TAA => "*", TAG => "*",
  CAT => "H", CAC => "H", CAA => "Q", CAG => "Q", AAT => "N", AAC => "N",
  AAA => "K", AAG => "K", GAT => "D", GAC => "D", GAA => "E", GAG => "E",
  TGT => "C", TGC => "C", TGA => "*", TGG => "W", CGT => "R", CGC => "R",
  CGA => "R", CGG => "R", AGT => "S", AGC => "S", AGA => "R", AGG => "R",
  GGT => "G", GGC => "G", GGA => "G", GGG => "G"};
is_deeply( $tCT, $rCT, "define codon table" );

# TESTING define_reverse_codon_table()
my $tRCT = define_reverse_codon_table($tCT);
my $rRCT = {A => [qw(GCC GCA GCG GCT)], C => [qw(TGT TGC)], D => [qw(GAC GAT)],
  E => [qw(GAG GAA)], F => [qw(TTT TTC)], G => [qw(GGA GGT GGG GGC)],
  H => [qw(CAC CAT)], I => [qw(ATC ATT ATA)], K => [qw(AAG AAA)],
  L => [qw(CTG CTC TTG CTT TTA CTA)], M => [qw(ATG)],
  N => [qw(AAC AAT)], P => [qw(CCG CCA CCC CCT)], Q => [qw(CAG CAA)],
  R => [qw(CGA AGA CGC AGG CGT CGG)], S => [qw(AGT AGC TCG TCA TCC TCT)],
  T => [qw(ACA ACC ACG ACT)], V => [qw(GTG GTT GTA GTC)], W => [qw(TGG)],
  Y => [qw(TAC TAT)], '*' => [qw(TGA TAA TAG)]};
is_deeply( $tRCT, $rRCT, "define reverse codon table" );									

# TESTING parse_RSCU_values()
my $tpRSCU = parse_RSCU_values("t/codon_tables/Saccharomyces_cerevisiae.rscu");
my $rpRSCU = {TTT => "0.19", TTC => "1.81", TTA => "0.49", TTG => "5.34",
  CTT => "0.02", CTC => "0.00", CTA => "0.15", CTG => "0.02", ATT => "1.26",
  ATC => "1.74", ATA => "0.00", ATG => "1.00", GTT => "2.07", GTC => "1.91",
  GTA => "0.00", GTG => "0.02", TCT => "3.26", TCC => "2.42", TCA => "0.08",
  TCG => "0.02", CCT => "0.21", CCC => "0.02", CCA => "3.77", CCG => "0.00",
  ACT => "1.83", ACC => "2.15", ACA => "0.00", ACG => "0.01", GCT => "3.09",
  GCC => "0.89", GCA => "0.03", GCG => "0.00", TAT => "0.06", TAC => "1.94",
  TAA => "1.00", TAG => "0.00", CAT => "0.32", CAC => "1.68", CAA => "1.98",
  CAG => "0.02", AAT => "0.06", AAC => "1.94", AAA => "0.16", AAG => "1.84",
  GAT => "0.70", GAC => "1.30", GAA => "1.98", GAG => "0.02", TGT => "1.80",
  TGC => "0.20", TGA => "0.00", TGG => "1.00", CGT => "0.63", CGC => "0.00",
  CGA => "0.00", CGG => "0.00", AGT => "0.06", AGC => "0.16", AGA => "5.37",
  AGG => "0.00", GGT => "3.92", GGC => "0.06", GGA => "0.00", GGG => "0.02"};
is_deeply( $tpRSCU, $rpRSCU, "parse RSCU values" );

# TESTING define_RSCU_values()
my $tRSCU = define_RSCU_values("Saccharomyces_cerevisiae", $tGD);
my $rRSCU = {TTT => "0.19", TTC => "1.81", TTA => "0.49", TTG => "5.34",
  CTT => "0.02", CTC => "0.00", CTA => "0.15", CTG => "0.02", ATT => "1.26",
  ATC => "1.74", ATA => "0.00", ATG => "1.00", GTT => "2.07", GTC => "1.91",
  GTA => "0.00", GTG => "0.02", TCT => "3.26", TCC => "2.42", TCA => "0.08",
  TCG => "0.02", CCT => "0.21", CCC => "0.02", CCA => "3.77", CCG => "0.00",
  ACT => "1.83", ACC => "2.15", ACA => "0.00", ACG => "0.01", GCT => "3.09",
  GCC => "0.89", GCA => "0.03", GCG => "0.00", TAT => "0.06", TAC => "1.94",
  TAA => "1.00", TAG => "0.00", CAT => "0.32", CAC => "1.68", CAA => "1.98",
  CAG => "0.02", AAT => "0.06", AAC => "1.94", AAA => "0.16", AAG => "1.84",
  GAT => "0.70", GAC => "1.30", GAA => "1.98", GAG => "0.02", TGT => "1.80",
  TGC => "0.20", TGA => "0.00", TGG => "1.00", CGT => "0.63", CGC => "0.00",
  CGA => "0.00", CGG => "0.00", AGT => "0.06", AGC => "0.16", AGA => "5.37",
  AGG => "0.00", GGT => "3.92", GGC => "0.06", GGA => "0.00", GGG => "0.02"};
is_deeply( $tRSCU, $rRSCU, "define RSCU values" );

# TESTING translate()
subtest "translation" => sub
{
  plan tests => 6;

  my $tpeptide = translate($orf, 1, $tCT);
  my $rpeptide = "MDRSWKQKLNRDTVKLTEVMTWRRPAAKWFYTLINANYLPPCPPDHQDHRQQQLPEPDQP";
     $rpeptide .= "EHQRPEQPHQAAPPDRLAAQAGPHLLLPPGDPPAREGPALPAGEGLEDHLSGQRPEEAG";
     $rpeptide .= "WRGHPDQRQDRLPAQGDQEGQGGPLHPDQGQDPAGGAEHSEHLRPQRPRRHLHQGHPRE";
     $rpeptide .= "AEGPHRSPHHHRRRPEHPPEQ*";
  is ($tpeptide, $rpeptide, "translate frame 1");

  my $teptide = translate($orf, 2, $tCT);
  my $reptide = "WTDLGSRS*TATP*S*PR**PGEDPPLNGFIL*LMLIICHHAHPTTKITGSNNYLSLISLN";
     $reptide .= "INGLNSPIKRHRLTDWLHKQDPTFCCLQETHLREKDRHYLRVKGWKTIFQANGLKKQAGV";
     $reptide .= "AILISDKIDFQPKVIKKDKEGHFILIKGKILQEELSILNIYAPNARAATFIKDTLVKLKA";
     $reptide .= "HIAPHTIIVGDLNTPLSS";
  is ($teptide, $reptide, "translate frame 2");

  my $tptide = translate($orf, 3, $tCT);
  my $rptide = "GQILEAEAEPRHREADRGDDLEKTRR*MVLYFN*C*LFATMPTRPPRSPAATTT*A*SA*TS";
     $rptide .= "TA*TAPSSGTA*PTGCTSRTPPSVASRRPTCARRTGTTCG*RAGRPSFRPTA*RSRLAWPS";
     $rptide .= "*SATRSTSSPR*SRRTRRATSS*SRARSCRRS*AF*TSTPPTPAPPPSSRTPS*S*RPTSL";
     $rptide .= "PTPSSSAT*TPP*AV";
  is ($tptide, $rptide, "translate frame 3");

  my $teditpep = translate($orf, -1, $tCT);
  my $reditpep = "SLLRGVFRSPTMMVWGAMWAFSFTRVSLMKVAARALGA*MFRMLSSSCRILPLIRMKWPS";
     $reditpep .= "LSFLITLGWKSILSLIRMATPACFFRPLA*KMVFQPFTRR*CRSFSRRWVSWRQQKVGS";
     $reditpep .= "CLCSQSVRRCRLMGLFRPLMFRLIRLR*LLLPVILVVGWAWWQIISIN*SIKPFSGGSS";
     $reditpep .= "PGHHLGQLHGVAVQLLLPRSVH";
  is ($teditpep, $reditpep, "translate frame -1");

  my $tditpep = translate($orf, -2, $tCT);
  my $rditpep = "HCSGGCSGRRR*WCGERCGPSASRGCP**RWRRGRWGRRCSECSAPPAGSCP*SG*SGPPC";
     $rditpep .= "PS*SPWAGSRSCR*SGWPRQPASSGRWPERWSSSPSPAGSAGPSRAGGSPGGNRRWGPAC";
     $rditpep .= "AASRSGGAA*WGCSGR*CSG*SGSGSCCCR*SWWSGGHGGK*LALIKV*NHLAAGLLQVI";
     $rditpep .= "TSVSFTVSRFSFCFQDLS";
  is ($tditpep, $rditpep, "translate frame -2");

  my $titpep = translate($orf, -3, $tCT);
  my $ritpep = "TAQGGVQVADDDGVGSDVGLQLHEGVLDEGGGAGVGGVDVQNAQLLLQDLALDQDEVALLVL";
     $ritpep .= "LDHLGLEVDLVADQDGHASLLLQAVGLKDGLPALHPQVVPVLLAQVGLLEATEGGVLLVQP";
     $ritpep .= "VGQAVPLDGAVQAVDVQADQAQVVVAAGDLGGRVGMVANN*H*LKYKTI*RRVFSRSSPRSA";
     $ritpep .= "SRCRGSASASKICP";
  is ($titpep, $ritpep, "translate frame -3");
};

# TESTING amb_transcription()
subtest "ambiguous nucleotide transcription" => sub
{
  plan tests => 2;

  my @tnopep = amb_transcription($shortamb, $tCT, undef);
  my $rnopep = [qw(AGGCTT ATGCGT ACGCAT ATGCAT ATGCTT AGGCAT ACGCTT AGGCGT
    ACGCGT)];
  is_deeply(\@tnopep, $rnopep, "ambiguous transcription, no peptide");

  my @twithpep = amb_transcription($shortamb, $tCT, MH);
  my $rwithpep = [qw(ATGCAT)];
  is_deeply(\@twithpep, $rwithpep, "ambiguous transcription with peptide");

};

# TESTING degcodon_to_aas()
subtest "ambiguous codon translation" => sub
{
  plan tests => 25;
  my $rDRCT = {GCN => [qw(A)], TGY => [qw(C)], GAY => [qw(D)], GAR => [qw(E)], TTY => [qw(F)],
    GGN => [qw(G)], CAY => [qw(H)], ATH => [qw(I)], AAR => [qw(K)], YTR => [qw(L)], CTY => [qw(L)],
    ATG => [qw(M)], AAY => [qw(N)], CCN => [qw(P)], CAR => [qw(Q)], MGR => [qw(R)], CGY => [qw(R)],
    WCY => [qw(S T)], TCR => [qw(S)], ACN => [qw(T)], GTN => [qw(V)], TGG => [qw(W)], TAY => [qw(Y)],
    TAR => [qw(*)], TGA => [qw(*)]};
  foreach my $r (keys %$rDRCT)
  {
    my @test = degcodon_to_aas($r, $tCT);
    is_deeply (\@test, $rDRCT->{$r}, "ambiguous translate $r");
  }
};

# TESTING amb_translation()
subtest "ambiguous nucleotide translation" => sub
{
  plan tests => 4;

  my @tpospeps = amb_translation($shortamb, $tCT, "s");
  @tpospeps = sort @tpospeps;
  my $rpospeps = [qw(*A* *AC *AF *AL *AS *AW *AY *CF *CI *CL *CM *CV *GF *GI *GL
    *GM *GV *RF *RI *RL *RM *RV *SI *SL *SM *SV ACF ACI ACL ACM ACV AGF AGI AGL
    AGM AGV ARF ARI ARL ARM ARV ASI ASL ASM ASV DA* DAC DAF DAL DAS DAW DAY EA*
    EAC EAF EAL EAS EAW EAY ECF ECI ECL ECM ECV EGF EGI EGL EGM EGV ERF ERI ERL
    ERM ERV ESI ESL ESM ESV GCF GCI GCL GCM GCV GGF GGI GGL GGM GGV GRF GRI GRL
    GRM GRV GSI GSL GSM GSV HA* HAC HAF HAL HAS HAW HAY ICF ICI ICL ICM ICV IGF
    IGI IGL IGM IGV IRF IRI IRL IRM IRV ISI ISL ISM ISV KA* KAC KAF KAL KAS KAW
    KAY KCF KCI KCL KCM KCV KGF KGI KGL KGM KGV KH KP KR KRF KRI KRL KRM KRV KSI
    KSL KSM KSV LCF LCI LCL LCM LCV LGF LGI LGL LGM LGV LRF LRI LRL LRM LRV LSI
    LSL LSM LSV MH ML MP MR NA* NAC NAF NAL NAS NAW NAY PCF PCI PCL PCM PCV PGF
    PGI PGL PGM PGV PRF PRI PRL PRM PRV PSI PSL PSM PSV QA* QAC QAF QAL QAS QAW
    QAY QCF QCI QCL QCM QCV QGF QGI QGL QGM QGV QRF QRI QRL QRM QRV QSI QSL QSM
    QSV RCF RCI RCL RCM RCV RGF RGI RGL RGM RGV RH RL RR RRF RRI RRL RRM RRV RSI
    RSL RSM RSV SCF SCI SCL SCM SCV SGF SGI SGL SGM SGV SRF SRI SRL SRM SRV SSI
    SSL SSM SSV TCF TCI TCL TCM TCV TGF TGI TGL TGM TGV TH TL TP TR TRF TRI TRL
    TRM TRV TSI TSL TSM TSV VCF VCI VCL VCM VCV VGF VGI VGL VGM VGV VRF VRI VRL
    VRM VRV VSI VSL VSM VSV YA* YAC YAF YAL YAS YAW YAY)];
  is_deeply(\@tpospeps, $rpospeps, "ambiguous nucleotide translation 6 frame");

  my @tpospepst = amb_translation($shortamb, $tCT, "t");
  @tpospepst = sort @tpospepst;
  my $rpospepst = [qw(*A* *AC *AF *AL *AS *AW *AY *CF *CI *CL *CM *CV *GF *GI
    *GL *GM *GV *RF *RI *RL *RM *RV ACF ACI ACL ACM ACV AGF AGI AGL AGM AGV ARF
    ARI ARL ARM ARV DA* DAC DAF DAL DAS DAW DAY EA* EAC EAF EAL EAS EAW EAY ECF
    ECI ECL ECM ECV EGF EGI EGL EGM EGV ERF ERI ERL ERM ERV GCF GCI GCL GCM GCV
    GGF GGI GGL GGM GGV GRF GRI GRL GRM GRV HA* HAC HAF HAL HAS HAW HAY ICF ICI
    ICL ICM ICV IGF IGI IGL IGM IGV IRF IRI IRL IRM IRV KA* KAC KAF KAL KAS KAW
    KAY KCF KCI KCL KCM KCV KGF KGI KGL KGM KGV KRF KRI KRL KRM KRV LCF LCI LCL
    LCM LCV LGF LGI LGL LGM LGV LRF LRI LRL LRM LRV MH ML MR NA* NAC NAF NAL NAS
    NAW NAY PCF PCI PCL PCM PCV PGF PGI PGL PGM PGV PRF PRI PRL PRM PRV QA* QAC
    QAF QAL QAS QAW QAY QCF QCI QCL QCM QCV QGF QGI QGL QGM QGV QRF QRI QRL QRM
    QRV RCF RCI RCL RCM RCV RGF RGI RGL RGM RGV RH RL RR RRF RRI RRL RRM RRV SCF
    SCI SCL SCM SCV SGF SGI SGL SGM SGV SRF SRI SRL SRM SRV TCF TCI TCL TCM TCV
    TGF TGI TGL TGM TGV TH TL TR TRF TRI TRL TRM TRV VCF VCI VCL VCM VCV VGF VGI
    VGL VGM VGV VRF VRI VRL VRM VRV YA* YAC YAF YAL YAS YAW YAY)];
  is_deeply(\@tpospepst, $rpospepst, "ambiguous nucleotide translation 3 frame");

  my @tpospeps1 = amb_translation($shortamb, $tCT, 1);
  @tpospeps1 = sort @tpospeps1;
  my $rpospeps1 = [qw(MH ML MR RH RL RR TH TL TR)];
  is_deeply(\@tpospeps1, $rpospeps1, "ambiguous nucleotide translation 1 frame");

  my @tpospeps3 = amb_translation($shortamb, $tCT, -3);
  @tpospeps3 = sort @tpospeps3;
  my $rpospeps3 = [qw(*CI *CL *CM *CV *RI *RL *RM *RV *SI *SL *SM *SV ACI ACL
    ACM ACV ARI ARL ARM ARV ASI ASL ASM ASV ECI ECL ECM ECV ERI ERL ERM ERV ESI
    ESL ESM ESV GCI GCL GCM GCV GRI GRL GRM GRV GSI GSL GSM GSV ICI ICL ICM ICV
    IRI IRL IRM IRV ISI ISL ISM ISV KCI KCL KCM KCV KRI KRL KRM KRV KSI KSL KSM
    KSV LCI LCL LCM LCV LRI LRL LRM LRV LSI LSL LSM LSV PCI PCL PCM PCV PRI PRL
    PRM PRV PSI PSL PSM PSV QCI QCL QCM QCV QRI QRL QRM QRV QSI QSL QSM QSV RCI
    RCL RCM RCV RRI RRL RRM RRV RSI RSL RSM RSV SCI SCL SCM SCV SRI SRL SRM SRV
    SSI SSL SSM SSV TCI TCL TCM TCV TRI TRL TRM TRV TSI TSL TSM TSV VCI VCL VCM
    VCV VRI VRL VRM VRV VSI VSL VSM VSV)];
  is_deeply(\@tpospeps3, $rpospeps3, "ambiguous nucleotide translation -3 frame");
};

# TESTING pattern_remover()
subtest "pattern removal" => sub
{
  plan tests => 2;

  my $tnewshort = pattern_remover("GACAGATCT", "CAGATCT", $tCT, $tRSCU);
  my $rnewshort = "GATAGATCT";
  is ($tnewshort, $rnewshort, "remove pattern scalar");

  my $tnewshorta = pattern_remover("GACAGATCT", ["AGA", "TCT"], $tCT, $tRSCU);
  my $rnewshorta = "GACCGTTCC";
  is ($tnewshorta, $rnewshorta, "remove pattern arrayref");
};

# TESTING pattern_aligner()
subtest "pattern alignment" => sub
{
  plan tests => 3;

  my $tnewaligned = pattern_aligner("GACAGATCT", "CCGGAGC", "DRS", $tCT);
  my $rnewaligned = "NNCCGGAGC";
  is ($tnewaligned, $rnewaligned, "align in frame 3");

	my $tnewaligned2 = pattern_aligner("GACAGATCT", "GACCGGA", "DRS", $tCT);
	my $rnewaligned2 = "GACCGGANN";
	is ($tnewaligned2, $rnewaligned2, "align in frame 1");
	
	my $tnewaligned3 = pattern_aligner("GACAGATCT", "ACAGATC", "DRS", $tCT);
	my $rnewaligned3 = "NACAGATCN";
	is ($tnewaligned3, $rnewaligned3, "algin in frame 2");
};

# TESTING pattern_adder()
my $tpattadd = pattern_adder("GACAGATCT", "NNCCGGAGC", $tCT);
my $rpattadd = "GACCGGAGC";
is($tpattadd, $rpattadd, "pattern adding");

# TESTING pattern_finder()
my @tstops = pattern_finder($orf, "*", 2, 2, $tCT);
my $rstops = [qw(8 13 15 18 19 32)];
is_deeply (\@tstops, $rstops, "pattern finding");

# TESTING define_aa_defaults()
my $tdefaults = define_aa_defaults($tCT, $tRSCU);
my $rdefaults = {A => "GCT", C => "TGT", D => "GAC", E => "GAA", F => "TTC",
  G => "GGT", H => "CAC", I => "ATC", K => "AAG", L => "TTG", M => "ATG",
  N => "AAC", P => "CCA", Q => "CAA", R => "AGA", S => "TCT", T => "ACC",
  V => "GTT", W => "TGG", Y => "TAC", "*" => "TAA"};
is_deeply($tdefaults, $rdefaults, "define amino acid defaults");

# TESTING reverse_translate()
my $trevorf = reverse_translate(translate($orf, 1, $tCT), $tdefaults);
my $rrevorf = "ATGGACAGATCTTGGAAGCAAAAGTTGAACAGAGACACCGTTAAGTTGACCGAAGTTATGACC";
   $rrevorf .= "TGGAGAAGACCAGCTGCTAAGTGGTTCTACACCTTGATCAACGCTAACTACTTGCCACCATG";
   $rrevorf .= "TCCACCAGACCACCAAGACCACAGACAACAACAATTGCCAGAACCAGACCAACCAGAACACC";
   $rrevorf .= "AAAGACCAGAACAACCACACCAAGCTGCTCCACCAGACAGATTGGCTGCTCAAGCTGGTCCA";
   $rrevorf .= "CACTTGTTGTTGCCACCAGGTGACCCACCAGCTAGAGAAGGTCCAGCTTTGCCAGCTGGTGA";
   $rrevorf .= "AGGTTTGGAAGACCACTTGTCTGGTCAAAGACCAGAAGAAGCTGGTTGGAGAGGTCACCCAG";
   $rrevorf .= "ACCAAAGACAAGACAGATTGCCAGCTCAAGGTGACCAAGAAGGTCAAGGTGGTCCATTGCAC";
   $rrevorf .= "CCAGACCAAGGTCAAGACCCAGCTGGTGGTGCTGAACACTCTGAACACTTGAGACCACAAAG";
   $rrevorf .= "ACCAAGAAGACACTTGCACCAAGGTCACCCAAGAGAAGCTGAAGGTCCACACAGATCTCCAC";
   $rrevorf .= "ACCACCACAGAAGAAGACCAGAACACCCACCAGAACAATAA";
is ($trevorf, $rrevorf, "reverse translate");

#TESTING codon_count()
subtest "codon count" => sub
{
  plan tests => 2;

  my $tcount = codon_count($orf, $tCT);
  my $rcount = {TTT => 1, TTA => 1, TTG => 2, CTT => 5, CTA => 3, CTG => 5,
    ATT => 1, ATG => 2, GTG => 2, TCT => 2, TCC => 1, TCA => 1, CCT => 15,
    CCC => 6, CCA => 9, CCG => 3, ACT => 1, ACC => 3, GCT => 5, GCC => 2,
    GCA => 6, GCG => 3, TAT => 2, CAT => 11, CAC => 7, CAA => 15, CAG => 6,
    AAT => 2, AAC => 1, AAA => 1, AAG => 3, GAT => 7, GAC => 6, GAA => 12,
    GAG => 4, TGC => 1, TGA => 1, TGG => 4, CGT => 3, CGC => 6, CGA => 5,
    CGG => 4, AGA => 3, GGT => 2, GGC => 4, GGA => 8, GGG => 3, ACG => 0,
    TAG => 0, GTC => 0, AGT => 0, TGT => 0, ATC => 0, AGC => 0, TAC => 0,
    ACA => 0, TCG => 0, AGG => 0, GTT => 0, ATA => 0, TTC => 0, CTC => 0,
    TAA => 0, GTA => 0};
  is_deeply($tcount, $rcount, "codon count starting blank");

  my $tcount2 = codon_count($orf, $tCT, $rcount);
  my %rcount2 = map {$_ => $tcount->{$_}+$tcount->{$_}} keys %$tcount;
  is_deeply($tcount2, \%rcount2, "codon count starting from a count");

};

#TESTING generate_RSCU_values()
my $tcount = codon_count($orf, $tCT);
my $tORSCU = generate_RSCU_values($tcount, $tCT);
my $rORSCU = {TTT => "2.00", TTA => 0.38, TTG => 0.75, CTT => 1.88,
  CTA => 1.12, CTG => 1.88, ATT => "3.00", ATG => "1.00", GTG => "4.00",
  TCT => "3.00", TCC => "1.50", TCA => "1.50", CCT => 1.82, CCC => 0.73,
  CCA => 1.09, CCG => 0.36, ACT => "1.00", ACC => "3.00", GCT => 1.25,
  GCC => "0.50", GCA => "1.50", GCG => 0.75, TAT => "2.00", CAT => 1.22,
  CAC => 0.78, CAA => 1.43, CAG => 0.57, AAT => 1.33, AAC => 0.67,
  AAA => "0.50", AAG => "1.50", GAT => 1.08, GAC => 0.92, GAA => "1.50",
  GAG => "0.50", TGC => "2.00", TGA => "3.00", TGG => "1.00", CGT => 0.86,
  CGC => 1.71, CGA => 1.43, CGG => 1.14, AGA => 0.86, GGT => 0.47,
  GGC => 0.94, GGA => 1.88, GGG => 0.71, TAG => "0.00", GTC => "0.00",
  AGT => "0.00", TGT => "0.00", ATC => "0.00", AGC => "0.00", TAC => "0.00",
  TCG => "0.00", ACA => "0.00", GTT => "0.00", AGG => "0.00", ATA => "0.00",
  TTC => "0.00", CTC => "0.00", TAA => "0.00", ACG => "0.00", GTA => "0.00"};
is_deeply($tORSCU, $rORSCU, "generate RSCU values");

#TESTING define_codon_percentages()
my $tpercents = define_codon_percentages($tCT, $tRSCU);
my $rpercents = {AAA => 0.08, AAC => 0.97, AAG => 0.92, AAT => 0.03, ACA => 0,
  ACC => 0.5375, ACG => 0.0025, ACT => 0.4575, AGA => 0.895,
  AGC => 0.0266666666666667, AGG => 0, AGT => 0.01, ATA => 0, ATC => 0.58,
  ATG => 1, ATT => 0.42, CAA => 0.99, CAC => 0.84, CAG => 0.01, CAT => 0.16,
  CCA => 0.9425, CCC => 0.005, CCG => 0, CCT => 0.0525, CGA => 0, CGC => 0,
  CGG => 0, CGT => 0.105, CTA => 0.025, CTC => 0, CTG => 0.00333333333333333,
  CTT => 0.00333333333333333, GAA => 0.99, GAC => 0.65, GAG => 0.01,
  GAT => 0.35, GCA => 0.0075, GCC => 0.2225, GCG => 0, GCT => 0.7725, GGA => 0,
  GGC => 0.015, GGG => 0.005, GGT => 0.98, GTA => 0, GTC => 0.4775,
  GTG => 0.005, GTT => 0.5175, TAA => 0.333333333333333, TAC => 0.97, TAG => 0,
  TAT => 0.03, TCA => 0.0133333333333333, TCC => 0.403333333333333,
  TCG => 0.00333333333333333, TCT => 0.543333333333333, TGA => 0, TGC => 0.1,
  TGG => 1, TGT => 0.9, TTA => 0.0816666666666667, TTC => 0.905, TTG => 0.89,
  TTT => 0.095};
is_deeply($tpercents, $rpercents, "define codon percentages");

#TESTING change_codons
subtest "change codons" => sub
{
  plan tests => 5;
  my $tccrand = change_codons($shortorf, $tCT, $tRSCU, 0);
  isnt($tccrand, $shortorf, "random codon replacement");

  my $tccopt = change_codons($shortorf, $tCT, $tRSCU, 1);
  my $rccopt = "ATGGACAGATCTTGGAAGCAAAAGTTGAACAGA";
  is($tccopt, $rccopt, "optimal codon replacement");

  my $tcclopt = change_codons($shortorf, $tCT, $tRSCU, 2);
  my $rcclopt = "ATGGATCGTTCCTGGAAACAAAAATTGAATAGA";
  is($tcclopt, $rcclopt, "less optimal codon replacement");

  my $tccmdiff = change_codons($shortorf, $tCT, $tRSCU, 3);
  my $rccmdiff = "ATGGATCGTAGCTGGAAACAAAAATTAAATAGA";
  is ($tccmdiff, $rccmdiff, "most different codon replacement");

  my $tccldiff = change_codons($shortorf, $tCT, $tRSCU, 4);
  my $rccldiff = "ATGGATAGATCCTGGAAGCAGAAGCTTAACCGA";
  is ($tccldiff, $rccldiff, "least different codon replacement");
};

#TESTING orf_finder()
my $torfmap = orf_finder($orf, $tCT);
my $rorfmap = [[1, 0, 199], [2, 34, 165], [-1, 11, 27], [-1, 39, 50]];
is_deeply($torfmap, $rorfmap, "ORF finding");

#TESTING is_ORF()
subtest "is ORF" => sub
{
  plan tests => 3;

  my $tisorf = is_ORF($orf, $tCT);
  my $tisntorf = is_ORF($notorf, $tCT);
  my $tisorf2 = is_ORF($nostop, $tCT);
  is($tisorf, 1, "true ORF");
  is($tisntorf, 0, "not an ORF");
  is($tisorf, 1, "stopless ORF");
}
