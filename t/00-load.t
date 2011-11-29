#!perl -T

use Test::More tests => 6;

BEGIN {
    use_ok( 'Bio::GeneDesign::Basic' ) || print "Bail out!\n";
    use_ok( 'Bio::GeneDesign::HTML' ) || print "Bail out!\n";
    use_ok( 'Bio::GeneDesign::SufTree' ) || print "Bail out!\n";
    use_ok( 'Bio::GeneDesign::Codons' ) || print "Bail out!\n";
    use_ok( 'Bio::GeneDesign::RestrictionEnzymes' ) || print "Bail out!\n";
    use_ok( 'Bio::GeneDesign::Random' ) || print "Bail out!\n";
}

diag( "Testing Bio::GeneDesign::Basic $Bio::GeneDesign::Basic::VERSION, Perl $], $^X" );
