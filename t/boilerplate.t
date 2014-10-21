#!perl -T

use 5.006;
use strict;
use warnings;
use Test::More tests => 14;

sub not_in_file_ok
{
  my ($filename, %regex) = @_;

  open(my $fh, '<', $filename) || die "couldn't open $filename for reading: $!";
  my $ref = do {local $/; <$fh>};
  close $fh;
  my @data = split(m{\n}msx, $ref);

  my %violated;

  foreach (my $line = @data)
  {
    while (my ($desc, $regex) = each %regex)
    {
      if ($line =~ $regex)
      {
        push @{$violated{$desc}||=[]}, $.;
      }
    }
  }
  close $fh;

  if (%violated)
  {
    fail("$filename contains boilerplate text");
    diag "$_ appears on lines @{$violated{$_}}" for keys %violated;
  }
  else
  {
    pass("$filename contains no boilerplate text");
  }
  return;
}

sub module_boilerplate_ok
{
  my ($module) = @_;
  not_in_file_ok($module =>
      'the great new $MODULENAME'   => qr/ - The great new /,
      'boilerplate description'     => qr/Quick summary of what the module/,
      'stub function definition'    => qr/function[12]/,
  );
  return;
}

TODO:
{
  local $TODO = "Need to replace the boilerplate text";

  not_in_file_ok(README =>
    "The README is used..."       => qr/The README is used/,
    "'version information here'"  => qr/to provide version information/,
  );

  not_in_file_ok(Changes =>
    "placeholder date/time"       => qr(Date/time)
  );

  module_boilerplate_ok('lib/Bio/GeneDesign.pm');
  module_boilerplate_ok('lib/Bio/GeneDesign/Basic.pm');
  module_boilerplate_ok('lib/Bio/GeneDesign/CodonJuggle.pm');
  module_boilerplate_ok('lib/Bio/GeneDesign/Codons.pm');
  module_boilerplate_ok('lib/Bio/GeneDesign/Graph.pm');
  module_boilerplate_ok('lib/Bio/GeneDesign/Oligo.pm');
  module_boilerplate_ok('lib/Bio/GeneDesign/Random.pm');
  module_boilerplate_ok('lib/Bio/GeneDesign/RestrictionEnzyme.pm');
  module_boilerplate_ok('lib/Bio/GeneDesign/RestrictionEnzymes.pm');
  module_boilerplate_ok('lib/Bio/GeneDesign/ReverseTranslate.pm');
  module_boilerplate_ok('lib/Bio/GeneDesign/PrefixTree.pm');
  module_boilerplate_ok('lib/Bio/GeneDesign/Exceptions.pm');
}

