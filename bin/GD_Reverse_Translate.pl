#!/usr/bin/env perl

use Bio::GeneDesign;
use File::Basename;
use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $VERSION = '5.52';
my $GDV = "GD_Reverse_Translate_$VERSION";
my $GDS = "_RT";

local $| = 1;

my %p = ();
GetOptions (
      'input=s'     => \$p{INPUT},
      'output=s'    => \$p{OUTPUT},
      'format=s'    => \$p{FORMAT},
      'algorithm=s' => \$p{ALGORITHM},
      'rscu=s'      => \$p{FILES},
      'organism=s'  => \$p{ORGS},
      'split'       => \$p{SPLIT},
      'help'        => \$p{HELP}
);

################################################################################
################################ SANITY  CHECK #################################
################################################################################
pod2usage(-verbose=>99, -sections=>"NAME|VERSION|DESCRIPTION|ARGUMENTS|USAGE")
  if ($p{HELP});

my $GD = Bio::GeneDesign->new();
  
#The input file must exist and be a format we care to read.
die "\n GDERROR: You must supply an input file.\n"
  if (! $p{INPUT});
my ($iterator, $filename, $suffix) = $GD->import_seqs($p{INPUT});

$p{FORMAT} = $p{FORMAT} || $suffix || "genbank";

#The output path must exist, and we'll need it to end with a slash
$p{OUTPUT} = $p{OUTPUT} || ".";
$p{OUTPUT} .= "/" if (substr($p{OUTPUT}, -1, 1) !~ /[\/]/);
die "\n GDERROR: $p{OUTPUT} does not exist.\n"
  if ($p{OUTPUT} && ! -e $p{OUTPUT});
  
#We must get a list of organisms or a set of rscu files
die "\n GDERROR: Neither an organism nor an RSCU table were supplied.\n"
  if (! $p{ORGS} && ! $p{FILES});

$p{ORGS}       = $p{ORGS}       ? $p{ORGS}       : q{};
$p{FILES}      = $p{FILES}      ? $p{FILES}      : q{};
$p{ALGORITHM}  = $p{ALGORITHM}  ? $p{ALGORITHM}  : "high";

################################################################################
################################# CONFIGURING ##################################
################################################################################
my @fileswritten;
my @seqstowrite;

my %works = ();
foreach my $org (split (",", $p{ORGS}))
{
  $works{$org} = {on => $org, path => undef};
}
foreach my $file ( split ( ",", $p{FILES} ) )
{
  $works{$file} = {on => basename($file), path => $file};
}

################################################################################
############################# REVERSE TRANSLATING ##############################
################################################################################
while ( my $obj = $iterator->next_seq() )
{
  foreach my $work (keys %works)
  {
    $GD->set_organism(
        -organism_name => $works{$work}->{on},
        -rscu_path     => $works{$work}->{path});
    
    my $newobj = $GD->reverse_translate(
        -peptide => $obj,
        -algorithm => $p{ALGORITHM});
        
    if ($p{SPLIT})
    {
      my $outputfilename = $obj->id . $GDS . "." . $p{FORMAT};
      $GD->export_seqs(
          -filename   => $outputfilename,
          -path       => $p{OUTPUT},
          -format     => $p{FORMAT},
          -sequences  => [$obj]);
      push @fileswritten, $outputfilename;
    }

    else
    {
      push @seqstowrite, $newobj;
    }
  }
}
if (scalar @seqstowrite)
{
  my $outputfilename = $filename . $GDS . "." . $p{FORMAT};
  $GD->export_seqs(
          -filename   => $outputfilename,
          -path       => $p{OUTPUT},
          -format     => $p{FORMAT},
          -sequences  => \@seqstowrite);
  push @fileswritten, $outputfilename;
}

print "\n";
print "Wrote " . $p{OUTPUT} . "$_\n" foreach @fileswritten;
print "\n";
print $GD->attitude() . " brought to you by $GDV\n\n";

exit;

__END__

=head1 NAME

  GD_Reverse_Translate.pl

=head1 VERSION

  Version 5.52

=head1 DESCRIPTION

  Given at least one protein sequence as input, the Reverse_Translate script
  generates synonymous nucleotide sequences using either a user-defined RSCU
  table or an RSCU table provided by GeneDesign
  
  If no algorithm is specified, the balanced algorithm will be used. These are
  the algorithms provided by default with GeneDesign; you can make your own; see
  developer docs.

  Output will be named according to the name of the input file, and will be
  tagged with the suffix _RT.
  
  Algorithms:
    high: The high algorithm replaces every codon in the input sequence with
        the most translationally optimal codon as specified by the input RSCU
        tables or known RSCU tables (if organism is specified). If the codon is
        already the ideal codon it is left alone.
    balanced: The balanced algorithm uses the rscu data to determine a
         likelihood of codon replacement.
         
=head1 USAGE

  -r OR -o must be provided. If both are given the table will be treated as
      another organism, named after the table's filename.
      
=head1 ARGUMENTS

Required arguments:

  -i,   --input : a file containing protein sequences.
  -org, --organism : an organism whose RSCU table can be found in the config
    directory, or several separated by commas
  AND/OR
  -r,   --rscu : path to an RSCU table generated by GD_Generate_RSCU_Table.pl,
    or several separated by commas
  -r OR -o must be provided. If both are given the table will be treated as
      another organism, named after the table's filename.

Optional arguments:

  -a,   --algorithm : protocol for replacing peptides; defaults to "balanced"
  -out, --output : path to an output directory; defaults to GeneDesign's default
  -f,   --format : default genbank
  -s,   --split : output all sequences as separate files
  -h,   --help : Display this message

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013, GeneDesign developers
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

* The names of Johns Hopkins, the Joint Genome Institute, the Lawrence Berkeley
National Laboratory, the Department of Energy, and the GeneDesign developers may
not be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE DEVELOPERS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=cut