#
# GeneDesign module for sequence segmentation
#

=head1 NAME

Bio::GeneDesign::RestrictionEnzymes

=head1 VERSION

Version 3.05

=head1 DESCRIPTION

GeneDesign functions for handling restriction enzymes

=head1 AUTHOR

Sarah Richardson <notadoctor@jhu.edu>

=cut

package Bio::GeneDesign::RestrictionEnzymes;
require Exporter;

use Bio::GeneDesign::Basic qw(:all);
use Bio::GeneDesign::Codons qw(:all);
use Bio::GeneDesign::RestrictionEnzyme;
use Bio::GeneDesign::SufTree;
use Perl6::Slurp;
use Carp;

use strict;
use warnings;
use warnings::unused;

use vars qw($VERSION @ISA @EXPORT_OK %EXPORT_TAGS);
$VERSION = 3.05;

@ISA = qw(Exporter);
@EXPORT_OK = qw(
  define_sites
  define_site_status
  siteseeker
  mutexclu
  first_base
  report_RE
  build_suffix_tree
  search_suffix_trees
  overhangP
  overhangA
  filter_by_stickiness
  filter_by_cleavagesite
  filter_by_overhangpalindromy
  filter_by_length
  filter_by_baseambiguity
  filter_by_inactivationtemp
  filter_by_incubationtemp
  filter_by_staractivity
  filter_by_cpgsensitivity
  filter_by_damsensitivity
  filter_by_dcmsensitivity
  filter_by_bufferactivity
  filter_by_vendor
  filter_by_price
  filter_by_sequence
  $IIPreg $IIAreg $IIBreg %REVENDORS
  $VERSION
);

%EXPORT_TAGS =  (all => \@EXPORT_OK);

our $IIPreg  = qr/   ([A-Z]*)   \^ ([A-Z]*)      /x;
our $IIAreg  = qr/\A \w+ \(([\-]*\d+) \/ ([\-]*\d+)\)\Z  /x;
our $IIBreg  = qr/\A\(([\-]*\d+) \/ ([\-]*\d+)\) \w+ \(([\-]*\d+) \/ ([\-]*\d+)\)\Z  /x;


our %REVENDORS = (B => "Invitrogen", C => "Minotech", E => "Stratagene", 
  F => "Thermo Scientific Fermentas", I => "SibEnzyme", J => "Nippon Gene Co.", 
  K => "Takara", M => "Roche Applied Science", N => "New England Biolabs", 
  O => "Toyobo Technologies", Q => "Molecular Biology Resources", 
  R => "Promega", S => "Sigma Aldrich", U => "Bangalore Genei", V => "Vivantis",
  X => "EURx", Y => "CinnaGen");

=head2 define_sites()

    Generates a hash reference where the keys are enzyme names and the values are
    Bio::GeneDesign::RestrictionEnzyme objects.

=cut

sub define_sites
{
  my ($file) = @_;
  my @data = split(/\n/, slurp($file));
  my $swit = $data[0] =~ /^name/ ? 1  :  0;
  shift @data if ($swit == 1);
  my %RES;
  foreach my $line (@data)
  {
    my ($name, $site, $temp, $inact, $buf1, $buf2, $buf3, $buf4, $bufu, $dam,
        $dcm, $cpg, $score, $star, $vendor) = split("\t", $line);
    my $recseq = $site;
    $recseq =~ s/\W*\d*//g;
    my %vendhsh = map {$_ => $REVENDORS{$_}} split ("", $vendor);
    my $buffhsh = {NEB1 => $buf1, NEB2 => $buf2, NEB3 => $buf3, NEB4 => $buf4, 
                   Other => $bufu};
    $star = "n" unless ($star eq "y");
    my $re = Bio::GeneDesign::RestrictionEnzyme->new(
                  -id => $name,
                  -cutseq => $site,
                  -recseq => $recseq,
                  -temp   => $temp,
                  -tempin => $inact,
                  -score  => $score,
                  -methdam => $dam,
                  -methdcm => $dcm,
                  -methcpg => $cpg,
                  -staract => $star,
                  -vendors => \%vendhsh,
                  -buffers => $buffhsh);
    $RES{$re->id} = $re;
  }
  #Make exclusion lists
  foreach my $re (values %RES)
  {
    my %excl;
    foreach my $ar (sort grep {$_->id ne $re->id} values %RES)
    {
      foreach my $arreg (@{$ar->regex})
      {
        $excl{$ar->id}++ if ($re->recseq =~ $arreg)
      }
      foreach my $rereg (@{$re->regex})
      {
        $excl{$ar->id}++ if ($ar->recseq =~ $rereg)
      }
    }
    my @skips = keys %excl;
    $re->exclude(\@skips);
  }
  return \%RES;
}

=head2 define_site_status()

  Generates a hash describing a restriction count of a nucleotide sequence.
  in: nucleotide sequence as a string
      reference to a hash containing enzyme names as keys and regular
        expressions in a reference to an array as values
        (SITE_REGEX from define_sites)
  out: reference to a hash where the keys are enzyme names and the value is a
        count of their occurence in the nucleotide sequence

=cut

sub define_site_status
{
  my($seq, $RES) = @_;
  my $SITE_STATUS = {};
  foreach my $re (values %$RES)
  {
    my $count = 0;
    foreach my $sit (@{$re->regex})
    {
      $count++ while ($seq =~ /(?=$sit)/ig);
    }
    $SITE_STATUS->{$re->id} = $count;
  }
  return $SITE_STATUS;    
}

=head2 siteseeker()

  Generates a hash describing the positions of a particular enzyme's recognition
  sites in a nucleotide sequence.
  in: nucleotide sequence as a string,
      Bio::GeneDesign::RestrictionEnzyme object
  out: reference to a hash where the keys are positions and the value is the
      recognition site from that position as a string.

=cut

sub siteseeker
{
  my ($seq, $enzyme) = @_;
  croak("$enzyme is not a Bio::GeneDesign::RestrictionEnzyme object\n")
    unless (ref $enzyme && $enzyme->isa("Bio::GeneDesign::RestrictionEnzyme") );
  my $total = {};
  foreach my $sit (@{$enzyme->regex})
  {
    while ($seq =~ /(?=$sit)/ig)
    {
      $total->{pos $seq} = substr($seq, pos $seq, $enzyme->length);
    }
  }
  return $total;
}

=head2 build_suffix_tree()

  Builds two suffix trees from a list of restriction enzyme recognition sites
  in: array of enzyme names,
      reference enzyme hash
  out: array of suffix trees, one for forward orientaion, one for reverse

  #NO UNIT TESTS

=cut

sub build_suffix_tree
{
  my ($list_ref, $CODON_TABLE, $xlationref) = @_;
  my ($fehsh, $rehsh) = ({}, {});
  my $ftree = new_aa Bio::GeneDesign::SufTree();
  my $rtree = new_aa Bio::GeneDesign::SufTree();
  foreach my $enzyme (@$list_ref)
  {
    my $site = $enzyme->recseq;
    #print "\n adding $site to tree\n";
    my $lagcheck = $site;
    $lagcheck = substr($lagcheck, 0, length($lagcheck)-1) while (substr($lagcheck, -1) eq "N");
    foreach my $peptide (amb_translation($site, $CODON_TABLE, 't', $xlationref))
    {
      $fehsh->{$peptide} = [] unless (exists $fehsh->{$peptide});
      push @{$$fehsh{$peptide}}, $enzyme;
    }
    my $etis = complement($site, 1);
    if ($etis ne $site && $lagcheck eq $site)
    {
      #print "\n adding $etis to tree\n";
      foreach my $peptide (amb_translation($etis, $CODON_TABLE, 't', $xlationref))
      {
        $rehsh->{$peptide} = [] unless (exists $rehsh->{$peptide});
        push @{$$rehsh{$peptide}}, $enzyme;
      }
    }
  }
  $ftree->add_aa_paths($fehsh);
  $rtree->add_aa_paths($rehsh);
  return [$ftree, $rtree];
}

=head2 search_suffix_trees()

  #NO UNIT TESTS

=cut

sub search_suffix_trees
{
  my ($treelist, $nucseq, $aaseq, $CODON_TABLE, $xlationref) = @_;
  my $results;
  my $flip = 0;
  foreach my $tree (@$treelist)
  {
  #print "\tSearching $flip tree...\n" if ($debug);
    $flip++;
    foreach my $rog ($tree->find_aa_paths($aaseq))
    {
      my $enzyme    = $$rog[0];
      my $enz       = $enzyme->id;
      my $recsite   = $enzyme->recseq();
         $recsite   = complement($recsite, 1) if ($flip > 1);
      my $nucstart  = $$rog[1] * 3;
      my $peptide   = $$rog[2];
      my $presence  = "p";
      my $ohang     = {};
      my ($mattersbit, $ohangstart, $ohangend, $fabric) = ("", 0, 0, "");
      my ($offset, $situ, $pproof, $sitelen) = (0, "", "", 0);

      #print "\tgot $enz @ $nucstart, $peptide\n" if ($debug);
    ##Figure Possible overhangs
      if ($enzyme->type() eq "b" || $enzyme->class() eq "IIB")
      {
        $situ = substr($nucseq, $nucstart, length($peptide)*3);
        $pproof = translate($situ, 1, $CODON_TABLE);
        $ohangstart = 0;
      }
      elsif ($enzyme->class() eq "IIP")
      {
        my ($lef, $rig) = (length($1), length($2)) if ($enzyme->cutseq() =~ $IIPreg);
        ($lef, $rig) = ($rig, $lef) if ($rig < $lef);
        $ohangstart = $enzyme->length - $rig + 1;
        $ohangend = $enzyme->length - $lef;
        $situ = substr($nucseq, $nucstart, length($peptide)*3);
        ($fabric, $offset) = pattern_aligner($situ, $recsite, $peptide, $CODON_TABLE, 1, $xlationref);
        $mattersbit = substr($situ, $ohangstart + $offset -1, $ohangend - $ohangstart + 1);
        $pproof = translate($situ, 1, $CODON_TABLE);
        $sitelen = $enzyme->length;
      }
      elsif ($enzyme->class() eq "IIA")
      {
        my ($lef, $rig) = ($1, $2)  if ($enzyme->cutseq() =~ $IIAreg);
        ($rig, $lef) = ($lef, $rig) if ($rig < $lef);
        $sitelen = $rig >= 0 ? $enzyme->length + $rig  : $enzyme->length;
        my $nuclen = length($peptide)*3;
        $nuclen++ while($nuclen % 3 != 0);
        $situ = substr($nucseq, $nucstart, $nuclen);
        ($fabric, $offset) = pattern_aligner($situ, $recsite, $peptide, $CODON_TABLE, 1, $xlationref);
        my $add;
        if ($flip == 1)
        {
          $ohangstart = $enzyme->length + $lef + 1;
          unless ($rig <= 0)
          {
            $add = $rig - (length($fabric) - ($offset + $enzyme->length));
            $add ++ while ($add % 3 != 0);
            $situ .= substr($nucseq, $nucstart + $nuclen, $add);
            $fabric .= "N" x $add;
          }

        }
        else
        {
          unless ($rig <= 0)
          {
            $add =  $rig - $offset;
            $add ++ while ($add % 3 != 0);
            $situ = substr($nucseq, $nucstart - $add, $add) . $situ;
            $fabric = "N" x $add . $fabric;
            $nucstart = $nucstart - $add;
            $ohangstart = $add - $rig + 1;
          }
          else
          {
            $ohangstart = $offset + abs($rig) + 1;
          }
        }
        $mattersbit = substr($nucseq, $nucstart + $ohangstart+1, $rig-$lef);
        $pproof = translate($situ, 1, $CODON_TABLE);
      }
      else
      {
        print "I don't recognize this type of enzyme: " . $enzyme->id;
      }
      unless ($enzyme->type() eq "b" || $enzyme->class() eq "IIB")
      {
        if ($fabric eq "0")
        {
          print "oh no bad fabric, $enz, $fabric, $peptide\n";
          next;
        }
        my $lenm = $mattersbit  ? length($mattersbit) : 0;
        my $matstart = $ohangstart + $offset - 1;
           $matstart-- while($matstart % 3 != 0);
        my $matend = $ohangstart + $offset + $lenm - 1;
           $matend++ while($matend % 3 != 2);
        my $matlen = $matend - $matstart + 1;
        my $peproof = substr($pproof, ($matstart/3), $matlen/3);
        my $what = substr($fabric, $matstart, $matlen);
      #print "\n\n$flip $enz, $peptide, fab: $fabric, situ: $situ, pproof: $pproof, " if ($debug);
      #print "ohangstart: $ohangstart, ohangend: $ohangend, matters: $mattersbit, " if ($debug);
      #print "offset: $offset, matstart: $matstart, matlen: $matlen, peproof: $peproof, what: $what\n" if ($debug);
        foreach my $swapseq (amb_transcription($what, $CODON_TABLE, $peproof))
        {
          substr($fabric, $matstart, $matlen) = $swapseq;
          my $tohang = substr($fabric, $ohangstart +  $offset - 1, $lenm);
        #print "\t\t\tsubbing $swapseq into $fabric\tconsidering $tohang\n" if ($debug);
          $$ohang{$tohang}++ if ($tohang);#if ($tohang ne complement($tohang, 1));
        }
      }

    ##Determine Presence
      $presence = "e" if ($situ =~ $enzyme->regex()->[$flip - 1]);
      #print "got $enz @ $nucstart, $presence, $peptide, $fabric\n" if ($presence eq "e" && $debug);
      unless ($enzyme->class() eq "IIB" && $presence eq "p")
      {
        my $ohangoffset = $ohangstart + $offset - 1;
        push @$results, [$enzyme, $nucstart, $presence, $sitelen, $ohang, $pproof, $ohangoffset, $mattersbit, $flip];
      }      
    }
  }
  return $results;
}


=head2 first_base

  #NO UNIT TESTS

=cut

sub first_base
{
  my ($relev, $index, $orient) = @_;
  my $term = (($index -($relev %3))%3);
  if ($orient eq "+")
  {
    return $index - $term;
  }
  else
  {
    return $term != 0  ?  $index + (3-$term)  :  $index;
  }
}

=head2 overhangP

  #NO UNIT TESTS

=cut

sub overhangP
{
  my ($dna, $enzyme) = @_;
  my ($ohangstart, $mattersbit) = (0, "");
  if ($enzyme->class eq "IIP")
  {
    my ($lef, $rig) = (length($1), length($2)) if ($enzyme->cutseq =~ $IIPreg);
    ($lef, $rig) = ($rig, $lef) if ($rig < $lef);
    $ohangstart = $lef + 1;
    $mattersbit = substr($dna, $ohangstart-1, $rig-$lef);
  }
  else
  {
    return undef;
  }
  return ($ohangstart, $mattersbit);
}

=head2 overhangA

  #NO UNIT TESTS

=cut

sub overhangA
{
  my ($dna, $context, $enzyme, $strand) = @_;
  my ($ohangstart, $mattersbit) = (0, "");
  if ($enzyme->class eq "IIA")
  {
    my ($lef, $rig) = ($1, $2) if ($enzyme->cutseq =~ $IIAreg);
    ($lef, $rig) = ($rig, $lef) if ($rig < $lef);
    if ($strand == 1)
    {
      $ohangstart = length($dna) + $lef + 1;
    }
    else
    {
      $ohangstart = length($context) - length($dna) - $rig + 1;
    }
    $mattersbit = substr($context, $ohangstart-1, $rig-$lef);
    $ohangstart = $strand == 1  ? length($dna) + $lef :   0 - ($rig);
  }
  return ($ohangstart, $mattersbit);
}

=head2 filter_by_stickiness

  #NO UNIT TESTS

=cut

sub filter_by_stickiness
{
  my ($RES, $hshref) = @_;
  foreach my $id (keys %$RES)
  {
    my $re = $RES->{$id};
    unless (exists $hshref->{$re->type})
    {
      next if ($re->onebpoverhang && exists $hshref->{1});
      delete $RES->{$id};
    }
  }
  return;
}

=head2 filter_by_cleavagesite

  #NO UNIT TESTS

=cut

sub filter_by_cleavagesite
{
  my ($RES, $hshref) = @_;
  foreach my $id (keys %$RES)
  {
    my $re = $RES->{$id};
    delete $RES->{$id} unless ( exists $hshref->{$re->class} );
  }
}

=head2 filter_by_overhangpalindromy

  #NO UNIT TESTS

=cut

sub filter_by_overhangpalindromy
{
  my ($RES, $hshref) = @_;
  foreach my $id (keys %$RES)
  {
    my $re = $RES->{$id};
    delete $RES->{$id} unless ( exists $hshref->{$re->palindromy} );
  }
}

=head2 filter_by_length

  #NO UNIT TESTS

=cut

sub filter_by_length
{
  my ($RES, $hshref) = @_;
  foreach my $id (keys %$RES)
  {
    my $re = $RES->{$id};
    delete $RES->{$id} unless ( exists $hshref->{length($re->recseq)} );
  }
}

=head2 filter_by_baseambiguity

  #NO UNIT TESTS

=cut

sub filter_by_baseambiguity
{
  my ($RES, $regex) = @_;
  foreach my $id (keys %$RES)
  {
    my $re = $RES->{$id};
    delete $RES->{$id} unless ( $re->recseq !~ $regex );
  }
}

=head2 filter_by_inactivationtemp

  #NO UNIT TESTS

=cut

sub filter_by_inactivationtemp
{
  my ($RES, $tempin) = @_;
  foreach my $id (keys %$RES)
  {
    my $re = $RES->{$id};
    unless ($re->tempin)
    { 
      delete $RES->{$id};
      next;
    }
    delete $RES->{$id} unless ( $re->tempin <= $tempin );
  }
}

=head2 filter_by_incubationtemp

  #NO UNIT TESTS

=cut

sub filter_by_incubationtemp
{
  my ($RES, $hshref) = @_;
  foreach my $id (keys %$RES)
  {
    my $re = $RES->{$id};
    delete $RES->{$id} unless ( exists($hshref->{$re->temp}) );
  }
}

=head2 filter_by_staractivity

  #NO UNIT TESTS

=cut

sub filter_by_staractivity
{
  my ($RES, $f) = @_;
  foreach my $id (keys %$RES)
  {
    my $re = $RES->{$id};
    delete $RES->{$id} unless ( $f eq $re->staract );
  }
}

=head2 filter_by_cpgsensitivity

  #NO UNIT TESTS

=cut

sub filter_by_cpgsensitivity
{
  my ($RES, $hshref) = @_;
  foreach my $id (keys %$RES)
  {
    my $re = $RES->{$id};
    delete $RES->{$id} unless ( exists($hshref->{$re->methcpg}) );
  }
}

=head2 filter_by_damsensitivity

  #NO UNIT TESTS

=cut

sub filter_by_damsensitivity
{
  my ($RES, $hshref) = @_;
  foreach my $id (keys %$RES)
  {
    my $re = $RES->{$id};
    delete $RES->{$id} unless ( exists($hshref->{$re->methdam}) );
  }
}

=head2 filter_by_dcmsensitivity

  #NO UNIT TESTS

=cut

sub filter_by_dcmsensitivity
{
  my ($RES, $hshref) = @_;
  foreach my $id (keys %$RES)
  {
    my $re = $RES->{$id};
    delete $RES->{$id} unless ( exists($hshref->{$re->methdcm}) );
  }
}

=head2 filter_by_bufferactivity

  #NO UNIT TESTS

=cut

sub filter_by_bufferactivity
{
  my ($RES, $hshref) = @_;
  foreach my $id (keys %$RES)
  {
    my $re = $RES->{$id};
    my $rebuff = $re->buffers;
    foreach my $buff (keys %$hshref)
    {
      my $val = $hshref->{$buff};
      delete $RES->{$id} unless ( exists($rebuff->{$buff}) && $rebuff->{$buff} eq $val );
    }
  }
}

=head2 filter_by_vendor

  #NO UNIT TESTS

=cut

sub filter_by_vendor
{
  my ($RES, $hshref) = @_;
  foreach my $id (keys %$RES)
  {
    my $re = $RES->{$id};
    my $revend = $re->vendors;
    foreach my $vend (keys %$hshref)
    {
      delete $RES->{$id} unless ( exists( $revend->{$vend} ) );
    }
  }
}

=head2 filter_by_price

  #NO UNIT TESTS

=cut

sub filter_by_price
{
  my ($RES, $p) = @_;
  foreach my $id (keys %$RES)
  {
    my $re = $RES->{$id};
    delete $RES->{$id} unless ( $re->score <= $p );
  }
}

=head2 filter_by_sequence

  #NO UNIT TESTS

=cut

sub filter_by_sequence
{
  my ($RES, $arrref, $flag) = @_;
  foreach my $id (keys %$RES)
  {
    my $re = $RES->{$id};
    foreach my $regex (@$arrref)
    {
      delete $RES->{$id} if ( $flag == 0 && $re->recseq =~ $regex );
      delete $RES->{$id} if ( $flag == 1 && $re->recseq !~ $regex );
    }
  }
}

=head2 mutexclu

  #NO UNIT TESTS

=cut

sub mutexclu
{
  my ($used_ref, $allsites_ref) = @_;
  #don't want second degree exclusion - only exclude if it's not a -2 site
  foreach my $c ( grep {$$used_ref{$_} != -2} keys %{$used_ref})
  {
    my $te = regres($$allsites_ref{$c});
    foreach my $d ( grep {$c ne $_} keys %{$allsites_ref})
    {
      my $ue = regres($$allsites_ref{$d});
      $$used_ref{$d} = -2  if ($$used_ref{$c} =~ $ue || $$allsites_ref{$d} =~ $te);
    }
  }
  return $used_ref;
}

1;

__END__

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011, GeneDesign developers
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the Johns Hopkins nor the
      names of the developers may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE DEVELOPERS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=cut
