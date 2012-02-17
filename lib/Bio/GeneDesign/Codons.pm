#
# GeneDesign module for codon analysis and manipulation
#

=head1 NAME

GeneDesign::Codons

=head1 VERSION

Version 3.05

=head1 DESCRIPTION

GeneDesign functions for codon analysis and manipulation

=head1 AUTHOR

Sarah Richardson <notadoctor@jhu.edu>.

=cut

package Bio::GeneDesign::Codons;
require Exporter;

use Bio::GeneDesign::Basic qw(regres %AACIDS $ambnt %NTIDES complement @NTS);
use List::Util qw(shuffle max first);
use Perl6::Slurp;
use Carp;

use strict;
use warnings;

use vars qw($VERSION @ISA @EXPORT_OK %EXPORT_TAGS);
$VERSION = 3.05;

@ISA = qw(Exporter);
@EXPORT_OK = qw(
  define_codons
  define_codon_table
  define_reverse_codon_table
  define_RSCU_values
  parse_RSCU_values
  translate
  amb_transcription
  reverse_translate
  degcodon_to_aas
  amb_translation
  pattern_finder
  pattern_remover
  pattern_adder
  pattern_aligner
  change_codons
  RSCU_filter
  is_ORF
  minimize_local_alignment_dp
  define_codon_percentages
  index_codon_percentages
  codon_count
  generate_RSCU_values
  define_aa_defaults
  orf_finder
  $strcodon
  $VERSION
);

%EXPORT_TAGS =  (
  all => \@EXPORT_OK,
  basic => [qw(define_codons define_codon_table define_RSCU_values translate
    is_ORF)]);

our $strcodon  = qr/[ATCG]{3}/;

#my @out =  map { $_ . " (" . $CODON_TABLE->{$_} . ") ". $RSCUVal->{$_} ."\n"}
#      sort {  $CODON_TABLE->{$a} cmp $CODON_TABLE->{$b}
#        ||  $RSCUVal->{$b} <=> $RSCUVal->{$a}}
#      keys %$RSCUVal;

=head2 define_codons

  Generates an array reference that contains every possible nucleotide codon
  out : list of codons (array reference)

=cut

sub define_codons
{
  my @codons;
  foreach my $a (@NTS)
  {
    foreach my $b (@NTS)
    {
      foreach my $c (@NTS)
      {
        push @codons, $a . $b . $c;
      }
    }
  }
  return \@codons;
}

=head2 define_codon_table()

  Brings into memory a codon table from a file in GeneDesign's configuration
   directory. The table is represented as a hashref, where each key is a three
   letter nucleotide codon and each value is a one letter amino acid residue.
  in: name of a .ct file and the GD config hashref
  out: codon table (hash reference)
  
=cut

sub define_codon_table
{
  my ($org, $GD) = @_;
  croak ("No organism provided\n") unless ($org);
  croak ("No configuration hashref provided\n") unless ($GD);
  my $CODON_TABLE = {};
  my $target = exists($GD->{CODONTABLES}->{$org}) ? $org  : "Standard";
  my $filepath = $GD->{codon_dir} . '/' . $target . ".ct";
  my @arr = slurp($filepath);
  my @codlines = grep {$_ !~ /\#/} @arr;
  foreach my $line (@codlines)
  {
    my ($cod, $rscu) = ($1, $2) if ($line =~ /\{([\w]{3})\} = ([\w\*]{1})/);
    $CODON_TABLE->{$cod} = $rscu;
  }
  my $codonaref = define_codons();
  foreach my $codon (@$codonaref)
  {
    croak("$target codon table missing definition for codon $codon\n")
      unless (exists $CODON_TABLE->{$codon});
  }
  return $CODON_TABLE;
}

=head2 define_reverse_codon_table()

  Takes a codon table hashref and reverses it such that each key is a one letter
   amino acid residue and each value is an array reference containing all of the
   codons that can code for that residue.
  in: codon table (hash reference)
  out: reverse codon table (hash reference)
  
=cut

sub define_reverse_codon_table
{
  my ($CODON_TABLE) = @_;
  my $REV_CODON_TABLE = {};
  foreach my $codon (keys %$CODON_TABLE)
  {
    my $aa = $CODON_TABLE->{$codon};
    $REV_CODON_TABLE->{$aa} = [] if ( ! exists $REV_CODON_TABLE->{$aa} );
    push @{$REV_CODON_TABLE->{$aa}}, $codon;
  }
  return $REV_CODON_TABLE;
}

=head2 define_RSCU_values()

  Brings into memory an RSCU table from a file in GeneDesign's configuration
   directory. The table is represented as a hashref, where each key is a three
   letter nucleotide codon and each value is a three digit float.
  in: path to a .rscu file and the GD config hashref
  out: RSCU value table (hash reference)
  
=cut

sub define_RSCU_values
{
  my ($org, $GD) = @_;
  croak ("No organism provided\n") unless ($org);
  croak ("No configuration hashref provided\n") unless ($GD);
  my $filepath = $GD->{codon_dir} . '/' . $org . ".rscu";
  my $RSCU_TABLE = parse_RSCU_values($filepath);
  return $RSCU_TABLE;
}

=head2 parse_RSCU_values

=cut

sub parse_RSCU_values
{
  my ($path) = @_;
  croak ("No path provided\n") unless ($path);
  my $RSCU_TABLE = {};
  my @arr = slurp($path);
  my @codlines = grep {$_ !~ /\#/} @arr;
  foreach my $line (@codlines)
  {
    my ($cod, $rscu) = ($1, $2) if ($line =~ /\{([\w]{3})\} = ([\d\.]+)/);
    $RSCU_TABLE->{$cod} = $rscu;
  }
  my $codonaref = define_codons();
  foreach my $codon (@$codonaref)
  {
    croak("RSCU table missing definition for codon $codon\n")
      unless (exists $RSCU_TABLE->{$codon});
  }
  return $RSCU_TABLE;
}

=head2 translate()

  takes a nucleotide sequence, a frame, and a codon table and returns that frame
  translated into amino acids.
  in: nucleotide sequence (string),
      switch for frame (+/-1, +/-2, or +/-3),
      codon table (hash reference)
  out: amino acid sequence (string)
  
=cut

sub translate
{
  my ($nucseq, $swit, $CODON_TABLE) = @_;
  croak ("translate needs an unambiguous nucleotide sequence\n")
    if ( ! $nucseq  ||  $nucseq =~ $ambnt);
  $nucseq = complement($nucseq, 1) if ($swit < 0);
  my $peptide = "";
  for (my $offset = abs($swit)-1; $offset < length($nucseq); $offset += 3)
  {
    my $codon = substr($nucseq, $offset, 3);
    $peptide .= $CODON_TABLE->{$codon} if (exists $CODON_TABLE->{$codon});
  }
  return $peptide;
}

=head2 reverse_translate()

  takes an amino acid sequence and a specific codon table and returns that frame
  translated into amino acids.  See gdRevTrans.cgi for use.
  in: nucleotide sequence (string),
      switch for frame (1, 2, or 3),
      codon table (hash reference)
  out: amino acid sequence (string)

=cut

sub reverse_translate
{
  my($aaseq, $codonhash) = @_;
  my $newseq = "";
  $newseq .= $codonhash->{$_} foreach (split('', $aaseq));
  return $newseq;
}

=head2 degcodon_to_aas()

  takes a codon that may be degenerate and a codon table and returns a list of
  all amino acids that codon could represent. If a hashref is provided with
  previous answers, it will run MUCH faster.
  in: codon (string),
      codon table (hash reference)
  out: amino acid list (vector)
  
=cut

sub degcodon_to_aas
{
  my ($codon, $CODON_TABLE, $xlationref) = @_;
  return if ( ! $codon  ||  length($codon) != 3);
  my (@answer, %temphash) = ((), ());
  if (exists $xlationref->{$codon})
  {
    return @{$xlationref->{$codon}};
  }
  elsif ($codon eq "NNN")
  {
    %temphash = map { $_ => 1} values %$CODON_TABLE;
    @answer = keys %temphash;
  }
  else
  {
    my $reg = regres($codon, 1);
    %temphash = map {$CODON_TABLE->{$_}  => 1} grep { $_ =~ $reg } keys %$CODON_TABLE;
    @answer = keys %temphash;
  }
  $xlationref->{$codon} = \@answer;
  return @answer;
}

=head2 amb_translation()

  takes a nucleotide that may be degenerate and a codon table and returns a list
  of all amino acid sequences that nucleotide sequence could be translated into.
  in: nucleotide sequence (string),
      codon table (hash reference),
      optional switch to force only a single frame of translation
      optional hashref of previous answers to speed processing
  out: amino acid sequence list (vector)
  
=cut

sub amb_translation
{
  my ($site, $CODON_TABLE, $swit, $xlationref) = @_;
  croak ("Bad frame argument\n") if ($swit ne "1" && $swit ne "2"
    && $swit ne "3"  && $swit ne "-1" && $swit ne "-2" && $swit ne "-3"
    && $swit ne "s" && $swit ne "t");
  if ($swit eq "s" || $swit eq "t")
  {
    my @frames = qw(1 2 3);
    push @frames, qw(-1 -2 -3) if ($swit eq "s");
    my (%RES);
    foreach my $s (@frames)
    {
      my @set = amb_translation($site, $CODON_TABLE, $s, $xlationref);
      $RES{$_}++ foreach(@set);
    }
    return keys %RES;
  }
  else
  {
    my (%RES, @SEED, @NEW);
    $site = complement($site, 1) if ($swit < 0);
    $site = 'N' x (abs($swit) - 1) . $site if (abs($swit) < 4);
    my $gothrough = 0;
    for (my $offset = 0; $offset < (length($site)); $offset +=3)
    {
      my $tempcodon = substr($site, $offset, 3);
      $tempcodon .= "N" while (length($tempcodon) % 3);
      if (!$swit)
      {
        $tempcodon .= 'N' while (length($tempcodon) < 3);
      }
      if ($gothrough == 0)
      {
        @SEED = degcodon_to_aas($tempcodon, $CODON_TABLE, $xlationref) ;
      }
      else
      {
        @NEW  = degcodon_to_aas($tempcodon, $CODON_TABLE, $xlationref);
        @SEED = combine(\@SEED, \@NEW);
      }
      $gothrough++;
    }
    $RES{$_}++ foreach(@SEED);
    return keys %RES;
  }
}

=head2 amb_transcription()

  takes a nucleotide sequence that may have ambiguous bases and returns a list
   of all possible non ambiguous nucleotide sequences it could represent
  in: nucleotide sequence (string),
      codon table (hash reference),
      reverse codon table (hash reference)
  out: nucleotide sequence list (vector)

=cut

sub amb_transcription
{
  # has test in t/02-codons.t
  my ($ntseq, $CODON_TABLE, $pepseq) = @_;
  my (@SEED, @NEW) = ((), ());
  my $offset = 0;
  if ( !$pepseq )
  {
    while ($offset < length($ntseq))
    {
      my $template = substr($ntseq, $offset, 3);
      my $regtemp = regres($template);
      if ($template !~ $ambnt)
      {
        @SEED = ( $template ) if ($offset == 0);
        @NEW  = ( $template ) if ($offset >  0);
      }
      else
      {
        my @TEMP = grep { $_ =~ $regtemp } keys %$CODON_TABLE;
        @SEED = @TEMP if ($offset == 0);
        @NEW  = @TEMP if ($offset >  0);
      }
      unless ($offset == 0)
      {
        @SEED = combine(\@SEED, \@NEW);
      }
      $offset += 3;
    }
  }
  else
  {
    my $REV_CODON = define_reverse_codon_table($CODON_TABLE);
    my @each = split("", $pepseq);
    while ($offset < scalar(@each))
    {
      my $codon = substr($ntseq, $offset*3, 3);
      my $peptide = $each[$offset];
      my @TEMP = grep {$_ =~ regres($codon, 1)} @{$REV_CODON->{$peptide}};
      @SEED = @TEMP if ($offset == 0);
      @NEW  = @TEMP if ($offset >  0);
      unless ($offset == 0)
      {
        @SEED = combine(\@SEED, \@NEW);
      }
      $offset++;
    }
  }
  my %hsh = map {$_ => 1 } @SEED;
  return keys %hsh;
#  return grep {translate($_, 1, $hashref) eq $pepseq} keys %SEED_TOTAL if ($pepseq);
}

=head2 combine

  Basically builds a list of tree nodes for the amb_trans* functions.
  in: 2 x peptide lists (array reference) out: combined list of peptides
  
=cut

sub combine
{
  my ($arr1ref, $arr2ref) = @_;
  my @arr3 = ();
  foreach my $do (@$arr1ref)
  {
    push @arr3, $do . $_ foreach (@$arr2ref)
  }
  return @arr3;
}

=head2 pattern_finder

=cut

sub pattern_finder
{
  # has test in t/02-codons.t
  my ($strand, $pattern, $swit, $frame, $CODON_TABLE) = @_;
  my @positions = ();
  if ($swit == 2)
  {
    return if (! $frame || ! $CODON_TABLE);
    $strand = translate($strand, $frame, $CODON_TABLE)
  }
  my $exp = regres($pattern, $swit);
  while ($strand =~ /(?=$exp)/ig)
  {
    push @positions, (pos $strand);
  }
  return @positions;
}

=head2 pattern_remover()

  takes a nucleotide sequence, a nucleotide "pattern" to be removed, and a few
  codon tables, and returns an edited nucleotide sequence that is missing the
  pattern (if possible).  Ranks codon replacements by RSCU differences to
  minimize expression damage.
  in: nucleotide sequence (string),
      nucleotide pattern (string),
      codon table (hash reference),
      RSCU value table (hash reference)
  out: nucleotide sequence (string) OR null
  
=cut

sub pattern_remover
{
  # has tests in t/02-codons.t
  my ($critseg, $pattern, $CODON_TABLE, $RSCU_TABLE) = @_;
  my @patterns;
  if (ref($pattern) eq "ARRAY")
  {
    @patterns = @{$pattern};
  }
  else
  {
    push @patterns, $pattern;
  }
  my $REV_CODON_TABLE = define_reverse_codon_table($CODON_TABLE);
  my %changes;
  # for each codon position, get RSCU difference possibilities
  for (my $offset = 0; $offset < length($critseg); $offset += 3)
  {
    my $codono = substr($critseg, $offset, 3);
    foreach my $codonp (  grep { $codono ne $_}
                @{ $REV_CODON_TABLE->{ $CODON_TABLE->{$codono} } } )
    {
      my $id = $offset . "." . $codono . "." . $codonp;
      $changes{$id}->{DRSCU} = abs($RSCU_TABLE->{$codonp} - $RSCU_TABLE->{$codono});
      $changes{$id}->{OFFSET} = $offset;
      $changes{$id}->{OLDCOD} = $codono;
      $changes{$id}->{NEWCOD} = $codonp;
    }
  }
  #Replace the least different codons one at a time, take first solution
  foreach my $id (sort {$changes{$a}->{DRSCU} <=> $changes{$b}->{DRSCU}} keys %changes)
  {
    my $copy = $critseg;
    substr($copy, $changes{$id}->{OFFSET}, 3) = $changes{$id}->{NEWCOD};
    next if ($copy =~ $patterns[0] || complement($copy, 1) =~ $patterns[0]);
    next if ($patterns[1] && ($copy =~ $patterns[1] || complement($copy, 1) =~ $patterns[1]));
    return $copy;
  }
  #Try pairwise combinations of codons, sorted by sum of RSCU difference
  my @singles = keys %changes;
  my @pairs;
  for my $x (0..scalar(@singles)-1)
  {
    for my $y ($x+1..scalar(@singles)-1)
    {
      if ($changes{$singles[$x]}->{OFFSET} != $changes{$singles[$y]}->{OFFSET})
      {
        push @pairs, [$singles[$x], $singles[$y]];
      }
    }
  }
  @pairs = sort {($changes{$a->[0]}->{DRSCU} + $changes{$a->[1]}->{DRSCU})
            <=>  ($changes{$b->[0]}->{DRSCU} + $changes{$b->[1]}->{DRSCU})} @pairs;
  foreach my $pair (@pairs)
  {
    my ($one, $two) = @$pair;
    my $copy = $critseg;
    substr($copy, $changes{$one}->{OFFSET}, 3) = $changes{$one}->{NEWCOD};
    substr($copy, $changes{$two}->{OFFSET}, 3) = $changes{$two}->{NEWCOD};
    next if ($copy =~ $patterns[0] || complement($copy, 1) =~ $patterns[0]);
    next if ($patterns[1] && ($copy =~ $patterns[1] || complement($copy, 1) =~ $patterns[1]));
    return $copy;
  }
  return 0;#$critseg;
}

=head2 pattern_adder()

  takes a nucleotide sequence, a nucleotide "pattern" to be interpolated, and
  the codon table, and returns an edited nucleotide sequence that contains the
  pattern (if possible).
  in: nucleotide sequence (string),
      nucleotide pattern (string),
      codon table (hash reference),
      RSCU value table (hash reference)
  out: nucleotide sequence (string) OR null
  
=cut

sub pattern_adder
{
  # has test in t/02-codons.t
  my ($oldpatt, $newpatt, $CODON_TABLE, $xlationref) = @_;
  #assume that critseg and pattern come in as complete codons
  # (i.e., have been run through pattern_aligner)
  my $REV_CODON_TABLE = define_reverse_codon_table($CODON_TABLE);
  my $copy = "";
  for (my $offset = 0; $offset < length($oldpatt); $offset += 3)
  {
    my $curcod = substr($oldpatt, $offset, 3);
    my $curtar = substr($newpatt, $offset, 3);
    foreach my $g (degcodon_to_aas($curcod, $CODON_TABLE, $xlationref))
    {
      $copy .= $curcod =~ regres($curtar)
           ? $curcod
           : first { compareseqs($curtar, $_) } @{$REV_CODON_TABLE->{$g}};
      # print "\t\tpatadd\t($curcod, $curtar)\t$copy<br>\n";
    }
  }
  # print "\t\tpatadd $copy from $oldpatt\n";
  return length($copy) == length($oldpatt)  ?  $copy  :  0;
}

=head2 is_ORF()

  takes a nucleotide sequence and a codon table and determines whether or not
  the nucleotide sequence is a simple ORF - that is, starts with an ATG codon
  and contains no stop codons in the first frame until the very end.
  in: nucleotide sequence (string),
      codon table (hash reference),
  out: 1 if ORF, 0 if not
  
=cut

sub is_ORF
{
  my ($nucseq, $CODON_TABLE) = @_;
  my @war2 = pattern_finder($nucseq, "*", 2, 1, $CODON_TABLE);
  return 0 if (scalar(@war2) > 1);
  return 0 if (scalar(@war2) &&
               ($war2[0] + 1) != length(translate($nucseq, 1, $CODON_TABLE)));
  return 1;
}

=head2 compareseqs

=cut

sub compareseqs
{
  my ($cur, $tar) = @_;
  return 1 if ($tar =~ regres($cur, 1) || $cur =~ regres($tar, 1));
  return 0;
}

=head2 pattern_aligner()

  takes a nucleotide sequence, a pattern, a peptide sequence, and a codon table
  and inserts Ns before the pattern until they align properly. This is so a
  pattern can be inserted out of frame.
  in: nucleotide sequence (string),
      nucleotide pattern (string),
      amino acid sequence (string),
      codon table (hash reference)
  out: nucleotide pattern (string)
  
=cut

sub pattern_aligner
{
  # has tests in t/02-codons.t
  my ($critseg, $pattern, $peptide, $CODON_TABLE, $swit, $xlationref) = @_;
  $swit = 0 if (!$swit);
  my $diff = length($critseg) - length($pattern);
  my ($newpatt, $nstring, $rounds, $offset, $check, $pcheck) = ("", "N" x $diff, 0, 0, "", "");
  #  print "seeking $pattern for $peptide from $critseg...\n";
  while ($check ne $peptide && $rounds <= $diff*2 + 1)
  {
    $newpatt = $rounds <= $diff
      ?  substr($nstring, 0, $rounds) . $pattern
      :  substr($nstring, 0, ($rounds-3)) . complement($pattern, 1);
    $newpatt .=  "N" while (length($newpatt) != length($critseg));
    #  print "\t$newpatt\n";
    my ($noff, $poff) = (0, 0);
    $check = "";
    while ($poff < length($peptide))
    {
      my @possibles = degcodon_to_aas( substr($newpatt, $noff, 3), $CODON_TABLE, $xlationref );
      #   print "\t\t@possibles\n";
      $check .= $_ foreach( grep { substr($peptide, $poff, 1) eq $_ } @possibles);
      $noff += 3;
      $poff ++;
    }
    $pcheck = translate(substr($critseg, $offset, length($peptide) * 3), 1, $CODON_TABLE);
    #      print "\t\t$check, $pcheck, $offset\n";
    $rounds++;
    $offset += 3 if ($rounds % 3 == 0);
#    $check = "" if ( $pcheck !~ $check);
  }
  $newpatt = "0" if ($check ne $peptide);
  # print "\t\tpataln $check, $pcheck, $rounds, $newpatt\n" if ($check ne $peptide);
  return $swit == 1  ?  ($newpatt, $rounds-1)  :  $newpatt;
}

=head2 change_codons()

  takes a nucleotide sequence and a few codon tables and tries to recode the
  nucleotide sequence to one of four algorithms, 0 random, 1 most optimal,
  2 next most optimal, 3 most different, 4 least different.
  in: nucleotide sequence (string),
      codon table (hash reference),
      reverse codon table (hash reference),
      RSCU value table (hash reference),
      algorithm number
      tag: try not to change the first two bases of the first codon in.
  out: nucleotide sequence (string)
  
=cut

sub change_codons
{
  # has tests in t/02-codons.t
  my ($oldseq, $CODON_TABLE, $RSCU_VALUES, $swit, $tag) = @_;
  my $REV_CODON_TABLE = define_reverse_codon_table($CODON_TABLE);
  my ($offset, $newcod, $curcod, $newseq, $aa) = (0, undef, undef, undef, undef);
  $tag = 0 if (! $tag);
  while ($offset < length($oldseq))
  {
    $curcod = substr($oldseq, $offset, 3);
    $newcod = $curcod;
    $aa = $CODON_TABLE->{$curcod};
    my @posarr = sort { $RSCU_VALUES->{$b} <=> $RSCU_VALUES->{$a} }
           grep { exists($RSCU_VALUES->{$_}) }
           @{$REV_CODON_TABLE->{$aa}};
    if (scalar(@posarr) != 1 && $aa ne '*')
    {
      if    ($swit == 0)  #Random
      {
        $newcod = first { $_ ne $curcod } shuffle @posarr;
      }
      elsif ($swit == 1)  #Optimal
      {
        $newcod = first {    1    } @posarr;
      }
      elsif ($swit == 2)  #Less Optimal
      {
        $newcod = first { $_ ne $curcod } @posarr;
      }
      elsif ($swit == 3)  #Most Different
      {
        my $lastbase = substr $curcod, 2, 1;
        my $frstbase = substr $curcod, 0, 1;
        if ($tag && $offset == 0)
        {
          $newcod = first {  substr($_, 0, 2) eq substr($curcod, 0, 2)
                  &&  substr($_, 2, 1) =~ sitver($lastbase, 1)}
                @posarr;
          if (!$newcod)
          {
            $newcod = first {  substr($_, 0, 2) eq substr($curcod, 0, 2)
                    &&  substr($_, 2, 1) ne $lastbase}
                  @posarr;
          }
          if (!$newcod)
          {
            die ("can't find a better cod for $curcod $aa<br>");
          }
        }
        elsif  (scalar(@posarr) == 2)
        {
          $newcod = first { $_ ne $curcod } @posarr;
        }
        elsif  (scalar(@posarr) <= 4)
        {
          $newcod = first { (substr $_, 2, 1) =~ sitver($lastbase, 1) } @posarr;
          if (!$newcod)
          {
            $newcod = first { (substr $_, 2) ne $lastbase } @posarr;
          }
        }
        elsif  (scalar(@posarr) > 4)
        {
          $newcod = first {  ((substr $_, 2) !~ sitver($lastbase, 0) )
                  && ((substr $_, 0, 1) ne $frstbase)}
                @posarr;
          if (!$newcod)
          {
            $newcod = first {  ((substr $_, 2) ne $lastbase )
                    && ((substr $_, 0, 1) ne $frstbase)}
                  @posarr;
          }
        }
      #  print "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$tag, $curcod, $newcod !<br><Br>";
      }
      elsif ($swit == 4)  #Least Different
      {
        my @sorarr = sort {abs($RSCU_VALUES->{$a} - $RSCU_VALUES->{$curcod}) <=> abs($RSCU_VALUES->{$b} - $RSCU_VALUES->{$curcod})} @posarr;
        $newcod = first { $_ ne $curcod } @sorarr;
        $newcod = $curcod if (abs($RSCU_VALUES->{$newcod} - $RSCU_VALUES->{$curcod}) > 1);
      }
    }
    $newseq .= $newcod;
    $offset += 3;
  }
  return $newseq;
}

=head2 sitver

  takes a base as a string and returns by request either a set of transitions
  or a set of transversions, used by change_codons()
  
=cut

sub sitver
{
  my ($base, $swit) = @_;
  $swit = 0 if (!$swit);
  if ($swit == 1) #return set of transversions
  {
    return $base =~ $NTIDES{Y}  ?  $NTIDES{R}  :  $NTIDES{Y};
  }
  else      #return set of transitions
  {
    return $base =~ $NTIDES{Y}  ?  $NTIDES{Y}  :  $NTIDES{R};
  }
}

=head2 RSCU_filter()

  Deletes anything from the RSCU_TABLE that doesn't meet minimum RSCU.
  in: rscu table (hash reference),
      minimum RSCU value (integer)
  out: rscu table (hash reference)
  #NO UNIT TEST
  
=cut

sub RSCU_filter
{
  my ($RSCU_TABLE, $min_value) = @_;
  my @bad_codons = grep { $RSCU_TABLE->{$_}  < $min_value }
           keys %$RSCU_TABLE;
  delete @$RSCU_TABLE{@bad_codons};
  return $RSCU_TABLE;
}

=head2 define_codon_percentages()

  Generates a hash.  KEYS: codons (string)
  VALUES: RSCU value over codon family size (float)
  in: codon table (hash reference),
      RSCU value table (hash reference)
  out: codon percentage table (hash)
  
=cut

sub define_codon_percentages
{
  # has a test in t/02-codons.t
  my ($CODON_TABLE, $RSCU_VALUES) = @_;
  my %AA_cod_count;
  $AA_cod_count{$CODON_TABLE->{$_}}++  foreach keys %$CODON_TABLE;
  my %CODON_PERC_TABLE =
    map { $_ => $RSCU_VALUES->{$_} / $AA_cod_count{$CODON_TABLE->{$_}} }
    keys %$CODON_TABLE;
  return \%CODON_PERC_TABLE;
}

=head2 index_codon_percentages()

  Generates two arrays for x and y values of a graph of codon percentage values.
  in: dna sequence (string),
      window size (integer),
      codon percentage table (hash reference)
  out: x values (array reference), y values (array reference)
  #NO UNIT TEST
  
=cut

sub index_codon_percentages
{
  my ($ntseq, $window, $cpthashref) = @_;
  my @xvalues; my @yvalues;
  my %CODON_PERCENTAGE_TABLE = %$cpthashref;
  my $index; my $sum;
  for (my $x = int($window * (3 / 2)) - 3;
          $x < (length($ntseq) - 3 * (int($window * (3 / 2)) - 3));
          $x += 3)
  {
    $sum = 0;
    for(my $y = $x; $y < 3*$window + $x; $y += 3)
    {
      $sum += $CODON_PERCENTAGE_TABLE{substr($ntseq, $y, 3)};
  #    $sum += $RSCU_TABLE{substr($nucseq, $y, 3)};
    }
    $sum = $sum / $window;
    $index = ($x / 3) + 1;
    push @xvalues, $index;
    push @yvalues, $sum;
  }
  return (\@xvalues, \@yvalues);
}

=head2 codon_count()

  takes a reference to an array of sequences and returns a hash with codons as
  keys and the number of times the codon occurs as a value.
  in: gene sequences (array reference)
  out: codon count (hash reference)
  
=cut

sub codon_count
{
  # has a test in t/02-codons.t
  my ($seq, $CODON_TABLE, $hashref) = @_;
  my %blank = map {$_ => 0} keys %$CODON_TABLE;
  my $codoncount = $hashref ? $hashref : \%blank;
  my $offset = 0;
  while ( $offset <= length($seq) - 3 )
  {
    my $codon = substr($seq, $offset, 3);
    if ($codon =~ $strcodon)
    {
      $codoncount->{$codon} ++;
    }
    else
    {
      $codoncount->{"XXX"} ++;
    }
    $offset += 3;
  }
  return $codoncount;
}

=head2 generate_RSCU_values()

  takes a hash reference with keys as codons and values as number of times
  those codons occur (it helps to use codon_count) and returns a hash with each
  codon and its RSCU value
  in: codon count (hash reference),
      reverse codon table (hash reference)
  out: RSCU values (hash reference)
  
=cut

sub generate_RSCU_values
{
  # has a test in t/02-codons.t
  my ($codon_count, $CODON_TABLE) = @_;
  my $REV_CODON_TABLE = define_reverse_codon_table($CODON_TABLE);
  my $RSCU_hash = {};
  foreach my $codon (sort grep {$_ ne "XXX"} keys %$codon_count)
  {
    my $x_j = 0;
    my $x = $codon_count->{$codon};
    my $family = $REV_CODON_TABLE->{$CODON_TABLE->{$codon}};
    my $family_size = scalar(@$family);
    $x_j += $codon_count->{$_} foreach (grep {exists $codon_count->{$_}} @$family);
    my $rscu = $x / ($x_j / $family_size) || 0.00;
    $RSCU_hash->{$codon} = sprintf("%.2f",  $rscu ) ;
  }
  return $RSCU_hash;
}

=head2 define_aa_defaults()

  Generates a hash.  KEYS: one letter amino acid code (string)
  VALUES: most highly expressed codon for that amino acid (string)
  in: reverse codon table (hash reference),
      RSCU value table (hash reference)
  out: amino acid default table (hash)
  
=cut

sub define_aa_defaults
{
  # has test in t/02-codons.t
  my ($CODON_TABLE, $RSCU_VALUES) = @_;
  my $REV_CODON_TABLE = define_reverse_codon_table($CODON_TABLE);
  my %aa_defaults = ();
  foreach my $aa (keys %AACIDS)
  {
    my $myrscu = -1;
    foreach my $codon (@{$REV_CODON_TABLE->{$aa}})
    {
      if ($RSCU_VALUES->{$codon} > $myrscu)
      {
        $aa_defaults{$aa} = $codon;
        $myrscu = $RSCU_VALUES->{$codon};
      }
    }
  }
  return \%aa_defaults;
}

=head2 orf_finder()

=cut

sub orf_finder
{
  my ($strand, $CODON_TABLE) = @_;
  my $answer = [];
  for my $frame (qw(1 2 3 -1 -2 -3))
  {
    my $strandaa = translate($strand, $frame, $CODON_TABLE);
    my $leng = length($strandaa);
    my $curpos = 0;
    my $orflength = 0;
    my $onnaorf = 0;
    my $orfstart = 0;
    while ($curpos <= $leng)
    {
      my $aa = substr($strandaa, $curpos, 1);
      if ($aa eq 'M' && $onnaorf eq '0')
      {
        $onnaorf = 1;
        $orfstart = $curpos;
      }
      if ($aa eq '*' || ($curpos == $leng && $onnaorf == 1))
      {
        $onnaorf= 0;
        push @$answer, [$frame, $orfstart, $orflength] if ($orflength >= .1*($leng));
        $orflength = 0;
      }
      $curpos++;
      $orflength++ if ($onnaorf == 1);
    }
  }
  return $answer;
}

=head2 minimize_local_alignment_dp()

  Repeatsmasher, by Dongwon Lee. A function that minimizes local alignment
  scores.
  in: gene sequence (string)
      codon table (hashref)
      RSCU table (hashref)
  out: new gene sequence (string)
  #NO UNIT TEST
  
=cut

sub minimize_local_alignment_dp
{
  my ($oldseq, $CODON_TABLE, $RSCU_VALUES) = @_;
  my $match = 5;
  my $transi = -3;
  my $transv = -4;
  my $score_threshold = $match*6; #count the scores only consecutive 6 nts
  my @s = ( [$match, $transv, $transi, $transv],
            [$transv, $match, $transv, $transi],
            [$transi, $transv, $match, $transv],
            [$transv, $transi, $transv, $match] );
  my %nt2idx = ("A"=>0, "C"=>1, "G"=>2, "T"=>3);

  my $REV_CODON_TABLE = define_reverse_codon_table($CODON_TABLE);

  #initial values
  my @optM = (0);
  my $optseq = "";

  #assumming that the sequence is in frame
  my ($offset, $cod, $aa) = (0, "", "");
  while ( $offset <= length($oldseq)-3 )
  {
    $cod = substr($oldseq, $offset, 3);
    $aa  = translate($cod, 1, $CODON_TABLE);
    my @posarr = sort { $RSCU_VALUES->{$b} <=> $RSCU_VALUES->{$a} }
                  @{$REV_CODON_TABLE->{$aa}};

    my @minM = ();
    my $min_seq = "";
    #assign an impossible large score
    my $min_score = $match*(length($oldseq)**2);

    if ($aa ne '*')
    {
      foreach my $newcod (@posarr)
      {
        my @prevM = @optM;
        my $prevseq = $optseq;

        foreach my $nt (split(//, $newcod))
        {
          my $currseq = $prevseq . $nt;
          my @currM = ();
          my $pos = 0;
          push @currM, 0;
          while($pos < length($currseq))
          {
            my $nt2 = substr($currseq, $pos, 1);
            my $nidx1 = $nt2idx{$nt};
            my $nidx2 = $nt2idx{$nt2};
            push @currM, max(0, $prevM[$pos]+$s[$nidx1][$nidx2]);
            $pos++;
          }
          @prevM = @currM;
          $prevseq = $currseq;
        }
        my $scoresum = 0;
        foreach my $i (@prevM)
        {
          if ($i >= $score_threshold)
          {
            $scoresum += $i;
          }
        }
        if ($min_score > $scoresum)
        {
          $min_score = $scoresum;
          $min_seq = $prevseq;
          @minM = @prevM;
        }
      }
    }
    else
    {
      $optseq = $optseq . $cod;
      last;
    }
    @optM = @minM;
    $optseq = $min_seq;
    $offset+=3;
  }
  return $optseq;
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
