#
# GeneDesign module for sequence segmentation
#

=head1 NAME

Bio::GeneDesign::RestrictionEnzyme

=head1 VERSION

Version 3.05

=head1 DESCRIPTION

GeneDesign object that represents a restriction enzyme

=head1 AUTHOR

Sarah Richardson <notadoctor@jhu.edu>

=cut

package Bio::GeneDesign::RestrictionEnzyme;

use Bio::GeneDesign::Basic qw(complement regres $ambnt);
use Switch;

use strict;

use base qw(Bio::Root::Root);

my $VERSION = 3.05;

my $IIPreg  = qr/   ([A-Z]*)   \^ ([A-Z]*)      /x;
my $IIAreg  = qr/\A \w+ \(([\-]*\d+) \/ ([\-]*\d+)\)\Z  /x;
my $IIBreg  = qr/\A\(([\-]*\d+) \/ ([\-]*\d+)\) \w+ \(([\-]*\d+) \/ ([\-]*\d+)\)\Z  /x;

=head1 CONSTRUCTOR METHODS

=head2 new

 Title   : new
 Function:
 Returns : a new Bio::GeneDesign::RestrictionEnzyme object
 Args    :

=cut

sub new
{
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($object, $id, $cutseq, $recseq, $temp, $tempin, $score, $methdam, 
      $methdcm, $methcpg, $vendors, $staract, $buffers, $start) =
     $self->_rearrange([qw(ENZYME
                           ID
                           CUTSEQ
                           RECSEQ
                           TEMP
                           TEMPIN
                           SCORE
                           METHDAM
                           METHDCM
                           METHCPG
                           VENDORS
                           STARACT
                           BUFFERS
                           START)], @args);

  if ($object)
  {
    $self->throw("object of class " . ref($object) . " does not implement ".
		    "Bio::GeneDesign::RestrictionEnzyme.")
		  unless $object->isa("Bio::GeneDesign::RestrictionEnzyme");
		$self = $object->clone();
  }
  else
  {
    $self->throw("No name or id defined") unless ($id);
    $self->{'id'} = $id;

    $self->throw("No recognition sequence defined") unless ($recseq);
    $self->{'recseq'} = $recseq;

    $self->throw("No cut sequence defined") unless ($cutseq);
    $self->{'cutseq'} = $cutseq;

    #Regular expression arrayref to use for enzyme searching
    #Should store as compiled regexes instead
    my $qescer = complement($recseq, 1);
    my $arr = [regres($recseq, 1)];
    push @$arr, regres($qescer, 1) if ($recseq ne $qescer);
    $self->regex($arr);

    my $sitelen = length($recseq);
    $self->{'length'} = $sitelen;

    #Enzyme Class and Palindromy

    my ($lef, $rig) = ("", "");
    if ($cutseq =~ $IIPreg)
    {
      $lef = length($1);
      $rig = length($2);
      $self->class("IIP");

      my $inlef = $lef;
      $inlef = length($recseq) - $inlef if ($inlef > (.5 * length($recseq)));
      my $mattersbit = substr($recseq, $inlef, length($recseq) - (2 * $inlef));
      if ($mattersbit =~ $ambnt && length($mattersbit) % 2 == 0)
      {
        $self->palindromy("pnon");
      }
      elsif ($mattersbit eq complement($mattersbit, 1))
      {
        $self->palindromy("pal");
      }
      else
      {
        $self->palindromy("nonpal");
      }
    }
    elsif ($cutseq =~ $IIBreg)
    {
      $lef = int($1);
      $rig = int($2);
      $self->class("IIB");
      $self->palindromy("pnon");
    }
    elsif ($cutseq =~ $IIAreg)
    {
      $lef = int($1);
      $rig = int($2);
      $self->class("IIA");
      $self->palindromy("pnon");
    }
    else
    {
      $self->class("unknown");
    }

    #Enzyme type
    my $type;
    if ($lef < $rig)
    {
      $type .= "5'";
    }
    elsif ($lef > $rig)
    {
      $type .= "3'";
    }
    elsif ($lef == $rig)
    {
      $type .= "b";
    }
    $self->onebpoverhang(1) if (abs($lef - $rig) == 1);
    $self->type($type);

    $temp && $self->temp($temp);
    $tempin && $self->tempin($tempin);

    $score && $self->score($score);

    $staract && $self->staract($staract);

    switch ($methdam)
    {
      case "b"    {$self->methdam("blocked")}
      case "i"    {$self->methdam("inhibited")}
      case "u"    {$self->methdam("unknown")}
      else        {$self->methdam("indifferent")} 
    }
    switch ($methdcm)
    {
      case "b"    {$self->methdcm("blocked")}
      case "i"    {$self->methdcm("inhibited")}
      case "u"    {$self->methdcm("unknown")}
      else        {$self->methdcm("indifferent")} 
    }
    switch ($methcpg)
    {
      case "b"    {$self->methcpg("blocked")}
      case "i"    {$self->methcpg("inhibited")}
      case "u"    {$self->methcpg("unknown")}
      else        {$self->methcpg("indifferent")}
    }  

    $vendors && $self->vendors($vendors);
    $buffers && $self->buffers($buffers);
  }  
  
  $start && $self->start($start);
  
  return $self;
}

=head1 USEFUL METHODS

=head2 clone

=cut

sub clone
{
   my ($self) = @_;
   my $copy;
   foreach my $key (keys %$self)
   {
     if (ref $self->{$key} eq "ARRAY")
     {
       $copy->{$key} = [@{$self->{$key}}];
     }
     elsif (ref $self->{$key} eq "HASH")
     {
       $copy->{$key} = {%{$self->{$key}}};
     }
     else
     {
      $copy->{$key} = $self->{$key};
     }
   }
   bless $copy, ref $self;
   return $copy;
}

=head2 siteseeker()

  Generates a hash describing the positions of a particular enzyme's recognition
  sites in a nucleotide sequence.
  in: nucleotide sequence as a string
  out: reference to a hash where the keys are positions and the value is the
      recognition site from that position as a string.

=cut

sub siteseeker
{
  my ($self, $seq) = @_;
  my $total = {};
  foreach my $sit (@{$self->regex})
  {
    while ($seq =~ /(?=$sit)/ig)
    {
      $total->{pos $seq} = substr($seq, pos $seq, $self->length);
    }
  }
  return $total;
}

=head2 display

=cut

sub display
{
  my ($self) = @_;
  my $staract = "*" if ($self->staract eq "y");
  my (@blocked, @inhibed) = ((), ());
  push @blocked, "cpg" if ($self->methcpg eq "blocked");
  push @blocked, "dam" if ($self->methdam eq "blocked");
  push @blocked, "dcm" if ($self->methdcm eq "blocked");
  push @inhibed, "cpg" if ($self->methcpg eq "inhibited");
  push @inhibed, "dam" if ($self->methdam eq "inhibited");
  push @inhibed, "dcm" if ($self->methdcm eq "inhibited");
  my $buffstr = undef;
  foreach (sort keys %{$self->buffers})
  {
    $buffstr .= "$_ (" . $self->buffers->{$_} . ") " if ($self->buffers->{$_});    
  }
  my $vendstr = join(", ", values %{$self->vendors});
  my $display = undef;
  my $inact = " (". $self->tempin . ")" if ($self->tempin);
  $display .= $self->id . "\t";
  $display .= $self->cutseq . $staract . "\t";
  $display .= $self->type . "\t";
  $display .= $self->start . "\t" if ($self->start);
  $display .= $self->temp . $inact . "\t";
  $display .= join(", ", @blocked) . "\t";
  $display .= join(", ", @inhibed) . "\t";
  $display .= $self->score . "\t";
  $display .= $buffstr . "\t";
  $display .= $vendstr . "\t";
  return $display;  
}

=head1 ACCESSOR METHODS

=head2 id

=cut

sub id
{
  my ($self) = @_;
  return $self->{'id'};
}

=head2 vendorlist

=cut

sub vendorlist
{
  my ($self, $value) = @_;
  if (defined $value)
  {
	  $self->{'vendorlist'} = $value;
  }
  return $self->{'vendorlist'};
}

=head2 star

=cut

sub star
{
  my ($self, $value) = @_;
  if (defined $value)
  {
	  $self->{'star'} = $value;
  }
  return $self->{'star'};
}

=head2 score

=cut

sub score
{
  my ($self, $value) = @_;
  if (defined $value)
  {
	  $self->{'score'} = $value;
  }
  return $self->{'score'};
}

=head2 length

=cut

sub length
{
  my ($self) = @_;
  return $self->{'length'};
}

=head2 methcpg

=cut

sub methcpg
{
  my ($self, $value) = @_;
  if (defined $value)
  {
	  $self->{'methcpg'} = $value;
  }
  return $self->{'methcpg'};
}

=head2 methdcm

=cut

sub methdcm
{
  my ($self, $value) = @_;
  if (defined $value)
  {
	  $self->{'methdcm'} = $value;
  }
  return $self->{'methdcm'};
}

=head2 methdam

=cut

sub methdam
{
  my ($self, $value) = @_;
  if (defined $value)
  {
	  $self->{'methdam'} = $value;
  }
  return $self->{'methdam'};
}

=head2 buffers

=cut

sub buffers
{
  my ($self, $value) = @_;
  if (defined $value)
  {
	  $self->{'buffers'} = $value;
  }
  return $self->{'buffers'};
}

=head2 vendors

=cut

sub vendors
{
  my ($self, $value) = @_;
  if (defined $value)
  {
	  $self->{'vendors'} = $value;
  }
  return $self->{'vendors'};
}

=head2 tempin

=cut

sub tempin
{
  my ($self, $value) = @_;
  if (defined $value)
  {
	  $self->{'tempin'} = $value;
  }
  return $self->{'tempin'};
}

=head2 temp

=cut

sub temp
{
  my ($self, $value) = @_;
  if (defined $value)
  {
	  $self->{'temp'} = $value;
  }
  return $self->{'temp'};
}

=head2 recseq

=cut

sub recseq
{
  my ($self) = @_;
  return $self->{'recseq'};
}

=head2 cutseq

=cut

sub cutseq
{
  my ($self) = @_;
  return $self->{'cutseq'};
}

=head2 regex

=cut

sub regex
{
  my ($self, $value) = @_;
  if (defined $value)
  {
	  $self->{'regex'} = $value;
  }
  return $self->{'regex'};
}

=head2 class

=cut

sub class
{
  my ($self, $value) = @_;
  if (defined $value)
  {
	  $self->{'class'} = $value;
  }
  return $self->{'class'};
}

=head2 type

=cut

sub type
{
  my ($self, $value) = @_;
  if (defined $value)
  {
	  $self->{'type'} = $value;
  }
  return $self->{'type'};
}

=head2 onebpoverhang

=cut

sub onebpoverhang
{
  my ($self, $value) = @_;
  if (defined $value)
  {
	  $self->{'onebpoverhang'} = $value;
  }
  return $self->{'onebpoverhang'};
}

=head2 exclude

=cut

sub exclude
{
  my ($self, $value) = @_;
  if (defined $value)
  {
	  $self->{'exclude'} = $value;
  }
  return $self->{'exclude'};
}

=head2 palindromy

=cut

sub palindromy
{
  my ($self, $value) = @_;
  if (defined $value)
  {
	  $self->{'palindromy'} = $value;
  }
  return $self->{'palindromy'};
}

=head2 staract

=cut

sub staract
{
  my ($self, $value) = @_;
  if (defined $value)
  {
	  $self->{'staract'} = $value;
  }
    return $self->{'staract'};
}

=head2 start

=cut

sub start
{
  my ($self, $value) = @_;
  if (defined $value)
  {
	  $self->{'start'} = $value;
  }
    return $self->{'start'};
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
