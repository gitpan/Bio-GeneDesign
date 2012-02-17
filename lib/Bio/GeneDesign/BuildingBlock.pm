#
# GeneDesign module for sequence segmentation
#

=head1 NAME

Bio::GeneDesign::BuildingBlock

=head1 VERSION

Version 3.05

=head1 DESCRIPTION

GeneDesign object that represents a segment of designed DNA

=head1 AUTHOR

Sarah Richardson <notadoctor@jhu.edu>.

=cut

package Bio::GeneDesign::BuildingBlock;

use strict;

use base qw(Bio::Root::Root);

my $VERSION = 3.05;

=head2 new

 Title   : new                                     
 Function:
 Returns : a new Bio::GeneDesign::BuildingBlock object
 Args    : 
       
=cut

sub new
{
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($seq, $start, $end, $id, $description, $enzyme5, $enzyme3, $comment, 
      $userlist) =
     $self->_rearrange([qw(SEQ
                           START
                           END
                           ID
                           DESCRIPTION
                           ENZYME5
                           ENZYME3
                           COMMENT
                           USERLIST)], @args); 
                           
  $self->throw("No sequence provided") unless ($seq);
  $self->seq($seq);
  
  $self->length(length($seq));
  
  unless (!($start && $end) || ($start && $end))
  {
    $self->throw("Must provide start and end");
  }
  $start && $self->start($start);
  $end && $self->end($end);
  
  $id && $self->id($id);
  
  unless ($description)
  {
    $description  = $self->length() . "bp ";
    $description .= "from $start to $end " if ($start && $end);
    $description .= "enzyme $enzyme5 " if ($enzyme5);
    $description .= "to enzyme $enzyme3 " if ($enzyme3);
    $description .= $comment if ($comment);
  }
  $description && $self->description($description);

  $enzyme5 && $self->enzyme5($enzyme5);
  $enzyme3 && $self->enzyme3($enzyme3);
  
  $comment && $self->comment($comment);
  
  $userlist && $self->userlist($userlist);
  
  return $self;
}


=head2 seq

=cut

sub seq
{
  my ($obj, $value) = @_;
  if (defined $value)
  {
	  $obj->{'seq'} = $value;
  }
  return $obj->{'seq'};
}

=head2 id

=cut

sub id
{
  my ($obj, $value) = @_;
  if (defined $value)
  {
	  $obj->{'id'} = $value;
  }
  return $obj->{'id'};
}

=head2 start

=cut

sub start
{
  my ($obj, $value) = @_;
  if (defined $value)
  {
	  $obj->{'start'} = $value;
  }
  return $obj->{'start'};
}

=head2 end

=cut

sub end
{
  my ($obj, $value) = @_;
  if (defined $value)
  {
	  $obj->{'end'} = $value;
  }
  return $obj->{'end'};
}

=head2 length

=cut

sub length
{
  my ($obj, $value) = @_;
  if (defined $value)
  {
	  $obj->{'length'} = $value;
  }
  return $obj->{'length'};
}

=head2 description

=cut

sub description
{
  my ($obj, $value) = @_;
  if (defined $value)
  {
	  $obj->{'description'} = $value;
  }
  return $obj->{'description'};
}

=head2 enzyme3

=cut

sub enzyme3
{
  my ($obj, $value) = @_;
  if (defined $value)
  {
	  $obj->{'enzyme3'} = $value;
  }
  return $obj->{'enzyme3'};
}

=head2 enzyme5

=cut

sub enzyme5
{
  my ($obj, $value) = @_;
  if (defined $value)
  {
	  $obj->{'enzyme5'} = $value;
  }
  return $obj->{'enzyme5'};
}

=head2 comment

=cut

sub comment
{
  my ($obj, $value) = @_;
  if (defined $value)
  {
	  $obj->{'comment'} = $value;
  }
  return $obj->{'comment'};
}

=head2 userlist

=cut

sub userlist
{
  my ($obj, $value) = @_;
  if (defined $value)
  {
	  $obj->{'userlist'} = $value;
  }
  return $obj->{'userlist'};
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
