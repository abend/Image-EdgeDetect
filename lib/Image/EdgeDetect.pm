package Image::EdgeDetect;

use strict;
use warnings;
use diagnostics;

use File::Path;
use List::AllUtils qw(:all);
use Image::Magick;
use Data::Dumper;
use Math::Trig;
use Math::Complex;
use POSIX qw(ceil floor);

BEGIN {
    use Exporter ();
    use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
    $VERSION     = '0.9.1';
    @ISA         = qw(Exporter);
    @EXPORT      = qw();
    @EXPORT_OK   = qw();
    %EXPORT_TAGS = ();
}


=head1 NAME

Image::EdgeDetect - An implementation of the Canny edge detection algorithm.

=head1 SYNOPSIS

  use Image::EdgeDetect;
  my $detector = Image::EdgeDetect->new();
  $detector->process($in, $out);

=head1 DESCRIPTION

Perform Canny edge detection on an input image, and write out the
resulting binary edge map image.

=head1 USAGE

  use Image::EdgeDetect;
  my $detector = Image::EdgeDetect->new();
  $detector->process($in, $out);

=head1 BUGS

You tell me.

=head1 AUTHOR

    Sasha Kovar
    CPAN ID: ABEND
    sasha-cpan@arcocene.org

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

Image::CornerDetect

=cut

our $GAUSSIAN_CUT_OFF = 0.005;
our $MAGNITUDE_SCALE = 10;
our $MAGNITUDE_LIMIT = 255;
# our $MAGNITUDE_SCALE = 1;
# our $MAGNITUDE_LIMIT = 1;
#our $MAGNITUDE_MAX = 256;# int($MAGNITUDE_SCALE * $MAGNITUDE_LIMIT);


=head2 my $detector = Image::EdgeDetect->new(\%args)

Create a new edge detector, passing in parameters to override the
defaults if desired.  Arguments are:

  low_threshold (default 2.5)
  high_threshold (default 7.5)
  kernel_radius (default 2.0)
  kernel_width (default 16)

=cut

sub new {
  my ($this, $args) = @_;
  my $class = ref($this) || $this;

  my %args = ref($args) eq 'HASH' ? %$args : ();

  my $self =
  {
   lowThreshold         => $args{low_threshold} || 2.5,
   highThreshold        => $args{high_threshold} || 7.5,
   gaussianKernelRadius => $args{kernel_radius} || 2.0,
   gaussianKernelWidth  => $args{kernel_width} || 16,

   height => undef,
   width => undef,
   picsize => undef,
   data => undef,
   magnitude => undef,
   sourceImage => undef,
  };

  bless $self, $class;

  return $self;
}

=head2 $out_img = $detector->process($in, [$out])

Perform edge detection on the input image, writing the edge map to the
output file.  Input image can be a filename or Image::Magick image.
Output is an optional filename to write the image to.  A
Image::Magick image object is returned.

=cut

sub process {
  my ($self, $in, $outFilename) = @_;

  my $image = $in;

  unless (ref($in) eq "Image::Magick") {
    $image = Image::Magick->new();
    my $status = $image->Read($in);
    die "read: $status\n" if $status;
  }

  $self->{sourceImage} = $image;

  $self->{width} = $image->Get('width');
  $self->{height} = $image->Get('height');
  $self->{picsize} = $self->{width} * $self->{height};

  $self->readLuminance();
  $self->computeGradients();

  my $low = int($self->{lowThreshold});
  my $high = int($self->{highThreshold});
  $self->performHysteresis($low, $high);
  $self->thresholdEdges();

  return $self->writeEdges($outFilename);
}

sub readLuminance {
  my ($self) = @_;

  @{$self->{data}} = map { $_ * 255 } $self->{sourceImage}->GetPixels(map=>'I',
                                                     height=>$self->{height},
                                                     width=>$self->{width},
                                                     normalize=>1);
  #say "data:\n",Dumper(@{$self->{data}});
}

sub normalize {
  my @data = @_;

  my $min = min @data;
  my $max = max @data;
  my $diff = $max - $min;

  if ($diff != 0) {
    my $scale = 256.0 / $diff;

  #say "min $min, max $max, scale $scale";

    for (my $i = 0; $i < @data; $i++) {
      my $before = $data[$i];
      $data[$i] = int(($data[$i] - $min) * $scale);
      #say "$before => ".$data[$i] if $before;
    }
  }

  return @data;
}

sub writeEdges {
  my ($self, $outfile) = @_;

  writeImage($self->{data},
             $self->{width}, $self->{height},
             $outfile);
}

sub writeImage {
  my ($data, $w, $h, $file) = @_;

  my $out = Image::Magick->new;
  $out->Set(size => $w.'x'.$h);
  $out->ReadImage('xc:black');

  for (my $y = 0; $y < $h; $y++) {
    for (my $x = 0; $x < $w; $x++) {
      my $offset = $x + $w * $y;
      my $val = $$data[$offset];
      #say "val $val at $x,$y ($offset)" if $val != 0;
      $val /= 256.0;
      $out->SetPixel(x => $x, y => $y, color => [$val, $val, $val]);
    }
  }

  $out->Write($file) if $file;

  return $out;
}

# NOTE: The elements of the method below (specifically the technique for
# non-maximal suppression and the technique for gradient computation)
# are derived from an implementation posted in the following forum (with the
# clear intent of others using the code):
#   http:# forum.java.sun.com/thread.jspa?threadID=546211&start=45&tstart=0
# My code effectively mimics the algorithm exhibited above.
# Since I don't know the providence of the code that was posted it is a
# possibility (though I think a very remote one) that this code violates
# someone's intellectual property rights. If this concerns you feel free to
# contact me for an alternative, though less efficient, implementation.

sub computeGradients {
  my ($self) = @_;

  my $kernelRadius = $self->{gaussianKernelRadius};
  my $kernelWidth = $self->{gaussianKernelWidth};

  my $width = $self->{width};
  my $height = $self->{height};

  # generate the gaussian convolution masks
  my @kernel;
  my @diffKernel;
  my $kwidth;
  for ($kwidth = 0; $kwidth < $kernelWidth; $kwidth++) {
    my $g1 = gaussian($kwidth, $kernelRadius);
    last if $g1 <= $GAUSSIAN_CUT_OFF && $kwidth >= 2;
    my $g2 = gaussian($kwidth - 0.5, $kernelRadius);
    my $g3 = gaussian($kwidth + 0.5, $kernelRadius);
    $kernel[$kwidth] = ($g1 + $g2 + $g3) / 3 / (2 * pi * $kernelRadius * $kernelRadius);
    $diffKernel[$kwidth] = $g3 - $g2;
  }

  my $initX = $kwidth - 1;
  my $maxX = $width - ($kwidth - 1);
  my $initY = $width * ($kwidth - 1);
  my $maxY = $width * ($height - ($kwidth - 1));

  my @data = @{$self->{data}};
  my $size = scalar(@data);
  my @xConv = (0)x$size;
  my @yConv = (0)x$size;

  # perform convolution in x and y directions
  for (my $x = $initX; $x < $maxX; $x++) {
    for (my $y = $initY; $y < $maxY; $y += $width) {
      my $index = $x + $y;
      my $sumX = $data[$index] * $kernel[0];
      my $sumY = $sumX;
      my $xOffset = 1;
      my $yOffset = $width;

      while ($xOffset < $kwidth) {
        $sumY += $kernel[$xOffset] * ($data[$index - $yOffset] + $data[$index + $yOffset]);
        $sumX += $kernel[$xOffset] * ($data[$index - $xOffset] + $data[$index + $xOffset]);
        $yOffset += $width;
        $xOffset++;
      }

      $xConv[$index] = $sumX;
      $yConv[$index] = $sumY;
    }
  }

  #writeImage([normalize(@xConv)], $width, $height, "xxxx-xconv.png");
  #writeImage([normalize(@yConv)], $width, $height, "xxxx-yconv.png");

  # for (my $i = 0; $i < $size; $i++) {
  #   say "$i: ".$yConv[$i] if $yConv[$i] > 0;
  # }

  my @xGradient = (0)x$size;

  for (my $x = $initX; $x < $maxX; $x++) {
    for (my $y = $initY; $y < $maxY; $y += $width) {
      my $sum = 0.0;
      my $index = $x + $y;
      for (my $i = 1; $i < $kwidth; $i++) {
        $sum += $diffKernel[$i] * ($yConv[$index - $i] - $yConv[$index + $i]);
      }

      #say "xgrad $index: $sum";
      $xGradient[$index] = $sum;
    }
  }

  my @yGradient = (0)x$size;

  for (my $x = $kwidth; $x < $width - $kwidth; $x++) {
    for (my $y = $initY; $y < $maxY; $y += $width) {
      my $sum = 0.0;
      my $index = $x + $y;
      my $yOffset = $width;
      for (my $i = 1; $i < $kwidth; $i++) {
        $sum += $diffKernel[$i] * ($xConv[$index - $yOffset] - $xConv[$index + $yOffset]);
        $yOffset += $width;
      }

      $yGradient[$index] = $sum;
    }
  }

  #writeImage([normalize(@xGradient)], $width, $height, "xxxx-xgrad.png");
  #writeImage([normalize(@yGradient)], $width, $height, "xxxx-ygrad.png");

  my @magnitude = (0)x$size;
  my ($gradmax, $gradmin) = (0, 0);

  $initX = $kwidth;
  $maxX = $width - $kwidth;
  $initY = $width * $kwidth;
  $maxY = $width * ($height - $kwidth);
  for (my $x = $initX; $x < $maxX; $x++) {
    for (my $y = $initY; $y < $maxY; $y += $width) {
      my $index = $x + $y;
      my $indexN = $index - $width;
      my $indexS = $index + $width;
      my $indexW = $index - 1;
      my $indexE = $index + 1;
      my $indexNW = $indexN - 1;
      my $indexNE = $indexN + 1;
      my $indexSW = $indexS - 1;
      my $indexSE = $indexS + 1;

      my $xGrad = $xGradient[$index];
      my $yGrad = $yGradient[$index];
      my $gradMag = hypot($xGrad, $yGrad);

      #say "grads at $x,$y: $xGrad $yGrad ($gradMag)";

      # perform non-maximal supression
      my $nMag  = hypot($xGradient[$indexN],  $yGradient[$indexN]);
      my $sMag  = hypot($xGradient[$indexS],  $yGradient[$indexS]);
      my $wMag  = hypot($xGradient[$indexW],  $yGradient[$indexW]);
      my $eMag  = hypot($xGradient[$indexE],  $yGradient[$indexE]);
      my $neMag = hypot($xGradient[$indexNE], $yGradient[$indexNE]);
      my $seMag = hypot($xGradient[$indexSE], $yGradient[$indexSE]);
      my $swMag = hypot($xGradient[$indexSW], $yGradient[$indexSW]);
      my $nwMag = hypot($xGradient[$indexNW], $yGradient[$indexNW]);
      my $tmp;

      #say "mags at $x,$y: ($xGrad $yGrad - $gradMag) $nMag $sMag $wMag $eMag / $neMag $seMag $swMag $nwMag";

      # An explanation of what's happening here, for those who want
      # to understand the source: This performs the "non-maximal
      # supression" phase of the Canny edge detection in which we
      # need to compare the gradient magnitude to that in the
      # direction of the gradient; only if the value is a local
      # maximum do we consider the point as an edge candidate.
      #
      # We need to break the comparison into a number of different
      # cases depending on the gradient direction so that the
      # appropriate values can be used. To avoid computing the
      # gradient direction, we use two simple comparisons: first we
      # check that the partial derivatives have the same sign (1)
      # and then we check which is larger (2). As a consequence, we
      # have reduced the problem to one of four identical cases that
      # each test the central gradient magnitude against the values at
      # two points with 'identical support'; what this means is that
      # the geometry required to accurately interpolate the magnitude
      # of gradient function at those points has an identical
      # geometry (upto right-angled-rotation/reflection).
      #
      # When comparing the central gradient to the two interpolated
      # values, we avoid performing any divisions by multiplying both
      # sides of each inequality by the greater of the two partial
      # derivatives. The common comparand is stored in a temporary
      # variable (3) and reused in the mirror case (4).
      if ($xGrad * $yGrad <= 0.0 # (1)
          ? abs($xGrad) >= abs($yGrad) # (2)
          ? ($tmp = abs($xGrad * $gradMag)) >= abs($yGrad * $neMag - ($xGrad + $yGrad) * $eMag) # (3)
          && $tmp > abs($yGrad * $swMag - ($xGrad + $yGrad) * $wMag) # (4)
          : ($tmp = abs($yGrad * $gradMag)) >= abs($xGrad * $neMag - ($yGrad + $xGrad) * $nMag) # (3)
          && $tmp > abs($xGrad * $swMag - ($yGrad + $xGrad) * $sMag) # (4)
          : abs($xGrad) >= abs($yGrad) # (2)
          ? ($tmp = abs($xGrad * $gradMag)) >= abs($yGrad * $seMag + ($xGrad - $yGrad) * $eMag) # (3)
          && $tmp > abs($yGrad * $nwMag + ($xGrad - $yGrad) * $wMag) # (4)
          : ($tmp = abs($yGrad * $gradMag)) >= abs($xGrad * $seMag + ($yGrad - $xGrad) * $sMag) # (3)
          && $tmp > abs($xGrad * $nwMag + ($yGrad - $xGrad) * $nMag) # (4)
         )
      {
        $gradmin = min($gradmin, $gradMag);
        $gradmax = max($gradmax, $gradMag);

        my $mag = int($gradMag * $MAGNITUDE_SCALE);
        $magnitude[$index] = min($mag, $MAGNITUDE_LIMIT);
        #say "setting to ".$magnitude[$index];
        #say "setting at $x,$y: ".$magnitude[$index];

        # NOTE: The orientation of the edge is not employed by this
        # implementation. It is a simple matter to compute it at
        # this point as: atan2(yGrad, xGrad);
      } else {
        $magnitude[$index] = 0;
      }
    }
  }

  #say "grad min $gradmin, max $gradmax";
  #say "mag from ".min(@magnitude). " to ".max(@magnitude);

  @{$self->{magnitude}} = @magnitude;

  #writeImage([normalize(@magnitude)], $width, $height, "xxxx-mag.png");

  #say "MAG:\n",Dumper(@magnitude);

}

sub hypot {
  my ($x, $y) = @_;
  return abs($x) + abs($y);
  #return sqrt($x*$x, $y*$y);
}

sub gaussian {
  my ($x, $sigma) = @_;
  return exp(-($x * $x) / (2 * $sigma * $sigma));
}

sub performHysteresis {
  my ($self, $low, $high) = @_;

  my $min = min @{$self->{data}};
  my $max = max @{$self->{data}};
  #say "hyster input range from $min to $max, thresholds $low to $high";
  #say "mag from ".min(@{$self->{magnitude}}). " to ".max(@{$self->{magnitude}});

  # NOTE: this implementation reuses the data array to store both
  # luminance data from the image, and edge intensity from the processing.
  # This is done for memory efficiency, other implementations may wish
  # to separate these functions.
  @{$self->{data}} = map { 0 } @{$self->{data}};

  my $offset = 0;
  for (my $y = 0; $y < $self->{height}; $y++) {
    for (my $x = 0; $x < $self->{width}; $x++) {
      #say $self->magnitude($offset) ." >= ". $high if $self->magnitude($offset);
      if ($self->data($offset) == 0 && $self->magnitude($offset) >= $high) {
        #no warnings 'recursion';
        #say("----- TOP -----");
        $self->follow($x, $y, $offset, $low);
      }
      $offset++;
    }
  }
}

sub follow {
  my ($self, $x1, $y1, $i1, $threshold) = @_;
  #say("$x1,$y1 $i1 - $threshold");

  $self->data($i1, $self->magnitude($i1));
  #say "set data $i1 = ". $self->magnitude($i1);

  my $x0 = $x1 == 0 ? $x1 : $x1 - 1;
  my $x2 = $x1 == $self->{width} - 1 ? $x1 : $x1 + 1;
  my $y0 = $y1 == 0 ? $y1 : $y1 - 1;
  my $y2 = $y1 == $self->{height} -1 ? $y1 : $y1 + 1;

  for (my $x = $x0; $x <= $x2; $x++) {
    for (my $y = $y0; $y <= $y2; $y++) {
      my $i2 = $x + $y * $self->{width};
      if (($y != $y1 || $x != $x1)
          && $self->data($i2) == 0
          && $self->magnitude($i2) >= $threshold) {
        no warnings 'recursion';
        $self->follow($x, $y, $i2, $threshold);
        return;
      }
    }
  }
}

sub thresholdEdges {
  my ($self) = @_;

  my $min = min $self->{data};
  my $max = max $self->{data};
  my $mid = ($min + $max) / 2;
  for (my $i = 0; $i < $self->{picsize}; $i++) {
    #say "data $i: ".$self->data($i) if $self->data($i);
    $self->data($i, $self->data($i) > $mid ? 1 : 0);
  }
}

# accessor/mutators
sub data {
  my ($self, $index, $value) = @_;
  ${$self->{data}}[$index] = $value if $value;
  return ${$self->{data}}[$index];
}

sub magnitude {
  my ($self, $index, $value) = @_;
  ${$self->{magnitude}}[$index] = $value if $value;
  return ${$self->{magnitude}}[$index];
}

1;

