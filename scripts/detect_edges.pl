#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use 5.010;

use File::Basename;
use Image::EdgeDetect;

use Carp qw(confess);
local $SIG{__DIE__} = \&Carp::confess;
local $SIG{__WARN__} = \&Carp::confess;


my $detector = Image::EdgeDetect->new({ kernel_radius => 1.5,
                                        kernel_width => 6,
                                        low_threshold => 32,
                                        high_threshold => 128,
                                      });

my $infile = $ARGV[0];
die "usage: $0 filename" unless $infile;
my ($name, $path, $suffix) = fileparse($infile, qr/\.[^.]*/);

$detector->process($infile, "$path/${name}-out$suffix");
