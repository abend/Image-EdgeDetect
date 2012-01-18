

use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME         => 'Image::EdgeDetect',
    VERSION_FROM => 'lib/Image/EdgeDetect.pm', # finds \$VERSION
    AUTHOR       => 'Sasha Kovar (sasha-cpan@arcocene.org)',
    ABSTRACT     => 'The Canny edge detection algorithm.',
    PREREQ_PM    => {
                     'Test::Simple' => 0.44,
                    },
);
