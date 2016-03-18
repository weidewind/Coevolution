#!/usr/bin/perl
use File::Spec;
use Cwd;

my $nullstr;
my $dir = File::Spec -> catdir(getcwd(),$nullstr);
my $file  = File::Spec -> catfile(($dir, "else", "also"), "protein");
print $file;
