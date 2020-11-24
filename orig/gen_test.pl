#!/usr/bin/perl

our $example=0;
require "orig/conf.pl";

my $code="require('blockmodels')\nset.seed(12)\n";

my $modelname=$ARGV[0];
open(my $c,"<orig/code/$modelname.R");
while(<$c>)
{
    $code.=$_;
}
close($c);

foreach my $key (keys(%config))
{
    $code =~ s/[^A-Z_]\K$key(?=[^A-Z_])/$config{$key}/g;
}

open(my $g,">tests/$modelname.R");
print $g $code;
close($g);

