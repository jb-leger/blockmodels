#!/usr/bin/perl

our $example=1;
require "orig/conf.pl";

my $modelname=$ARGV[0];

open(my $f,"<orig/man/$modelname.Rd");

my $man_content="";

while(<$f>)
{
    $man_content.=$_;
}
close($f);



$man_content.="\\examples{\\dontrun{\n\n";
open(my $c,"<orig/code/$modelname.R");
while(<$c>)
{
    s/%/\\%/g;
    $man_content.=$_;
}
close($c);
$man_content.="}}\n";

foreach my $key (keys(%config))
{
    $man_content =~ s/[^A-Z_]\K$key(?=[^A-Z_])/$config{$key}/g;
}

open(my $g,">man/$modelname.Rd");
print $g $man_content;
close($g);

