#!/usr/bin/perl

require "orig_man/conf.pl";

my $man_name=$ARGV[0];

open(my $f,"<orig_man/$man_name");

my $man_content="";

while(<$f>)
{
    $man_content.=$_;
}
close($f);

foreach my $key (keys(%config))
{
    $man_content =~ s/[^A-Z_]\K$key(?=[^A-Z_])/$config{$key}/g;
}

open(my $g,">man/$man_name");
print $g $man_content;
close($g);

