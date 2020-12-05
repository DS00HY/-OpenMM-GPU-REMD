use strict;


my $infile = shift;
my $outfile = shift;
my $start_t = shift;
my $dt = shift;
my $start_step = shift;

my $deltat=0.002;
my $ds = int($dt/$deltat);


my $ii = 0;

#my $start_t=150000;
#my $dt= 100 ;
#my $start_step = 5000;
#my $ds = 500 ;

open (IF, "<$infile") || die "open infile $infile error\n";
open (OF, ">$outfile") || die "open outfile $outfile error\n";
while (<IF>)
 {
  my $line =$_;
  if ($line=~m/t=/)
    {
      my $tt =sprintf ("%12.5f",  $start_t    + $dt * $ii) ;
      my $ss =sprintf ("%d", $start_step + $ds * $ii) ;
      printf "t= $tt step= $ss\n";
      $line=~s/t= .+ step= ([0-9]+)/t= $tt step= $ss/g;
      #$line=~s/step= ([0-9]+)/step= $ss/g;
      $ii++;
    }
  print OF $line;
 }
close(IF);
close(OF);
