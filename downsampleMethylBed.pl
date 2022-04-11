#!/usr/bin/env perl
# Format: chr1	1141933	1141935	1	16

use strict;
use List::Util qw(shuffle);
use Getopt::Long;
use File::Basename qw/ basename /;

$::USAGE = "downsampleMethylBed.pl --coverageLevels 1E-5,2E-5,5E-5 --fracTotalFieldsFrom0 3,4 -ncpgs_genome 29175448 ".
  "fn1.hg38.bedgraph fn2.hg38.bedgraph ...\n".
  "(IF YOU WANT THE FRAC FIELD INTERPRETED AS AN INTEGER NUMBER METHYLATED, SPECIFY -3 instead of 3)\n";

my $DEFAULT_METHYL_FRAC_FIELDNUM_FROM_0 = 3;
my $DEFAULT_METHYL_TOTAL_FIELDNUM_FROM_0 = 4;

my @DEFAULT_FRACS = (0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1);
@DEFAULT_FRACS = (1E-5, 2E-5, 5E-5, 1E-4); #, 2E-4, 5E-4, 1E-3);
@DEFAULT_FRACS = (0.5,1);
my $DEFAULT_NCPGS_GENOME = 29175448; # hg38 from cgcontext 

# Get parameters
my $fracTotalFieldsFrom0 = join(",",$DEFAULT_METHYL_FRAC_FIELDNUM_FROM_0,$DEFAULT_METHYL_TOTAL_FIELDNUM_FROM_0);
my $coverageLevels = join(",",@DEFAULT_FRACS);
my $ncpgs_genome = $DEFAULT_NCPGS_GENOME;
GetOptions ('fracTotalFieldsFrom0=s', \$fracTotalFieldsFrom0, 'coverageLevels=s', \$coverageLevels, 'ncpgsGenome=i',\$ncpgs_genome) || die "$::USAGE\n";

# Input params
die "$::USAGE\n" unless (@ARGV>=1);
my (@fns) = @ARGV;


foreach my $fn (@fns)
  {

    my @coverages = split(",",$coverageLevels);
    die "Incorrect format for coverageLevels: $coverageLevels\n"
      unless (scalar(grep {/^[0-9\.E\-]+$/} @coverages)==scalar(@coverages));

    
    # Read in orig
    my $allreadcpgs = readFn($fn,$fracTotalFieldsFrom0);
    # Get coverage stats
    my $fulln = scalar(@$allreadcpgs);
    my ($cvg_full) = readCpgStats($ncpgs_genome, @$allreadcpgs);
    
    CVG: foreach my $cvg_requested (@coverages)
      {

        # We want the "coverage" to be interpreted as mean read coverage
        my $frac = ($cvg_requested/$cvg_full);
        print STDERR sprintf("%0.2e\t%0.2e\t%0.2e\tfrac requesed, coverage requested, coverage full\n",$frac,$cvg_requested,$cvg_full);

        if ($frac>1.5)
          {
            print STDERR sprintf("Coverage %0.2e is >1.5* actual coverage. Skipping\n");
            next CVG;
          }
        elsif ($frac>=0.95)
          {
            $frac = 1.0;
          }
        else
          {
          }

        
        my $fracn = int($fulln*$frac);
        
        print STDERR "\tShuffling...\n";
        my @readcpgs_shuf = shuffle(@$allreadcpgs);
        print STDERR "\tDone shuffling...\n";

        my @downsampled_readcpgs = @readcpgs_shuf[0..($fracn-1)];
        my ($cvg) = readCpgStats($ncpgs_genome, @downsampled_readcpgs);

        my $outfn = $fn;
        my $fracstr = sprintf("cvg%1.6e_frac%1.6e_%dof%d",$cvg_requested,$frac,$fracn,$fulln);
        if ($outfn =~ /\.?hg\d+.bed[a-z]*$/)
          {
            $outfn =~ s/\.?hg(\d+).bed[a-z]*$/-${fracstr}.hg$1.bedgraph/g;
          }
        elsif ($outfn =~ /.bed[a-z]*$/)
          {
            $outfn =~ s/.bed[a-z]*$/-${fracstr}.bedgraph/g;
          }
        else
          {
            $outfn .= "-${fracstr}.bedgraph";
          }
        print STDERR "\tWriting to ${outfn}\n";
        printReadCpgs($outfn,@downsampled_readcpgs);
      }
  }







# - - - - - - Functions

sub readCpgStats
  {
    my ($ntotal, @readcpgs) = @_;

    my $count_by_cpgkey = {};
    
    my $total = 0;
    foreach my $readcpg (@readcpgs)
      {
        my ($linenum,$chr,$s,$e,$m) = @$readcpg;
        my $key = cpg_to_key($linenum,$chr,$s,$e);

        $count_by_cpgkey->{$key}++;
        $total++;
      }
    my $nkeys=scalar(keys(%$count_by_cpgkey));
#    my $NTOTAL = 453143; # This is total for 450k!
    my $cvg = $total/$ntotal;
    print STDERR sprintf("\tcpgs/total/cvg\t%d\t%d\t%0.2f\n",$nkeys,$total,$cvg);

    return($cvg);
  }
  
sub readFn
  {
    my ($fn,$fracTotalFieldsFrom0) = @_;
    my @all = ();

    die "Incorrect format for --flds: $fracTotalFieldsFrom0\n"
      unless ($fracTotalFieldsFrom0 =~ /^\-?\d+,\d+$/);
    my ($frac_fld_from_0, $total_fld_from_0) = split(",",$fracTotalFieldsFrom0);

    die "Couldn't read $fn\n" unless (open(F,$fn));
    my $nlines = 0;
    LINE: while (my $line = <F>) {
      chomp $line;
      my @f = split(/\t/,$line);

      next LINE if (($f[0]=~/random/)||($f[0]=~/chrUn/));
      
      my @readCpgs = getReadCpgs($nlines,$frac_fld_from_0, $total_fld_from_0,@f);

      print STDERR sprintf("On line %d\n",$nlines) unless ($nlines%50000);

      push(@all, @readCpgs);
      $nlines++;
    }

    print STDERR sprintf("Found %d readcpgs in %d cpgs\n",scalar(@all),$nlines);
    return \@all;
  }

  sub printReadCpgs
    {
      my ($outfn,@readcpgs) = @_;

      die "Can't write to $outfn\n" unless (open(OUTF,">$outfn"));
      
      # Tricky because we have to group ones that are at the same position
      my $counts_by_cpgkey = {};

      foreach my $readcpg (@readcpgs)
        {
          my ($linenum,$chr,$s,$e,$m) = @$readcpg;
          my $key = cpg_to_key($linenum,$chr,$s,$e);

          my $counts = $counts_by_cpgkey->{$key};
          $counts = [0,0] unless ($counts);
          @{$counts}[1]++;
          @{$counts}[0]++ if ($m);
          $counts_by_cpgkey->{$key} = $counts;
        }

      foreach my $key (sort(keys(%$counts_by_cpgkey)))
        {
          my $counts = $counts_by_cpgkey->{$key};
          my ($m,$t) = @$counts;
#          print STDERR "Found key $key $m/$t\n";

          my ($linenum,$chr,$st,$end) = key_to_cpg($key);
          print OUTF join("\t",$chr,$st,$end,sprintf("%0.3f",$m/$t),sprintf("%d",$t))."\n";
        }

      close(OUTF);
    }

    sub cpg_to_key
      {
        my ($linenum,$chr,$st,$end) = @_;

        return join("__",sprintf("%08d",$linenum),$chr,$st,$end); # Up to 99M lines
      }

    sub key_to_cpg
      {
        my ($key) = @_;

        my ($linenum,$chr,$st,$end) = split(/\_\_/,$key);
        return ($linenum,$chr,$st,$end);
      }


  sub getReadCpgs
    {
      my ($linenum, $frac_fld_from_0, $total_fld_from_0, @f) = @_;

      my $use_meth_count = 0;
      if ($frac_fld_from_0<0)
        {
          $use_meth_count = 1;
          $frac_fld_from_0 *= -1;
        }
      #print STDERR "frac_fld=$frac_fld_from_0\n";
      
      my $chr = $f[0];
      my $s = $f[1];
      my $e = $f[2];
      my $m = $f[$frac_fld_from_0];
      my $t = $f[$total_fld_from_0];

      # Megalodon puts 0-100 rather than 0-1.  Not sure how to detect this generally
      if (!$use_meth_count)
        {
          $m=($m/100) if ($frac_fld_from_0 == 10);
        }

      # Sanity check.
      if ($use_meth_count)
        {
          die sprintf("meth count not an integer.  Did you specify the correct fields in --fracTotalFieldsFrom0?\n%s\n%s\n",
                      join(", ",@f),$::USAGE)
            unless ($m=~/^\d+$/);
        }
      else
        {
          die sprintf("meth not beween 0-1.  Did you specify the correct fields in --fracTotalFieldsFrom0?\n%s\n%s\n",
                      join(", ",@f),$::USAGE)
            unless (($m>=0)&&($m<=1));
        }

      # Transform to counts
      if (!$use_meth_count)
        {
          $m = int($m*$t);
        }
      


      
      my $u = $t-$m;

      my @out = ();
      for (my $i=0; $i<$m; $i++)
        {
          push(@out,[$linenum,$chr,$s,$e,1]);
        }
      for (my $i=0; $i<$u; $i++)
        {
          push(@out,[$linenum,$chr,$s,$e,0]);
        }

      return @out;
    }
