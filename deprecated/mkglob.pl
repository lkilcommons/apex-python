#!/usr/bin/perl -w

#          Make a file with global interpolation tables for geomagnetic coordinates
#          conversions; usage:
#
#            mkglob [-h] [-a altmx] [-d date1:date2:...] [-r nvert] [file]
#
#          where
#            -a altmx = Maximum required altitude (km).  The actual maximum altitude is
#                       higher; it is the next grid point above 'altmx' above which
#                       interpolation is not possible; i.e. a fatal error.
#            -d dates = One or more dates for which tables are to be created, each
#                       expressed as year and fraction (yyyy.fraction), separated by a
#                       colon and in chronological order.  The earliest may not precede
#                       1900, the first IGRF epoch as defined in magfld.f.
#            -h       = Help option: print this help message and quit.
#            -r nvert = Spatial resolution parameter must be at least 2 but 30 or more is
#                       recommended.  Resolution increases with increasing nvert until
#                       about 100 where accuracy degrades close to the poles in east-west
#                       gradients.  nvert greater than 1000 degrades accuracy at all
#                       latitudes. More details are given in ggrid.f.
#            file     = Name of the interpolation tables file to be created.
#
#          When not specified via command arguments (-a -d -r file), 'altmx' is $altmxd, 'dates'
#          (epochs) are $datds,
#          'nvert' is $nvertd, and 'file' is '$filed'.
#
#          A pre-existing 'file' is saved as 'file.old' before creating the new version;
#          this is reported on stderr as well as the time required to create the new version.
#
#          Once created, 'file' is suitable for conversion from geographic coordinates to
#          Apex, Modified Apex and Quasi-Dipole coordinates and back using entries in
#          apxntrp.f as exemplified in xglob.f.
#
#          Suitable table design must address the adequate maximum grid altitude
#          constraint, above which it is impossible to interpolate, and may also wish to
#          avoid non-fatal 'dates' warnings issued by apxntrp.f when there is more than one
#          year difference in a single epoch 'file', more than five years between epochs,
#          before the first epoch in a multiple epoch 'file' or after the last epoch.
#
#          SEE ALSO:
#          apxntrp.f - interpolation routines                           - $apxntr
#          ggrid.f   - routine to get a grid for interpolation tables   - $ggrid
#          magfld.f  - IGRF definition routines                         - $magfld
#          mkglob.f  - program to create global interpolation tables    - $mkglob
#          xglob     - driver script for xglob.f                        - $xglobs
#          xglob.f   - example program of apxntrp.f and prepared tables - $xglob
#
#          HISTORY:
#          May 2004:  Written for flexible assigment of variables formerly hard-wired in
#                     apex/src/mkglob.f which are now command arguments
#
#          Configuration variables:
#          Environment variable definition is prerequisite: 'source envapex' where envapex
#          is in this directory (apex/bin)
use Cwd;                                        # Current working directory module for cwd()

$altmxd = 1000;                                 # default maximum required altitude
$cdir   = cwd();                                # current directory
$datds  = '1965.0';                             # default date (epoch) list string
for ($y =  1970; $y < (gmtime (time()))[5] + 1904.9999; $y += 5) {
   $datds .= " $y.0";
}
$filed  =  "apxntrp_grid2011";                  # default write file name
$mkge   = "./mkglob";                       # Compiled executable to do the heavy lifting
$nvertd = 40;                                   # default nvert
$|      = 1;                                    # hot pipes
$snmt   = ($0 =~ m#/*([^/]+)$#o) ? $1 : $0;     # script name tail


#          Parse command arguments to get options
@datl   = ();                                    # date (epoch) list
$hopt   = 0;                                     # help option
$optc   = 'hH';                                  # option letters
$i      = -1;
while ($i < $#ARGV) {
   $i++;
   if ($ARGV[$i] =~ /^-([$optc]+)$/) {
      for $char (split //, $1) {
	 $hopt = 1 if $char eq 'h' || $char eq 'H';
      }

   } elsif ($ARGV[$i] =~ /^-a(.*)$/) {
      $altmx = $1;
      if (! $altmx) {
	 $i++;
	 die "$snmt missing -a option value, maximum required altitude\n" if $i > $#ARGV;
	 $altmx = $ARGV[$i];
      }
      die "$snmt: maximum altitude '$altmx' may not be less than zero\n" if $altmx !~ /^\d+\.?\d*$/;

   } elsif ($ARGV[$i] =~ /^-d(.*)$/) {
      $dats = $1;
      if (! $dats) {
	 $i++;
	 die "$snmt missing -d option value, colon delimited date list\n" if $i > $#ARGV;
	 $dats = $ARGV[$i];
      }
      @datl = split ':', $dats;
      for (@datl) {
	 next if /^(\d{4}|\d{4}\.|\d{4}\.\d+)$/;
	 die "$snmt: date ($_) must be formatted as year and fraction: yyyy.fraction\n";
      }

   } elsif ($ARGV[$i] =~ /^-r(.*)$/) {
      $nvert = $1;
      if (! $nvert) {
	 $i++;
	 die "$snmt missing -r option value, nvert\n" if $i > $#ARGV;
	 $nvert = $ARGV[$i];
      }
      die "$snmt: -r option value ($nvert) not a positive integer\n" if $nvert !~ /^\d+$/;
      die "$snmt: -r option value ($nvert) must be at least 2\n"     if $nvert < 2;

   } else {
      $file = $ARGV[$i];
      $file = "$cdir/$file" if $file !~ /^\//;
   }
}
prthelp() if $hopt;

$altmx = $altmxd if ! $altmx;
$file  = $filed  if ! $file;
$nvert = $nvertd if ! $nvert;
$dats  = (@datl) ? join (' ', @datl) : $datds;

#          Rename a pre-existing $file, since they're not easy to create
if (-e $file) {
   $filo = "$file.old";
   warn "$snmt: Existing file '$file' is being saved as '$filo'\n";
   rename ($file, $filo) == 1 or die "$snmt: trouble renaming $file: $!\n";
}


#            Build command string: mkglob.exe nvert altmx apxntrp_file epoch [epoch2 ...]
$altmx = $altmxd if ! $altmx;
$nvert = $nvertd if ! $nvert;
$cmd   = "$mkge $nvert $altmx $file ";
$cmd  .= (@datl) ? join (' ', @datl) : $datds;
@tbeg  = times;
$retc  = 0xfff & system ($cmd);
@tend  = times;
if ($retc) {
   if ($retc == 0x255) {
      die "$snmt: '$cmd' failed: $!\n";
   } elsif ($retc > 0x80) {
      $retc >>= 8;
      die "$snmt: '$cmd' ran with non-zero exit status $retc\n";
   } else {
      $emsg = "$snmt: '$cmd' ran with";
      if ($retc & 0x80) {
	 $retc &= ~0x80;
	 $emsg .= " coredump from";
      }
      die "$emsg signal $retc\n";
   }
}
printf STDERR "execution times (sec): %.2f (user) and %.2f (system)\n", $tend[2]-$tbeg[2], $tend[3]-$tbeg[3];
exit 0;

#=================================================================================

sub prthelp {
   print "Make a file with global interpolation tables for geomagnetic coordinates
conversions; usage:

  mkglob [-h] [-a altmx] [-d date1:date2:...] [-r nvert] [file]

where
  -a altmx = Maximum required altitude (km).  The actual maximum altitude is
	     higher; it is the next grid point above 'altmx' above which
	     interpolation is not possible; i.e. a fatal error.
  -d dates = One or more dates for which tables are to be created, each
	     expressed as year and fraction (yyyy.fraction), separated by a
	     colon and in chronological order.  The earliest may not precede
	     1900, the first IGRF epoch as defined in magfld.f.
  -h       = Help option: print this help message and quit.
  -r nvert = Spatial resolution parameter must be at least 2 but 30 or more is
	     recommended.  Resolution increases with increasing nvert until
	     about 100 where accuracy degrades close to the poles in east-west
	     gradients.  nvert greater than 1000 degrades accuracy at all
	     latitudes. More details are given in ggrid.f.
  file     = Name of the interpolation tables file to be created.

When not specified via command arguments (-a -d -r file), 'altmx' is $altmxd, 'dates'
(epochs) are $datds,
'nvert' is $nvertd, and 'file' is '$filed'.

A pre-existing 'file' is saved as 'file.old'.

Once created, 'file' is suitable for conversion from geographic coordinates to
Apex, Modified Apex and Quasi-Dipole coordinates and back using entries in
apxntrp.f as exemplified in xglob.f.

Suitable table design must address the adequate maximum grid altitude
constraint, above which it is impossible to interpolate, and may also wish to
avoid non-fatal 'dates' warnings issued by apxntrp.f when there is more than one
year difference in a single epoch 'file', more than five years between epochs,
before the first epoch in a multiple epoch 'file' or after the last epoch.

SEE ALSO:
apxntrp.f - interpolation routines                           - $apxntr
ggrid.f   - routine to get a grid for interpolation tables   - $ggrid
magfld.f  - IGRF definition routines                         - $magfld
mkglob.f  - program to create global interpolation tables    - $mkglob
xglob     - driver script for xglob.f                        - $xglobs
xglob.f   - example program of apxntrp.f and prepared tables - $xglob
";
  exit;
}
