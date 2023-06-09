#!/bin/sh
#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/scripts/pwi2xsf.sh
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

# set locales to C
LANG=C 
LC_ALL=C
export LANG LC_ALL

#
# pwi2xsf.sh: PW-input to XSF converison
#
# Usage: pwi2xsf [-r] pw-input-file
#
# Written by Tone Kokalj on Tue May  8 20:43:44 2001
#            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

if [ "$#" -lt 1 ]; then
    echo "
Usage: pwi2xsf.sh [-r] pw-input

Option for PWscf version < 1.2:

-r ... for pw.x < 1.2 one must spefify the ityp->nat conversion, and
       the corresponding data are writen to file nuclei.charges. The
       -r flag deletes this file as to force new ityp->nat specification.
"
    exit 1
fi

r=0
if [ "$1" = "-r" ]; then
    r=1
    shift
fi


#######################################
if test "x`type tee`" = "x"; then
    # no tee cmd; make a function-substitute
    tee() {
	cat - > $1
	cat $1
    }
fi

if test "x`type readlink`" = "x"; then
    # no readlink cmd; make a function-substitute
    readlink() {
	echo `ls -l $1 | awk '{print $1}'`
    }
fi


pathname() {
    file=`type -p $1`
    if test $? -gt 0; then
	file=`which $1`
	if test $? -gt 0; then
	    # give-up
	    file=$1
	fi
    fi
    echo $file
}


pathdir() {
    file=`pathname $1`
    
    while test -h $file; do  
	file=`readlink $file`
    done

    dir=`dirname $file`
    ( cd $dir; pwd )
}

if test -z $XCRYSDEN_TOPDIR; then
    # XCRYSDEN_TOPDIR does not exists, guess it from the process
    script_dir=`pathdir $0`
    export XCRYSDEN_TOPDIR=`(cd $script_dir/..; pwd)`
fi

if test -f $XCRYSDEN_TOPDIR/scripts/pwLib_old.sh ; then
    . $XCRYSDEN_TOPDIR/scripts/pwLib_old.sh
    load_old_lib=1
else
    load_old_lib=0
fi
#######################################

#
# check if we have OLD or NEW PW.X input format
#
new_format1=`grep 'ATOMIC_POSITIONS' $1`
new_format2=`grep -i '&system' $1`

if [ "$new_format1" != ""  -a  "$new_format2" != "" ]; then
    # we have NEW PW.X input format
    #
cat $1 | awk 'BEGIN {RS=",";} {print $0}' | awk '
BEGIN {
  calculation="";
  num_of_images="";
  nml_end=0;
  nml_end_string="";
}

toupper($0) ~ /&SYSTEM/          { print; }

/=/ { 
  if ( toupper($1) ~ /^IBRAV($|=)|^CELLDM\([1-6]\)($|=)|^NAT($|=)|^A($|=)|^B($|=)|^C($|=)|^COSAB($|=)|^COSAC($|=)|^COSBC($|=)/ ) { print; } 
  
  if ( toupper($1) ~ /^CALCULATION($|=)/ ) { calculation=toupper($0); }

  if ( toupper($1) ~ /^NUM_OF_IMAGES($|=)/ ) { num_of_images=toupper($0); }
}

/ATOMIC_POSITIONS|CELL_PARAMETERS/ {
  if ( !nml_end) {
     # first finish the namelist
     nml_end=1;
     if (calculation != "")   print calculation;
     if (num_of_images != "") print num_of_images;
     print nml_end_string;
  }
  # now print the current record
  print_line=1;
  print toupper($0);   
  next;
}


toupper($0) ~ /&END|^\/|^ +\// { 
  nml_end_string=$0;
}

/a*/ {
  if ( print_line == 1 ) {
    print $0;
  }
}'> pw.$$
    PWI2XSF=pwi2xsf
else 
    # we have OLD PW.X input format	

    if test $load_old_lib = 0; then
	echo "
ERROR: cannot convert to XSF, because loading of pwLib_old.sh failed
"
	exit 1
    fi

    pwNucleiCharges $1 /dev/null

    cat $1 | awk 'BEGIN {RS=",";} {print}' | awk '
BEGIN {
    end=0;
}
toupper($0) ~ /&INPUT|CELLDM|NAT|LTAUCRY/ { print; }
toupper($0) ~ /IBRAV/ { 
    print;
    split($0,a,"=");
    split(a[1],b,",");
    ibrav = b[1];
}
toupper($0) ~ /&END|^\/|^ \// { end=1; }
/a*/ {
    if ( end == 1 ) print;
}' > pw.$$
    PWI2XSF=pwi2xsf_old
fi
#
# execute $PWI2XSF fortran program and print the XSF file
#
if test -f $XCRYSDEN_TOPDIR/bin/$PWI2XSF ; then
    $XCRYSDEN_TOPDIR/bin/$PWI2XSF < pw.$$ | tee pwi2xsf.xsf_out
elif test -f $XCRYSDEN_LIB_BINDIR/$PWI2XSF ; then
    $XCRYSDEN_LIB_BINDIR/$PWI2XSF < pw.$$ | tee pwi2xsf.xsf_out
else
    $PWI2XSF < pw.$$ | tee pwi2xsf.xsf_out
fi
#rm -f pw.$$

if [ "$r" -eq 1 ]; then
    if [ -f nuclei.charges ]; then rm nuclei.charges; fi
fi

for file in pw.$$
do
    if test -f $file; then rm -f $file; fi
done
     
exit 0
