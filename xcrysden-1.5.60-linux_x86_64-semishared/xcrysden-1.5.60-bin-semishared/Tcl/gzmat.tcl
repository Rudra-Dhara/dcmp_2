#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/gzmat.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

# ------------------------------------------------------------------------
# loads structure from Gaussian input file (requires babel)
# ------------------------------------------------------------------------
proc gzmat {filedir {viewmol_exists {}}} {
    global xcMisc env system

    if { ! [info exists xcMisc(babel)] } {
	if [winfo exists .title] { 
	    destroy .title 
	}
	ErrorDialog "--gzmat option requires the definition of xcMisc(babel) in the ~/.xcrysden/custom-definitions"
	exit
    } else {
	set env(BABEL) $xcMisc(babel)
    }

    # OpenBabel seems not to need BABEL_DIR variable !!!

    #if { ! [info exists xcMisc(babel_dir)] && ! [info exists env(BABEL_DIR)] } {
    #	if [winfo exists .title] { 
    #	    destroy .title 
    #	}
    #	ErrorDialog "neither xcMisc(babel_dir) nor environmental BABEL_DIR variables defined. Can't run BABEL."
    #	exit
    #}
	
    # xcMisc(babel_dir) has priority over env(BABEL_DIR)

    if { [info exists xcMisc(babel_dir)] } {
	set env(BABEL_DIR) $xcMisc(babel_dir)
    }

    set head [file rootname [file tail $filedir]]
    puts stderr "Executing babel ..."
    xcCatchExecReturnRedirectStdErr sh $system(TOPDIR)/scripts/gzmat2xsf.sh $filedir > $system(SCRDIR)/$head.xsf
    
    if { ![file exists $system(SCRDIR)/$head.xsf] \
	     || [file size $system(SCRDIR)/$head.xsf] == 0 } {
	if [winfo exists .title] { 
	    destroy .title 
	}
	ErrorDialog "an error occured while executing BABEL program."
	exit
    }

    if { $viewmol_exists == "" } {
	ViewMol .
    }
    xsfOpen $system(SCRDIR)/$head.xsf
}
	

proc gzmat_menu {openwhat} {
    global system
    set filedir [tk_getOpenFile -defaultextension .xsf \
		     -filetypes { 
			 { {All Files}  {.*}  }
		     } -initialdir $system(PWD) -title $openwhat]
    if { $filedir == "" } {
	return
    }
    gzmat $filedir viewmol_exists
}