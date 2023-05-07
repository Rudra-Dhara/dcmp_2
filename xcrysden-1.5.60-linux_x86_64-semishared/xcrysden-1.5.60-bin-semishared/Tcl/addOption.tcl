#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/addOption.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

#
#
# addOption --unknown /path/converterProgram -xterm
proc addOption {option converterProgram description {args {}}} {
    global addOption

    lappend addOption(optionList)     $option
    lappend addOption(converterList)  $converterProgram
    lappend addOption(description)    $description
    #
    # here some preprocessing of the args
    # ... insert HERE ...
    lappend addOption(args)           $args
    
    lappend addOption(usage) [format "  %s  <file>\n%s\n" $option $description]
}


proc addOption:parse {option filedir} {
    global addOption system

    if { ! [info exists addOption(optionList)] } {
	#wm withdraw .title
	if { [winfo exists .title] } {
	    destroy .tilte
	}
	ErrorDialog "bad option \"$option\"" 
	XCrySDenUsage
    }
    
    set head [file rootname [file tail $filedir]]
    set ind [lsearch $addOption(optionList) $option]
    if { $ind >= 0 } {
	set converter [lindex $addOption(converterList) $ind]	
	#
	# take care of args (tey to be done)
	if { [eval xcCatchExec $converter $filedir > $system(SCRDIR)/$head.xsf] } {
	    exit 1
	}
	ViewMol .
	xsfOpen $system(SCRDIR)/$head.xsf
    } else {
	ErrorDialog "ERROR: Bad option: $option" 
	XCrySDenUsage
    }
}
	
 
proc addOption:register {} {
    global system env

    if { [file exists $system(TOPDIR)/addOptions] } {
	source $system(TOPDIR)/addOptions
    }
}


proc addOption:printCustomUsage {} {
    global addOption

    if { ! [info exists addOption(optionList)] } {
	return
    }

    puts stderr "-----------------------"
    puts stderr "   Custom User options"
    puts stderr "-----------------------\n"
    foreach usage $addOption(usage) {
	puts stderr $usage
    }
    puts stderr {}
}


proc addOption:hardcoded {converter {filedir {}} {openwhat {}} {viewmol_exists {}}} {
    global system

    if { $filedir == {} } {
	set filedir [tk_getOpenFile -defaultextension .xsf \
			 -filetypes { 
			     { {All Files}  {.*}  }
			 } -initialdir $system(PWD) -title $openwhat]
	if { $filedir == "" } {
	    return
	}
	set viewmol_exists 1
    }
    
    set head [file rootname [file tail $filedir]]
    if { [eval xcCatchExec $converter $filedir > $system(SCRDIR)/$head.xsf] } {
	exit 1
    }
    if { $viewmol_exists == "" } {
	ViewMol .
    }
    xsfOpen $system(SCRDIR)/$head.xsf
}
