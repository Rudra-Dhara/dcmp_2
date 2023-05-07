#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/wnOpenSFile.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################


proc wnOpenDirectory {} {
    global system

    set dir [tk_chooseDirectory \
	    -initialdir $system(PWD) \
	    -title      "Choose WIEN2k case directory" \
	    -mustexist 1]

    if { $dir == "" } {
	return ""
    }

    if { ![file isdirectory $dir] } {
	ErrorDialog "file \"$file\" does not exists !!!"
	return ""
    }
    
    return $dir
}


proc wnOpenSFile {{file 0}} {
    global fileselect system periodic species geng xcMisc wn

    ##########################################
    # in WIEN2k there it is always 3D system #
    ##########################################
    set periodic(dim) 3
    set species crystal
 
    if { $file == 0} {
	fileselect "Open WIEN Struct File" 
	if { $fileselect(path) != "" } {
	    set file $fileselect(path)
	} else {
	    puts stderr "WARNING:: \fileselect(path) = \"\""
	    flush stderr
	    return
	}    
    }

    if { ![file exists $file] } {
	ErrorDialog "file \"$file\" does not exists !!!"
	return
    }
    #
    # get the lattice type from struct-file
    #
    set fID [open $file r]
    gets $fID; # title line
    set wn(lattice_type) [lindex [gets $fID] 0]
    close $fID

    #
    # update the title of "."
    #
    set xcMisc(titlefile) $file
    wm title . "XCrySDen: [file tail $file]"
    file copy -force $file $system(SCRDIR)/xc_wienstruct.$system(PID) 

    #########################################################
    # get the filehead out of filename (case.struct --> case)
    set pwd [pwd]
    cd [file dirname $file]
    set filehead [file tail [FileHead $file]]
    if { [catch {exec $system(FORDIR)/str2xcr $filehead}] } {
	ErrorDialog "error while executing \"str2xcr\" program"
	return
    }
    # now copy $filehead.xcr to $system(SCRDIR)/xc_str2xcr.$$
    file rename -force ${filehead}.xcr $system(SCRDIR)/xc_str2xcr.$system(PID)
    cd $pwd

    #
    # WIEN2k struct file is in BOHRs, thatwhy xc_str2xcr.$$ is in BOHRs
    #
    set geng(M3_ARGUMENT) [GetGengM3Arg BOHR 95]
    xcAppendState wien
    
    #
    # determine periodic(igroup)
    #
    cd $system(SCRDIR)
    xcCatchExecReturn $system(BINDIR)/gengeom $geng(M1_INFO) 1 $geng(M3_ARGUMENT) \
	    1  1 1 1 $system(SCRDIR)/xc_gengeom.$system(PID) \
	    $system(SCRDIR)/xc_str2xcr.$system(PID)
    set fileID [open "$system(SCRDIR)/xc_gengeom.$system(PID)" r]
    GetDimGroup periodic(dim) periodic(igroup) $fileID
    close $fileID

    ResetDispModes
    DisplayDefaultMode
    CellMode
    OpenStruct .mesa $system(SCRDIR)/xc_struc.$system(PID)
    # append 'render' to XCState(state) if "render" is not defines yet
    xcAppendState render
    xcUpdateState

    return
}


proc wnOpenRenderDensity {} {
    
    set filedir  [wnOpenDirectory]        
    if { $filedir == "" } {
	return
    }
    set filehead [file tail $filedir]
    set file     $filedir/$filehead.struct
    wnOpenSFile $file
    if { ! [file exists $file] } {
	return
    }
    wnDensity   $filedir    
}


proc wnOpenCalcAndRenderDensity {} {

    set filedir  [wnOpenDirectory]        
    if { $filedir == "" } {
	return
    }
    set filehead [file tail $filedir]
    set file     $filedir/$filehead.struct

    wnOpenSFile       $file
    wnDetComOpt       $filedir
    wnDensity2D_or_3D $filedir
}


proc wnOpenKPath {} {

    set filedir   [wnOpenDirectory]        
    if { $filedir == "" } {
	return
    }

    wnDetComOpt  $filedir
    wnKPath      $filedir
}


proc wnOpenFS {} {
    
    set filedir  [wnOpenDirectory]       
    if { $filedir == "" } {
	return
    }
    
    wnDetComOpt  $filedir
    wnFSInit     $filedir
}