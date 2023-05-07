#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/wnFS.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

#
# EXECUTED BY BUTTON:  "Generate k-mesh"
#
proc wnGenKMesh {} {
    global wn

    #
    # order of questions in KGEN:
    #
    # 1. "NUMBER OF K-POINTS IN WHOLE CELL:"
    # 2. "Shift of k-mesh allowed. Do you want to shift:" 

    ###########
    cd $wn(dir)
    ###########

    set input    "$wn(fs_nkp)\n"
    append input "0\n"; # shift of k-mesh is not allowed

    xcDebug -debug "KGEN INPUT: $input"

    WriteFile xc_kgen.inp $input w
    catch {exec x kgen < xc_kgen.inp > xc_kgen.out} exit_status

    set out [ReadFile xc_kgen.out]
    dialog [WidgetName] "Notification" \
	"K-MESH GENERATED !!!.\nOutput of kgen:\n\n $out\n\nExit status:\n$exit_status" \
	info 0 Done
}


# not used anymore
# #
# # EXECUTED BY BUTTON: "Change to unit 5 in *in1*"
# #
# proc wnTo5In1 in1 {
#     global wn system
# 
#     ###########
#     cd $wn(dir)
#     ###########
# 
#     #
#     # KLIST file is:   wn(filehead).klist
#     #
#     set klistf $wn(filehead).klist
# 
#     if ![file exists $klistf] {
# 	tk_dialog [WidgetName] ERROR \
# 	    "Please generate the k-mesh first !!!" \
# 	    error 0 OK
# 	return
#     }
# 
#     if ![info exists wn(fs_in1_copy_index)] {
# 	set wn(fs_in1_copy_index) 0
# 	set wn(in1_orig) ${in1}_orig.0
#     } else {
# 	incr wn(fs_in1_copy_index)
#     }
# 
#     file copy -force -- $in1 ${in1}_orig.$wn(fs_in1_copy_index)
# 
#     # 5090 FORMAT(20X,I1,2f10.1)
#     if [catch {exec $system(awk) \
# 	    -v unit=5 -f $system(AWKDIR)/in1.awk ${in1}_orig.0 > ${in1}_tmp}] {
# 	tk_dialog [WidgetName] ERROR \
# 	    "ERROR while executing \"in1.awk\" program" \
# 	    error 0 OK
# 	return
#     }
#     exec cat ${in1}_tmp $klistf > $in1
#     exec rm -f ${in1}_tmp
# 
#     tk_messageBox -message "Notification:\n\nUnit was set to 5 and k-mesh was inserted !!!\n" \
# 	-type ok -icon info
# }


#
# EXECUTED BY BUTTON:   "Render Fermi Surface"
#
proc wnFSGo {outkgen {spin {}}} {
    global wn system

    SetWatchCursor
    update

    xcDebug -stderr "DEBUG -- args: $outkgen $spin"
    ###########
    cd $wn(dir)
    ###########

    set outkgen [file tail $outkgen]

    # try outkgen
    # try $wn(dir)/$wn(filehead).fs_outputkgen
    # try $wn(dir)/$wn(filehead).outputkgen !!! DANGEROUS
    if { ! [file exists $outkgen] } {
	tk_dialog [WidgetName] ERROR \
	    "Please generate the k-mesh first !!!" error 0 OK
	ResetCursor; return
    }
    
    set ok 0
    foreach line [split [ReadFile $outkgen] \n] {
	if [string match "*k-vectors:*" $line] {
	    set ok 1
	    break
	}
    }
    
    if { ! $ok } {
	ErrorDialog "Please generate the k-mesh first !!!"
	ResetCursor
	return	
    }

    #
    # read output1 file by wn_readbands prog !!!
    #
    append deffile "6,'$wn(filehead).outputba$spin','unknown','formatted',0\n"
    append deffile "7,'$wn(filehead).output1$spin' ,'old',    'formatted',0\n"
    append deffile "8,'$wn(filehead).outputbw$spin','unknown','formatted',0\n"
    
    set wn($spin,fs_bandfile)      $wn(filehead).outputba$spin
    set wn($spin,fs_bandwidthfile) $wn(filehead).outputbw$spin
    
    WriteFile band.def $deffile w
    if { [catch {exec $system(BINDIR)/wn_readbands band.def} error_msg] } {
	ErrorDialog "while executing wn_readbands program !!!" $error_msg
	ResetCursor
	return
    }
    
    #
    # execute program wn_readbakgen
    #
    set deffile {}
    append deffile "7, '$wn(filehead).outputba$spin','unknown','formatted',0\n"
    append deffile "8, '$outkgen'                   ,'old'    ,'formatted',0\n"
    append deffile "10,'$wn(filehead).outputfs$spin','unknown','formatted',0\n"

    WriteFile bakgen.def $deffile w
    if { [xcCatchExec $system(BINDIR)/wn_readbakgen bakgen.def] } {
	ResetCursor
	return
    }
    set wn($spin,fs_fsfile) $wn(filehead).outputfs$spin    
    file rename -force -- $wn($spin,fs_fsfile) $system(SCRDIR)/$wn(filehead)
    set wn($spin,fs_fsfile) $system(SCRDIR)/$wn(filehead)/$wn($spin,fs_fsfile)

    ResetCursor
    update

    ##############################################################
    #
    # BAND SELECTION
    #
    set wn(fs_Efermi) 0.0
    catch {set wn(fs_Efermi) \
	       [exec grep :FER $wn(filehead).output2$spin | \
		    tail -1 | awk "{print \$NF}"]}
    OneEntryToplevel [WidgetName] "Fermi Energy" "Ferm Energy" \
	"Specify the Fermi Energy:" 15 wn(fs_Efermi) float 300 20

    #
    # display the bandwidths in a Text widget
    #
    set text [FS_displayBandWidths $wn($spin,fs_bandwidthfile) $spin]    
    wm geometry $text +0-0
    raise $text 

    #
    # make a band-width graph !!!
    #
    set xlabel "Band Widths" 
    if { $spin != {} } {
	append xlabel " (spin type: $spin)"
    }    
    GraphInit
    grapher_BARGraph $wn($spin,fs_bandwidthfile) \
	    -Xtitle $xlabel \
	    -Ytitle "E / Ry" \
	    -Yline  $wn(fs_Efermi) \
	    -Yline_text Ef
    set graph [Grapher BARGraph]
    wm geometry $graph -0-0
    raise $graph

    global grafdata
    set gID [CurrentGrapherID]
    set wn($spin,nbands) $grafdata($gID,N_point,1)

    for {set i 1} {$i <= $wn($spin,nbands)} {incr i} {
	set wn($spin,$i,minE) $grafdata($gID,$i,1,1)
	set wn($spin,$i,maxE) $grafdata($gID,$i,2,1)
    }

    #
    # select bands window
    #

    set t [xcToplevel [WidgetName] "Select bands" "Select Bands" .fs_init 0 0 1]
    wm geometry $t -0+0
    raise $t    
    tkwait visibility $t

    label $t.l \
	-text "Select bands for Fermi Surface drawing:" \
	-relief ridge -bd 2
    pack $t.l -side top -expand 1 -fill x -padx 2m -pady 3m \
	-ipadx 2m -ipady 2m

    #
    # we should make a scrolled window
    #
    # CANVAS & SCROLLBAR in CANVAS
    set scroll_frame [frame $t.f -relief sunken -bd 1]
    pack $scroll_frame -side top -expand true -fill y -padx 5
    set c [canvas $scroll_frame.canv \
	       -yscrollcommand [list $scroll_frame.yscroll set]]
    set scb [scrollbar $scroll_frame.yscroll \
		 -orient vertical -command [list $c yview]]
    pack $scb -side right -fill y
    pack $c -side left -fill both -expand true

    bind $c <4> { %W yview scroll -5 units }
    bind $c <5> { %W yview scroll +5 units }

    # create FRAME to hold all checkbuttons
    set f [frame $c.f -bd 0]
    $c create window 0 0 -anchor nw -window $f -tags frame
 
    # CHECKBUTTONS
    for {set i 1} {$i <= $wn($spin,nbands)} {incr i} {
	set wn($spin,$i,band_selected) 0
	set cb [checkbutton [WidgetName $f] -text "Band number: $i" \
		    -variable wn($spin,$i,band_selected) -relief ridge -bd 2]
	pack $cb -side top -padx 2m -pady 1m -fill x -expand 1
    }

    # make correct DISPLAY
    set child [lindex [pack slaves $f] 0]       
    update
    #tkwait visibility $child
    set width  [winfo reqwidth $c]
    set height [winfo reqheight $f]
    if { $wn($spin,nbands) < 8 } {
    	 $c config -width $width -height $height 
    } else {
    	 $c config \
	     -width $width -height [expr 8*($height / $wn($spin,nbands))] \
	     -scrollregion "0 0 $width $height"
    }

    set b [button [WidgetName $t] -text "Selected" \
	       -command [list wnToFS $t $spin]]
    pack $b -side top -expand 1 -fill x -padx 2m -pady 3m    

    if { $spin != {} } {
	set l [label [WidgetName $t] -text " Spin type: $spin" \
		   -relief ridge -bd 4 -anchor w]
	pack $l -side top -expand 1 -fill x -padx 1m -pady 1m \
	    -ipadx 2m -ipady 2m
    }
}

# not used anymore
# proc wnFSResetIn1 in1 {
#     global wn
#     
#     cd $wn(dir)
#     if { ! [info exists wn(in1_orig)] } {
# 	set wn(in1_orig) ${in1}_orig.0
#     }
# 
#     if { ! [file exists $wn(in1_orig)] } {
# 	tk_messageBox \
# 	    -message "Warning:\n\nCan't reset back to unit 4 !!!\n" \
# 	    -type ok -icon warning
# 	return
#     } else {
# 	exec cp $wn(in1_orig) $in1
# 	
# 	tk_messageBox \
# 	    -message "Notification:\n\nThe $in1 has been reseted !!!\n" \
# 	    -type ok -icon info
# 	return
#     }
# }


proc wnToFS {t {spin {}}} {
    global wn fs xcMisc

    CancelProc $t

    set fs(titlefile)      $wn(filehead)
    set fs($spin,nbands)   $wn($spin,nbands)
    set fs(Efermi)         $wn(fs_Efermi)
    for {set i 1} {$i <= $fs($spin,nbands)} {incr i} {
	set fs($spin,$i,band_selected) $wn($spin,$i,band_selected)
	set fs($spin,$i,minE)          $wn($spin,$i,minE) 
	set fs($spin,$i,maxE)          $wn($spin,$i,maxE)
	set fs($spin,$i,isolevel)      $fs(Efermi)
    }

    # lets read (band)XSF file; xc_readbandXSF return info structure
    # which look like:
    #
    # 1 { grid_index 3D grid_ident grid_nband 1 {subgrid0_ident}}
    set sinfo [xc_readbandXSF $wn($spin,fs_fsfile)]
    xcDebug -debug "DEBUG: FS-sinfo: $sinfo"
    # parse $sinfo
    set slist [lindex $sinfo 1]
    set fs($spin,grid_index)    [lindex $slist 0]
    set fs($spin,grid_subindex) [expr [lindex $slist 4] - 1]
    if { [lindex $slist 3] != $wn($spin,nbands) } {
	dialog [WidgetName] "ERROR" \
	    "Mismatch occured while reading FermiSurface File: $wn($spin,fs_fsfile)" error 0 Done
    }
    FS_GoFermi $spin
}

