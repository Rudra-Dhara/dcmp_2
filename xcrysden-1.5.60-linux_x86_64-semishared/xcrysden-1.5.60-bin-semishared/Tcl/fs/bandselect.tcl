# 
# FS_bandSelect --
#
# Select a bands for Fermi-surface plotting among the all the bands from
# small toplevel window. Also a abnd-widths are ploted for orientation
#
proc FS_bandSelect {spin} {
    global fs

    # specify fermi energy

    OneEntryToplevel [WidgetName] "Fermi Energy" "Ferm Energy" \
	"Specify the Fermi Energy:" 15 fs(Efermi) float 300 20

    #
    # display the bandwidths in a Text widget
    #
    set text [FS_displayBandWidths $fs($spin,bandwidthfile)]
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
    grapher_BARGraph $fs($spin,bandwidthfile) \
	-Xtitle     $xlabel \
	-Ytitle     "E" \
	-Yline      $fs(Efermi) \
	-Yline_text Ef
    set graph [Grapher BARGraph]
    wm geometry $graph +0+0
    raise $graph

    #
    # select bands to plot Fermi surface
    #
    global grafdata
    set gID [CurrentGrapherID]
    set fs($spin,nbands) $grafdata($gID,N_point,1)
    set t [xcToplevel [WidgetName] "Select bands" "Select Bands" . 0 0 0]
    raise $t
    
    for {set i 1} {$i <= $fs($spin,nbands)} {incr i} {
	set fs($spin,$i,minE) $grafdata($gID,$i,1,1)
	set fs($spin,$i,maxE) $grafdata($gID,$i,2,1)
    }

    label $t.l \
	-text "Select bands for Fermi Surface drawing:" \
	-relief ridge -bd 2
    pack $t.l -side top -expand 1 -fill x -padx 2m -pady 3m \
	-ipadx 2m -ipady 2m

    
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

    # create FRAME to hold all checkbuttons
    set f [frame $c.f -bd 0]
    $c create window 0 0 -anchor nw -window $f -tags frame
 
    # CHECKBUTTONS
    for {set i 1} {$i <= $fs($spin,nbands)} {incr i} {
	set fs($spin,$i,band_selected) 0
	set cb [checkbutton [WidgetName $f] -text "Band number: $i" \
		    -variable fs($spin,$i,band_selected) -relief ridge -bd 2]
	pack $cb -side top -padx 2m -pady 1m -fill x -expand 1
    }

    # make correct DISPLAY
    set child [lindex [pack slaves $f] 0]       
    tkwait visibility $child
    set width  [winfo width $c]
    set height [winfo height $f]
    if { $fs($spin,nbands) < 8 } {
	$c config -width $width -height $height 
    } else {
	$c config \
	    -width $width -height [expr 8*($height / $fs($spin,nbands))] \
	    -scrollregion "0 0 $width $height"
    }
    
    #
    # press the "Selected" button when done
    #
    set b [button [WidgetName $t] -text "Selected" \
	       -command [list FS_bandSelect:selected $t $spin]]
    pack $b -side top -expand 1 -fill x -padx 2m -pady 3m    
    
    if { $spin != {} } {
	set l [label [WidgetName $t] -text " Spin type: $spin" \
		   -relief ridge -bd 4 -anchor w]
	pack $l -side top -expand 1 -fill x -padx 1m -pady 1m \
	    -ipadx 2m -ipady 2m
    }
    tkwait variable fs($spin,selected)
}
proc FS_bandSelect:selected {t spin} {
    global fs

    set fs($spin,selected) 1
    CancelProc $t
}