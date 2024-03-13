#This is a vmd tcl script to print out coordinates for different grids for each of the leaflets as a function of simulation time to compute the time evolution of curvature out of it 
mol new ../P.gro
set selp [atomselect top "name P"]
set memb [atomselect top "resname DOPC DOPE DOPS DOPA" frame 0]
set com [measure center $memb]
set zcom [lindex $com 2]
set P_upper [atomselect top "name P and z < $zcom"]
set P_lower [atomselect top "name P and z > $zcom"]
set P_upper_resid [$P_upper get resid]
set P_lower_resid [$P_lower get resid]
set R [list $P_upper_resid]
set X [molinfo top get a]
set xint 10
set xgrid [expr {round($X/$xint)}]
set Y [molinfo top get b]
set yint 10
set ygrid [expr {round($Y/$yint)}]
set xyzmin [lindex [measure minmax $selp] 0]
set xmin [lindex $xyzmin 0]
set ymin [lindex $xyzmin 1]

animate read xtc ../production-0ns-260ns-P.xtc waitfor all
#animate delete beg 0 end 0

set k [molinfo top get numframes]

#for {set r 1} {$r<24} {incr r} {
#set f3 [open "resid-$r.dat" w]
#for {set f 500} {$f<$k} {incr f} {
#set sel [atomselect top "resid $r and noh" frame $f]
#set O [measure center $sel]
#puts $f3 $O
#}
#puts $r
#close $f3
#}

for {set k1 0} {$k1 < $k} {incr k1} {
set all [atomselect top "all" frame $k1]
set p1 [atomselect top "name P" frame $k1]
set frame0 [atomselect top "name P" frame 0]
$p1 frame $k1
$frame0 frame 0
$all move [measure fit $p1 $frame0]
}
animate delete beg 0 end 0
#set x1 $xmin
#set y1 $ymin
for {set p1 0} {$p1<=$xgrid} {incr p1} {
set x1 [expr {$xmin + $p1*$xint}]
set p11 [expr {$p1+1}]
set x2 [expr {$xmin + $p11*$xint}]
for {set p2 0} {$p2<=$ygrid} {incr p2} {
set y1 [expr {$ymin + $p2*$yint}]
set p22 [expr {$p2+1}]
set y2 [expr {$ymin + $p22*$yint}]
set f1 [open "P-coordinate-upper-grid-$p1-$p2.dat" w]
set f2 [open "P-coordinate-lower-grid-$p1-$p2.dat" w]
for {set f 50} {$f<$k} {incr f} {
set sel [atomselect top "resid $P_upper_resid and name P and x > $x1 and x < $x2 and y > $y1 and y < $y2" frame $f]
set M1 [$sel get resid]
#puts $M1
set l1 [llength $M1]
if {$l1>0} {
set pcoor [measure center $sel]
} else {set pcoor [list 0 0 0]}
puts $f1 $pcoor
$sel delete
set sel [atomselect top "resid $P_lower_resid and name P and x > $x1 and x < $x2 and y > $y1 and y < $y2" frame $f]
set M1 [$sel get resid]
set l1 [llength $M1]
if {$l1>0} {
set pcoor [measure center $sel]
} else {set pcoor [list 0 0 0]}
puts $f2 $pcoor
$sel delete
#set y1 $y2
}
close $f1
close $f2
#set x1 $x2
puts "$p1 $p2"
}
#puts $p1
}

exit

