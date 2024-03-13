#This is a vmd tcl script to compute d1 (index pair 4,299) and d2 (index pair 3165,3460) as described in Fig-1 for the trajectories production-0-69ns-P.xtc and production-69-1000ns-2-P.xtc  
mol new P.gro
animate read xtc production-0-69ns-P.xtc waitfor all
animate read xtc production-69-1000ns-2-P.xtc waitfor all
animate delete beg 0 end 0
set nf [molinfo top get numframes]
set f [open "vectors-for-angle.dat" w]
for {set k 0} {$k < $nf} {incr k} {
set sel [atomselect top "index 4" frame $k]
set h1a [measure center $sel]
set sel [atomselect top "index 299" frame $k]
set h1b [measure center $sel]
set h1 [vecnorm [vecsub $h1b $h1a]]

set sel [atomselect top "index 3165" frame $k]
set h2a [measure center $sel]
set sel [atomselect top "index 3460" frame $k]
set h2b [measure center $sel]
set h2 [vecnorm [vecsub $h2b $h2a]]

puts $f "$h1 $h2"
puts $k
}
close $f
exit
