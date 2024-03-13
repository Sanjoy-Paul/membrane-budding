#This is a vmd tcl script to compute distance between index 0-62, 461-523, 922-984 and 1383-1445 for the trajectory prod-0-10mus-now-prot-center.xtc
mol new solu-memb.gro
animate read xtc prod-0-10mus-now-prot-center.xtc waitfor all
set nf [molinfo top get numframes]
set f [open "distances.dat" w]
for {set kf 0} {$kf < $nf} {incr kf} {
set h1 [atomselect top "index 0 to 62" frame $kf]
set h2 [atomselect top "index 461 to 523" frame $kf]
set h3 [atomselect top "index 922 to 984" frame $kf]
set h4 [atomselect top "index 1383 to 1445" frame $kf]
set c1 [measure center $h1]
set c2 [measure center $h2]
set c3 [measure center $h3]
set c4 [measure center $h4]
set d1 [vecdist $c1 $c2]
set d2 [vecdist $c1 $c3]
set d3 [vecdist $c1 $c4]
set d4 [vecdist $c2 $c3]
set d5 [vecdist $c2 $c4]
set d6 [vecdist $c3 $c4]
puts $f "$d1 $d2 $d3 $d4 $d5 $d6"
$h1 delete 
$h2 delete 
$h3 delete 
$h4 delete 
puts $kf
}
close $f
exit
