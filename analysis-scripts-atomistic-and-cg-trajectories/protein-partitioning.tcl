#This is a vmd tcl script to compute partition coefficient as described in eqn. 1  
mol new solu-memb.gro
animate read xtc step7_production_now.xtc waitfor all
#animate read xtc step7_production_80_160ns_now.xtc waitfor all
animate delete beg 0 end 0
set r 5
set f1 [open "hydrophilic-atoms-dimer-inside-lipid-0-160ns.dat" w]
set nf [molinfo top get numframes]
set memb [atomselect top "resname DOPC or resname DOPA or resname DOPE or resname DOPS"]
set hydrophilic_ref [atomselect top "protein and resname LYS ASP GLN GLY SER"]
set num_up_ref [$hydrophilic_ref num]
set hydrophobic_ref [atomselect top "protein and not resname LYS ASP GLN GLY SER"]
set num_down_ref [$hydrophobic_ref num]
set com [measure center $memb]
for {set k 0} {$k < $nf} {incr k} {
set zcom [lindex $com 2]
set P_upper [atomselect top "name P and z > $zcom" frame $k]
set P_lower [atomselect top "name P and z < $zcom" frame $k]

$P_upper frame $k
$P_lower frame $k

set z_lower [lindex [measure center $P_lower] 2]
set z_upper [lindex [measure center $P_upper] 2]
set hydrophilic [atomselect top "protein and resname LYS ASP GLN GLY SER and z > $z_lower and z > $z_upper" frame $k]
set hydrophobic [atomselect top "protein and not resname LYS ASP GLN GLY SER and z > $z_lower and z < $z_upper" frame $k]

set num_up_1 [$hydrophilic num]
set num_down_1 [$hydrophobic num]

set num_up [expr {$num_up_1/$num_up_ref}]
set num_down [expr {$num_down_1/$num_down_ref}]

puts $f1 "$k $num_up_1 $num_up_ref $num_down_1 $num_down_ref"
puts $k
}
close $f1
exit
