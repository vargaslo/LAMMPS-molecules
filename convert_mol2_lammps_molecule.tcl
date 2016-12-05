# This script is run from within VMD in the usual way:
# source convert_mol2_lammps_molecule.tcl
# 
# From the command line:
# vmd -e convert_mol2_lammps_molecule.tcl


package require topotools

# load molecule into vmd -- must contain bond information
set infile pentane_eh.mol2
mol new $infile

# set atom types manually
[atomselect top "serial>5"] set type HC
[atomselect top "serial 1 5"] set type CH3
[atomselect top "serial 2 3 4"] set type CT
topo retypebonds

# Amber-ii atoms have zero charge
[atomselect top all] set charge 0

# find 1-3 and 1-4 neighbors
topo guessangles
topo guessdihedrals


#------------------------------------------------------------------------------

proc molec_file {outfile} {

variable infile

set fileID [open $outfile "w"]

# header information
puts $fileID  "# LAMMPS Molecule File created from $infile\n"
puts $fileID  "[format %4d [topo numatoms]] atoms"
puts $fileID  "[format %4d [topo numbonds]] bonds"
puts $fileID  "[format %4d [topo numangles]] angles"
puts $fileID  "[format %4d [topo numdihedrals]] dihedrals"
puts $fileID  "[format %4d [topo numimpropers]] impropers"

# Coordinates
set i 0
puts $fileID  "\nCoords\n"
foreach at [[atomselect top all] get {x y z}] {
    incr i 1
    lassign $at x y z
    puts $fileID  [format "%4d %10.6f %10.6f %10.6f" $i $x $y $z]
}

# Atom types
set i 0
set types [topo atomtypenames]
puts $fileID  "\nTypes\n"
foreach type [[atomselect top all] get type] {
    incr i
    set atomtype [expr "[lsearch $types $type] + 1"]
    puts $fileID  [format "%4d %4d  # %s" $i $atomtype $type]
}

# Charges
set i 0
puts $fileID  "\nCharges\n"
foreach ch [[atomselect top all] get charge] {
    incr i
    puts $fileID  [format "%4d %10.4f" $i $ch]
}

# Bonds
set i 0
set types [topo bondtypenames]
puts $fileID  "\nBonds\n"
foreach b [topo getbondlist type] {
    incr i
    lassign $b ind1 ind2 type
    set bondtype [expr "[lsearch $types $type] + 1"]
    set ser1 [expr "$ind1+1"]
    set ser2 [expr "$ind2+1"]
    puts $fileID  [format "%4d %4d %4d %4d  # %s" $i $bondtype $ser1 $ser2 $type]
}

# Angles
set i 0
set types [topo angletypenames]
puts $fileID  "\nAngles\n"
foreach a [topo getanglelist] {
    incr i
    lassign $a type ind1 ind2 ind3
    set angletype [expr "[lsearch $types $type] + 1"]
    set ser1 [expr "$ind1+1"]
    set ser2 [expr "$ind2+1"]
    set ser3 [expr "$ind3+1"]
    puts $fileID  [format "%4d %4d %4d %4d %4d  # %s" $i $angletype $ser1 $ser2 $ser3 $type]
}

# Dihedrals
set i 0
set types [topo dihedraltypenames]
puts $fileID  "\nDihedrals\n"
foreach d [topo getdihedrallist] {
    incr i
    lassign $d type ind1 ind2 ind3 ind4
    set dihedraltype [expr "[lsearch $types $type] + 1"]
    set ser1 [expr "$ind1+1"]
    set ser2 [expr "$ind2+1"]
    set ser3 [expr "$ind3+1"]
    set ser4 [expr "$ind4+1"]
    puts $fileID  [format "%4d %4d %4d %4d %4d %4d  # %s" $i $dihedraltype $ser1 $ser2 $ser3 $ser4 $type]
}

# Impropers
set i 0
puts $fileID  "\nImpropers\n"
set impropers [list ]
foreach atom_ind [[atomselect top all] get index] {
    set atom_type [[atomselect top "index $atom_ind"] get type]
    set bonded_atoms [[atomselect top "index $atom_ind"] getbonds]
    if {[llength [lindex $bonded_atoms 0]]==3} {
        set bonded_atomtypes [[atomselect top "index [lindex $bonded_atoms 0]"] get type]
        set sortorder [lsort -dictionary -indices $bonded_atomtypes]
                
        lappend bonded_atomtypes_ordered [lindex $bonded_atomtypes [lindex $sortorder 0]]
        lappend bonded_atomtypes_ordered [lindex $bonded_atomtypes [lindex $sortorder 1]]
        lappend bonded_atomtypes_ordered [lindex $bonded_atomtypes [lindex $sortorder 2]]
        
        # identify improper type (or assign new type if needed)
        set improper_type [expr {1 + [lsearch impropers $bonded_atomtypes_ordered]}]
        if {$improper_type==0} {
            lappend impropers $bonded_atomtypes_ordered
            set improper_type [expr {1 + [lsearch impropers $bonded_atomtypes_ordered]}]
        }

        incr i          
        set ind1 $atom_ind
        set ind2 [lindex $bonded_atoms 0 [lindex $sortorder 0]]
        set ind3 [lindex $bonded_atoms 0 [lindex $sortorder 1]]
        set ind4 [lindex $bonded_atoms 0 [lindex $sortorder 2]]

        lappend type $atom_type
        lappend type [lindex $bonded_atomtypes_ordered 0]
        lappend type [lindex $bonded_atomtypes_ordered 1]
        lappend type [lindex $bonded_atomtypes_ordered 2]
        
        set ser1 [expr "$ind1+1"]
        set ser2 [expr "$ind2+1"]
        set ser3 [expr "$ind3+1"]
        set ser4 [expr "$ind4+1"]
        puts $fileID [format "%4d %4d %4d %4d %4d %4d  # %s" $i $improper_type $ser1 $ser2 $ser3 $ser4 $type]
    }
}

puts "\n\nMolecule file written: $outfile"
puts "\nManually check for impropers!\n\n"

}


molec_file mymolecule.molec
