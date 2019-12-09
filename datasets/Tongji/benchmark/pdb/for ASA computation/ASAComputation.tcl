source _sasa_res.tcl
set fr [open "training1.txt" r]
set file_data [read $fr]
set data [split $file_data "\n"]
set ind 0
set radius {0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4}
foreach line $data {
	set contents [split $line ":"]
	set complex [lindex $contents 0]
	set chain [lindex $contents 1]
	mol load pdb "$complex.pdb"
	#set command "\"protein and chain $chain or nucleic\" 1.4 1 1 $complex\_$chain\_RNA"
	#puts $command
	foreach radii $radius {
		getAllResSASA  "protein and chain $chain" $radii 1 1 "$complex\_$chain\_$radii"
	}
	puts $ind
	set ind [expr $ind+1]
}
close $fr

set fr [open "testing.txt" r]
set file_data [read $fr]
set data [split $file_data "\n"]
foreach line $data {
	set contents [split $line ":"]
	set complex [lindex $contents 0]
	set chain [lindex $contents 1]
	mol load pdb "$complex.pdb"
	#set command "\"protein and chain $chain or nucleic\" 1.4 1 1 $complex\_$chain\_RNA"
	#puts $command
	foreach radii $radius {
		getAllResSASA  "protein and chain $chain" $radii 1 1 "$complex\_$chain\_$radii"
	}
	puts $ind
	set ind [expr $ind+1]
}
close $fr
