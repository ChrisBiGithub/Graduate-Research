## 
## Get list of residues (use 'residue' and not 'resid' so we don't get 
## duplicate residues from unusual PDB files..) 
## 
set allsel [atomselect top all] 
set residlist [lsort -unique [$allsel get residue]] 
## 
## Make an atom selection, set the User field with the SASA value for 
## the selected atom 
## 
foreach r $residlist { 
  set sel [atomselect top "residue $r"] 
  set rsasa [measure sasa 1.4 $allsel -restrict $sel] 
  $sel set user $rsasa 
  $sel delete 
  set residueName [$r get resname]
  set residueID [$r get resid]
  puts "residue $residueID-$residueName, sasa: $rsasa" 
} 
## 
## change the "color by" and "trajectory" tab settings to color by SASA 
## 
mol modcolor 0 [molinfo top] User 
mol colupdate 0 [molinfo top] 1 
mol scaleminmax [molinfo top] 0 auto 