#!/bin/sh
old_mol_name='OLD_MOL'
start_frame_path='STARTFRAMEPATH'
start_frame='STARTFRAMENAME'
mol_name='NEWMOL'
repn='REPN'
end_string="ENDSTRING"

##get tleap/cpptraj template
cp TEMPDIR/temp_getpdb.cpptraj getpdb.cpptraj
if [ "${repn}" = "unbound" ]; then
	cp TEMPDIR/temp_tleap.unbound.in tleap.in
	#edit cpptraj template
	sed -i "s/NAME/${mol_name}/g" getpdb.cpptraj
	#edit tleap template
	sed -i "s/OLDMOL/${old_mol_name}/g" tleap.in
	sed -i "s/MOL_NAME/${mol_name}/g" tleap.in
	
else
	read -r line < ../../../../${start_frame_path}
	XXX=$(echo "$line" | awk '{print $2}')
	YYY=$(echo "$line" | awk '{print $3}')
	ZZZ=$(echo "$line" | awk '{print $4}')
	#get starting frame
	cp ../../../../${start_frame_path} .
	#get target molecule PDB, shared atoms should have same XYZ coords as starting frame initial comp
	cp ../../../00-compounds/${mol_name}_cut.${repn}.pdb ${mol_name}.pdb
	cp TEMPDIR/temp_tleap.in tleap.in
	#edit cpptraj template
	sed -i "s/NAME/hNKCC1_${mol_name}/g" getpdb.cpptraj
	#fix and polish target molecule PDB for insertion
	sed -i "s/${old_mol_name}/${mol_name}/g" ${mol_name}.pdb
	sed -i "s/ X /   /g" ${mol_name}.pdb
	sed -i "s/END/TER/g" ${mol_name}.pdb
	#sed -i '1d' ${mol_name}.pdb
	#insert target molecule
	line_number=$(grep -n "$end_string" ${start_frame} | cut -d ':' -f 1)
	echo "New molecule inserted at line $line_number"
	sed -i "${line_number}r ${mol_name}.pdb" ${start_frame}
	#get full name of start frame
	file=$(find . -type f -name ${start_frame} -print -quit)
	full_name=$(basename "$file")
	#edit tleap template
	sed -i "s/OLDMOL/${old_mol_name}/g" tleap.in
	sed -i "s/MOL_NAME/${mol_name}/g" tleap.in
	sed -i "s/START_FRAME/${full_name}/g" tleap.in
	#sed -i "s/START_FRAME/${start_frame}/g" tleap.in
	sed -i "s/XXX/${XXX}/g" tleap.in
	sed -i "s/YYY/${YYY}/g" tleap.in
	sed -i "s/ZZZ/${ZZZ}/g" tleap.in

fi

