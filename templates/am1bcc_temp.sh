#!/bin/bash -e
mol_name=NEWCOMP
chrg=CHARGE

antechamber -i "../../00-compounds/$mol_name.pdb" -fi pdb -fo mol2 -o "$mol_name.mol2" -c bcc -nc $chrg -rn $mol_name -at gaff2
antechamber -i "sqm.pdb" -fi pdb -fo mol2 -o "$mol_name.mol2" -c bcc -nc $chrg -rn $mol_name -at gaff2

parmchk2 -i $mol_name.mol2 -f mol2 -o $mol_name.frcmod -s gaff2

exit 0

