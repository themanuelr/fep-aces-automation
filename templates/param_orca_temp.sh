#!/bin/bash

mol_names=("NEWCOMP")
start_struc=../00-compounds/pre_${newcomp}.mol2

for newcomp in "${mol_names[@]}";do
	echo Start $newcomp 01-opt
	echo $(date)
 	mkdir ${newcomp}
 	mkdir ${newcomp}/01-opt
 	cp ../../../../PROPS/FEP_inputs/AUTOMATIZATION/orca_prm/orca_opt_fileprep.py ${newcomp}/01-opt
 	cd ${newcomp}/01-opt;python orca_opt_fileprep.py ${newcomp};orca opt_${newcomp}.inp > opt_${newcomp}.out; cd ../../

	echo Start $newcomp 02-esp	       
       	echo $(date)
 	mkdir ${newcomp}/02-esp
 	cp ../../../../PROPS/FEP_inputs/AUTOMATIZATION/orca_prm/orca_esp_fileprep.py ${newcomp}/02-esp
 	cd ${newcomp}/02-esp;python orca_esp_fileprep.py ${newcomp}; orca esp_${newcomp}.inp > esp_${newcomp}.out ;cd ../../
	
        echo Start $newcomp 03-resp
        echo $(date)
 	mkdir ${newcomp}/03-resp
 	cp ../../../../PROPS/FEP_inputs/AUTOMATIZATION/orca_prm/orca_resp_fileprep.py ${newcomp}/03-resp
 	cd ${newcomp}/03-resp;python orca_resp_fileprep.py ${newcomp}
 	charge=$(grep "Sum of atomic charges:" ../02-esp/esp_${newcomp}.out | awk '{print $5}')
 	charge_int=$(echo "$charge/1" | bc)
 	pyresp_gen.py -i esp_${newcomp}.vpot -f1 resp1.in -f2 resp2.in -p chg -q $charge_int
 	py_resp.py -O -i resp1.in -o resp1.out -t q1 -e esp_${newcomp}.vpot -s resp1.esp
 	py_resp.py -O -i resp2.in -o resp2.out -t q2 -q q1 -e esp_${newcomp}.vpot -s resp2.esp
 	cd ../../

        echo Start $newcomp 04-antechamber 
        echo $(date)
	mkdir ${newcomp}/04-antechamber
	cp ../../../../PROPS/FEP_inputs/AUTOMATIZATION/orca_prm/orca_antechamber_fileprep.py ${newcomp}/04-antechamber
	cp ../../../../PROPS/FEP_inputs/AUTOMATIZATION/orca_prm/orca_antechamber_badchrg.py ${newcomp}/04-antechamber
	cd ${newcomp}/04-antechamber; python orca_antechamber_fileprep.py ${newcomp};
        charge=$(grep "Sum of atomic charges:" ../02-esp/esp_${newcomp}.out | awk '{print $5}')
        charge_int=$(echo "$charge/1" | bc)
	antechamber -i opt_${newcomp}.pdb -fi pdb -o badchrg_${newcomp}.mol2 -fo mol2 -c gas -nc ${charge_int} -rn ${newcomp} -at gaff2
	python orca_antechamber_badchrg.py ${newcomp}
	parmchk2 -i ${newcomp}.mol2 -f mol2 -o ${newcomp}.frcmod
	cp ${newcomp}.mol2 ../
	cp ${newcomp}.frcmod ../
	echo "_" in mol2:
	grep _ ${newcomp}.mol2
	cd ../../
done	

