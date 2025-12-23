import os
import glob
import numpy as np

template_dir='templates'

def campaign_dir_setup(campaign_name='FEP_campaign'):
    DIRS=['00-compounds','01-mol_param','02-leap','03-md','04-analysis']
    os.makedirs(campaign_name)
    for dir in DIRS:
        os.makedirs(f'{campaign_name}/{dir}')
        if dir=='00-compounds':
            os.makedirs(f'{campaign_name}/{dir}/rep_struc')
    print('Copy your 3 starting structures into the path:')
    PWD=os.getcwd()
    print(f'{PWD}/{campaign_name}/00-compounds/rep_struc')
    print('Naming convention of starting representative structures is as follows:')
    print('OLDCOMP_repN_fNNNN.pdb')
    print(' OLDCOMP is the 3-letter code for the old compound')
    print(' repN is rep1, rep2 or rep3')
    print(' fNNNN is the frame number that was extracted from the OLDCOMP equilibrium MD')
    return

def extract_oldcomp_pdb(oldcomp,campaign_name='FEP_campaign',repN='rep1'):
    os.system(f'grep {oldcomp} {campaign_name}/00-compounds/rep_struc/{oldcomp}_{repN}*.pdb > {campaign_name}/00-compounds/{oldcomp}.pdb')

def get_newcomp_list(oldcomp,campaign_name='FEP_campaign'):
    PWD=f'{os.getcwd()}/{campaign_name}/00-compounds/'
    mol2_files=glob.glob(os.path.join(PWD,'*.mol2'))
    mol2_files=[file.split('/')[-1].split('_')[-1].split('.')[0] for file in mol2_files]
    print('This is the list of new compounds that will be used going forward. \nPlease recreate this variable if you wish to work with a smaller subset of the new compounds')
    print(mol2_files)
    return mol2_files

def highlight_newatoms(oldcomp,newcomp,PWD):
    #Outputs a list of the unique atoms in the original/oldcomp compound and the new compound
    #Even if atoms have the same atom type, they will be considered unique if their coordinates are not equal between compounds
    oldcomp_pdb=np.genfromtxt(f'{PWD}/{oldcomp}.pdb',dtype=str,skip_footer=1)
    newcomp_pdb=np.genfromtxt(f'{PWD}/{newcomp}.pdb',dtype=str,skip_footer=1,skip_header=1)
    ori_atomT=[]
    ori_coords=[]
    ori_dict=dict()
    newc_atomT=[]
    newc_coords=[]
    newc_dict=dict()
    for atom in oldcomp_pdb:
        ori_atomT.append(atom[2])
        ori_coords.append(atom[5:8])
        ori_dict[ori_atomT[-1]]=ori_coords[-1]
    for atom in newcomp_pdb:
        newc_atomT.append(atom[2])
        newc_coords.append(atom[6:9])
        newc_dict[newc_atomT[-1]]=newc_coords[-1]
        def get_atom_list(A_atomT,B_atomT,A_atomD,B_atomD):
            diff_atoms=[]
            for new_atom in A_atomT:
                if new_atom not in B_atomT:
                    diff_atoms.append(new_atom)
                elif new_atom in B_atomT:
                    are_equal=np.array_equal(B_atomD[new_atom],A_atomD[new_atom])
                    if not are_equal:
                        diff_atoms.append(new_atom)
            return diff_atoms
    new_atoms_SC=get_atom_list(newc_atomT,ori_atomT,newc_dict,ori_dict)
    ori_atoms_SC=get_atom_list(ori_atomT,newc_atomT,ori_dict,newc_dict)
    return new_atoms_SC,ori_atoms_SC

def generate_cutpdb(oldcomp,newcomp,campaign_name='FEP_campaign',reps=['rep1','rep2','rep3']):
    #Takes the output of highlight_newatoms() and creates a PDB file, for each replica, which includes exclusively the common core atoms between oldcomp and newcomp
    #Also, creates a logfile with the softcore atoms for each compound (this will be read for input file creation)
    ##Possible duplication issue - if re run will continue to add the softcore atoms to the logfile
    ##Interaction unknown when this file is read for input file creation
    PWD=f'{campaign_name}/00-compounds'
    new_atoms_SC,ori_atoms_SC=highlight_newatoms(oldcomp,newcomp,PWD)
    for rep in reps:
        fname=glob.glob(f'{PWD}/rep_struc/{oldcomp}_{rep}_f*.pdb')[0]
        with open(fname,'r') as pdb_file:
            with open(f'{PWD}/{newcomp}_cut.{rep}.pdb','w') as output_file:
                for atom in pdb_file:
                    if oldcomp in atom:
                        if atom.split()[2] not in ori_atoms_SC:
                            output_file.write(atom)
    log_file=f'{PWD}/SC_atoms.log'
    def fix_list(atom_list):
        fixed=''
        for atom in atom_list:
            fixed=fixed+f'{atom},'
        return fixed
    ori_amber=fix_list(ori_atoms_SC)
    new_amber=fix_list(new_atoms_SC)
    with open(log_file,'a') as log:
        output=f'_____\n{oldcomp} to {newcomp}\n{oldcomp} SC atoms {ori_amber}\n{newcomp} SC atoms {new_amber}\n_____\n'
        log.write(output)
    print(f'Cut PDBs for {newcomp} have been created')

def setup_param(param_type,newcomp_list,campaign_name='FEP_campaign',charge=0,execute_param=False):
    #This function will execute the appropiate function to setup or setup+execute compound parametrization, depending on user choice
    #orca takes about 8 hours, whilst am1bcc takes seconds-minutes; which is why orca is always executed off-script
    #am1bcc is executed sequentially, but the function can be told to not execute am1bcc, 
    ##only to write the script, allowing the user to run multiple instances in parallel through command line
    if param_type=='orca':
        setup_param_orca(newcomp_list,campaign_name)
    elif param_type=='am1bcc':
        setup_param_am1bcc(newcomp_list,campaign_name,charge,execute_param)


def setup_param_orca(newcomp_list,campaign_name):
    #Copies the relevant template script and changes necessary variables to parametrize using orca
    print(f'Go to:\n{campaign_name}/01-mol_param/\nand execute the script:')
    for newcomp in newcomp_list:
        os.system(f'cp {template_dir}/param_orca_temp.sh {campaign_name}/01-mol_param/param_orca_{newcomp}.sh')
        os.system(f"sed -i 's/NEWCOMP/{newcomp}/g' {campaign_name}/01-mol_param/param_orca_{newcomp}.sh")
        print(f'./param_orca_{newcomp}.sh')

def setup_param_am1bcc(newcomp_list,campaign_name,charge,execute_param):
    #creates a directory per compound, copies the relevant template script and changes the necessary variables to parametrize using am1bcc
    if not execute_param:print(f'Go to:\n{campaign_name}/01-mol_param/COMPOUND\nand execute the script:')
    for newcomp in newcomp_list:
        os.makedirs(f'{campaign_name}/01-mol_param/{newcomp}')
        os.system(f'cp {template_dir}/am1bcc_temp.sh {campaign_name}/01-mol_param/{newcomp}/am1bcc_{newcomp}.sh')
        os.system(f"sed -i 's/NEWCOMP/{newcomp}/g' {campaign_name}/01-mol_param/{newcomp}/am1bcc_{newcomp}.sh")
        os.system(f"sed -i 's/CHARGE/{charge}/g' {campaign_name}/01-mol_param/{newcomp}/am1bcc_{newcomp}.sh")
        if execute_param:
            os.system(f'cd {campaign_name}/01-mol_param/{newcomp}/;./am1bcc_{newcomp}.sh > am1bcc_out.log ;cd ../../../')
            print(f'Compound {newcomp} parametrized')
        else:
            print(f'./am1bcc_{newcomp}.sh')

def get_oldcomp_endstring(oldcomp,replica,campaign_name='FEP_campaign'):
    #returns the line which contains both the resname of oldcomp and TER flag
    #this should be the line after which newcomp will be inserted
    if replica=='unbound':replica='rep1'
    fname=glob.glob(f'{campaign_name}/00-compounds/rep_struc/{oldcomp}_{replica}_*.pdb')[0]
    with open(fname,'r') as f:
        for line in f:
            if 'TER' in line:
                if oldcomp in line:
                    return line

def tleap_setup_organize_dirs(newcomp_list,oldcomp,campaign_name='FEP_campaign',replicas=['rep1','rep2','rep3','unbound']):
    #will create a directory for each compound, and a subdirectory for each replica
    #will also copy every necessary file to be able to run tleap and create .top and .rst
    #will copy and edit get_files.sh from templates, which will prepare the tleap.in file and polish all necessary files
    for newcomp in newcomp_list:
        os.makedirs(f'{campaign_name}/02-leap/{newcomp}')
        for rep in replicas:
            print(f'Setting up: {newcomp}-{rep}')
            endstring=get_oldcomp_endstring(oldcomp,rep)[:-1]
            if rep !='unbound':startframe_path=glob.glob(f'{campaign_name}/00-compounds/rep_struc/{oldcomp}_{rep}_*.pdb')[0]
            else:startframe_path=glob.glob(f'{campaign_name}/00-compounds/rep_struc/{oldcomp}_{replicas[0]}_*.pdb')[0]
            dirname=f'{campaign_name}/02-leap/{newcomp}/{rep}'
            os.makedirs(dirname)
            os.system(f'cp {template_dir}/get_files_temp.sh {dirname}/get_files.sh')
            os.system(f'sed -i "s/OLD_MOL/{oldcomp}/g" {dirname}/get_files.sh')
            os.system(f'sed -i "s=STARTFRAMEPATH={startframe_path}=g" {dirname}/get_files.sh')
            startframe_name=startframe_path.split('/')[-1]
            os.system(f'sed -i "s=STARTFRAMENAME={startframe_name}=g" {dirname}/get_files.sh')
            os.system(f'sed -i "s/NEWMOL/{newcomp}/g" {dirname}/get_files.sh')
            os.system(f'sed -i "s/REPN/{rep}/g" {dirname}/get_files.sh')
            os.system(f'sed -i "s/ENDSTRING/{endstring}/g" {dirname}/get_files.sh')
            os.system(f'sed -i "s=TEMPDIR={template_dir}=g" {dirname}/get_files.sh')
            os.system(f'cd {dirname};./get_files.sh; cd ../../../../')

def tleap_execute(newcomp_list,oldcomp,campaign_name='FEP_campaing',replicas=['rep1','rep2','rep3','unbound']):
    #will go into each compound-replica directory and execute the tleap file and the extract pdb cpptraj file
    #in the end, each compound-replica directory will contain the .top, .rst, and .pdb files
    for newcomp in newcomp_list:
        for rep in replicas:
            print(f'Executing {newcomp}-{rep}')
            dirname=f'{campaign_name}/02-leap/{newcomp}/{rep}'
            if rep=='unbound':
                os.system(f'cd {dirname}; grep {oldcomp} ../rep1/hNKCC1_{newcomp}.pdb > hNKCC1_{newcomp}.pdb; cd ../../../../')
                os.system(f'cd {dirname}; grep {newcomp} ../rep1/hNKCC1_{newcomp}.pdb >> hNKCC1_{newcomp}.pdb; cd ../../../../')
            os.system(f'cd {dirname}; tleap -f tleap.in > tleap.log; cd ../../../../')
            os.system(f'cd {dirname}; cpptraj getpdb.cpptraj > getpdb.log ; cd ../../../../')

def mdfep_organize_dirs(newcomp_list,campaign_name,replicas=['rep1','rep2','rep3','unbound']):
    #will create a directory for every compound, then a subdirectory for both aces and stdr FEP runs, and within each of these a subdirectory for each replica
    #will then populate the replica subdirectories with the relevant .top and .rst files
    #will also rename all but unbound replica restart to 00-* 
    for newcomp in newcomp_list:
        newcomp_dir=f'{campaign_name}/03-md/{newcomp}'
        os.makedirs(newcomp_dir)
        os.makedirs(f'{newcomp_dir}/stdr')
        os.makedirs(f'{newcomp_dir}/aces')
        for rep in replicas:
            os.makedirs(f'{newcomp_dir}/stdr/{rep}')
            os.makedirs(f'{newcomp_dir}/aces/{rep}')
            if rep=='unbound':
                fname=newcomp
                os.makedirs(f'{newcomp_dir}/stdr/{rep}/mdsolv')
            else:
                fname=f'hNKCC1_{newcomp}'
            os.system(f'cp {campaign_name}/02-leap/{newcomp}/{rep}/{fname}.* {newcomp_dir}/stdr/{rep}')
            os.system(f'cp {campaign_name}/02-leap/{newcomp}/{rep}/{fname}.top {newcomp_dir}/aces/{rep}')
            if rep!='unbound':
                os.system(f'mv {newcomp_dir}/stdr/{rep}/{fname}.rst {newcomp_dir}/stdr/{rep}/00-{fname}.md.rst')
            
def extract_props(PROPS,fep_type):
    #extract the values from the PROPS dictionary and turn it into individual variables
    NTPR=PROPS['NTPR']
    NSTLIM=PROPS['NSTLIM']
    DT=PROPS['DT']
    NTC=PROPS['NTC']
    GTI_BAT_SC=PROPS['GTI_BAT_SC']
    if fep_type=='aces':
        NUMEXCHG=PROPS['NUMEXCHG']
        return NTPR,NSTLIM,DT,NUMEXCHG,NTC,GTI_BAT_SC
    else:
        return NTPR,NSTLIM,DT,NTC,GTI_BAT_SC

def extract_FEP_masks(FEP_masks):
    #extract the values from the FEP_masks dictionary and turn it into individual variables
    ini_rname=FEP_masks['ini_rname']
    tgt_rname=FEP_masks['tgt_rname']
    ini_scmask=FEP_masks['ini_scmask']
    tgt_scmask=FEP_masks['tgt_scmask']
    return ini_rname,tgt_rname,ini_scmask,tgt_scmask

def extract_clambdas(gap):
    #return a list of lambdas, if gap is a list, returns gap. If gap is just a float, representing the gap between lamdba windows, return a list of values between 0 and 1, separated by a gap
    if isinstance(gap,float):
        clambdas=np.arange(0,1+gap,gap)
    elif isinstance(gap,list):
        clambdas=gap
    return clambdas

def produce_input_file(fep_type,gap,out_pwd,PROPS,FEP_masks,skip=[],gap_ROUND=2,dividedWindow=False):
    #Write all necessary FEP md.in files
    temp_min=f'{template_dir}/temp-FEP{fep_type}-min.in'
    temp_md=f'{template_dir}/temp-FEP{fep_type}-md.in'
    clambdas=extract_clambdas(gap)
    all_lambda=str([round(i,gap_ROUND) for i in clambdas]).replace(' ','').replace('[','').replace(']','')
    if fep_type=='stdr':NTPR,NSTLIM,DT,NTC,GTI_BAT_SC=extract_props(PROPS,fep_type)
    if fep_type=='aces':NTPR,NSTLIM,DT,NUMEXCHG,NTC,GTI_BAT_SC=extract_props(PROPS,fep_type)
    ini_rname,tgt_rname,ini_scmask,tgt_scmask=extract_FEP_masks(FEP_masks)
    
    for ndx,clambda in enumerate(clambdas):
        num=ndx+1
        if num in skip:
            continue
        os.system(f'cp {temp_md} {out_pwd}/{num:02d}-md.in')
        if fep_type=='stdr':os.system(f'cp {temp_min} {out_pwd}/{num:02d}-min.in')
    
        os.system(f'sed -i "s/NTPR/{NTPR}/g" {out_pwd}/{num:02d}-md.in')
        os.system(f'sed -i "s/NSTLIM/{NSTLIM}/g" {out_pwd}/{num:02d}-md.in')
        os.system(f'sed -i "s/DT/{DT}/g" {out_pwd}/{num:02d}-md.in')
        if fep_type=='aces':os.system(f'sed -i "s/NUMEXCHG/{NUMEXCHG}/g" {out_pwd}/{num:02d}-md.in')
        os.system(f'sed -i "s/NTC/{NTC}/g" {out_pwd}/{num:02d}-md.in')
        os.system(f'sed -i "s/GTI_BAT_SC/{GTI_BAT_SC}/g" {out_pwd}/{num:02d}-md.in')
        os.system(f'sed -i "s/CLAMBDA/{clambda:.2f}/g" {out_pwd}/{num:02d}-md.in')
        os.system(f'sed -i "s/TIMASK1/:{ini_rname}/g" {out_pwd}/{num:02d}-md.in')
        os.system(f'sed -i "s/TIMASK2/:{tgt_rname}/g" {out_pwd}/{num:02d}-md.in')
        os.system(f'sed -i "s/SCMASK1/:{ini_rname}@{ini_scmask}/g" {out_pwd}/{num:02d}-md.in')
        os.system(f'sed -i "s/SCMASK2/:{tgt_rname}@{tgt_scmask}/g" {out_pwd}/{num:02d}-md.in')
        os.system(f'sed -i "s/N_WINDOWS/{len(clambdas)}/g" {out_pwd}/{num:02d}-md.in')
        os.system(f'sed -i "s/ALL_WINDOWS/{all_lambda}/g" {out_pwd}/{num:02d}-md.in')

        if fep_type=='stdr':
            os.system(f'sed -i "s/NTC/{NTC}/g" {out_pwd}/{num:02d}-min.in')
            os.system(f'sed -i "s/GTI_BAT_SC/{GTI_BAT_SC}/g" {out_pwd}/{num:02d}-min.in')       
            os.system(f'sed -i "s/CLAMBDA/{clambda:.2f}/g" {out_pwd}/{num:02d}-min.in')
            os.system(f'sed -i "s/TIMASK1/:{ini_rname}/g" {out_pwd}/{num:02d}-min.in')
            os.system(f'sed -i "s/TIMASK2/:{tgt_rname}/g" {out_pwd}/{num:02d}-min.in')
            os.system(f'sed -i "s/SCMASK1/:{ini_rname}@{ini_scmask}/g" {out_pwd}/{num:02d}-min.in')
            os.system(f'sed -i "s/SCMASK2/:{tgt_rname}@{tgt_scmask}/g" {out_pwd}/{num:02d}-min.in')
    if out_pwd.split('/')[-1]=='unbound' and fep_type=='stdr':
        temp_1min=f'{template_dir}/temp-FEP{fep_type}-01-min.in'
        temp_2md=f'{template_dir}/temp-FEP{fep_type}-02-md.in'
        temp_3md=f'{template_dir}/temp-FEP{fep_type}-03-md.in'
        temp_4md=f'{template_dir}/temp-FEP{fep_type}-04-md.in'
        os.system(f'cp {temp_1min} {out_pwd}/mdsolv/01-min.in')
        os.system(f'cp {temp_2md} {out_pwd}/mdsolv/02-md.in')
        os.system(f'cp {temp_3md} {out_pwd}/mdsolv/03-md.in')
        os.system(f'cp {temp_4md} {out_pwd}/mdsolv/04-md.in')
        for file in [f'{out_pwd}/mdsolv/01-min.in',f'{out_pwd}/mdsolv/02-md.in',f'{out_pwd}/mdsolv/03-md.in',f'{out_pwd}/mdsolv/04-md.in']:
            os.system(f'sed -i "s/NTC/{NTC}/g" {file}')
            os.system(f'sed -i "s/DT/{DT}/g" {file}')
            os.system(f'sed -i "s/TIMASK1/:{ini_rname}/g" {file}')
            os.system(f'sed -i "s/TIMASK2/:{tgt_rname}/g" {file}')
            os.system(f'sed -i "s/SCMASK1/:{ini_rname}@{ini_scmask}/g" {file}')
            os.system(f'sed -i "s/SCMASK2/:{tgt_rname}@{tgt_scmask}/g" {file}')

def create_groupfile(comp,fep_type,campaign_name,REPS=['rep1','rep2','rep3'],divWindow=False):
    #Write the groupfile for the aces section of the transformation
    mdins=glob.glob(f'{campaign_name}/03-md/{comp}/{fep_type}/rep1/*md.in')
    num_list=[item.split('/')[-1].split('-')[0] for item in mdins]
    num_list.sort()
    
    if REPS[0]=='unbound':
        fname=f'{comp}'
    else:
        fname=f'hNKCC1_{comp}'
    TOP=f'{fname}.top'
    TYPE=f'FEP{fep_type}'
    sNAME=f'{fname}.{TYPE}'
    rst00_name=f'{fname}.FEPstdr'
    loops=divWindow+1
    for rep in REPS:
        for win_divi in range(loops):
            NAME=f'{sNAME}.{win_divi}.md'
            fname=f'{campaign_name}/03-md/{comp}/{fep_type}/{rep}/groupfile_windowpart_{win_divi}'
            with open(fname,'w') as groupfile:
                for num in num_list:
                    if win_divi==0:
                        line=f'-O -p {TOP} -i {num}-md.in -c ../../stdr/{rep}/{num}-{rst00_name}.md.rst -r {num}-{NAME}.rst -x {num}-{NAME}.crd -o {num}-{NAME}.out -inf {num}-{NAME}.mdinfo \n'
                    elif win_divi==1:
                        prev_NAME=f'{sNAME}.0.md'
                        line=f'-O -p {TOP} -i {num}-md.in -c {num}-{prev_NAME}.rst -r {num}-{NAME}.rst -x {num}-{NAME}.crd -o {num}-{NAME}.out -inf {num}-{NAME}.mdinfo \n'
                    groupfile.write(line)

def get_FEP_masks(comp,old_comp,campaign_name):
    ############################################
    #Redundant code, if it returns individual variables for each mask, we can remove the function
    #extract_FEP_masks()
    ############################################
    #read the softcore atoms logfile and turn it into a dictionary
    fname=f'{campaign_name}/00-compounds/SC_atoms.log'
    with open(fname,'r') as logfile:
        check=False
        for line in logfile:
            if check:
                if line.startswith(old_comp):
                    old_comp_SC=line.split(' ')[-1][:-2]
                elif line.startswith(comp):
                    comp_SC=line.split(' ')[-1][:-2]
                    check=False
                    break
            if comp in line and old_comp in line:
                check=True
    FEP_masks=dict()
    FEP_masks['ini_rname']=old_comp
    FEP_masks['tgt_rname']=comp
    FEP_masks['ini_scmask']=old_comp_SC
    FEP_masks['tgt_scmask']=comp_SC
    return FEP_masks            

def prepare_FEP_files(comp,old_comp,campaign_name,lambda_gap,PROPS_stdr,PROPS_acesR,PROPS_acesUB,gap_ROUND=3):
    #combins all functions relating to FEP input files, and generates all files necessary, for each replica
    FEP_masks=get_FEP_masks(comp,old_comp,campaign_name)

    fep_type='stdr'
    main_PWD=f'{campaign_name}/03-md/{comp}/{fep_type}'
    R1=f'{main_PWD}/rep1'
    R2=f'{main_PWD}/rep2'
    R3=f'{main_PWD}/rep3'
    UB=f'{main_PWD}/unbound'
    OUT_PWD_stdr=[R1,R2,R3,UB]
    for out_pwd in OUT_PWD_stdr:
        produce_input_file(fep_type,lambda_gap,out_pwd,PROPS_stdr,FEP_masks,skip=[],gap_ROUND=gap_ROUND)

    fep_type='aces'
    main_PWD=f'{campaign_name}/03-md/{comp}/{fep_type}'
    R1=f'{main_PWD}/rep1'
    R2=f'{main_PWD}/rep2'
    R3=f'{main_PWD}/rep3'
    UB=f'{main_PWD}/unbound'
    OUT_PWD_acesR=[R1,R2,R3]
    for out_pwd in OUT_PWD_acesR:
        produce_input_file(fep_type,lambda_gap,out_pwd,PROPS_acesR,FEP_masks,skip=[],gap_ROUND=gap_ROUND)
    OUT_PWD_acesUB=[UB]
    for out_pwd in OUT_PWD_acesUB:
        produce_input_file(fep_type,lambda_gap,out_pwd,PROPS_acesUB,FEP_masks,skip=[],gap_ROUND=gap_ROUND)
    create_groupfile(comp,'aces',campaign_name,divWindow=False)
    create_groupfile(comp,'aces',campaign_name,REPS=['unbound'],divWindow=False)

def distribute_submission_jobs(newcomp_list,campaign_name,windows=25,groupfile='groupfile_windowpart_0',replicas=['rep1','rep2','rep3','unbound']):
    fep_types=['stdr','aces']
    for fep_type in fep_types:
            for newcomp in newcomp_list:
                for rep in replicas:
                    last_char=rep[-1]
                    fname=f'{campaign_name}/03-md/{newcomp}/{fep_type}/{rep}/job_FEP{fep_type}_1.sh'
                    os.system(f'cp {template_dir}/temp-job_FEP{fep_type}_1.sh {fname}')
                    os.system(f'sed -i "s/COMP/{newcomp}/g" {fname}')
                    os.system(f'sed -i "s/REPN/{last_char}/g" {fname}')
                    os.system(f'sed -i "s/NREP/{rep}/g" {fname}')
                    os.system(f'sed -i "s/WIN/{windows}/g" {fname}')
                    if fep_type=='aces':os.system(f'sed -i "s/GROUPFILE/{groupfile}/g" {fname}')
                    elif fep_type=='stdr':
                        if rep=='unbound':
                            os.system(f'sed -i "s/PREFIX//g" {fname}')
                            os.system(f'cp {template_dir}/temp-job_mdsolv_1.sh {campaign_name}/03-md/{newcomp}/{fep_type}/{rep}/mdsolv/job_mdsolv_1.sh')
                            os.system(f'sed -i "s/COMP/{newcomp}/g" {campaign_name}/03-md/{newcomp}/{fep_type}/{rep}/mdsolv/job_mdsolv_1.sh')
                        elif rep!='unbound':
                            os.system(f'sed -i "s/PREFIX/hNKCC1_/g" {fname}')
                if fep_type=='aces':
                    print(f'Copy the complete 03-md/{newcomp} directory into the HPC')
                    print(f'In the HPC, go to:\n {newcomp}/stdr/unbound/mdsolv\n and submit: \nqsub job_mdsolv_1.sh')
