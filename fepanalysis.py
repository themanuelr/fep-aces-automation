import numpy as np
import os
import glob
import matplotlib.pyplot as plt
import scipy as sp
import pandas as pd
Z=18

#from rdkit import Chem
#from rdkit.Chem import AllChem
#from rdkit.Chem import Draw
#from rdkit.Chem import rdDepictor
#from rdkit.Chem import rdFMCS
#from rdkit.Chem.Draw import IPythonConsole
#from rdkit.Chem.Draw import rdMolDraw2D
#from cairosvg import svg2png
#rdkit part of the script is not fully implemented, because of how AmberTools environment is set up

#overall figure saving is broken - needs directory setup and figure routing

#DDG processing and calc##############################################################################################################
def process_lines(lines):
    #takes in the lines from the *out file and extracts dU/dlambda
    #this is done one lambda window at a time
    #also reads in lambda from the file
    #calculates the mean dU/dlambda and the STD
    dU_dl=[]
    Lambda=None
    c=0
    for line in lines:
        c+=1
        if 'Total dU/dl' in line:
            ls_line=line.split(' ')
            ls_line=list(filter(None,ls_line))
            #dU_dl.append(float(ls_line[6]))
            try:dU_dl.append(float(ls_line[6]))
            except:
                try:dU_dl.append(float(ls_line[5].split(':')[-1]))
                except:
                    print(ls_line)
                    print(c)
            while Lambda is None:
                Lambda=float(ls_line[2])
    mean_dU_dl=np.mean(dU_dl)
    std_dU_dl=np.std(dU_dl)
    return Lambda,mean_dU_dl,std_dU_dl,dU_dl

def get_files(comp,rep,fep_type,campaign_name='FEP_campaign'):
    #creates a list of every md.out file of a specific replica
    files=os.listdir(f'{campaign_name}/03-md/{comp}/{fep_type}/{rep}/')
    files=[f'{campaign_name}/03-md/{comp}/{fep_type}/{rep}/{file}' for file in files if file.endswith('md.out')]
    files.sort()
    return files
            
def remove_extremes(data):
    #removes first and last data point - since they should always be zero
    #function might be currently unused
    return data[1:-1]

def calc_running_mean(data,window=0.03):
    #calculates running mean on a window of len(data)*window
    window=int(window*len(data))
    return np.convolve(data,np.ones(window)/window,mode='valid')

def get_data(comp,rep,fep_type,campaign_name='FEP_campaign'):
    #collects all dU/dlambda data for all 25 windows and unifies it
    LAMBDAS=[]
    MEAN_dU_dl=[]
    STD_dU_dl=[]
    ALL_dU_dl=[]
    files=get_files(comp,rep,fep_type,campaign_name)
    for fname in files:
        with open(fname,'r') as f: lines=f.readlines()
        Lambda,mean_dpoints,std_dpoints,all_du_dl=process_lines(lines)
        ALL_dU_dl.append(np.array(all_du_dl))
        LAMBDAS.append(Lambda)
        MEAN_dU_dl.append(mean_dpoints)
        STD_dU_dl.append(std_dpoints)
    return LAMBDAS,MEAN_dU_dl,STD_dU_dl,ALL_dU_dl

def calc_dG(comp,rep,fep_type,LAMBDAS=None,MEAN_dU_dl=None):
    #calculates deltaG by getting the area under the curve dU/dlambda using trapezoids. does not require mean_dU_dl to have been previously calculated.
    if LAMBDAS is None or MEAN_dU_dl is None:
        LAMBDAS,MEAN_dU_dl,_,_=get_data(comp,rep,fep_type)
        #LAMBDAS,MEAN_dU_dl=remove_extremes(LAMBDAS),remove_extremes(MEAN_dU_dl)
        LAMBDAS,MEAN_dU_dl=np.array(LAMBDAS),np.array(MEAN_dU_dl)
    dG=sp.integrate.trapezoid(MEAN_dU_dl,LAMBDAS)
    return dG

def remove_cutoff(data,cutoff):
    #removes a certain chunk of data at the beginning of solution
    #used for simulation validation
    chunk=int(len(data)*cutoff)
    return data[chunk:]
#########################################################################################################################################

#DDG plotting and outputting#############################################################################################################
def plot_curve(comp,rep,fep_type,save=False):
    #plots the dU/dlamda curve of a single replica
    LAMBDAS,MEAN_dU_dl,STD_dU_dl,_=get_data(comp,rep,fep_type)
    #LAMBDAS=remove_extremes(LAMBDAS)
    #MEAN_dU_dl=remove_extremes(MEAN_dU_dl)
    #STD_dU_dl=remove_extremes(STD_dU_dl)
    fig,ax=plt.subplots(figsize=(8,6))
    ax.errorbar(LAMBDAS,MEAN_dU_dl,yerr=STD_dU_dl,fmt='o',c='k',capsize=3,ls='-',lw=1)
    ax.set_xlabel('λ',fontsize=Z)
    ax.set_ylabel('∂V/∂λ',fontsize=Z)
    ax.tick_params(axis='both',labelsize=Z-2)
    ax.set_xticks(np.arange(0,1.1,0.1))
    if save:
        plt.savefig('dUdl-vs-l.png',dpi=300)
    plt.show()
    plt.close
def plot_curves(comp,REPS,fep_type,save=False):
    #plots the dU/dlambda curve for all replicas, including unbound
    curve_clrs={'rep1':'k','rep2':'dimgray','rep3':'darkgray','unbound':'lightskyblue'}
    LAMBDAS,MEAN_dU_dl,STD_dU_dl=dict(),dict(),dict()
    fig,ax=plt.subplots(figsize=(6,4))
    for rep in REPS:
        LAMBDAS[rep],MEAN_dU_dl[rep],STD_dU_dl[rep],_=get_data(comp,rep,fep_type)
        ax.errorbar(LAMBDAS[rep],MEAN_dU_dl[rep],yerr=STD_dU_dl[rep],fmt='o',c=curve_clrs[rep],capsize=3,ls='-',ms=2,lw=1,label=rep)
    ax.set_xlabel('λ',fontsize=Z)
    ax.set_ylabel('∂V/∂λ',fontsize=Z)
    ax.tick_params(axis='both',labelsize=Z-2)
    ax.set_xticks(np.arange(0,1.1,0.1))
    ax.legend()
    if save:
        plt.savefig('dUdl-vs-l.png',dpi=300)
    plt.show()
    plt.close

def plot_timeseries(comp,rep,fep_type,cutoff=None,save=False):
    #Plots the time evolution of the dU/dlambda for every lambda window
    LAMBDAS,MEAN_dU_dl,STD_dU_dl,ALL_dU_dl=get_data(comp,rep,fep_type)
    LAMBDAS,MEAN_dU_dl,STD_dU_dl,ALL_dU_dl=np.array(LAMBDAS),np.array(MEAN_dU_dl),np.array(STD_dU_dl),np.array(ALL_dU_dl)
    fig,axs=plt.subplots(len(LAMBDAS),1,figsize=(8,1*len(LAMBDAS)))
    for ndx,ax in enumerate(axs):
        X=np.linspace(0,1,len(ALL_dU_dl[ndx]))
        if cutoff is not None:
            if ndx==0 or ndx==(len(LAMBDAS)-1):
                ax.vlines(cutoff,-0.05,0.05,ls='--',lw=1,color='r')
            ax.vlines(cutoff,np.min(ALL_dU_dl[ndx]),np.max(ALL_dU_dl[ndx]),ls='--',lw=1,color='r')
        
        ax.plot(X,ALL_dU_dl[ndx],c='k',zorder=11)
        #ax.plot(XRA,RA,c='r',ls='-',lw=1)
        ax.hlines(MEAN_dU_dl[ndx],X[0],X[-1],ls='-',lw=1.5,color='royalblue')
        ax.hlines([MEAN_dU_dl[ndx]+STD_dU_dl[ndx],MEAN_dU_dl[ndx]-STD_dU_dl[ndx]],X[0],X[-1],ls='--',lw=1,color='royalblue')
        ax.fill_between(X,MEAN_dU_dl[ndx]+STD_dU_dl[ndx],MEAN_dU_dl[ndx]-STD_dU_dl[ndx],color='royalblue',alpha=0.3)

        ax.set_ylabel(f'λ={LAMBDAS[ndx]:.2f}\n∂V/∂λ')
        ax.set_xlim(X[0],X[-1])
        if ndx!=len(LAMBDAS)-1:
            ax.set_xticklabels([])
            #ax.set_xticks([])
        if ndx==len(LAMBDAS)-1:
            ax.set_xlabel('Simulation Time',fontsize=Z)
            ax.tick_params(axis='x',labelsize=Z-2)
    if save:
        plt.savefig('timeseries.png')
    plt.tight_layout()
    plt.show()
    plt.close()

def plot_convergence(comp,rep,fep_type,cutoff=None,plot_curve_chunks=False,save=False):
    #plots the convergence of dG for a specific replica
    #calculates the dG for progressively longer chunks of simulation:
    #   first 10%, first 20%, and so on.
    #an initial chunk taken as equilibration can be removed from calculation by using cutoff
    LAMBDAS,_,_,ALL_dU_dl=get_data(comp,rep,fep_type)
    LAMBDAS=np.array(LAMBDAS)
    original_length=len(ALL_dU_dl[0])
    if cutoff is not None:
        try:
            ALL_dU_dl=np.array([remove_cutoff(data,cutoff) for data in ALL_dU_dl])
            title=f'Convergence of ΔG - removed initial {cutoff*100}%'
        except: title='Convergence of ΔG - full trajectory'
    else: title='Convergence of ΔG - full trajectory'
    sim_chunks=np.linspace(0,1,11)[1:]
    deltaG={}
    for chunk_perc in sim_chunks:
        chunk=int(original_length*chunk_perc)
        means=[np.mean(data[:chunk]) for data in ALL_dU_dl]
        #dG=calc_dG(comp,rep,fep_type,remove_extremes(LAMBDAS),remove_extremes(means))
        dG=calc_dG(comp,rep,fep_type,LAMBDAS,means)
       
        deltaG[chunk_perc]=dG
        if plot_curve_chunks:
            stds=[np.std(data[:chunk]) for data in ALL_dU_dl]
            fig,ax=plt.subplots(figsize=(8,6))
            ax.errorbar(LAMBDAS,means,yerr=stds,fmt='o',c='k',capsize=3,ls='-',lw=1)
            ax.set_xlabel('λ',fontsize=Z)
            ax.set_ylabel('∂V/∂λ',fontsize=Z)
            ax.tick_params(axis='both',labelsize=Z-2)
            ax.set_xticks(np.arange(0,1.1,0.1))
            ax.set_title(f'First {chunk} steps')
            plt.show()
            plt.close()
    dGs=[deltaG[chunk] for chunk in sim_chunks]
    fig,ax=plt.subplots(figsize=(6,4))
    ax.plot(sim_chunks,dGs,ls='-',marker='o',c='k')
    #add text to the points
    for i,txt in enumerate(dGs):
        ax.annotate(f'{txt:.2f}',(sim_chunks[i],dGs[i]),textcoords='offset points',xytext=(0,10),ha='center')
    ax.set_title(f'{title}',fontsize=Z)
    ax.set_xlabel('Fraction of Simulation',fontsize=Z)
    ax.set_ylabel('ΔG (kcal/mol)',fontsize=Z)
    ax.tick_params(axis='both',labelsize=Z-4)
    ax.set_xticks(np.arange(0,1.1,0.1))
    if save:
        if cutoff is None:
            cutoff=0
        plt.savefig(f'convergence_cut_{cutoff}.png',dpi=300)
    plt.show()
    plt.close


#def get_FEPcompound_grid(names,save=False,core_align=False,campaign_name='FEP_campaign'):
#    #draws the compounds based on their pre_COMP.mol2 structure
#    mol2_files=[f'{campaign_name}/00-compounds/pre_{comp}.mol2' for comp in names]
#    mols=[Chem.rdmolfiles.MolFromMol2File(mol) for mol in mol2_files]
#    for mol in mols:AllChem.Compute2DCoords(mol)
#    rdDepictor.SetPreferCoordGen(True)
#
#    if core_align:
#        core_MCS=rdFMCS.FindMCS(mols)
#        core=Chem.MolFromSmarts(core_MCS.smartsString)
#        rdDepictor.Compute2DCoords(core)
#        for mol in mols:
#            _=rdDepictor.GenerateDepictionMatching2DStructure(mol,core)
#    if len(mols)<5:MPR=len(mols)
#    else:MPR=5
#    img=Draw.MolsToGridImage(mols,molsPerRow=MPR,subImgSize=(300,300),legends=names,useSVG=True)
#   
#    if save:
#        name=f'{names[0][0]}_compounds'
#        svg2png(bytestring=img.data,write_to=f'{name}.png')
#
#    return img
    
def get_ddG_table(COMPS,REPS,fep_type='aces'):
    #calculates all ddG, creates and formats a table for easy understanding
    DG_OUT=dict()
    for comp in COMPS:
        DG_OUT[comp]=dict()
        for rep in REPS:
            dG=calc_dG(comp,rep,fep_type)
            DG_OUT[comp][rep]=float(dG)
    df=pd.DataFrame(DG_OUT)
    for column in df.columns:
        df[column]=df[column]-df.loc['unbound',column]
    df=df.drop('unbound')
    mean_row=[]
    for column in df.columns:
        mean=round(df[column].mean(),1)
        std=round(df[column].std(),1)
        if mean>0:mean=f'+{mean}'
        else:mean=str(mean)
        mean_std=f'{mean} ({std})'
        mean_row.append(mean_std)
    df.loc['mean']=mean_row
    df=df.sort_index()
    def format_row(val):
        if val>0:val=f'+{val}'
        else:val=str(val)
        return val
    for row in df.index:
        if row=='mean':continue
        df.loc[row]=df.loc[row].apply(lambda x: round(x,1))
        df.loc[row]=df.loc[row].apply(lambda x: format_row(x))

        #df.style.format({row:"{0:+g}"})
    print(df)
    return df
def get_all_plots(COMPS,REPS,save=False,fep_type='aces',cutoff=0.05):
    #plots curves and convergence for all compounds and all replicas
    for comp in COMPS:
        print(comp)
        plot_curves(comp,REPS,fep_type,save=save)
        for rep in REPS:
            print(comp,rep)
            #plot_curve(comp,rep,fep_type,save=save)
            plot_convergence(comp,rep,fep_type,cutoff=cutoff,save=save)
            print('-'*20)
#########################################################################################################################################

#Traj conversion#########################################################################################################################
def create_comp_dir(comp,campaign_name='FEP_campaign'):
    #create a directory for the compound in the 04-analysis directory, if one doesn't exist
    dirname=f'{campaign_name}/04-analysis/{comp}'
    if not os.path.exists(dirname):
        os.mkdir(dirname)

def create_traj_dir(comp,rep,fep_type,campaign_name='FEP_campaign'):
    #create a directory where to store the xtc trajectories
    MD_PWD=f'{campaign_name}/03-md/{comp}/{fep_type}/{rep}'
    create_comp_dir(comp,campaign_name)
    trajDIR=f'{campaign_name}/04-analysis/{comp}/{comp}_{rep}_traj'
    if not os.path.exists(trajDIR):
        os.mkdir(trajDIR)
        TOP_PWD=MD_PWD+f'/hNKCC1_{comp}.top'
        os.system(f'cp {TOP_PWD} {trajDIR}')

def get_crd_numbers(comp,rep,fep_type,campaign_name='FEP_campaign'):
    #make a list of all the crd filenames and numbers
    MD_PWD=f'{campaign_name}/03-md/{comp}/{fep_type}/{rep}'
    crd_files=glob.glob(MD_PWD+'/*md.crd')
    crd_files.sort()
    crd_numbers=[crd.split('/')[-1].split('.')[0].split('-')[0] for crd in crd_files]
    crd_files=np.array(crd_files)
    return crd_numbers,crd_files

def check_pre_converted(crd_nums,comp,rep,campaign_name='FEP_campaign'):
    #check and create a mask of which trajectories have already been converted
    trajDIR=f'{campaign_name}/04-analysis/{comp}/{comp}_{rep}_traj'
    crd_nums=np.array(crd_nums)
    crd_mask=list()
    for i in range(len(crd_nums)):
        if os.path.exists(f'{trajDIR}/{crd_nums[i]}_{comp}.xtc'):
            crd_mask.append(False)
            print(f'{trajDIR}/{crd_nums[i]}_{comp}.xtc already exists')
        else:
            crd_mask.append(True)
    crd_mask=np.array(crd_mask)
    return crd_mask

def write_converter(comp,rep,num,crd_file,fast=False,campaign_name='FEP_campaign'):
    #write the cpptraj that will be ran to convert each trajectory file
    trajDIR=f'{campaign_name}/04-analysis/{comp}/{comp}_{rep}_traj'
    with open(f'{trajDIR}/conv_{num}_{comp}.cpptraj', 'w') as fh:
        fh.write(f'parm {trajDIR}/hNKCC1_{comp}.top\n')
        if fast:
            fh.write(f'trajin {crd_file} 1 last 10\n')
        else:
            fh.write(f'trajin {crd_file}\n')
        fh.write('autoimage\n')
        fh.write(f'trajout {trajDIR}/{num}_{comp}.xtc\n')

def convert_traj(comp,rep,fep_type,fast,campaign_name='FEP_campaign'):
    #use all previous functions to create the list of trajectories to convert, 
    #   check which ones have been converted already, then write the converter file and convert the remaining trajectories
    trajDIR=f'{campaign_name}/04-analysis/{comp}/{comp}_{rep}_traj'
    create_traj_dir(comp,rep,fep_type,campaign_name)
    crd_nums,crd_files=get_crd_numbers(comp,rep,fep_type,campaign_name)
    to_convert_mask=check_pre_converted(crd_nums,comp,rep,campaign_name)
    conv_files=list()
    for ndx,num in enumerate(crd_nums):
        if to_convert_mask[ndx]:
            write_converter(comp,rep,num,crd_files[ndx],fast,campaign_name)
            conv_files.append(f'{trajDIR}/conv_{num}_{comp}.cpptraj')
    #execute command in shell
    for conv_file in conv_files:
        os.system(f'cpptraj -i {conv_file} > {conv_file[:-7]}log')
        os.system(f'rm {conv_file}')
        os.system(f'rm {conv_file[:-7]}log')
        print(f'Converted {conv_file}')
#########################################################################################################################################

#RMSD calculation########################################################################################################################
def create_rmsd_dir(comp,rep,campaign_name='FEP_campaign'):
    #creates a directory in 04-analysis where to store the RMSD results
    RMSD_DIR=f'{campaign_name}/04-analysis/{comp}/{comp}_{rep}_rmsd'
    create_comp_dir(comp,campaign_name)
    if not os.path.exists(RMSD_DIR):
        os.mkdir(RMSD_DIR)

def execute_cpptraj_rmsd(comp,rep,campaign_name='FEP_campaign'):
    #executes the command to run cpptraj and calculate the rmsd
    #requires AmberTools
    RMSD_DIR=f'{campaign_name}/04-analysis/{comp}/{comp}_{rep}_rmsd'
    cmd=f'cd {RMSD_DIR};cpptraj -i rmsd_{comp}.cpptraj > rmsd_{comp}.log;cd ../../../../'
    os.system(cmd)
    print(f'RMSD calculated for {comp}_{rep}')

def calc_rmsd(comp,rep,fep_type,comp_mask,campaign_name='FEP_campaign'):
    #calls the functions to set up the directories, writes the cpptraj input files, and then executes it to calculate RMSD
    create_rmsd_dir(comp,rep)
    RMSD_DIR=f'{campaign_name}/04-analysis/{comp}/{comp}_{rep}_rmsd'
    trajDIR=f'{campaign_name}/04-analysis/{comp}/{comp}_{rep}_traj'
    MD_PWD=f'{campaign_name}/03-md_inputs/{comp}/{fep_type}/{rep}'
    with open(f'{RMSD_DIR}/rmsd_{comp}.cpptraj', 'w') as fh:   
        fh.write(f'parm ../{trajDIR.split("/")[-1]}/hNKCC1_{comp}.top\n')
        #fh.write(f'reference {MD_PWD}/00*rst \n')
        fh.write(f'trajin ../{trajDIR.split("/")[-1]}/*xtc \n')
        fh.write('rms rms1 :7-34,36-64,81-118,126-146,147-173,206-230,236-263,319-352,371-388,391-418@N,CA,C,O&!@H= out rms_back_h first\n')
        fh.write(f'rms rms2 :{comp}&@{comp_mask} nofit out rms_{comp} first\n')
        fh.write(f'rms rms3 :{comp}&@{comp_mask} nofit out rms_{comp}_prev previous\n')
        #fh.write(f'rms rms2 :{comp}&!@H= nofit out rms_{comp} first\n')
        #fh.write(f'rms rms3 :{comp}&!@H= nofit out rms_{comp}_prev previous\n')
        fh.write('go\n')
    execute_cpptraj_rmsd(comp,rep)

def runave(points,winperc=0.03):
    #similar to calc_running_mean()
    #must check the difference and if we can delete one function
    output=[]
    window=int(len(points)*winperc)
    for i in range(len(points[:-window])):
        avpoint=np.average(points[i:i+window])
        output.append(avpoint)
    return output

def rmsd_selector(comp,rep,rmsd_type,campaign_name='FEP_campaign'):
    #uses keywords to select the correct file for RMSD and correct plot title
    RMSD_DIR=f'{campaign_name}/04-analysis/{comp}/{comp}_{rep}_rmsd'   
    data_dict={
        'prot_bb':np.loadtxt(f"{RMSD_DIR}/rms_back_h"),
        'lig':np.loadtxt(f"{RMSD_DIR}/rms_{comp}"),
        'lig_prev':np.loadtxt(f"{RMSD_DIR}/rms_{comp}_prev")
    }
    title_dict={
        'prot_bb':f"Protein backbone RMSD - {comp} bound ",
        'lig':f"{comp} RMSD -nofit",
        'lig_prev':f"{comp} RMSD prev -nofit"
    }
    return data_dict[rmsd_type],title_dict[rmsd_type]
def plot_rmsd(comp,rep,rmsd_type,factor=10):
    #plots the rmsd for a specific replica
    Z=14
    data,title=rmsd_selector(comp,rep,rmsd_type)
    title=title+f' {rep}'

    time=data[:,0]
    time=time/factor
    rmsd=data[:,1]
    ra_rmsd=runave(rmsd)
    ra_time=time[:len(ra_rmsd)]

    fig=plt.figure(figsize=(4,3))
    ax=plt.axes()

    ax.plot(time,rmsd,c='royalblue',alpha=0.5)
    ax.plot(ra_time,ra_rmsd,c='royalblue',lw=2)

    ax.set_xlim(time[0],time[-1])

    ax.set_title(title,fontsize=Z)
    ax.set_xlabel("Time (ns)",fontsize=Z)
    ax.set_ylabel("RMSD (Å)",fontsize=Z)
    ax.tick_params(axis='both',labelsize=Z)
    plt.tight_layout()

    plt.show()
    plt.close
    
