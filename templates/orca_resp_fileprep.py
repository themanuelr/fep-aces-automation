import argparse
import numpy as np
from PyAstronomy import pyasl

parser=argparse.ArgumentParser()
parser.add_argument('compound',type=str)
args=parser.parse_args()
comp=args.compound

def create_resp_inp(comp):
    elems=np.genfromtxt(f'../01-opt/opt_{comp}.xyz',skip_header=2,dtype=str)[:,0]
    an=pyasl.AtomicNo()
    atomNums=[an.getAtomicNo(elem) for elem in elems]
    fname=f'../02-esp/esp_{comp}.scfp.esp_{comp}.vpot'
    vpot_lines=[]
    with open(fname,'r') as vpot_file:
        for ndx,line in enumerate(vpot_file):
            if ndx==0:
                Natoms,gridpoints=line.split()
                new_line=f'   {Natoms}{gridpoints}\n'
                vpot_lines.append(new_line)
            elif ndx<=int(Natoms):
                old_line=line.strip('\n')
                if len(elems[ndx-1])>1:new_line=f'{old_line}{atomNums[ndx-1]:>10}{elems[ndx-1]:>3}\n'
                elif len(elems[ndx-1])==1:new_line=f'{old_line}{atomNums[ndx-1]:>10}{elems[ndx-1]:>2}\n'
                vpot_lines.append(new_line)
            else:
                vpot_lines.append(line)
    with open(f'esp_{comp}.vpot','w') as inp_file:
        for line in vpot_lines:
            inp_file.write(line)

create_resp_inp(comp)
