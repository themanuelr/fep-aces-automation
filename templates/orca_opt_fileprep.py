import argparse

parser=argparse.ArgumentParser()
parser.add_argument('compound',type=str)
args=parser.parse_args()
comp=args.compound

def create_opt_inp(comp,mol_fmt='mol2'):
    if mol_fmt!='mol2':return 'error- fmt must be mol2'
    fname=f'../../../00-compounds/pre_{comp}.{mol_fmt}'
    elem_coords=[]
    with open(fname,'r') as ini_comp:
        check=False
        for line in ini_comp:
            if '@<TRIPOS>BOND' in line:
                break
            if check:
                splitted=line.split()
#                elem=splitted[5][0]
                elem=splitted[5]
                if '.' in elem: elem=elem.split('.')[0]
                coords=splitted[2:5]
                formatted_coords = [f"{coord.split('.')[0]}.{coord.split('.')[-1]:0<6}" for coord in coords]
                elem_coords.append([elem,formatted_coords])
            if '@<TRIPOS>ATOM' in line:
                check=True
    with open(f'opt_{comp}.inp','w') as inp_file:
        inp_file.write('! B3LYP 6-31g* opt\n\n')
        inp_file.write('*xyz 0 1\n')
        for elem_coord in elem_coords:
            atom,(x,y,z)=elem_coord
            formatted_string=f'{atom:>2}{x:>17}{y:>16}{z:>16}\n'
            inp_file.write(formatted_string)
        inp_file.write('*')

create_opt_inp(comp)
