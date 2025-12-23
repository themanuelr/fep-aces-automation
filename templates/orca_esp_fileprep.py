import argparse

parser=argparse.ArgumentParser()
parser.add_argument('compound',type=str)
args=parser.parse_args()
comp=args.compound

def create_esp_inp(comp):
    fname=f'../01-opt/opt_{comp}.xyz'
    elem_coords=[]
    with open(fname,'r') as opt_comp:
        for line in opt_comp:
            if len(line.split())==4 and not line.startswith('Coordinates'):
                elem_coords.append(line)
    with open(f'esp_{comp}.inp','w') as inp_file:
        inp_file.write('! HF 6-31g* chelpg\n\n')
        inp_file.write('*xyz 0 1\n')
        for elem_coord in elem_coords:
            inp_file.write(elem_coord)
        inp_file.write('*')

create_esp_inp(comp)