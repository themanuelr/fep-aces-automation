import argparse

parser=argparse.ArgumentParser()
parser.add_argument('compound',type=str)
args=parser.parse_args()
comp=args.compound

def create_antechamber_inp(comp):
    with open(f'../../../00-compounds/{comp}.pdb', 'r') as file:
        lines_pdb = file.readlines()
    with open(f'../01-opt/opt_{comp}.xyz', 'r') as file:
        lines_opt = file.readlines()

    coordinates_opt = []
    for line in lines_opt:
        parts = line.split()
        if len(parts) == 4 and not line.startswith('Coordinates'):
            _, x, y, z = parts
            coordinates_opt.append((float(x), float(y), float(z)))

    new_lines = []
    coord_idx = 0
    for line in lines_pdb:
        if line.startswith("ATOM"):
            parts = line.split()
            if len(parts) >= 11:
                new_x = f"{coordinates_opt[coord_idx][0]:8.3f}"
                new_y = f"{coordinates_opt[coord_idx][1]:8.3f}"
                new_z = f"{coordinates_opt[coord_idx][2]:8.3f}"
                new_line = f"{line[:30]}{new_x:>8}{new_y:>8}{new_z:>8}{line[54:]}"
                new_lines.append(new_line)
                coord_idx += 1
        else:
            new_lines.append(line)

    # Write the new contents to C.pdb
    with open(f'opt_{comp}.pdb', 'w') as file:
        file.writelines(new_lines)

create_antechamber_inp(comp)