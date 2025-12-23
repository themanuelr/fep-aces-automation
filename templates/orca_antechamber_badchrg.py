import argparse

HAND_EDIT_TOR_N=True

parser=argparse.ArgumentParser()
parser.add_argument('compound',type=str)
args=parser.parse_args()
comp=args.compound

resp2_file_path = '../03-resp/resp2.out'
mol2_file_path = f'badchrg_{comp}.mol2'
output_file_path = f'{comp}.mol2'

q_opt_values = []
with open(resp2_file_path, 'r') as file:
    lines = file.readlines()
    for line in lines:
        if 'q(opt)' in line:
            # Assuming the values start immediately after the line containing 'q(opt)'
            for data_line in lines[lines.index(line)+1:]:
                if data_line.strip() == "":  # Assuming empty line means end of table
                    break
                parts = data_line.split()
                q_opt_values.append(parts[-3])  # Assuming q(opt) is the last column

# Step 2: Read and modify badchrg_COMP.mol2
with open(mol2_file_path, 'r') as file:
    lines = file.readlines()

new_lines = []
inside_atom_section = False
q_opt_index = 0

for line in lines:
    if line.startswith('@<TRIPOS>ATOM'):
        inside_atom_section = True
        new_lines.append(line)
    elif line.startswith('@<TRIPOS>BOND'):
        inside_atom_section = False
        new_lines.append(line)
    elif inside_atom_section and line.strip():
        if HAND_EDIT_TOR_N:
            if 'ne' in line:
                line=line.replace('ne','n2')
        parts = line.split()
        if len(parts) >= 9:  # Ensure it is a valid atom line
            if float(q_opt_values[q_opt_index])<0:
                q_opt=str(q_opt_values[q_opt_index])
            elif float(q_opt_values[q_opt_index])>0:
                q_opt=f' {q_opt_values[q_opt_index]}'
            new_line = f"{line[:-10]}{q_opt}\n"
            q_opt_index += 1
            new_lines.append(new_line)
        else:
            new_lines.append(line)
    else:
        new_lines.append(line)



# Step 3: Write the new contents to a new file
with open(output_file_path, 'w') as file:
    file.writelines(new_lines)
