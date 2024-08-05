import shlex
import json
import re
import pymol2

aa_codes ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q','ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y','ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A','GLY':'G', 'PRO':'P', 'CYS':'C', 'SEP':'[SEP]','TPO':'[TPO]', 'LPD':'P', 'DVA':'V', 'DAL':'A', 'B3L':'L', 'B3A':'A', 'UXQ':'X', 'ALY':'K'}
one_letter ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q','ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y','ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A','GLY':'G', 'PRO':'P', 'CYS':'C', 'SEP':'S','TPO':'T', 'LPD':'P', 'DVA':'V', 'DAL':'A', 'B3L':'L', 'B3A':'A', 'UXQ':'X', 'ALY':'K'}


def read_database(path):
    ids = []
    with open(path) as f:
        for line in f:
            ids.extend(line.split(','))
    return [id.lower() for id in ids]

def get_database(database, pymol):
    for id in database:
        pymol.cmd.reinitialize()
        pymol.cmd.fetch(id)
    pymol.cmd.reinitialize()

def clean_resi(resi):
    idx = re.sub(r"\D", "", resi)
    if idx != '': return int(idx)
    return None

def clean_resi_list(start_resi):
    return [clean_resi(resi) for resi in start_resi if re.sub(r"\D", "", resi).isdigit()]

def parse_cif_for_chain_info(pdb_code):
    """
    Parses CIF file to extract chain information, focusing on chains under 20 residues in length,
    or the shortest chain if none are under 20 residues. Returns PDB ID, chain ID, protein name, and length.

    Args:
    file_path (str): Path to the CIF file.

    Returns:
    tuple: (PDB ID, chain ID, protein name, chain length)
    """
    file_path = './' + pdb_code + '.cif'

    # Storing chain information and protein names
    chain_info = {}
    protein_names = {}
    accession_code = {}

    with open(file_path, 'r') as file:
        #lines = file.readlines()
        current_loop = None
        multiline = False
        headers = []
        loop_id = None
        length = {}

        for line in file:
            # Identify the beginning of a loop
            if line.startswith('loop_'):
                current_loop = 'loop'
                headers = []
                continue

            # Identify the end of a loop
            if line.startswith('#'):
                current_loop = None
                continue

            # Process headers in the loop
            if current_loop and line.startswith('_'):
                headers.append(line.strip())
                loop_id = line.split('.')[0]
                continue

            # Process values in the loop
            if current_loop and line.strip() and not line.startswith('_') and not line.startswith(';'):
                if loop_id in ['_struct_ref_seq','_struct_ref','_entity']: # Speed optimization
                    values = []
                    values.extend(shlex.split(line))

                    if len(values) < len(headers):
                        while len(values) < len(headers):
                            print(len(headers),len(values),values)
                            value = ""
                            for continued_line in file:
                                print(loop_id,continued_line.strip())
                                if len(value) == 0 and ';' not in continued_line:
                                    values.extend(shlex.split(continued_line.strip()))
                                    print("END MULTILINE 2",continued_line.strip())
                                    break
                                if len(value) > 0 and ';' in continued_line or '#' in continued_line:
                                    print("END MULTILINE 1",value)
                                    values.append(value)
                                    break
                                else:
                                    value += continued_line.strip(';').strip()

                    if '_struct_ref_seq.pdbx_strand_id' in headers:
                        chain_id_index = headers.index('_struct_ref_seq.pdbx_strand_id')
                        beg_index = headers.index('_struct_ref_seq.pdbx_auth_seq_align_beg')
                        end_index = headers.index('_struct_ref_seq.pdbx_auth_seq_align_end')
                        ref_id_index = headers.index('_struct_ref_seq.ref_id')
                        ref_id = values[ref_id_index]
                        chain_id = values[chain_id_index]

                        if chain_id not in length: length[chain_id] = 0
                        beg = clean_resi(values[beg_index])
                        end = clean_resi(values[end_index])
                        if beg != None and end != None: length[chain_id] += (end - beg + 1)

                        chain_info[chain_id] = {'ref_id': ref_id, 'length': length[chain_id]}

                    elif '_struct_ref.pdbx_db_accession' in headers:
                        ref_id_index = headers.index('_struct_ref.id')
                        acc_index = headers.index('_struct_ref.pdbx_db_accession')
                        if len(values) <= acc_index: continue
                        ref_id = values[ref_id_index]
                        acc = values[acc_index]
                        accession_code[ref_id] = acc

                    elif '_entity.pdbx_description' in headers:
                        ref_id_index = headers.index('_entity.id')
                        name_index = headers.index('_entity.pdbx_description')
                        if len(values) <= name_index: continue
                        ref_id = values[ref_id_index]
                        name = values[name_index]
                        protein_names[ref_id] = name


        # Cross-referencing to get full chain details
        chains_under_20 = []
        shortest_chain = None
        shortest_length = float('inf')

        for chain_id, info in chain_info.items():
            ref_id = info['ref_id']
            length = info['length']
            name = protein_names.get(ref_id, 'Unknown')
            acc = accession_code.get(ref_id, 'Unknown')

            entry = {'pdb':pdb_code,'chain':chain_id, 'name':name, 'length':length, 'accession_code':acc}

            if length <= 30:
                chains_under_20.append(entry)

            if length < shortest_length:
                shortest_chain = entry
                shortest_length = length

        # Decide what to return based on the findings
        if chains_under_20:
            return chains_under_20
        else:
            return [shortest_chain]

def parse_cif_for_uniprot_info(pdb_code):
    file_path = './' + pdb_code + '.cif'

    # Storing chain information and protein names
    chain_info = {}
    protein_names = {}
    accession_info = {}

    with open(file_path, 'r') as file:
        #lines = file.readlines()
        current_loop = None
        multiline = False
        headers = []
        uniprot_id = None

        for line in file:
            # Identify the beginning of a loop
            if line.startswith('loop_'):
                current_loop = 'loop'
                headers = []
                continue

            # Identify the end of a loop
            if line.startswith('#'):
                current_loop = None
                continue

            # Process headers in the loop
            if current_loop and line.startswith('_'):
                headers.append(line.strip().split('.')[1])
                loop_id = line.split('.')[0]
                continue

            # Process values in the loop
            if current_loop and line.strip() and not line.startswith('_') and not line.startswith(';'):
                values = []
                values.extend(shlex.split(line))

                if len(values) < len(headers):
                    while len(values) < len(headers):
                        print(len(headers),len(values),values)
                        value = ""
                        for continued_line in file:
                            print(loop_id,continued_line.strip())
                            if len(value) == 0 and ';' not in continued_line:
                                values.extend(shlex.split(continued_line.strip()))
                                print("END MULTILINE 2",continued_line.strip())
                                break
                            if len(value) > 0 and ';' in continued_line or '#' in continued_line:
                                print("END MULTILINE 1",value)
                                values.append(value)
                                break
                            else:
                                value += continued_line.strip(';').strip()

                if loop_id == '_struct_ref':
                    db_name_index = headers.index('db_name')
                    db_code_index = headers.index('db_code')
                    ref_id_index = headers.index('id ')
                    db_name = values[db_name_index]
                    db_code = values[db_code_index]
                    entity_id = values[ref_id_index]
                    if accession_info[ref_id] is None: accession_info[ref_id] = []
                    
                    accession_info[ref_id].append({'db':db_name,'code':db_code})

                if loop_id == '_struct_ref_seq':
                    chain_id_index = headers.index('pdbx_strand_id')
                    ref_id_index = headers.index('ref_id')
                    ref_id = values[ref_id_index]
                    chain_id = values[chain_id_index]

                    chain_info[chain_id] = {'ref_id': ref_id, 'chain': chain_id}


        # Collect data for chain and databases
        for chain in chain_info:
            acc = accession_info[chain['ref_id']]

def analyze_cis_trans_isomerization(pdb_id, pymol):
    prolines = []

    pymol.cmd.reinitialize()
    pymol.cmd.fetch(pdb_id)
    chains = pymol.cmd.get_chains()
    chain_ids = []
    for chain in chains:
        # Count the number of residues in each chain
        residue_count = pymol.cmd.count_atoms(f"chain {chain} and name CA")
        if residue_count < 30: chain_ids.append(chain)


    for chain_id in chain_ids:
        # Select the chain
        chain_selection = f"{pdb_id} and chain {chain_id}"

        # Find the starting residue index for the chain
        start_resi = []
        cmd.iterate(f"{chain_selection} and present", "start_resi.append(resi)", space={'start_resi': start_resi})
        if not start_resi:
            continue  # Skip if no residues found

        # Filter out non-numeric elements
        start_resi = clean_resi_list(start_resi)
        if not start_resi:
            continue  # Skip if no numeric residues found

        start_resi = int(min(start_resi, key=int)) #get minimum resi in list

        # Find proline residues in the chain
        proline_indices = []
        pymol.cmd.iterate(f"{chain_selection} and present and resn PRO and name CA", 
                    "proline_indices.append(resi)",
                    space={'proline_indices': proline_indices})
        proline_indices = clean_resi_list(proline_indices)
        proline_indices = [resi for resi in proline_indices if resi != start_resi]

        # Find preceeding residues
        tmp_sequence = {}
        pymol.cmd.iterate(f"{chain_selection} and present and name CA", "sequence[resi] = resn", space={'sequence': tmp_sequence})
        sequence = {}
        for resi in tmp_sequence: 
            sequence[clean_resi(resi)] = tmp_sequence[resi] # Convert to int

        # Analyze each XP bond
        if len(proline_indices) == 0: continue

        # Analyze each XP bond
        for index in proline_indices:
            context = ""
            for c_idx in range(index-5,index+5):
                aa = X
                if c_idx >= 0 and c_idx in sequence:
                    aa = one_letter[sequence[c_idx]]
                context += aa

            prev_residue_index = index - 1
            if prev_residue_index >= 0 and prev_residue_index in sequence:
                prev_residue_type = sequence[prev_residue_index]

                # Define the selection for the peptide bond
                peptide_bond_selection = f"{chain_selection} and resi {prev_residue_index}-{index}"

                peptide_bond_selection_1 = f"{chain_selection} and resi {prev_residue_index}"
                peptide_bond_selection_2 = f"{chain_selection} and resi {index}"

                try:
                    # Measure the omega angle
                    omega_angle = pymol.cmd.get_dihedral(f"{peptide_bond_selection_1} and name CA", 
                                               f"{peptide_bond_selection_1} and name C", 
                                               f"{peptide_bond_selection_2} and name N", 
                                               f"{peptide_bond_selection_2} and name CA")

                    if omega_angle < -90: omega_angle += 360

                    entry = {'chain':chain_id, 'resi':index, 'prev':aa_codes[prev_residue_type], 'context':context, 'omega':omega_angle}
                    prolines.append(entry)
                except:
                    print(pdb_id,chain_id,index,'error')
                
    return prolines

def group_data_by_resi(data):
    grouped_data = {}
    for item in data:
        key = (item['resi'], item['prev'])
        if key not in grouped_data:
            grouped_data[key] = {'resi': item['resi'], 'prev': item['prev'], 'context': item['context'], 'omegas': [], 'chains':[]}
        grouped_data[key]['omegas'].append(item['omega'])
        grouped_data[key]['chains'].append(item['chain'])

    # Convert the grouped data into a list of dictionaries
    return [value for value in grouped_data.values()]


def process():

    data = {}
    database = {}

    with pymol2.PyMOL() as pymol:
        ids = read_database('../pdb_lists/231129_sub30res_peptides.txt')[0:20]

        #get_database(ids, pymol)

        for id in ids:
            print(id)
            #chain_information = parse_cif_for_chain_info(id)

            #database[id] = {'info':chain_information,'data':[]}

            dat = analyze_cis_trans_isomerization(id, pymol)
            database[id] = group_data_by_resi(dat)

        with open('../database.json','w') as json_file:
            json.dump(database, json_file)

        totalcomplexes = 0
        totalchains = 0
        prolines = 0
        cisprolines = 0
        transprolines = 0
        prolinecontainingpeptides = 0
        cisprolinecontainingpeptides = 0

    with open('../pdb_angles.txt','w') as f:
        for id in database:
            entry = database[id]
            if len(entry) > 0:
                prolinecontainingpeptides += 1

                for proline in entry:
                    out = id + '_' + str(proline['resi']) + ','
                    prolines += 1

                    cis = False
                    trans = False

                    for omega in proline['omegas']:
                        out += str(omega) + ','

                        if abs(omega) < 45: cis = True
                        elif abs(omega) > 135 and abs(omega) < 225: trans = True
                    
                    if cis: cisprolines += 1
                    if trans: transprolines += 1
                
                    f.write(out.strip(',') + '\n')

    print('IDs:',len(ids))
    print('Total Complexes Analyzed:   ',len(database))
    print('Total Prolines Found:       ',prolines)
    print('Proline Containing Peptides:',prolinecontainingpeptides)
    print('Total Trans Prolines:       ',transprolines)
    print('Total Cis Prolines:         ',cisprolines)

process()