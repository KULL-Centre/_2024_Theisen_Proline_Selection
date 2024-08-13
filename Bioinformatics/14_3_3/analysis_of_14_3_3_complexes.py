# Python script for analysis of 14-3-3 complexes
# Should be executed in pymol command line
# Frederik Friis Theisen 2024, University of Copenhagen

import re
import statistics

# Mapping of amino acid codes to single-letter codes
ONE_LETTER_STANDARD = {
    'VAL': 'V', 'ILE': 'I', 'LEU': 'L', 'GLU': 'E', 'GLN': 'Q', 'ASP': 'D', 'ASN': 'N', 
    'HIS': 'H', 'TRP': 'W', 'PHE': 'F', 'TYR': 'Y', 'ARG': 'R', 'LYS': 'K', 'SER': 'S', 
    'THR': 'T', 'MET': 'M', 'ALA': 'A', 'GLY': 'G', 'PRO': 'P', 'CYS': 'C' }

ONE_LETTER = {
    'VAL': 'V', 'ILE': 'I', 'LEU': 'L', 'GLU': 'E', 'GLN': 'Q', 'ASP': 'D', 'ASN': 'N', 
    'HIS': 'H', 'TRP': 'W', 'PHE': 'F', 'TYR': 'Y', 'ARG': 'R', 'LYS': 'K', 'SER': 'S', 
    'THR': 'T', 'MET': 'M', 'ALA': 'A', 'GLY': 'G', 'PRO': 'P', 'CYS': 'C', 'SEP': 'S', 
    'TPO': 'T', 'LPD': 'P', 'DVA': 'V', 'DAL': 'A', 'B3L': 'L', 'B3A': 'A', 'UXQ': 'X', 
    'ALY': 'K', 'ACE':'', 'ACY':'','NH2':'','TSA':'X', 'BSE':'S', 'CSO':'C', 'PIP':'X',
    'B3Q': 'X'
}

DATABASE_FILE = "./14_3_3_ComplexPDBIDs.txt"  # File with PDB ids
SKIP_IDS = {'7ZMW', '7ZMU', '6TCH', '3MVH', '5OMA', '1G6G', '7MFF', '1IB1', '6Q0K', '6G6X', '8A65','3EFZ','7QI1','7Q16'} # PDB IDs have some disqualifying property
ISOFORMS = ['alpha','alpha/beta','beta/alpha','theta','zeta','delta','sigma','eta','epsilon','gamma']
GROUP_BY_MOTIF = True
PHOS_BUFFER = [2, 2]
DICTIONARYKEYS = []
SmallMoleculeMassIgnoreList = []

ALIGNMODEL = ""

def clamp(n, smallest, largest):
    return max(smallest, min(n, largest))

def setup_smallmolecule_clearlist():
    with open("smallmol_ignorelist.txt") as f:
        for line in f:
            dat = line.split()
            mw = float(dat[-1])
            name = " ".join(dat[:-1])

            SmallMoleculeMassIgnoreList.append(mw)
            SmallMoleculeMassIgnoreList.append(name)

def write_groups(peptides, filename=""):
    with open(f"{filename}_groups.tsv", "w") as f:
        header = ['MotifSequence']
        for key in DICTIONARYKEYS:
            header.append(key)
        f.write('\t'.join(header) + '\n')

        for pepid in peptides:
            peptide = peptides[pepid]
            
            grouplen = len(peptide['PDB'])
            
            for i in range(grouplen):
                line = [pepid]
                for key in DICTIONARYKEYS:
                    line.append(str(peptide[key][i])) 

                f.write(f"{'\t'.join(line)}\n")

def align_structured(groups, parentgroup = ""):
    global ALIGNMODEL

    struc_selection = f'sel_{parentgroup}_1433'
    pep_selection =   f'sel_{parentgroup}_pep'
    cmd.select(struc_selection, 'None')
    cmd.select(pep_selection, 'None')

    models = []
    for group_seq,group in groups.items():
        group_mdls = []
        for i in range(len(group['PDB'])):
            pdb_id = group['PDB'][i]
            chain = group['CHAIN'][i]
            resi = group['PHOSRESI'][i]
            contains_smallmol = group['SMALLMOL'][i] != ""
            print(pdb_id)
            cmd.fetch(pdb_id)
            cmd.do('remove solvent')
            cmd.select('PROBE', f'{pdb_id} and chain {chain} and resi {resi}')
            cmd.do('select PROBE, PROBE extend 7')

            selection_radius = 2
            radius_stepsize = 1

            # Try to increase selection radius until two chains are selected
            while True:
                cmd.select('MONOMER', f'{pdb_id} within {selection_radius} of PROBE')
                stored = {'chains':[]}
                cmd.iterate('MONOMER','chains.append(chain)',space=stored)

                print(f"{selection_radius} {len(set(stored['chains']))}")
                if len(set(stored['chains'])) == 2: # Two chains selected, exit loop
                    break
                if len(set(stored['chains'])) > 2: # More than two chains selected, go back, reduce radius increase, try again
                    selection_radius -= radius_stepsize
                    radius_stepsize /= 2

                selection_radius += radius_stepsize

                if radius_stepsize < 0.05 or selection_radius > 10: break # Give up...

            # Complete monomer selection
            cmd.do('select MONOMER, bymolecule MONOMER')
            mdl = f"{pdb_id}_{len([n for n in models if pdb_id in n]) + 1}_{parentgroup}" # Assign name to specific monomer
            if contains_smallmol: mdl += '_sm'
            cmd.extract(mdl, 'MONOMER')
            models.append(mdl)
            group_mdls.append(mdl)

            cmd.select(f'{struc_selection}',f'{struc_selection} ({mdl} and not chain {chain})')
            cmd.select(f'{pep_selection}',f'{pep_selection} ({mdl} and chain {chain})')

            if len(group_mdls) > 1: cmd.disable(mdl) # Only show one of each peptide model by default

            cmd.delete(pdb_id)

        cmd.group(f'{group_seq}_{parentgroup}',' '.join(group_mdls)) # Group peptide groups

    cmd.group(parentgroup,' '.join([f'{g}_{parentgroup}' for g in list(groups.keys())])) # Group all peptide group in cis or trans

    if ALIGNMODEL == "": ALIGNMODEL = models[1] # To align both cis and trans group to the same model
    
    for mdl in models:
        cmd.do(f'align {mdl}, {ALIGNMODEL}')

    cmd.center()

def group_by_peptide(peptides):
    output = {}
    
    for peptide in peptides:
        logged = False
        pepseq = peptide['SEQUENCE_PEPTIDE_RAW']

        if GROUP_BY_MOTIF:
            start = clamp(peptide['PHOSRESI_LOCAL'] - PHOS_BUFFER[0], 0, len(pepseq))
            end =   clamp(peptide['PHOSRESI_LOCAL'] + PHOS_BUFFER[1] + 1, 0, len(pepseq))
            pepseq = pepseq[start:end]

        for seq, details in output.items():
            if pepseq in seq:
                for key, value in peptide.items():
                    details[key].append(value)
                logged = True
                break
            elif seq in pepseq:
                new_details = output.pop(seq)
                for key, value in peptide.items():
                    new_details[key].append(value)
                output[pepseq] = new_details
                logged = True
                break
        
        if not logged:
            output[pepseq] = {key: [value] for key, value in peptide.items()}
    
    return output

def process_pdb_id(pdb_id, chains_and_phosres):
    ENTRIES = []
    print('Processing:',chains_and_phosres)
    # Loop chain ids in PDB file. dict[chainid] = resi
    for chain_and_phosres in chains_and_phosres:
        print('CHAIN, RESI:',chain_and_phosres)
        ENTRY = {}
        CHAIN = chain_and_phosres[0]
        PHOSRESI = int(chain_and_phosres[1])

        ENTRY['PDB'] = pdb_id
        ENTRY['CHAIN'] = CHAIN
        ENTRY['PHOSRESI'] = PHOSRESI
        ENTRY['ISOFORM'],ENTRY['SMALLMOL'] = get_isoform_and_substance_info_from_cif(pdb_id)
        ENTRY['MESSAGE'] = ""
        ENTRY['REMOVE'] = True

        cmd.select('PHOSRES','chain ' + CHAIN + " and resi " + str(PHOSRESI)) #select SEP or TPO
        cmd.select('PEPTIDE','bychain PHOSRES') #Complete chain selection
        cmd.select('PEPTIDE_RESOLVED', 'PEPTIDE and present') #discard non resolved residues

        # Get sequence data of theoretical and resolved peptide
        RESIDUES  = {'data':[]}
        RESIDUES_RESOLVED = {'data':[]}
        cmd.iterate('PEPTIDE and n. CA','data.append([resi,resn])',space=RESIDUES)
        cmd.iterate('PEPTIDE_RESOLVED and n. CA','data.append([resi,resn])',space=RESIDUES_RESOLVED)

        # Filter duplicated residues and identify phosres
        PHOS_RES_TYPE = ""
        LOCAL_PHOSRESIDX = None
        plus_two_is_proline = False
        local_residx = 0
        RESI_MOTIF = []
        sequencedata = []
        sequence_raw = ""
        sequence_annotated = ""
        motifsequence = ""
        occupancy = [False,False,False]
        for residue in RESIDUES['data']: #loop through [[resi,resn],...]
            residx = int(residue[0])
            resname = residue[1]

            # Check if resi is already in sequence
            if residue in sequencedata: continue # Skip if already in sequence data

            # We have iterated to PHOSRESI, now find RESI for phosres +1 and +2 (sometimes not numerically consecutive ???)
            if occupancy[0] and len(RESI_MOTIF) < 2:
                RESI_MOTIF.append(residx)
                motifsequence += ONE_LETTER[resname]
                if residue in RESIDUES_RESOLVED['data']:
                    occupancy[len(RESI_MOTIF)] = True # Position is resolved
                if len(RESI_MOTIF) == 2 and resname == 'PRO': plus_two_is_proline = True # +2 position is proline

            #Check if residue is phosres
            if residx == PHOSRESI: # This should be the phosres we are looking for

                if residue in RESIDUES_RESOLVED['data']: 
                    occupancy[0] = True
                    LOCAL_PHOSRESIDX = local_residx

                if resname == "SEP": 
                    PHOS_RES_TYPE = "S"
                    phosresi_local = len(sequencedata)
                elif resname == "TPO": 
                    PHOS_RES_TYPE = "T"
                    phosresi_local = len(sequencedata)

            sequencedata.append(residue)
            if resname in ONE_LETTER_STANDARD: 
                sequence_annotated += ONE_LETTER_STANDARD[resname] 
            else: 
                sequence_annotated += '[' + resname + ']'

            local_residx += 1

        ENTRY['PHOS_RES_TYPE'] = PHOS_RES_TYPE
        ENTRY['SEQUENCE_PEPTIDE_ANNOTATED'] = sequence_annotated
        ENTRY['SEQUENCE_PEPTIDE_RAW'] = ''.join([ONE_LETTER[aa[1]] for aa in sequencedata])
        ENTRY['SEQUENCE_MOTIF'] = motifsequence
        ENTRY['PHOSRESI_LOCAL'] = LOCAL_PHOSRESIDX

        if PHOS_RES_TYPE == "":   # No phosphorylated residue resolved in chain
            ENTRY['MESSAGE'] = "No phosphorylated residue found"
        
        elif len(ENTRY['SEQUENCE_PEPTIDE_RAW']) < ENTRY['PHOSRESI_LOCAL'] + 1 + 2: # Peptide too short for analysis
            ENTRY['MESSAGE'] = "Peptide did not extend +2 residues from phosres"
        
        elif False in occupancy: # Peptide must contain PHOSRES + 2 or more residues
            ENTRY['MESSAGE'] = "Phosres +2 not resolved"

        else: #Phosresi and occupancy necessary to calculate omega angle
            # Get omega angle of +1 to +2 peptide bond
            cmd.select('atom1','chain ' + CHAIN + " and (resi " + str(RESI_MOTIF[0]) + " and n. CA)")
            cmd.select('atom2','chain ' + CHAIN + " and (resi " + str(RESI_MOTIF[0]) + " and n. C)")
            cmd.select('atom3','chain ' + CHAIN + " and (resi " + str(RESI_MOTIF[1]) + " and n. N)")
            cmd.select('atom4','chain ' + CHAIN + " and (resi " + str(RESI_MOTIF[1]) + " and n. CA)")
            omega = cmd.get_dihedral('first atom1','first atom2','first atom3','first atom4')

            if omega < -90: omega += 360 # Fix +-180 degree switch, without messing up cis angles

            ENTRY['OMEGA'] = omega
            ENTRY['REMOVE'] = False
        
        ENTRIES.append(ENTRY)

    return ENTRIES

def get_isoform_and_substance_info_from_cif(pdb_id):
    isoform = ""
    small_mol_name = ""
    small_mol_mw = 0

    with open(pdb_id + ".cif") as f:
        lines = f.readlines()
        for line in lines:
            if "14-3-3" in line:
                for form in ISOFORMS:
                    if form in line:
                        if form not in isoform:
                            isoform += form + ','
                if isoform != "": break
        
        if isoform == "":
            isoform = "unknown"
        else: isoform = isoform[:-1]

        # Detect substances
        _detectsm = False
        for line in lines:
            if "_entity.details" in line: _detectsm = True
            elif _detectsm and line[0] == '#': break
            elif _detectsm:
                molecule = re.findall("(?:\'.*?\'|\S)+", line)
                if len(molecule) < 5: continue
                if molecule[1] == 'non-polymer' and molecule[2] == 'syn':

                    try:
                        molname = molecule[3].strip('\'')
                        mw_sm = float(molecule[4])
                        if mw_sm < 90 or mw_sm in SmallMoleculeMassIgnoreList or molname in SmallMoleculeMassIgnoreList: continue # Substance not banned
                        elif mw_sm > small_mol_mw: # Return the largest small molecule found
                            small_mol_name = molname
                            small_mol_mw = mw_sm
                    except:
                        print(molecule)

    return isoform, small_mol_name

def main():
    global DICTIONARYKEYS
    output = {}
    peptides = []
    incompatiblecomplexes = []

    setup_smallmolecule_clearlist()

    with open(DATABASE_FILE) as f:
        filedata = f.readlines()[0]
        pdb_ids = filedata.split(',')

    for pdb_id in pdb_ids:
        if pdb_id in SKIP_IDS:
            continue
        # Fetch and process data for each pdb_id
        cmd.reinitialize()
        cmd.fetch(pdb_id)
        cmd.do('remove solvent')
        # Select all phosphorylated residues of correct types
        cmd.select('PHOS','present and (resn SEP or resn TPO)')

        # Identify chains with phos residues and index of residues
        chains_and_phosres = {'data':[]}
        cmd.iterate('PHOS and n. CA','data.append([chain,resi])',space=chains_and_phosres)

        BOUND_MOTIFS = process_pdb_id(pdb_id, chains_and_phosres['data'])

        if len(BOUND_MOTIFS) > 0: 
            peptides.extend([pep for pep in BOUND_MOTIFS if not pep['REMOVE']])
            incompatiblecomplexes.extend([pep for pep in BOUND_MOTIFS if pep['REMOVE']])

        #if len(peptides) > 60: break

    DICTIONARYKEYS = list(peptides[0].keys())

    cis_peptides = [pep for pep in peptides if pep['OMEGA'] < 30]
    trans_peptides = [pep for pep in peptides if pep['OMEGA'] > 150]

    cis_groups = group_by_peptide(cis_peptides)
    trans_groups = group_by_peptide(trans_peptides)

    write_groups(cis_groups,'cis')
    write_groups(trans_groups,'trans')

    cmd.reinitialize()
    align_structured(cis_groups,'cis')
    align_structured(trans_groups,'trans')

    print()
    print("14-3-3 ANALYSIS SCRIPT COMPLETED")
    print(f"TOTAL CIS PEPTIDES COUNTED:     {len(cis_peptides)}")
    print(f"    UNIQUE PEPTIDES:            {len(cis_groups)}")
    print(f"TOTAL TRANS PEPTIDES COUNTED:   {len(trans_peptides)}")
    print(f"    UNIQUE PEPTIDES:            {len(trans_groups)}")
    print()
    if GROUP_BY_MOTIF: print(f"PEPTIDE GROUPING RANGE: PHOSRES -{PHOS_BUFFER[0]} to +{PHOS_BUFFER[1]}")

    with open('./incomp_complexes.txt','w') as f:
        for pep in incompatiblecomplexes:
            f.write(f"{pep['PDB']} {pep['CHAIN']} {pep['MESSAGE']}\n")

main()
