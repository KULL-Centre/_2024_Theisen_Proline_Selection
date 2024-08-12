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
SKIP_IDS = {'7ZMW', '7ZMU', '6TCH', '3MVH', '5OMA', '1G6G', '7MFF', '1IB1', '6Q0K', '6G6X', '8A65','3EFZ'}
ISOFORMS = ['alpha','alpha/beta','beta/alpha','theta','zeta','delta','sigma','eta','epsilon','gamma']
SMALL_MOLECULE_MASS_IGNORE_LIST = set()
GROUP_BY_MOTIF = False
PHOS_BUFFER = [0, 2]
DICTIONARYKEYS = []

def clamp(n, smallest, largest):
    return max(smallest, min(n, largest))

def write_groups(peptides, filename=""):
    with open(f"{filename}_groups.txt", "w") as f:
        header = ['MotifSequence']
        for key in DICTIONARYKEYS:
            header.append(key)
        f.write(','.join(header) + '\n')

        for pepid in peptides:
            peptide = peptides[pepid]
            
            grouplen = len(peptide['PDB'])
            
            for i in range(grouplen):
                line = [pepid]
                for key in DICTIONARYKEYS:
                    line.append(str(peptide[key][i])) 

                f.write(f"{','.join(line)}\n")

# def align_structured(peptides):
#     models = []
#     for id in (v[0] for v in peptides):
#         cmd.fetch(id)
#         cmd.select('PHOS', f'{id} and (resn SEP or resn TPO) and visible')
#         residues = {'data': []}
#         cmd.iterate('PHOS and n. CA', 'data.append([resi, chain])', space=residues)
#         unique_chains = set()
#         for res in residues['data']:
#             if res[1] not in unique_chains:
#                 unique_chains.add(res[1])
#                 cmd.select('probe', f'{id} and chain {res[1]} and resi {res[0]}')
#                 cmd.do('select probe, probe extend 7')
#                 cmd.select('monomer', f'{id} within 7 of probe')
#                 cmd.do('select monomer, bymolecule monomer')
#                 mdl = f"{id}_{len(models) + 1}"
#                 cmd.extract(mdl, 'monomer')
#                 models.append(mdl)
#         cmd.delete(id)
#     cmd.do('remove solvent')
#     for mdl in models:
#         cmd.do(f'align {mdl}, {models[4]}')

def group_by_peptide(peptides):
    output = {}
    
    for peptide in peptides:
        logged = False
        pepseq = peptide['SEQUENCE_PEPTIDE_RAW']

        print(pepseq)

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
        ENTRY['ISOFORM'] = get_isoform_from_cif(pdb_id)

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
        ENTRY['SEQUENCE_PEPTIDE_RAW'] = ''.join([ONE_LETTER[aa[1]] for aa in sequencedata]) #data sequence_annotated.replace('[TPO]','T').replace('[SEP]','S')''
        ENTRY['SEQUENCE_MOTIF'] = motifsequence
        ENTRY['PHOSRESI_LOCAL'] = LOCAL_PHOSRESIDX

        if PHOS_RES_TYPE == "": continue # No phosphorylated residue resolved in chain, continue to next chain
        if False in occupancy: continue # Peptide must contain PHOSRES + 2 or more residues

        # Get omega angle of +1 to +2 peptide bond
        cmd.select('atom1','chain ' + CHAIN + " and (resi " + str(RESI_MOTIF[0]) + " and n. CA)")
        cmd.select('atom2','chain ' + CHAIN + " and (resi " + str(RESI_MOTIF[0]) + " and n. C)")
        cmd.select('atom3','chain ' + CHAIN + " and (resi " + str(RESI_MOTIF[1]) + " and n. N)")
        cmd.select('atom4','chain ' + CHAIN + " and (resi " + str(RESI_MOTIF[1]) + " and n. CA)")
        omega = cmd.get_dihedral('first atom1','first atom2','first atom3','first atom4')

        if omega < -90: omega += 360 # Fix +-180 degree switch, without messing up cis angles

        ENTRY['OMEGA'] = omega
        ENTRIES.append(ENTRY)

    return ENTRIES

def get_isoform_from_cif(pdb_id):
    isoform = ""

    with open(ID + ".cif") as f:
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

    return isoform

def main():
    global DICTIONARYKEYS
    output = {}
    peptides = []
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

        if len(BOUND_MOTIFS) > 0: peptides.extend(BOUND_MOTIFS)

        #if len(peptides) > 40: break

    DICTIONARYKEYS = list(peptides[0].keys())

    cis_peptides = [pep for pep in peptides if pep['OMEGA'] < 30]
    trans_peptides = [pep for pep in peptides if pep['OMEGA'] > 150]

    cis_groups = group_by_peptide(cis_peptides)
    trans_groups = group_by_peptide(trans_peptides)

    write_groups(cis_groups,'cis')
    write_groups(trans_groups,'trans')

    print()
    print("14-3-3 ANALYSIS SCRIPT COMPLETED")
    print(f"TOTAL CIS PEPTIDES COUNTED:     {len(cis_peptides)}")
    print(f"    UNIQUE PEPTIDES:            {len(cis_groups)}")
    print(f"TOTAL TRANS PEPTIDES COUNTED:   {len(trans_peptides)}")
    print(f"    UNIQUE PEPTIDES:            {len(trans_groups)}")
    print()
    if GROUP_BY_MOTIF: print(f"PEPTIDE GROUPING RANGE: PHOSRES -{PHOS_BUFFER[0]} to +{PHOS_BUFFER[1]}")

main()


