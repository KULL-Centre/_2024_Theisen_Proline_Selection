FILE = 'cis'

def levenshtein_distance(s1, s2):
    if len(s1) < len(s2):
        return levenshtein_distance(s2, s1)

    if s2 in s1:
    	return 0

    if len(s2) == 0: # Should not be possible
        print("Should not be possible",s1,s2)
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row

    return previous_row[-1]

def count_unique_strings(amino_acid_strings):
    unique_strings = set(amino_acid_strings)  # Remove exact duplicates first

    # Prepare to track non-unique strings
    non_unique_strings = set()

    for s1 in unique_strings:
        for s2 in unique_strings:
            if s1 != s2 and s1 in s2:
                non_unique_strings.add(s1)
                break

    # Filter out non-unique strings to get the final list of unique strings
    final_unique_strings = unique_strings - non_unique_strings

    return len(final_unique_strings), final_unique_strings

def count_lenient_unique_strings(amino_acid_strings, tolerance=2):
    unique_strings = set(amino_acid_strings)  # Remove exact duplicates first

    if len(unique_strings) == 1:
    	return 1, unique_strings

    non_unique_strings = set()

    for s1 in unique_strings:
        for s2 in unique_strings:
            if s1 != s2 and levenshtein_distance(s1, s2) <= int(0.1 * len(s1)) + 0.5:
                non_unique_strings.add(s1)
                break

    # Filter out non-unique strings to get the final list of unique strings
    final_unique_strings = unique_strings - non_unique_strings

    return len(final_unique_strings), final_unique_strings

def retain_one_similar_string(amino_acid_strings, tolerance=0.2):
    retained_strings = []  # This will store the strings we decide to keep

    for s1 in amino_acid_strings:
        keep_s1 = True
        for s2 in retained_strings:
            if levenshtein_distance(s1, s2) <= tolerance * len(s1):
                keep_s1 = False
                break
        if keep_s1:
            retained_strings.append(s1)

    return len(retained_strings), retained_strings

def script():
	unique = 0
	peptides = []
	groups = []
	unique_ids = []


	with open(FILE + '_groups.txt') as f:
		group = []
		ids = []

		for line in f:
			line = line.strip()
			if line == "" and len(group) > 0:

				n, p = retain_one_similar_string(group)
				unique += n
				peptides.extend(p)
				group_ids = []

				for i in range(len(group)):
					if group[i] in p: group_ids.append(ids[i])

				groups.append([len(groups) + 1, n, p,  group_ids])

			elif "GROUP" in line:
				group = []
				ids = []
			elif line != "":
				group.append(line.split('|')[1])
				ids.append(line.split('|')[0])


	print(unique)


	with open(FILE + '_peptides_unique.txt','w') as f:
		for p in peptides:
			f.write(p + '\n')

	with open(FILE + '_peptides_unique_groups.txt','w') as f:
		for g in groups:
			out = "["
			for v in g:
				out += str(v) + ","
			f.write(out[0:-1] + ']\n')

	with open(FILE + '_unique_ids.txt', 'w') as f:
		for g in groups:
			ids = g[3]
			out = "#"+ str(g[0]) + " "
			for id in ids:
				out += id.upper() + " "
			f.write(out + '\n')

FILE = 'cis'

script()

print("CIS DONE")

FILE = 'trans'

script()

print("TRANS DONE")

print("SCRIPT DONE")

