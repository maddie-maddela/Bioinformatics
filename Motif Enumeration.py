#The probability of existing a 9-length sub-string among all 9-length sub-strings = (0.25)**9

#The number of occurrences of a 9-length sub-string in a string having length 1000 = (1000−9+1)∗(0.25)**9

#If the number of such string becomes 500, then the number of occurrences would be = 500∗(1000−9+1)∗(0.25)**9

k_mer_length = 9
probability_single_string = 1 / (4 ** k_mer_length)
number_of_strings = 500
length_of_each_string = 1000
expected_occurrences = probability_single_string * number_of_strings * (length_of_each_string-k_mer_length+1)
print(expected_occurrences) #1.8920898


import itertools

def combination(k):
    return (''.join(p) for p in itertools.product('ATCG', repeat=k))

def hamming_distance(pattern, seq):
    return sum(c1 != c2 for c1, c2 in zip(pattern, seq))

def window(s, k):
    for i in range(1 + len(s) - k):
        yield s[i:i+k]

def motif_enumeration(k, d, DNA):
    pattern = set()
    for combo in combination(k):
        if all(any(hamming_distance(combo, pat) <= d 
                for pat in window(string, k)) for string in DNA):
            pattern.add(combo)
    return pattern

     
DNA = ['CAATAATAGTCGTCGCAATTGGCTC' 'CAATCCTTCGTGCAGGGGGTGGTTC' 'TTCTATAGTCACCGGCTCGCCCTCG' 'TACCCTAACCCTTCGGACCAAAGGA' 'CTGTCTTGGTCTTCGGAGGTCAGTA' 'CCTTCCATCGGTATCTATACTCATG']
print(*motif_enumeration(5, 1, DNA))



# from itertools import product

# def hamming_distance(str1, str2):
#     return sum(c1 != c2 for c1, c2 in zip(str1, str2))

# def generate_kmers(k):
#     return [''.join(p) for p in product('ACGT', repeat=k)]

# def motif_enumeration(dna, k, d):
#     patterns = set()

#     for pattern in generate_kmers(k):
#         found_in_all = True
#         for sequence in dna:
#             found_in_sequence = False
#             for i in range(len(sequence) - k + 1):
#                 kmer = sequence[i:i+k]
#                 if hamming_distance(pattern, kmer) <= d:
#                     found_in_sequence = True
#                     break
#             if not found_in_sequence:
#                 found_in_all = False
#                 break
#         if found_in_all:
#             patterns.add(pattern)

#     return patterns

# # Example usage:          
# dna_sequences = ['TGACCGAATGTATGGCGGGGTAAAC' 'CTGCTCGGCGATCTTCTTTAGCACG' 'CGAATTTAGGCAATAGCCTACGGTG' 'AGCCTCGGTGTCAAAGGACAATGCT' 'CGGGTCGGAGCGTGCGGTCGTGAAC' 'ATGGTGACTATCGGTCGGTGGTCAG']#["ATCGTACG", "TACGTGAC", "CGTACGTG", "GTACGTAC"]
# k_size = 5
# mismatches_allowed = 1

# result = motif_enumeration(dna_sequences, k_size, mismatches_allowed)
# print("Motif Enumeration Result:", *result)
