from collections import defaultdict

# Define the encoding for each base
base_encoding = {
    "A": "0*01",
    "C": "11*0",
    "T": "10*0",
    "G": "0*11",
    "N": "****",  # N for unknown base
}

def encode_sequence(sequence):
    """
    :param sequence: DNA Sequence (string, e.g., "ACTG")
    :return: Sequence after encoding (list of strings)
    """
    encoded_sequence = []
    for base in sequence:
        if base in base_encoding:
            encoded_sequence.append(base_encoding[base])
        else:
            encoded_sequence.append("????") # Placeholder for unknown base
    return encoded_sequence


def generate_encoded_seeds(read, seed_length):
    """
    :param read: Gene read data (string)
    :param read:(string)
    :param seed_length:  (int)
    :return: (list of list of strings)
    """
    seeds = []
    for i in range(len(read) - seed_length + 1):
        seed = read[i:i + seed_length]
        encoded_seed = encode_sequence(seed)
        seeds.append(encoded_seed)
    return seeds


def tcam_lookup_with_encoding(encoded_seed, encoded_reference, max_hamming_distance):
    """
    :param encoded_seed: (list of strings)
    :param encoded_reference: (list of strings)
    :param max_hamming_distance: (int)
    :return: 匹配位置列表 (list of int)
    """
    matches = []
    seed_length = len(encoded_seed)
    for i in range(len(encoded_reference) - seed_length + 1):
        reference_window = encoded_reference[i:i + seed_length]
        mismatches = 0
        for seed_base, ref_base in zip(encoded_seed, reference_window):
            if not match_with_wildcard(seed_base, ref_base):
                mismatches += 1
                if mismatches > max_hamming_distance:
                    break
        if mismatches <= max_hamming_distance:
            matches.append(i)
    return matches

def match_with_wildcard(encoded_base, encoded_reference):
    """
    Check if the encoded base matches the encoded reference.
    :param encoded_base: (string)
    :param encoded_reference: (string)
    :return:(bool)
    """
    for b, r in zip(encoded_base, encoded_reference):
        if b != "*" and r != "*" and b != r:
            return False
    return True



def voting(matches, locality_size):
    """
    Count votes for each locality.
    :param matches: (list of int)
    :param locality_size: (int)
    :return: (dict of int -> int)
    """
    vote_counts = defaultdict(int)
    for match in matches:
        locality_start = (match // locality_size) * locality_size
        vote_counts[locality_start] += 1
    return vote_counts


def filtering(vote_counts, vote_threshold):
    """
    filter the candidates with the given vote threshold.
    :param vote_counts: (dict of int -> int)
    :param vote_threshold:(int)
    :return: (list of int)
    """
    candidates = [position for position, count in vote_counts.items() if count >= vote_threshold]
    return candidates


def smith_waterman(query, reference):
    """
    Perform Smith-Waterman algorithm to find the best alignment score.
    :param query:(string)
    :param reference:(string)
    :return:(int)
    """
    m, n = len(query), len(reference)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    max_score = 0

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = dp[i - 1][j - 1] + (2 if query[i - 1] == reference[j - 1] else -1)
            delete = dp[i - 1][j] - 1
            insert = dp[i][j - 1] - 1
            dp[i][j] = max(0, match, delete, insert)
            max_score = max(max_score, dp[i][j])

    return max_score


def fsva_with_encoding(read, reference, seed_length, max_hamming_distance, locality_size, vote_threshold):
    """
    :param read:(string)
    :param reference:(string)
    :param seed_length:(int)
    :param max_hamming_distance:(int)
    :param locality_size:int)
    :param vote_threshold:(int)
    :return:(list of int)
    """
    # Encode the read and reference
    encoded_read = encode_sequence(read)
    encoded_reference = encode_sequence(reference)
    
    # Parse the read into fixed-length seeds and encode them
    encoded_seeds = generate_encoded_seeds(read, seed_length)
    
    # TCAM lookup
    all_matches = []
    for encoded_seed in encoded_seeds:
        matches = tcam_lookup_with_encoding(encoded_seed, encoded_reference, max_hamming_distance)
        all_matches.extend(matches)
    
    # Voting and filtering
    vote_counts = voting(all_matches, locality_size)
    candidates = filtering(vote_counts, vote_threshold)
    
    return candidates



read = "ACGTACGT"
reference = "ACGTACGTACGTACGT"
seed_length = 4
max_hamming_distance = 1
locality_size = 4
vote_threshold = 2

candidates = fsva_with_encoding(read, reference, seed_length, max_hamming_distance, locality_size, vote_threshold)
print("候選位置:", candidates)

