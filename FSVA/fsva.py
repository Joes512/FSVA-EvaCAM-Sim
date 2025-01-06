from collections import defaultdict
from tqdm import tqdm

def generate_seeds(read, seed_length):
    """
    將讀取序列分解為固定長度的種子。
    :param read: 待處理的基因讀取數據 (string)
    :param seed_length: 種子長度 (int)
    :return: 種子列表 (list of strings)
    """
    seeds = []
    for i in range(len(read) - seed_length + 1):
        seeds.append(read[i:i + seed_length])
    return seeds


def tcam_lookup(seed, reference, max_hamming_distance):
    """
    模擬 TCAM 查找過程，支持近似匹配。
    :param seed: 待查找的種子 (string)
    :param reference: 參考基因組 (string)
    :param max_hamming_distance: 最大 Hamming 距離 (int)
    :return: 匹配位置列表 (list of int)
    """
    matches = []
    for i in range(len(reference) - len(seed) + 1):
        window = reference[i:i + len(seed)]
        # 計算 Hamming 距離
        mismatches = sum(1 for a, b in zip(seed, window) if a != b)
        if mismatches <= max_hamming_distance:
            matches.append(i)
    return matches


def voting(matches, locality_size):
    """
    根據種子匹配結果統計投票數。
    :param matches: 匹配位置列表 (list of int)
    :param locality_size: 投票區域大小 (int)
    :return: 投票計數字典 (dict of int -> int)
    """
    vote_counts = defaultdict(int)
    for match in matches:
        locality_start = (match // locality_size) * locality_size
        vote_counts[locality_start] += 1
    return vote_counts


def filtering(vote_counts, vote_threshold):
    """
    根據投票數過濾候選位置。
    :param vote_counts: 投票計數字典 (dict of int -> int)
    :param vote_threshold: 投票數門檻 (int)
    :return: 候選位置列表 (list of int)
    """
    candidates = [position for position, count in vote_counts.items() if count >= vote_threshold]
    return candidates


def smith_waterman(query, reference):
    """
    使用 Smith-Waterman 演算法執行局部比對。
    :param query: 待比對的序列 (string)
    :param reference: 參考基因組 (string)
    :return: 最佳比對分數 (int)
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



def fsva(read, reference, seed_length, max_hamming_distance, locality_size, vote_threshold):
    """
    完整的 Fast Seed-and-Vote Algorithm (FSVA)。
    :param read: 基因讀取數據 (string)
    :param reference: 參考基因組 (string)
    :param seed_length: 種子長度 (int)
    :param max_hamming_distance: 最大 Hamming 距離 (int)
    :param locality_size: 投票區域大小 (int)
    :param vote_threshold: 投票數門檻 (int)
    :return: 候選位置和其對應的匹配分數 (list of tuple)
    """
    seeds = generate_seeds(read, seed_length)
    all_matches = []
    for seed in seeds:
        matches = tcam_lookup(seed, reference, max_hamming_distance)
        all_matches.extend(matches)
    
    vote_counts = voting(all_matches, locality_size)
    candidates = filtering(vote_counts, vote_threshold)
    
    results = []
    for candidate in tqdm(candidates):
        query_segment = read
        reference_segment = reference[candidate:candidate + len(read)]
        score = smith_waterman(query_segment, reference_segment)
        results.append((candidate, score))
    
    return results

# 主程式執行
with open('reference.txt', 'r', encoding='utf-8') as ref_file:
    reference = ref_file.read().replace('\n', '')
    # print(reference)
    

with open('971_200.txt', 'r', encoding='utf-8') as read_file:
    read = read_file.read().replace('\n', '')

seed_length = 4
max_hamming_distance = 1
locality_size = 4
vote_threshold = 2

# 執行 FSVA
results = fsva(read, reference, seed_length, max_hamming_distance, locality_size, vote_threshold)

# 根據分數排序，取分數最高的三個結果
top_results = sorted(results, key=lambda x: x[1], reverse=True)[:3]

# 輸出結果
print("分數最高的三個匹配位置與分數:")
for position, score in top_results:
    print(f"位置: {position}, 分數: {score}")