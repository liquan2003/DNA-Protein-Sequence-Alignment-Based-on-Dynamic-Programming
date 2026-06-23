from alignment import Alignment
from scoring_matrices import BLOSUM62, PAM250
import time

def main():
    aligner = Alignment(match_score=2, mismatch_penalty=-1, gap_penalty=-2)
    
    test_cases = [
        ("ACGT", "AGT"),
        ("ATGCT", "AGCT"),
        ("GGTTGACTA", "TGTTACGG"),
        ("AATGCT", "AGCTT")
    ]
    
    print("=== Needleman-Wunsch 全局比对 ===")
    for seq1, seq2 in test_cases:
        print(f"\n序列1: {seq1}")
        print(f"序列2: {seq2}")
        
        start = time.time()
        result = aligner.needleman_wunsch(seq1, seq2)
        elapsed = time.time() - start
        
        align1, align2 = result['alignment']
        print(f"得分: {result['score']}")
        print(f"比对结果:")
        print(f"  {align1}")
        print(f"  {align2}")
        print(f"耗时: {elapsed:.6f}s")
    
    print("\n=== Smith-Waterman 局部比对 ===")
    for seq1, seq2 in test_cases:
        print(f"\n序列1: {seq1}")
        print(f"序列2: {seq2}")
        
        start = time.time()
        result = aligner.smith_waterman(seq1, seq2)
        elapsed = time.time() - start
        
        align1, align2 = result['alignment']
        print(f"得分: {result['score']}")
        print(f"最大得分位置: {result['max_position']}")
        print(f"比对结果:")
        print(f"  {align1}")
        print(f"  {align2}")
        print(f"耗时: {elapsed:.6f}s")

if __name__ == "__main__":
    main()