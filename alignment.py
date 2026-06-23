import numpy as np

class Alignment:
    def __init__(self, match_score=2, mismatch_penalty=-1, gap_penalty=-2):
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty
        self.gap_penalty = gap_penalty
    
    def _score(self, a, b):
        if a == b:
            return self.match_score
        return self.mismatch_penalty
    
    def needleman_wunsch(self, seq1, seq2):
        m, n = len(seq1), len(seq2)
        dp = np.zeros((m + 1, n + 1), dtype=int)
        
        for i in range(1, m + 1):
            dp[i, 0] = i * self.gap_penalty
        for j in range(1, n + 1):
            dp[0, j] = j * self.gap_penalty
        
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                match = dp[i-1, j-1] + self._score(seq1[i-1], seq2[j-1])
                delete = dp[i-1, j] + self.gap_penalty
                insert = dp[i, j-1] + self.gap_penalty
                dp[i, j] = max(match, delete, insert)
        
        align1, align2 = '', ''
        i, j = m, n
        while i > 0 or j > 0:
            if i > 0 and j > 0 and dp[i, j] == dp[i-1, j-1] + self._score(seq1[i-1], seq2[j-1]):
                align1 = seq1[i-1] + align1
                align2 = seq2[j-1] + align2
                i -= 1
                j -= 1
            elif i > 0 and dp[i, j] == dp[i-1, j] + self.gap_penalty:
                align1 = seq1[i-1] + align1
                align2 = '-' + align2
                i -= 1
            else:
                align1 = '-' + align1
                align2 = seq2[j-1] + align2
                j -= 1
        
        return {
            'score': dp[m, n],
            'alignment': (align1, align2),
            'matrix': dp
        }
    
    def smith_waterman(self, seq1, seq2):
        m, n = len(seq1), len(seq2)
        dp = np.zeros((m + 1, n + 1), dtype=int)
        
        max_score = 0
        max_pos = (0, 0)
        
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                match = dp[i-1, j-1] + self._score(seq1[i-1], seq2[j-1])
                delete = dp[i-1, j] + self.gap_penalty
                insert = dp[i, j-1] + self.gap_penalty
                dp[i, j] = max(0, match, delete, insert)
                
                if dp[i, j] > max_score:
                    max_score = dp[i, j]
                    max_pos = (i, j)
        
        align1, align2 = '', ''
        i, j = max_pos
        
        while dp[i, j] > 0:
            if i > 0 and j > 0 and dp[i, j] == dp[i-1, j-1] + self._score(seq1[i-1], seq2[j-1]):
                align1 = seq1[i-1] + align1
                align2 = seq2[j-1] + align2
                i -= 1
                j -= 1
            elif i > 0 and dp[i, j] == dp[i-1, j] + self.gap_penalty:
                align1 = seq1[i-1] + align1
                align2 = '-' + align2
                i -= 1
            else:
                align1 = '-' + align1
                align2 = seq2[j-1] + align2
                j -= 1
        
        return {
            'score': max_score,
            'alignment': (align1, align2),
            'matrix': dp,
            'max_position': max_pos
        }