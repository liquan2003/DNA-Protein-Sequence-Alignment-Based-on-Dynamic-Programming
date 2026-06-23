# DNA/蛋白质序列比对算法

基于动态规划的序列比对算法实现，包含Needleman-Wunsch全局比对和Smith-Waterman局部比对。

## 功能特性

- ✅ Needleman-Wunsch全局序列比对
- ✅ Smith-Waterman局部序列比对
- ✅ 支持自定义评分参数
- ✅ 包含BLOSUM62和PAM250评分矩阵
- ✅ 详细的比对结果输出

## 安装要求

- Python 3.7+
- NumPy

## 快速开始

### 安装依赖

```bash
pip install numpy
```

### 运行示例

```bash
cd src
python main.py
```

### 使用方法

```python
from alignment import Alignment

# 创建比对器实例
aligner = Alignment(
    match_score=2,      # 匹配得分
    mismatch_penalty=-1, # 错配惩罚
    gap_penalty=-2      # 空位惩罚
)

# Needleman-Wunsch全局比对
seq1 = "ACGT"
seq2 = "AGT"
result = aligner.needleman_wunsch(seq1, seq2)
print(f"全局比对得分: {result['score']}")
print(f"比对结果: {result['alignment']}")

# Smith-Waterman局部比对
result = aligner.smith_waterman(seq1, seq2)
print(f"局部比对得分: {result['score']}")
print(f"比对结果: {result['alignment']}")
```

## 项目结构

```
.
├── src/
│   ├── alignment.py          # 核心比对算法
│   ├── scoring_matrices.py   # 评分矩阵定义
│   ├── main.py               # 主程序入口
│   └── __init__.py           # 模块初始化
├── report.md                 # 项目报告
├── poster.md                 # 海报内容
└── README.md                 # 项目说明
```

## 算法复杂度

| 算法             | 时间复杂度 | 空间复杂度 |
| ---------------- | ---------- | ---------- |
| Needleman-Wunsch | O(mn)      | O(mn)      |
| Smith-Waterman   | O(mn)      | O(mn)      |

## 参考文献

1. Needleman, S. B., & Wunsch, C. D. (1970). A general method applicable to the search for similarities in the amino acid sequence of two proteins.
2. Smith, T. F., & Waterman, M. S. (1981). Identification of common molecular subsequences.

## 作者

**李泉** - 2026年6月
