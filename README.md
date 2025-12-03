# 🧬 Day 09: Differential Expression Analysis

## 📋 Overview

This project performs differential expression analysis on gene expression data to identify upregulated and downregulated genes between control and treatment groups.

## 🛠️ Technologies Used

- Python 🐍
- Pandas 📊
- NumPy 🔢

## 📂 Files

- `differential_expression_analysis.ipynb` - Main analysis notebook
- `dea_data.csv` - Input gene expression data
- `differential_expression_results.csv` - Output results

## 🔬 Analysis Steps

### 1️⃣ Data Loading

```python
import pandas as pd
import numpy as np
import math

df = pd.read_csv('dea_data.csv')
```

Imports necessary libraries and loads gene expression data from CSV file.

### 2️⃣ Group Separation

```python
control_cols = [col for col in df.columns if 'Control' in col]
treatment_cols = [col for col in df.columns if 'Treatment' in col]
```

Separates columns into control and treatment groups for comparison.

### 3️⃣ Mean Calculation

```python
df['Control Mean'] = df[control_cols].mean(axis=1)
df['Treatment Mean'] = df[treatment_cols].mean(axis=1)
```

Calculates average expression levels for both control and treatment groups.

### 4️⃣ Log2 Fold Change

```python
epsilon = 1e-10
df['log2FC'] = np.log2((df['Treatment Mean'] + epsilon) / (df['Control Mean'] + epsilon))
```

Computes log2 fold change to measure expression differences between groups. Uses epsilon (1e-10) to avoid division by zero errors.

### 5️⃣ Gene Classification

```python
df['Status'] = 'Not Significant'
df.loc[df['log2FC'] >= 1, 'Status'] = 'Upregulated'
df.loc[df['log2FC'] <= -1, 'Status'] = 'Downregulated'
```

Categorizes genes based on fold change threshold:

- **Upregulated** 📈: log2FC ≥ 1 (2-fold increase)
- **Downregulated** 📉: log2FC ≤ -1 (2-fold decrease)
- **Not Significant** ➖: |log2FC| < 1

### 6️⃣ Results Display

```python
result_df = df[['Gene', 'Control Mean', 'Treatment Mean', 'log2FC', 'Status']]
print(result_df)
```

Creates a clean DataFrame with essential columns and displays formatted results.

### 7️⃣ Summary Statistics

```python
print(f"Total genes analyzed: {len(result_df)}")
print(f"Upregulated genes: {sum(result_df['Status'] == 'Upregulated')}")
print(f"Downregulated genes: {sum(result_df['Status'] == 'Downregulated')}")
```

Provides counts of upregulated, downregulated, and non-significant genes.

### 8️⃣ Top Differentially Expressed Genes

```python
top_upregulated = result_df[result_df['Status'] == 'Upregulated'].sort_values('log2FC', ascending=False).head(10)
top_downregulated = result_df[result_df['Status'] == 'Downregulated'].sort_values('log2FC', ascending=True).head(10)
```

Identifies and displays the top 10 most upregulated and downregulated genes.

### 9️⃣ Export Results

```python
result_df.to_csv('differential_expression_results.csv', index=False)
```

Saves final results to CSV file for further analysis.

## 📊 Output

The analysis generates:

- ✅ Complete gene expression results table
- 📈 Summary statistics
- 🏆 Top differentially expressed genes
- 💾 CSV export of results

## 🚀 How to Run

1. Ensure `dea_data.csv` is in the same directory
2. Open `differential_expression_analysis.ipynb` in Jupyter/VS Code
3. Run all cells sequentially
4. Check `differential_expression_results.csv` for results

## 📝 Notes

- Log2 fold change threshold of ±1 represents 2-fold change in expression
- Small epsilon value prevents mathematical errors with zero values
- Higher |log2FC| values indicate stronger differential expression

---

💡 **Tip**: Adjust the fold change threshold in the classification step to modify sensitivity of gene categorization.
