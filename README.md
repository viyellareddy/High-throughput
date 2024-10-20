

## **Title:** Highthroughput assignments

## Assignment 1: miRNA Cancer Expression Analysis using Pearson Correlation

### **Overview:**
This project analyzes miRNA expression data across various cancer types using Pearson correlation. The analysis aims to understand the similarity between different cancer types based on their miRNA expression profiles. Two datasets, referred to as Matrix 1 and Matrix 2, are used to compute pairwise correlations between cancer types. Heatmaps are generated to visualize these correlations, and a final Pearson correlation is computed to compare the two datasets.

### **Dataset:**
- **Matrix 1 and Matrix 2:** 
  These matrices contain RPKM (Reads Per Kilobase Million) values representing miRNA expression levels across 12 cancer types for a specific patient group. Each matrix is a 12x12 table, with each column corresponding to a different cancer type.

### **Objectives:**
1. Compute a Pearson correlation matrix (12x12) for the cancer types in both Matrix 1 and Matrix 2.
2. Visualize the correlation matrices as heatmaps to reveal the similarity between the cancer types.
3. Compute the Pearson correlation between the two matrices to compare their overall similarity.

### **Key Steps:**
1. **Data Preparation:** Load and clean the miRNA expression datasets (`matrix1.txt` and `matrix2.txt`), ensuring proper formatting for analysis.
2. **Correlation Computation:** 
   - Generate 12x12 Pearson correlation matrices to evaluate pairwise relationships between cancer types in Matrix 1 and Matrix 2.
3. **Heatmap Visualization:** 
   - Use Seaborn to create heatmaps of the correlation matrices, highlighting the similarities or differences between cancer types.
4. **Comparison Between Matrices:**
   - Calculate the Pearson correlation between the two correlation matrices (Matrix 1 and Matrix 2) to measure how similar these datasets are.

### **Requirements:**
- Python 3.x
- Libraries:
  - `pandas`
  - `numpy`
  - `seaborn`
  - `scipy`

To install the required libraries, use the following command:
```bash
pip install pandas numpy seaborn scipy
```

### **Running the Analysis:**
1. Clone the repository:
   ```bash
   git clone https://github.com/your-username/repo-name.git
   ```
2. Navigate to the project directory:
   ```bash
   cd repo-name
   ```
3. Open the Jupyter Notebook (`miRNA_cancer_analysis.ipynb`) to execute the analysis step-by-step.

### **Project Structure:**
```
├── matrix1.txt                   # Matrix 1 containing miRNA expression data
├── matrix2.txt                   # Matrix 2 containing miRNA expression data
├── miRNA_cancer_analysis.ipynb    # Jupyter notebook with the code for the analysis
├── README.md                     # Project description and instructions
```

### **Results:**
- **Correlation Matrices:** Generated 12x12 Pearson correlation matrices for both Matrix 1 and Matrix 2, representing the relationships between cancer types based on miRNA expression.
- **Heatmaps:** Created heatmaps to visualize the similarity between cancer types in each matrix.
- **Comparison:** The Pearson correlation between the two correlation matrices is **0.2976**, indicating a moderate similarity between the miRNA expression patterns in the two datasets.


### **Conclusion:**
The moderate Pearson correlation (**0.2976**) between the two miRNA expression matrices suggests that there are some similarities in the cancer types' miRNA expression patterns, but they are not highly correlated. This may indicate biological differences in the datasets or varying levels of miRNA expression across the cancers.



# Assignment: 2 Transcript Half-Life Calculation and Functional Enrichment Analysis

## Objective

This project aims to calculate the transcript half-lives for yeast genes using a 60-minute time series dataset. The analysis involves three sets of time series data to calculate the average half-life for each transcript, followed by the identification of genes with very high (top 10%) and very low (bottom 10%) half-lives. Finally, we perform a functional enrichment analysis on these identified genes.


## Analysis Steps

### Step 1: Load and Prepare the Data

We load the yeast decay time course data using pandas and extract the necessary columns for each time course replicate.

```python
import pandas as pd
import numpy as np

data = pd.read_csv('DecayTimecourse.txt', sep='\t')
timecourse_1 = pd.DataFrame(data[data.columns[1:10]])
timecourse_2 = pd.DataFrame(data[data.columns[10:19]])
timecourse_3 = pd.DataFrame(data[data.columns[19:28]])
```

### Step 2: Log Transformation and Regression Model

To compute the decay rate, we apply a logarithmic transformation and fit a regression model (Ordinary Least Squares) to estimate the decay rate (slope) for each transcript.

```python
log_timecourse_1 = np.log(timecourse_1)
import statsmodels.api as sm

results1 = {}
for i in range(len(log_timecourse_1.index)):
    if log_timecourse_1.iloc[i].isnull().any():
        x = log_timecourse_1.columns[~np.isnan(log_timecourse_1.iloc[i])]
        y = log_timecourse_1.iloc[i].dropna().values.reshape(-1, 1)
        x = sm.add_constant(x)
        model = sm.OLS(y, x).fit()
        slope = model.params[1]
        results1[log_timecourse_1.index[i]] = slope
```

### Step 3: Calculate Half-Lives

Using the slope obtained from the regression model, we calculate the half-life for each gene.

```python
half_life1 = {}
for i, slope in results1.items():
    half_life = np.log(2) / slope
    half_life1[i] = half_life
```

Repeat this process for the other two time courses to get `half_life2` and `half_life3`.

### Step 4: Calculate the Average Half-Life

Compute the average half-life for each transcript across the three replicates.

```python
average_half_life = np.mean([half_life1, half_life2, half_life3], axis=0)
```

### Step 5: Identify Top and Bottom 10% of Genes

Sort the genes based on their average half-lives and identify the top 10% and bottom 10%.

```python
df_sorted = df4.sort_values('average_half_life')
top_genes = df_sorted.head(int(len(df_sorted) * 0.1))
bottom_genes = df_sorted.tail(int(len(df_sorted) * 0.1))
```

### Step 6: Save Results

Save the top and bottom genes to a text file.

```python
top_genes.to_csv('top_genes.txt', index=False)
bottom_genes.to_csv('bottom_genes.txt', index=False)
```

## Functional Enrichment Analysis

To perform functional enrichment analysis using Gene Ontology (GO):

1. Open [g:Profiler](https://biit.cs.ut.ee/gprofiler/).
2. Load the `top_genes.txt` and `bottom_genes.txt` data to get the Gene Ontology of the respective genes.

You will receive two files containing the gene ontology of the top and bottom genes.

## Results

The results of the half-life calculations and the average values are stored in the `average_half_life` variable. The identified top and bottom genes are saved in `top_genes.txt` and `bottom_genes.txt`, respectively.

