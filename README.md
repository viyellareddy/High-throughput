

## **Title:** miRNA Cancer Expression Analysis using Pearson Correlation

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

### **Heatmaps:**
1. **Heatmap of Matrix 1:**
   ![Matrix 1 Heatmap](path/to/matrix1_heatmap.png)
   
2. **Heatmap of Matrix 2:**
   ![Matrix 2 Heatmap](path/to/matrix2_heatmap.png)

### **Conclusion:**
The moderate Pearson correlation (**0.2976**) between the two miRNA expression matrices suggests that there are some similarities in the cancer types' miRNA expression patterns, but they are not highly correlated. This may indicate biological differences in the datasets or varying levels of miRNA expression across the cancers.

