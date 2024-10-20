

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


# Assignment 3: Operon Prediction

## Objective

The objective of this assignment is to predict operons, defined as the longest contiguous multi-gene transcriptional units, using PTT files for various genomes. Additionally, we will predict operons in a crop microbiome derived from the Hoatzin.

### Genomes to Analyze

1. **Escherichia coli K12**
2. **Bacillus subtilis**
3. **Halobacterium**
4. **Synechocystis**

### Additional Analysis

- Predict operons in the crop microbiome from Hoatzin (IMG ID: 2088090036) using the attached GFF file. The first column comprises the contig, while the fourth and fifth columns represent the gene start and stop positions.

## Data Files

- **E_coli_K12_MG1655.ptt**
- **B_subtilis_168.ptt**
- **Halobacterium_NRC1.ptt**
- **Synechocystis_PCC6803_uid159873.ptt**
- **2088090036.gff** (for Hoatzin microbiome)

## Code Implementation

The following Python code performs the operon prediction using the provided datasets.

### Step 1: Load PTT Data

```python
import pandas as pd

# Load PTT data for Escherichia coli K12
data1 = pd.read_csv('E_coli_K12_MG1655.ptt', sep='\t', skiprows=[0, 1])
```

### Step 2: Define Operon Prediction Function

The function `predict_operons(data)` identifies operons based on the defined criteria.

```python
def predict_operons(data):
    operons = []
    current_operon = []
    for _, row in data.iterrows():
        if current_operon == []:
            current_operon.append(row)
        elif row['Strand'] == current_operon[-1]['Strand'] and int(row['Location'].split('..')[0]) - int(current_operon[-1]['Location'].split('..')[1]) < 50:
            current_operon.append(row)
        else:
            operons.append(current_operon)
            current_operon = [row]
    operons.append(current_operon)
    return operons
```

### Step 3: Predict Operons for Each Genome

This section demonstrates how to predict operons for each genome.

```python
# Predict operons for Escherichia coli K12
operons = predict_operons(data1)
for i, operon in enumerate(operons):
    print(f"Operon {i+1}: {', '.join(gene['Gene'] for gene in operon)}")
```

Repeat the above steps for **Bacillus subtilis**, **Halobacterium**, and **Synechocystis** by loading the respective PTT files and calling the `predict_operons` function.

### Step 4: Load GFF Data for Hoatzin Crop Microbiome

```python
# Load GFF data for Hoatzin crop microbiome
data = pd.read_csv('2088090036.gff', sep='\t', header=None, comment='#')
data.columns = ['Contig', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Attributes']
# Sort data by contig and start position
data = data.sort_values(['Contig', 'Start'])
```

### Step 5: Predict Operons in Hoatzin Data

```python
# Extract gene names from attributes column
data['Gene'] = data['Attributes'].apply(lambda x: x.split(';')[0].split('=')[1])

# Predict operons for Hoatzin data
operons = predict_operons(data)
# Print predicted operons
for i, operon in enumerate(operons):
    print(f"Operon {i+1}: {', '.join(row['Gene'] for row in operon)}")
```

## Results

The results of the operon predictions for each genome will be printed to the console, listing the genes within each identified operon.



# Assignment 4: Motif Finding with Position Weight Matrix (PWM)

## Overview
This project focuses on motif finding using a Position Weight Matrix (PWM) to analyze the binding sites of the transcription factor argR. The goal is to compute the frequency and weight matrices from a provided counts matrix, and then use these matrices to identify potential binding sites in the upstream regulatory regions of genes. The results include the top 30 gene IDs that show the highest similarity to the argR binding motif.

## Requirements
- Python 3.x
- Pandas
- NumPy

## Installation
1. Clone the repository:
   ```bash
   git clone <repository_url>
   cd <repository_name>
   ```

2. Install the required packages:
   ```bash
   pip install pandas numpy
   ```

## Files
- `argR-counts-matrix.txt`: Contains the counts matrix representing the frequency of each base at specified positions derived from 27 documented binding sites.
- `E_coli_K12_MG1655.400_50`: Contains upstream regulatory regions of genes in FASTA format, including Gene IDs and corresponding sequences.
- `motif_finding.py`: The main Python script that performs the motif finding analysis.

## Usage
1. Ensure that the input files (`argR-counts-matrix.txt` and `E_coli_K12_MG1655.400_50`) are in the same directory as `motif_finding.py`.

2. Run the script:
   ```bash
   python motif_finding.py
   ```

3. The output will display the top 30 gene IDs with their corresponding scores, indicating the strength of binding site similarity to the argR motif.

## Code Explanation
1. **Loading Data**: The counts matrix is loaded and processed to compute the frequency matrix \( F(b, j) \) with pseudocounts added to avoid logarithmic issues.

2. **Frequency Matrix Calculation**:
   - A pseudocount of +1 is added to each count to calculate the adjusted frequency matrix \( F'(b, j) \).
   - The background frequency is assumed to be 0.25 for all bases.

3. **Weight Matrix Computation**: The weight matrix is calculated using the log-odds ratio of the frequency matrix to the background frequency.

4. **Scanning for Binding Sites**:
   - The upstream sequences of genes are read, and each sequence is scanned for potential argR binding motifs.
   - The highest motif scores are recorded for each gene.

5. **Results Output**: The top 30 gene IDs are sorted and displayed based on their motif scores.

## Conclusion
This project provides insights into the binding behavior of the argR transcription factor, aiding in the understanding of gene regulation mechanisms. The identified binding sites can serve as a foundation for further experimental validation and biological studies.


# Assignment 5: Protein Interaction Network Analysis

## Overview

This assignment involves analyzing the human protein-protein interaction (PPI) network to calculate key network metrics such as degree, clustering coefficient, and shortest path lengths between nodes (proteins). We will also visualize the degree distribution to assess the network's scale-free structure and compare path length distributions between two sets of proteins using statistical tests.

## Files Required

1. **`Human-PPI.txt`**: A text file containing the edge list of the protein interaction network.
2. **`Protein-list1.txt`**: A text file containing a list of proteins for the first set.
3. **`Protein-list2.txt`**: A text file containing a list of proteins for the second set.

## Dependencies

The script requires the following Python packages:

- `numpy`
- `pandas`
- `networkx`
- `matplotlib`
- `scipy`

To install these packages, use the following command:

```bash
pip install numpy pandas networkx matplotlib scipy
```

1. **Network Metrics Calculation**: Computes the degree, clustering coefficient, and plots the degree distribution.
2. **Shortest Path Length Analysis**: Calculates the shortest path lengths between nodes in the two provided protein lists and compares the distributions.

## Code Explanation

### Part 1: Network Metrics Calculation

1. **Importing Libraries**:
    ```python
    import os
    import networkx as nx
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    ```

   - Imports necessary libraries for data manipulation, network analysis, and visualization.

2. **Loading the Edge List**:
    ```python
    network_df = pd.read_csv('Human-PPI.txt', sep='\s+')
    ```
   - Loads the edge list of the PPI network into a pandas DataFrame.

3. **Creating the Graph**:
    ```python
    graph = nx.Graph()
    for _, row in network_df.iterrows():
        graph.add_node(row['OFFICIAL_SYMBOL_A'])
        graph.add_node(row['OFFICIAL_SYMBOL_B'])
        graph.add_edge(row['OFFICIAL_SYMBOL_A'], row['OFFICIAL_SYMBOL_B'])
    ```

   - Initializes an undirected graph and adds nodes and edges based on the edge list.

4. **Calculating Degree and Clustering Coefficient**:
    ```python
    degrees = nx.degree(graph)
    clustering = nx.clustering(graph)
    avg_clustering = nx.average_clustering(graph)
    ```

   - Computes the degree and clustering coefficient for each node and the average clustering coefficient of the network.

5. **Degree Distribution Plot**:
    ```python
    degree_values = np.array(list(degrees.values()))
    degree_freq = np.bincount(degree_values)
    plt.loglog(np.arange(len(degree_freq)), degree_freq, 'o')
    ```

   - Extracts degree values and their frequency, then plots the degree distribution on a log-log scale to visualize the scale-free structure.

### Part 2: Shortest Path Length Analysis

1. **Loading Protein Lists**:
    ```python
    proteins1 = pd.read_csv('Protein-list1.txt')
    proteins2 = pd.read_csv('Protein-list2.txt')
    ```

   - Loads the two protein lists into pandas DataFrames and renames columns for clarity.

2. **Calculating Shortest Path Lengths**:
    ```python
    subgraph1_path_lengths = []
    for i in range(len(proteins1)-1):
        for j in range(i+1, len(proteins1)):
            p1 = proteins1['protein A'][i]
            p2 = proteins1['protein A'][j]
            if p1 in list(graph.nodes()) and p2 in list(graph.nodes()):
                subgraph1_path_lengths.append(len(nx.shortest_path(graph,p1,p2)))
    ```

   - Loops through each pair of proteins in the first list to compute the shortest path lengths. A similar process is followed for the second protein list.

3. **Statistical Comparison of Path Length Distributions**:
    ```python
    from scipy.stats import ranksums
    wilcoxon_test = ranksums(subgraph1_path_lengths, subgraph2_path_lengths)
    ```

   - Performs a Wilcoxon rank-sum test to compare the path length distributions between the two protein sets.

## Outputs

- The average clustering coefficient of the network is printed.
- A log-log plot of the degree distribution is displayed.
- The result of the Wilcoxon rank-sum test is printed, indicating whether there is a statistically significant difference between the path lengths of the two protein sets.

## Conclusion

This analysis provides insights into the structure of the human protein interaction network, its degree distribution, and the comparative analysis of shortest path lengths between different protein sets. The results can inform further studies in protein interactions and their biological significance.

