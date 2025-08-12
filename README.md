# Identification of Essential Proteins

This repository implements the method described in the paper **"Identification of Essential Proteins Based on Improved HITS Algorithm."**

---

## Overview

Essential proteins are crucial for the survival and reproduction of organisms. Their loss often results in lethality or infertility, making them indispensable for cellular life. Compared to non-essential proteins, essential proteins tend to be more conserved throughout biological evolution. These proteins serve as promising targets for developing new drugs and vaccines aimed at treating and preventing diseases.

Advancements in high-throughput technologies—such as yeast two-hybrid systems and mass spectrometry—have generated extensive protein–protein interaction (PPI) datasets. These datasets enable the study of essential proteins at the network level.

---

## Implimentation stages

1. **Construct a bidirectional Protein-Protein Interaction (PPI) network** using the Gavin_PPI dataset.  
2. **Assign weights to network edges** based on:  
   - **Topological features:** Edge Clustering Coefficient of the graph nodes.  
   - **Biological information:** Pearson Correlation Coefficient (PCC) computed from gene expression data.  
3. **Apply the Hypertext Induced Topic Search (HITS) algorithm** to identify essential proteins (key nodes) within the PPI network.  
4. **Evaluate the method** using statistical metrics including:  
   - Sensitivity (SN)  
   - Specificity (SP)  
   - Positive Predictive Value (PPV)  
   - F-measure (F)  
   - Accuracy (ACC)  

---

## Datasets

- Saccharomyces cerevisiae PPI network from the **Gavin database**.  
- Gene expression data of *Saccharomyces cerevisiae* from the **GEO database**.  
- Known essential protein list from the **Saccharomyces Genome Database**.

