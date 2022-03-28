# Identification-of-Essential-Proteins
<b>Implement the paper "Identification of Essential Proteins Based on Improved HITS Algorithm".</b>

Essential proteins are vital for the life and reproduction of organisms, 
and play a very important role in maintaining cellular life. If the 
destruction of a certain protein would cause lethality or infertility, it can 
be said that it is essential to an organism, that is, the organism could not 
survive without it. Compared with non-essential proteins, essential 
proteins are more likely to remain in biological evolution. For example, 
essential proteins are excellent targets for the development of new 
potential drugs and vaccines for the treatment and prevention of 
diseases.
With the development of high-throughput technologies, such as yeast 
two-hybrid system and mass spectrometry analysis various protein–
protein interaction (PPI) data are available, which facilitates the studies 
of essential proteins at the network level.

<b>Task:</b>
<br/>
1. Build a Protein-Protein Interaction bidirectional Network using the
Gavin_PPI dataset.
2. Use the biological information and network topological features to 
weigh the edges separately:<br/>
  I. Use Edge Clustering Coefficient on nodes of graph for 
  topological features. <br/>
  II. Use Pearson Correlation Coefficient (PCC) on gene 
  expression data for biological information.
3. Use Hypertext Induced Topic Search (HITS) algorithm to identify 
essential nodes (Proteins).
4. Use several statistical measures, namely sensitivity (SN), 
specificity (SP), positive predictive value (PPV), F-measure (F), 
and accuracy (ACC) to evaluate your method.

<b>Datasets: </b>
<br/>
➢ Set of Saccharomyces cerevisiae PPI network from Gavin 
database.<br/>
➢ The gene expression data of Saccharomyces cerevisiae from GEO 
database.<br/>
➢ The list of known essential proteins collected from Saccharomyces 
Genome Database.<br/>
