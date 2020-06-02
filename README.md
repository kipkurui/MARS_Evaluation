# MARS Evaluation

MARS_Evaluation is a repository that host Jupyter notebooks and code necessary to reproduce the evaluations of [MARS](www.bioinf.ict.ru.ac.za) and the [MARSTools](https://github.com/kipkurui/MARS_Suite). The MARSTools should have been forked from the link provided. 

However, analysis for [Chapter 3](Chapter3) ships with its code.

First, you will need to install the requirements from the requirements files. We use conda environment (.yml), but a list is also provided for Virtual env (requirements.txt). There are two types evaluations here:

### Additional_Data
Contains additional data necessary for a comprehensive analysis. See the README file in the folder for more details. 

### Chapter3
Folder contains Ipython notebooks and Data for the systematic comparative evaluation of the motif assessment techniques.  A README in the folder provides more details of the folder. 
    * Chapter_3.ipynb: main code for the chapter
    * Run_chip_seq_analysis.ipynb
    * Run_PBM_analysis.ipynb

### MARS_Evaluation
MARS_Evaluation the central repository for the complete breakdown of the analysis undertaken for the MARSTools.

#### Data
Files contain a list of TFs used in the analysis with the  TF names mapped to their TF IDS. 

#### Example_data
Contains data utilised for the work through the examples in the Ipython notebook. 