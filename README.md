# Auto-VTNA: User-friendly tool for carrying out automatic VTNA
Auto-VTNA is a Python package that provides automatic variable time normalisation analysis (VTNA), offering several advantages over traditional methods:
- Manual and automatic VTNA in a more time efficient manner;
- VTNA with several normalised reaction species;
- Visualisation of overlay score across different reaction orders;
- Quantification of uncertainty;
- Improved accessibility via a graphical user interface which can be downloaded from https://drive.google.com/file/d/1ovT8l8WbG6TWcwYMoXxbiizgpe542vOC/view?usp=drive_link (version 2, recommended) or https://drive.google.com/file/d/1p6gPdmGTt2o32ueVD2Jnts269cob1koo/view?usp=drive_link (version 1). 

## Background
VTNA is a modern kinetic analysis method pioneered by Jordi Bures based on normalising the time axis of experiments with respect to the concentration profile of reactants, catalysts or additives whose intial concentration is varied between experiments raised to a selected order value. The correct order value is the one that gives concentration profile overlay: https://onlinelibrary.wiley.com/doi/abs/10.1002/anie.201609757. VTNA is traditionally carried out in Excel, for example with an Excel spreadsheet tool which can be downloaded from the supplementary information in a review by Jordi Bures and Christian Nielsen on Visual Kinetic Analysis: https://pubs.rsc.org/en/content/articlelanding/2019/sc/c8sc04698k. 
The kinetic data required for conventional VTNA analysis can be collected by setting up and monitoring reactions via concentration measurements for a standard set of initial concentrations (Exp. 1 in Table 1). Next, different excess experiments in which one initial concentration is altered at a time are set up and monitored (Exp. 2-4 in Table 1). The time axis of the standard experiment and the relevant different experiment is then normalised with respect to the concentration profile of the relevant reaction species raised to its order value using numerical integration via the trapezoid rule (Exp. 1 and 2 for A, Exp. 1 and 3 for B, Exp. 1 and 4 for cat). The correct reaction order is the one that causes the concentration profiles to overlay.

Table 1: Standard VTNA experiment design.  
| Experiment | [A]<sub>0</sub> | [B]<sub>0</sub> | [cat]<sub>0</sub> | 
|-----------------|-----------------|-----------------|-----------------|
| Exp 1: Standard | 2.0 M   | 2.0 M  | 0.01 M    | 
| Exp 2: Different excess in A  | 1.0 M   | 2.0 M  | 0.01 M  | 
| Exp 3: Different excess in B  | 2.0 M | 1.0 M    | 0.01 M   | 
| Exp 4: Different excess in P   | 2.0 M | 2.0 M    | 0.005 M   | 

## What is Auto-VTNA?
Auto-VTNA utilises a robust method to quantify the degree of concentration profile overlay based on total monotonic or ordinary polynomial fitting of all normalised reaction progress profiles. This computational overlay score enables VTNA to be carried out automatically via an algorithm which identifies the reaction order values that optimise the overlay score. Importantly, the automatic VTNA algorithm can identify the reaction order of several reaction species in one calculation if the initial concentration of each reaction species is varied at least once in the kinetic dataset (automatic total VTNA). For example, the reaction orders in A, B and cat can be identified in one calculation by selecting all experiments and both 'A', 'B' and 'cat' as normalised reaction species, and does not require only one initial concentration to be varied at a time. 
The conventional approach of identifying one reaction order value at a time can also be performed automatically by selecting the experiments in which only the initial concentration of the relevant reaction species is altered (automatic sequential VTNA). Auto-VTNA also offers significantly improved analysis of the uncertainty of the calculated order value(s) as well as visualisation of the inputted kinetic data and results from automatic or ordinary VTNA calculations. A more detailed description of the Python package can be found in this pre-print: https://chemrxiv.org/engage/chemrxiv/article-details/65fddc2d66c1381729948bb2 and the corresponding ESI. 

## The Automatic VTNA Calculator
To give users access to the advantages of Auto-VTNA without the requirement of coding experience, a graphical user interface (GUI) named the Automatic VTNA Calculator was developed using PySimpleGUI. This file can be found as Auto_VTNA_GUI.py in the Python package or simply downloaded as an executable application (currently available from https://drive.google.com/file/d/1p6gPdmGTt2o32ueVD2Jnts269cob1koo/view?usp=sharing). As the GUI is designed to help users access data visualization, conventional and automatic VTNA functionality and other features of Auto_VTNA in only a few mouse clicks, it is recommended even for users with coding experience. 

## Dependencies of the Auto-VTNA python package:
- pandas==2.2.2
- numpy==1.26.4
- matplotlib==3.9.0
- num2words==0.5.13
- mplcursors==0.5.3
- scipy==1.13.1
- polyfit==1.0
- PySimpleGUI==5.0.5
- openpyxl==3.1.3

Note that these dependencies also work for several older versions of packages, but the versions above have been shown to work.

## User manuals:

For user manuals for the python package and graphical user interface: See attachments in respository files. 
Also, check out the YouTube tutorial on using the Automatic VTNA Calculator: https://youtu.be/FEDwvUY1Dl4?si=zHkratsrWyWnk0tW 

