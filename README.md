# Auto-VTNA: User-friendly tool for carrying out automatic VTNA
Auto-VTNA is Python package that provides automatic variable time normalisation analysis (VTNA), offering several advantages to the traditional methods, such as:
- Manual and automatic VTNA in a more time efficient manner.
- VTNA with several normalised reaction species
- Visualisation of overlay score across different reaction orders
- Quantification of error analysis
- Improved accessibility via a graphical user interface (can be downloaded from https://drive.google.com/file/d/1PeoumECdGHQeyxp1p1m_wrZpIlOLTKXP/view?usp=sharing)

## Background
VTNA is a modern kinetic analysis method pioneered by Jordi Bures based on normalising the time axis of experiments with respect to the concentration profile of reactants, catalysts or additives whose intial concentration is varied between experiments raised to a selected order value. The correct order value is the one that gives concentration profile overlay: https://onlinelibrary.wiley.com/doi/abs/10.1002/anie.201609757. VTNA is traditionally carried out in Excel, for example with an Excel spreadsheet tool which can be downloaded from the supplementary information in a review by Jordi Bures and Christian Nielsen on Visual Kinetic Analysis: https://pubs.rsc.org/en/content/articlelanding/2019/sc/c8sc04698k. 
The kinetic data required for conventional VTNA analysis can be collected by setting up and monitoring reactions via concentration measurements for a standard set of initial concentrations (Exp. 1 in Table 1). Next, different excess experimentsin which one initial concentration is altered at the time are set up and monitored (Exp. 2-4 in Table 1). The time axis of the standard experiment and the relevant different experiment is then normalised with respect to the concentration profile of the relevant reaction species raised to its order value using numerical integration via the trapezoid rule (Exp. 1 and 2 for A, Exp. 1 and 3 for B, Exp. 1 and 4 for cat). The correct reaction order is the one that causes the concentration profiles to overlay.

Table 1: Standard VTNA experiment design.  
| Experiment | [A]<sub>0</sub> | [B]<sub>0</sub> | [cat]<sub>0</sub> | 
|-----------------|-----------------|-----------------|-----------------|
| Exp 1: Standard | 2.0 M   | 2.0 M  | 0.01 M    | 
| Exp 2: Different excess in A  | 1.0 M   | 2.0 M  | 0.01 M  | 
| Exp 3: Different excess in B  | 2.0 M | 1.0 M    | 0.01 M   | 
| Exp 4: Different excess in P   | 2.0 M | 2.0 M    | 0.005 M   | 

## What is Auto-VTNA?
Auto-VTNA utilises a robust method to quantify the degree of concentration profile overlay based on total monotonic or ordinary polynomial fitting of all normalised reaction progress profiles. This computational overlay score enables VTNA to be carried out automatically via an algorithm which identifies the reaction order values that otpimise the overlay score. Importantly, the automatic VTNA algorithm can identify the reaction order of several reaction species in one calculation if the initial concentration of each reaction species is varied at least once in the kinetic dataset (automatic total VTNA). For example, the reaction orders in A, B and cat can be identified in one calculation by selecting all experiments and both 'A', 'B' and 'cat' as normalised reaction species, and does not require only one initial concentration to be varied at the time. 
The conventional approach of identifying one reaction order value at the time can also be performed automatically by selecting the experiments in which only the initial concentration of the relevant reaction species is altered (automatic sequential VTNA). Auto-VTNA also offers significantly improved analysis of the uncertainty of the calculated order value(s) as well as visualisation of the inputed kinetic data and results from automatic or ordinary VTNA calculations. A more detailed description of the Python package can be found in this pre-print: https://chemrxiv.org/engage/chemrxiv/article-details/65fddc2d66c1381729948bb2 and the corresponding ESI. 

## The Automatic VTNA Calculator
To give users access to the advantages of Auto-VTNA without the requirement of coding experience, a graphical user interface (GUI) was developed using PySimpleGui named the Automatic VTNA Calculator. This file can be found as Auto_VTNA_GUI.py in the Python package or simply downloaded as an executible application (currently from https://drive.google.com/file/d/1p6gPdmGTt2o32ueVD2Jnts269cob1koo/view?usp=sharing). As the GUI is written to help the user access the data visualisation, conventional and automatic VTNA functionality and other features of Auto_VTNA in only a few mouse clicks, it is recommended also for users with coding experience. 

### NB: soon to be provided:
- Youtube tutorial on how to use the Automatic VTNA Calculator.
- Jupyter notebook tutorial files to illustrate the use of the Auto_VTNA Python package on simulated kinetic data.
- Next version of the package and GUI with some bugs fixed.
