# Auto-VTNA: User-friendly tool for carrying out automatic VTNA analysis. 
Auto-VTNA is Python package that provides automatic variable time normalisation analysis (VTNA), offering several advantages to the traditional methods, such as:
- Manual and automatic VTNA in a more time efficient manner.
- VTNA with several normalised reaction species
- Visualisation of overlay score across different reaction orders
- Quantification of error analysis
- Improved accessibility via a graphical user interface (can be downloaded from https://drive.google.com/file/d/1Rf0Hbi5fWXm2KILTKx5eOIQMjuHAvbPp/view?usp=drive_link)

## Background
VTNA is a modern kinetic analysis method pioneered by Jordi Bures based on normalising the time axis of experiments with respect to the concentration profile of reactants, catalysts or additives whose intial concentration is varied between experiments raised to a selected order value. The correct order value is the one that gives concentration profile overlay: https://onlinelibrary.wiley.com/doi/abs/10.1002/anie.201609757. VTNA is traditionally carried out in Excel, for example with an Excel spreadsheet tool which can be downloaded from the supplementary information in a review by Jordi Bures and Christian Nielsen on Visual Kinetic Analysis: https://pubs.rsc.org/en/content/articlelanding/2019/sc/c8sc04698k.

## What is Auto-VTNA?
Auto-VTNA utilises a robust method to quantify the degree of concentration profile overlay based on total monotonic or ordinary polynomial fitting of all normalised concentration profiles. This enables the automatic determination of the reaction order value that maximises concentration profile overlay via an efficient algorithm capable of determining the order value of one or several reaction species simultaneously, provided kinetic data in which the initial concentrations of the relevant reactants, additives or catalysts have been altered between experiments. Auto-VTNA also offers significantly improved analysis of the uncertainty of the calculated order value(s) as well as visualisation of the inputed kinetic data and results from automatic or ordinary VTNA calculations. A more detailed description of the Python package can be found in this pre-print: https://chemrxiv.org/engage/chemrxiv/article-details/65fddc2d66c1381729948bb2 and the corresponding ESI. 

## The Automatic VTNA Calculator
To give users access to the advantages of Auto-VTNA without the requirement of coding experience, a graphical user interface (GUI) was developed using PySimpleGui named the Automatic VTNA Calculator. This file can be found as Auto_VTNA_GUI.py in the Python package or simply downloaded as an executible application (currently from https://drive.google.com/file/d/10Tp8gTcDM4zp7xC7g6l0A6xrT56eidsN/view?usp=drive_link). As the GUI is written to help the user access the data visualisation, conventional and automatic VTNA functionality and other features of Auto_VTNA in only a few mouse clicks, it is recommended also for users with coding experience. 

### NB: soon to be provided:
- Youtube tutorial on how to use the Automatic VTNA Calculator.
- Jupyter notebook tutorial files to illustrate the use of the Auto_VTNA Python package on simulated kinetic data.
- Next version of the package and GUI with some bugs fixed.
