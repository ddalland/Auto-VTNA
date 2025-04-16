# Auto-VTNA: User-friendly tool for carrying out automatic VTNA
Auto-VTNA is a Python package that provides automatic variable time normalisation analysis (VTNA), offering several advantages over traditional methods:
- Manual and automatic VTNA in a more time efficient manner;
- VTNA with several normalised reaction species;
- Visualisation of overlay score across different reaction orders;
- Quantification of uncertainty;
- Improved accessibility via a graphical user interface which can be downloaded from https://drive.google.com/file/d/11XVQQUvcE7aQkmxGXVwMmGLOyNtIpg8J/view?usp=sharing (version 3, RECOMMENDED)
- Older versions of the GUI:  https://drive.google.com/file/d/1wsXq9gAJyrCMphbnrqYXViKLNYKChxDn/view?usp=drive_link (version 2) or https://drive.google.com/file/d/1p6gPdmGTt2o32ueVD2Jnts269cob1koo/view?usp=drive_link (version 1).
- Mac version: https://drive.google.com/file/d/1f7olytkpr0IVIT7xbB41zup4rweVwIWc/view?usp=sharing
- Unfortunately, PySimpleGUI now requires that users pay for a commercial licence for using it. So the GUI is only accessible by paying for this, unless you got a hobbyist licence before it became unavailable. Hopefully, I will be able to remake a GUI with a different package in the future. However, the Python package and use of Auto-VTNA in Python is unaffected by this. 
Youtube tutorial on the first version of the Automatic VTNA Calculator GUI: https://youtu.be/FEDwvUY1Dl4?si=zHkratsrWyWnk0tW
Follow-up Youtube video on the second version of the GUI: https://youtu.be/6JAytwnDCRE

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

## Version 1.1.1:
On 24.09.2024, a new version of Auto-VTNA (1.1.1) was finalized and made available on PyPi. This version mostly contains improvements to the GUI, which are incorporated into its second version (see Google drive link above). The improvements include: 

- Cleaner layout, especially in the settings menus, with column elements ensuring that input boxes, etc., are aligned. Improved error handling. If the imported data is faulty (e.g., containing negative concentration values or inconsistent column headings), the user will be notified immediately. Calculations that would crash the GUI are stopped with an appropriate warning.
- Improved range mode data cropping. The user can now specify how the density of the kinetic data should be reduced, i.e., with respect to the time or one of the concentration axes. The option for automatically removing negative values, as well as shifting the time axis if the first data point(s) have been cropped, is included.
- The "Plot all profiles for selected species" option in the Inspect Kinetic Data window has also been updated to enable the simultaneous visualization of all concentration profiles for more than one selected reaction species.

Additionally, the class Same_excess has been added to the package. This class inputs two experiment sheets where one experiment has initial concentrations lower than the other by the same amount in each reactant. The time axis of the most dilute experiment is then shifted horizontally to "touch" the profile of the more concentrated experiment. The modified data can be visualized using the method plot_same_excess, which enables the user to adjust the horizontal shift as desired to generate a good same-excess plot.

## User manuals:

For user manuals for the python package and graphical user interface: See attachments in respository files. 
Also, check out the YouTube tutorial on using the Automatic VTNA Calculator: https://youtu.be/FEDwvUY1Dl4?si=zHkratsrWyWnk0tW as well as the video outlining the new features of the second version of the GUI executible: https://youtu.be/6JAytwnDCRE

