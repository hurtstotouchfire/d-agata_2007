Source
======

Individual Based Model (IBM) and a Differential Equation Model (DEM) examining the impact of treatment duration on antibiotic resistance in hospitals. Model from: [D'Agata, P. Magal, D. Olivier, S. Ruan, G.F. Webb, Journal of Theoretical Biology 249 (2007)](http://www.math.u-bordeaux1.fr/~pmagal100p/Nosocomial%20infection/nosocomial.htm)

DEM Instructions
================

In order to run the DEM model, please copy all the files into a directory called for example "C:\DEM-model\". Then by using MATLAB, you run the program "SimulationDEM.m". Note that in SimulationDEM.m you can modify all the parameters of the DEM. 

IBM Instructions
================

In order to run the IBM model, please copy all the files into a directory called for example "C:\IBM-model\". Then by using MATLAB, you run the program "SimulationIBM.m".

Note that in SimulationIBMmodel.m you can modify all the parameters of the IBM. You can also modify the directory where the figure created by the program is saved. You may for example create a directory C:\IBM-model\figures. But then in the file SimulationIBM.m you must replace the line 14 by 

    filename='C:\IBM-model\figures\name-of-the-figure.fig';


 
