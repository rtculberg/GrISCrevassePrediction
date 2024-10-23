# GrISCrevassePrediction
Workflow for predicting the extent of crevasses in Greenland's ice slabs using a logistic regression model. 

These scripts were used to generate the results published in Culberg, Lai & MacKie (202X) "von Mises stress a robust predictor of ice slab fracture in Greenland". The instructions below describe the library and data dependencies, as well as the order in which scripts should be run to reproduce the results.    

Library Dependencies:  
These scripts are mix of MATLAB scripts and Jupyter Notebooks written in Python. The MATLAB scripts were developed under version R2022b Update 7 on a 64-bit Windows system. The Mapping Toolbox and Image Processing Toolbox are required. The Python Jupyter Notebooks were developed using Python 3.12.2 in a conda virtual environment. The environment.yml provides a list of the other library dependencies for easy install. (Note that environment.yml install may not work correctly if you are not installing to a Windows 64-bit machine.) 

To fully reproduce the paper results, scripts should be run in the order they are provided (a-j). Details of each script's functionality and data dependencies are given below:   

### a_GenerateRGF.ipynb   
This script generates random Gaussian fields with exponential or gaussian covariance structures and a 5km isotropic correlation length.  
- Data Dependencies: None    
- Package Dependencies:     
	- matplotlib    
	- numpy   
	- pandas   
	- scikit-gstat   
	- GStatSim     

### b_GenerateStrainMaps.m     
This script generates strain rate maps for the western and eastern training areas where crevasse observations are available. It is specifically designed to generate an ensemble of strain rate maps from different velocity data sources with different smoothing kernels and smoothing length scales in supporting of the optimization experiments shown in Supplementary Figure S1.    
- Data Dependencies:   
	- X component of surface velocity (available in accompanying dataset)   
	- Y component of surface velocity (available in accompanying dataset)     
	- Ice thickness (available in accompanying dataset)   
- Package Dependencies: None (MATLAB)            

### c_ParameterOptimization.ipynb        
This script runs the parameter optimization experiments described in Section 2.3 and which are used to produced Supplementary Figure S1. It also runs the regional cross-validation of the optimal model described in Section 2.4 and produces the example confusion matrix shown in Supplementary Figure S3.     
- Data Dependencies:   
	- Surface velocities (available in accompanying dataset)    
	- Strain rates over the training area (produced by running Script B)    
	- Surface temperatures (available in accompanying dataset)    
	- Crevasse observations (available in accompanying dataset)   
	- Test set definitions (available in accompanying dataset)   
	- K-fold cross validation fold definitions (available in accompanying dataset)   
- Package Dependencies:    
	- numpy    
	- scipy   
	- pandas   
	- matplotlib   
	- scikit-learn    
	- rasterio    

### d_CalculateFractureOrientation.m     
This script calculates the dominant orientation of fractures in the northwest training areas from individually segmented crevasses within velocity grid cells, as well as the orientation of the maximum principal stress in that grid cell. It outputs a csv file with crevasse orientation and stress orientation + 90 degrees for each grid cell with crevassing. These results are used to plot Figure 2c. Included in ./data/FractureOrientation.csv is the original version of this output text file.     
- Data Dependencies:     
	- GeoTIFF files with segmented crevasses in northwest Greenland (available in accompanying dataset)   
	- Strain rate maps over northwest Greenland (produced by Script B)   
	- Annual surface temperatures over northwest Greenland (available in accompanying dataset)     
- Package Dependencies: None (MATLAB)        

### e_PlotDonuts.ipynb      
This script reproduces the individual panels in Figure 2.  
- Data Dependencies:     
	- Surface velocities (available in accompanying dataset)    
	- Strain rates over the training area (produced by running Script B)    
	- Surface temperatures (available in accompanying dataset)    
	- Crevasse observations (available in accompanying dataset)   
	- Fracture orientation observations (generated by Script D)
- Package Dependencies:    
	- numpy    
	- scipy   
	- pandas   
	- matplotlib   
	- scikit-learn    
	- rasterio   
	- mpl-scatter-density     

### f_GeneratePerturbedStrainMaps.m     
This script generate a series of 20 strain rate maps over the training area from 20 velocity maps produced by perturbing the original velocity observations by the RGFs scaled by the velocity error.   
- Data Dependencies:   
	- X and Y surface velocity components (available in accompanying dataset)    
	- X and Y velocity error components (available in accompanying dataset)    
	- Random gaussian fields (can be generated using Script A)  
- Package Dependencies: None (MATLAB)       

### g_TensileStrengthPlots.ipynb       
This script runs the 2 x 120 model ensembles that quantify uncertainty in ice slab tensile strength and produces Figure 4. 
- Data Dependencies:      
	- Strain rate maps for the training area produced from perturbed velocities (generated by Script F)        
	- Surface temperatures from MAR, RACMO, and CARRA (available in accompanying dataset)   
	- Crevasse observations (available in accompanying dataset)   
	- Test region definitions (available in accompanying dataset)   
- Package Dependencies:    
	- numpy    
	- scipy   
	- pandas   
	- matplotlib   
	- scikit-learn    
	- rasterio    

### h_GenerateStrainMaps_AllGrIS.m        
This script generate a series of 20 strain rate maps over the entire ice sheet from 20 velocity maps produced by perturbing the original velocity observations by the RGFs scaled by the velocity error.   
- Data Dependencies:   
	- X and Y surface velocity components (available in accompanying dataset)    
	- X and Y velocity error components (available in accompanying dataset)    
	- Random gaussian fields (can be generated using Script A)  
- Package Dependencies: None (MATLAB)     

### i_AllGreenlandStress.ipynb       
This script generates the data final fracture probability map shown in Figure 4 and saves the results to a text file. 
- Data Dependencies:      
	- Strain rate maps for the entire ice sheet produced from perturbed velocities (generated by Script H)        
	- Surface temperatures from MAR, RACMO, and CARRA (available in accompanying dataset)   
	- Crevasse observations (available in accompanying dataset)   
	- Test region definitions (available in accompanying dataset)
	- Optimized tensile strengths estimates for each model in the ensemble
- Package Dependencies:    
	- numpy    
	- scipy   
	- pandas   
	- matplotlib     
	- rasterio      
	
### j_txt2tif.m   
This script takes the text file of fracture probabilities and generates a GeoTIFF for visualization.       
- Data Dependencies:   
	- Text files of fracture probabilities generated by Script I.    
- Package Dependencies:        
	- numpy 
	- scipy 
	- pandas  
	- matplotlib 
	- rasterio      