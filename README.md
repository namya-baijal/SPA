# South Pole-Aitken High-Resolution Oblique Numerical Simulations 

## About

GitHub repository for SPA impact crater ejecta analysis using iSALE3D hydrocode.

The main objective of this project is to determine the distribution of impact ejecta during the excavation stage of the impact. The model parameters for the analysis are listed in the table below. Two evaluation scenarios are analysed: 50 K/km (hotter) with 10 km/s impact velocity and 10 K/km (colder) with 15 km/s impact velocity. Data from the first 400 steps of the simulations is used for the analysis.

![Screenshot from 2021-12-30 01-20-51](https://user-images.githubusercontent.com/57348392/147739919-f094713d-97d8-4b1d-867e-7ad381a88a1c.png)

The filtered tracer data outputs from the full workflow is available [here](https://drive.google.com/drive/folders/17YVryG2ypDe_iZ0UzQi58W-ngpHmO8pp?usp=sharing) and can be used to generate final visualizations. Check access with iSALE3D developers for simulation data.

## Running Python Notebooks

### Full Workflow (full_workflow_no_plots.ipynb)

This notebook filters the simulation data file and generates an output csv file for all ejecta tracers that land on the moon. This csv file contains all relevant ballistic parameters for analysis / plotting. Note: for large simulation data files, the compute required by this notebook is large.

### Plotting (Plotting.ipynb)

Takes in the mirror csv file and generates the ejecta distribution on mollweide projection. Tracers are filtered by ejection velocity and target and impactor
material.

![material_map_v10_mollweide](https://user-images.githubusercontent.com/57348392/149660945-32fdf00c-eafe-46a1-b3a8-1fe5ee1eb1cf.png)

![velocity_map](https://user-images.githubusercontent.com/57348392/149660965-21ba0297-4524-42a3-8d1f-1ade65aa43d4.png)

### Layer Thickness Distribution (binning_histogram.ipynb)

This notebook creates 2x2 degree equal surface area bins. The layer thickness is determined by computing the histogram of landed ejecta volume and dividing by the surface area of the bins. The plots are filtered by depth of origin and shock pressure of ejecta particles.



![v15_troctolite](https://user-images.githubusercontent.com/57348392/149661040-40835266-3aca-4bcf-9417-77300a87571a.png)


### Radial Profile (radial_profile.ipynb)

This notebook creates ejecta thickness as a function of depth of origin along a 45 degree radial line from the SPA basin center with inverted stratigraphy. This notebook filters tracers with azimuth between 40-50 and bins these by radial distance. The histogram is divided by surface area as before to produce the layer thickness data.

![45_azimuth](https://user-images.githubusercontent.com/57348392/149661052-6bd1bddc-a233-446e-ab98-75e07b76ebc0.png)

### Find Location (find_loc.ipynb)

This notebook calculates the great circle distance and bearing between any two latitude longitude coordinates on a sphere using the haversine formula.



