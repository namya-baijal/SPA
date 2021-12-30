# South Pole-Aitken High-Resolution Oblique Numerical Simulations 

## About

GitHub repository for SPA impact crater ejecta analysis using iSALE3D hydrocode.

The main objective of this project is to determine the distribution of impact ejecta during the excavation stage of the impact. The model parameters for the analysis are listed in the table below. Two evaluation scenarios are analyzed: 50 K/km (hotter) with 10 km/s impact velocity and 10 K/km (colder) with 15 km/s impact velocity. Data from the first 400 steps of the simulations is used for the analysis.

![Screenshot from 2021-12-30 01-20-51](https://user-images.githubusercontent.com/57348392/147739919-f094713d-97d8-4b1d-867e-7ad381a88a1c.png)


The simulation data is available at [insert simulation data link here]. The filtered tracer data outputs from the full workflow is available at [insert csv file link here] and can be used to generate final visualizations.

## Running Python Notebooks

### Full Workflow (full_workflow_no_plots.ipynb)

This notebook filters the simulation data file and generates an output csv file for all ejecta tracers that land on the moon. This csv file contains all relevant ballistic parameters for analysis / plotting. Note: for large simulation data files, the compute required by this notebook is large.

### Plotting (Plotting.ipynb)

Takes in the mirror csv file and generates the ejecta distribution on mollweide projection. Tracers are filtered by ejection velocity and target and impactor
material.

![p1](https://user-images.githubusercontent.com/57348392/147739872-cce30c95-eccf-4c4e-82fa-8868822099ff.png)

![p2](https://user-images.githubusercontent.com/57348392/147739856-c223032d-88c2-41a8-bedb-80ee3b17415f.png)

### Layer Thickness Distribution (binning_histogram.ipynb)

This notebook creates 2x2 degree equal surface area bins. The layer thickness is determined by computing the histogram of landed ejecta volume and dividing by the surface area of the bins. The plots are filtered by depth of origin and shock pressure of ejecta particles.

![v15_troctolite](https://user-images.githubusercontent.com/57348392/147740488-fbd596f7-2791-4063-b8c6-fd624409760b.png)

### Radial Profile (radial_profile.ipynb)

This notebook creates ejecta thickness as a function of depth of origin along a 45 degree radial line from the SPA basin center with inverted stratigraphy. This notebook filters tracers with azimuth between 40-50 and bins these by radial distance. The histogram is divided by surface area as before to produce the layer thickness data.

![radial_plot_45_bearing_v15](https://user-images.githubusercontent.com/57348392/147741051-c093b199-d786-4faa-ac35-6d1faeb913bf.png)


### Find Location (find_loc.ipynb)

This notebook calculates the great circle distance and bearing between any two latitude longitude coordinates on a sphere using the haversine formula.



