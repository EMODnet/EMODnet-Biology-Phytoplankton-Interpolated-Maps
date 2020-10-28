# EMODnet-Biology-Phytoplankton-Interpolation


This directory provides the code for the reading of data files and performing spatial interpolation using `DIVAnd` software tool.

* [Heatmap](https://github.com/gher-ulg/EMODnet-Biology-Interpolated-Maps/blob/master/analysis/interp_presence_absence.ipynb) [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/gher-ulg/EMODnet-Biology-Interpolated-Maps/master?filepath=analysis%2Finterp_presence_absence.ipynb)
* [Neural network](https://github.com/gher-ulg/EMODnet-Biology-Interpolated-Maps/blob/master/analysis/emodnet_bio_DIVAndNN.ipynb) (binder does only provide 2GB of memory which is insufficient for this notebook)


## Directory structure

```
{{EMODnet-Biology-Interpolated-Maps}}/
├── analysis
├── data/
├── docs/
├── product/
│   ├── figures/
│   └── netCDF/
└── scripts/
```

* **analysis** - Jupyter notebooks used to perform the data analysis, create the figures and the `netCDF` files.
* **data** - contains a text file with the URLs of the datafiles.
* **docs** - Rendered reports
* **product** - Output product files: `netCDF` containing the gridded, probability fields and the corresponding figures in `PNG` format.
* **scripts** - Reusable code: functions employed in the Jupyter notebooks.

## Data

Data files have been produced by Deltares (Luuk van der Heijden, Willem Stolte).      
They consist of `CSV` files containing the dates, coordinates and occurences (presence of absence) of 200 species in the North Sea.      
The content of a data file looks like this (example for `Gymnodinium-1995-2020.csv`):
```bash
"abbr","date_year","genus","date_xUTM_yUTM","date","xUTM","yUTM","season","eventID","wint_year","occurs","gridnr","middleXgrid","middleYgrid"
"dome-phytoplankton",2010,"Gymnodinium","2010-01-11_477850.085718031_5756212.06783869",2010-01-11,477850.085718031,5756212.06783869,"winter",NA,NA,0,201,472500,5752500
"dome-phytoplankton",2010,"Gymnodinium","2010-01-11_560048.806065514_5746827.70084791",2010-01-11,560048.806065514,5746827.70084791,"winter",NA,NA,0,207,562500,5752500
"dome-phytoplankton",2010,"Gymnodinium","2010-01-11_563052.641129395_5744324.75671705",2010-01-11,563052.641129395,5744324.75671705,"winter",NA,NA,0,174,562500,5737500
"dome-phytoplankton",2010,"Gymnodinium","2010-01-11_742467.610873757_5940711.86346217",2010-01-11,742467.610873757,5940711.86346217,"winter",NA,NA,0,648,742500,5947500
```

## Analysis

This directory contains the notebooks for the preparation and analysis of the data.

* `scripts/emodnet_bio2020.jl`: compute the probability map of phytoplankton in the North Sea using `DIVAnd` and a neural network.
* `heatmap_testcase.ipynb`: several simple test cases to illustrate how the probability is derived from heatmaps.
* `interp_presence_absence.ipynb`: documented application showing the interpolation of a single species: _Biddulphia Sinensis_.
* `interp_presence_absence_prod.ipynb`

* `plot_data_histograms.ipynb`
* `plot_validationscore_L.ipynb`


The environment variable `DATADIR` should contain the path where the results should be written to.

## Citation

Please cite this product as:
*A. Barth, Willem Stolte, C. Troupin & Luuk van der Heijden (2020). Probability maps
for different phytoplankton species in the North Sea.*
