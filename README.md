# bycatch-data-improves-SDM
This repository is an implementaion of the research work *Bycatch information complements the understanding of spatial distribution for commercially-important fish species*.

The repository contains:
- Source code

### About Raw Data
The fishery data for hairtail that support the findings of this study are available on request from the corresponding author. The data are not publicly available due to privacy or ethical restrictions.
The environment data used in this study is public available.  

1. The bathymetric (Bathy) data was downloaded from the GEBCO website (https://www.gebco.net/data_and_products/gridded_bathymetry_data/). 
2. The chlorophyll-a (Chla) data was downloaded from the NASA website (https://oceandata.sci.gsfc.nasa.gov/opendap/MODISA/L3SMI/contents.html). 
3. The five monthly satellite environmental data, including sea surface temperature (SST), sea surface salinity (SSS), sea bottom temperature (SBT), sea surface height (SSH), and ocean mixed layer thickness (MLD) were downloaded from the CMEMS website (https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description). The data from the CMEMS website are available for users who have registered a free account.

### Repository structure
I first categorized the code into two groups, which is stored as two folders: 
- environment: this contains the code for compiling the environment data
- Integrate_bycatch: this folder contains the code for reproducing the results in the study.

Other Scripts with **"Function"** and **"Plot"** at the beginning are the function and needed in the scripts in the folder "integrate_bycatch" and are for ploting


# Step by Step Analysis
1. Download data from the Dryad or from the authors, and save it with the name "Fishtactics_SDM_data_check". 
2. Download code and (1) save at the same path as "Fishtactics_SDM_data_check" and (2) named it "Fishtactics_SDM_code_check"
3. Remember to change the contents saved in the variables **"wd"** and **"Root_dir"** in the scripts to your saving path. 
4. Run the scripts in the folder, integrate_bycatch, one by one in order (From 1st to 13th scripts).

# Contact
If you find any bugs or have any questions about the implementation, pelease contact us via elthina02017@gmail.com
