# [DataBio](https://www.databio.eu/en/)
## Change detection for SAR and optical Data
### [Citation](http://www2.imm.dtu.dk/pubdb/views/publication_details.php?id=7131)
### The usage of the Standalone software and individual codes are free, however attribution to the publication is required.

Pirzamanbein B., Nielsen A. A. (2019). Standalone software for for detecting changes in {SAR} and optical images. ESA Conference on Big Data from Space (BiDS'19), DOi: 10.2760/848593.

Download the paper from (http://www2.imm.dtu.dk/pubdb/views/publication_details.php?id=7131)


For downloading software MADChange and WISHATChange based on the source code in the src folder check:

Software: [Download](https://behnaz.pirzamanbin.name/PostDoc/)

The page contains:

1. MADChange
   * MADChange_app_Linux
   * MADChange_app_Windows
   * MADChange_app_Mac
   * MADChange_commandline_Linux
   * MADChange_commandline_Windows
   * MADChange_commandline_Mac
2. WISHARTChange
   * WISHARTChange_app_Linux
   * WISHARTChange_app_Windows
   * WISHARTChange_app_Mac
   * WISHARTChange_commandline_Linux
   * WISHARTChange_commandline_Windows
   * WISHARTChange_commandline_Mac
3. Test Data for WISHARTChange
4. Test Data for MADChange

Each of the folders in 1. and 2. contains
* Packaged_Installer
* Users

`Packaged_Installer` folder contains the `MATLAB Runtime` which you `must` install before using any of the software versions.

or

you can download the full MATLAB Runtime from
[MATLAB Runtime](https://se.mathworks.com/products/compiler/matlab-runtime.html)

In folder `Test Data for WISHARTChange` a time series of SAR data is available for diagonal dual polarization (VV and VH) in two formats: GeoTIFF and ENVI with header files. The data is acquired from Google Earth Engine (GEE), for more detail see the ReadMe file.

`Test Data for MADChange` contains two time points of optical images with 4 bands (B2,B3,B4,B8) from GEE in two formats, GeoTIFF and ENVI with a header file.

For more specific detail, see the `Software_description.pdf`.
