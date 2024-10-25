**Data available at** [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13988116.svg)](https://doi.org/10.5281/zenodo.13988116)

# NSSM

This repository contains example code and data accompanying the paper: "Nonlinear sound-sheet microscopy: imaging opaque organs at the capillary and cellular scale".

## NSSM Verasonics Imaging Sequence
The file `NSSM_Sequence.m` contains a Verasonics imaging sequence that allows
volumetric nonlinear sound-sheet microscopy (NSSM). 

While running the sequence, two orthogonal sound sheets will be displayed, and
the active displayed sound sheet can be changed using a slider in the VSX GUI. 

After running the sequence, the full volume can be reconstructed using the
script `MSSM_Sequence_Postprocess.m`. This will reconstruct all acquired sound 
sheets, and compound the two arrays. 

## Example Data
Two example datasets are available on Zenodo:
10.5281/zenodo.13988116. Download the `.mat` files
from Zenodo and store them in the `data` folder. 


### Wells with E. Coli 
Representative dataset for paper figure 2D,E. The file
`RFData_planes_figure_2E.mat` contains the RFData and the required parameters to
reconstruct two orthogonal sound sheets, in SSM and NSSM mode. The file
`demo_reconstruct_orthogonal_NSSM_images.m` shows how to reconstruct the sound sheet data. 

The file `beamformed_volume_figure_2E.mat` contains beamformed SSM and NSSM
volumes as displayed in figure 2G. The file `demo_navigate_3D_NSSM_data.m` can
be used to view the volumes. 

### mARG expression in orthopic tumors 
Representative dataset for paper figure 3B. The file
`RFData_planes_figure_3B.mat` contains the RFData and the required parameters to
reconstruct two orthogonal sound sheets, in SSM and NSSM mode. The file
`demo_reconstruct_orthogonal_NSSM_images.m` shows how to reconstruct the sound sheet data. 

The file `beamformed_volume_figure_3B.mat` contains beamformed SSM and NSSM
volumes as displayed in figure 3C. The file `demo_navigate_3D_NSSM_data.m` can
be used to view the volumes. 

## `nssm` code package

The functions to reconstruct the sound sheet data are packaged in the `+nssm`
namespace. There are three inner namespaces:
- `+recon` contains all functions required to reconstruct sound sheet RFData to
  images. 
- `+sequence` contains some helper functions for real time beamforming in the
  NSSM Verasonics imaging sequennce. 
- `+utils` contains some utility functions (e.g. for rearranging RFData).

All functions are documented and help can be viewed using e.g.:
```
help nssm.recon.reconDAS
``` 

## License
[CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
