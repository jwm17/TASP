# TASP
Welcome to Tracing Antarctic Sediment Provenance (TASP)!

**For full guidance on how to set up and use the code, please refer to the User Guide.**

This collection of MATLAB functions takes the output of an Antarctic ice sheet model simulation and produces a 
map of marine geochemical sediment provenance around the continent. Data coverage currently spans West 
Antarctica and adjacent parts of East Antarctica, and TASP is currently configured to estimate neodymium (Nd) 
isotope compositions. To produce the estimate, alongside the ice sheet model NetCDF file, TASP uses an estimated
map of subglacial epsilon Nd, plus ocean reanalysis data and an estimate of sub-ice shelf melt rates. It also 
requires data on seafloor surface sediment epsilon Nd values and coordinates of recently active volcanoes.

Future development aims to adapt the code to use the outputs of different ice sheet models, as well as
different sediment provenance proxies.
