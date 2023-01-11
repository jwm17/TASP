function callprov
% This function uses the output from an ice sheet model to estimate
% subglacial erosion and predict the epsilon Nd values of debris at the ice
% sheet margin. It contains the user-defined variables and options.
% The following functions are required:
% prov.m             - this reads in files and variables and calls the 
%                      other subfunctions.
% terrestrial_part.m - this performs the terrestrial (i.e. sub-ice)
%                      provenance tracing.
% surf_currents.m    - this estimates detritus transport in surface
%                      currents as an IRD proxy.
% bot_currents.m     - this estimates bottom current transport of detritus.
% grav_flows.m       - this uses bathymetry to route submarine gravity
%                      flows.
% modern_analysis.m  - this performs some analyses, makes plots and makes a  
%                      file storing the best fitting method to observations 
%                      across the domain, if using a modern ice sheet.
% mappingplotter.m   - called within the above subfunctions. A short function
%                      to plot maps using the MATLAB mapping toolbox.
% The script uses functions from Antarctic Mapping Tools: 
% https://www.mathworks.com/matlabcentral/fileexchange/47638-antarctic-mapping-tools
% Greene, C.A. et al. (2017). Antarctic Mapping Tools for Matlab. 
% Computers & Geosciences, 104, 151–57, doi:10.1016/j.cageo.2016.08.003.

% Author: Jim Marschalek, Imperial College London (November 2022).
% Email: jwm17@ic.ac.uk.

clc;
close all;
clear all;

runname = 'test'; %Name of this run. Used in figure file names.

%Select directory and add relevant paths. Can adjust as necessary.
PC = 0;% (0=personal, 1=HPC, other=Imp.)

if PC == 0
    cd 'C:\Users\Jim\OneDrive - Imperial College London\MATLAB\Modelling';
    addpath 'C:\Users\Jim\OneDrive - Imperial College London\MATLAB';
    W = 4; %Workers in parfor loop
elseif PC == 1
    cd '/rds/general/user/jwm17/home/MATLAB/provenance/modern/21oct22';
    addpath '/rds/general/user/jwm17/home/MATLAB/provenance/scripts';
    addpath '/rds/general/user/jwm17/home/MATLAB/provenance/data_files';
    addpath '/rds/general/user/jwm17/home/MATLAB/provenance/scripts/Ant_Mapping_Tools';
    W = 16; %Workers in parfor loop

else
    cd 'C:\Users\jwm17\OneDrive - Imperial College London\MATLAB\Modelling';
    addpath 'C:\Users\jwm17\OneDrive - Imperial College London\MATLAB';
    addpath 'C:\Users\jwm17\OneDrive - Imperial College London\MATLAB\Ant_Mapping_Tools';
    W = 4; %Workers in parfor loop
end

%% Input Parameters and Options
tuning = 0; %Option to reduce output when tuning parameters (1 = on, 0 = off).
shelfmelt = 1; %Option to set ocean velocities to ice flow velocities to simulate sub ice shelf melt (1 = on, 0 = off).
biasWAIS = 0; %Option when collapsed WAIS to bias quarryrate around MBL (1 = on, 0 = off).
undershelf = 0; %Option to track provenance under ice shelves (1 = on, 0 = off). Only used if shelfmelt = 0.
remshelf = 0; %Option to remove ice shelves from the provenance estimate (1 = on, 0 = off).
invdistweight = 0; %Option to use inverse distance weighting interpolation from ice sheet margin instead of full method (1 = on, 0 = off).
plot_figs = 1; %Option to plot and save figures (1 = on, 0 = off).
sed_dist_figs = 0; %option to plot figures showing sediment distribution (1 = on, 0 = off).
sourcetrace_figs = 0; %Option to plot figures tracing the source of debris near a given site (1 = on, 0 = off).
nc1 = 'Modern_control_deconto_fort.92.nc';%'DP16_fort.92.nc';% 'G21_T2_GRL_upload.nc'; % Ice sheet model file name. 
model = 1; %1 = P&D, 2 = PISM. For reading outputs from different ISMs.
palaeo = 0; %Is the ISM run modern (0) or not (1)?
forcebergs = 1; %Force Wilkes Subglacial Basin ocean to mirror modern ice velocities (1 = on, 0 = off).
method_file = 'SedTransMethCalibration_t1_11mar.mat'; % Method file used - only needed if palaeo = 1, otherwise a new one is created.
t = 1; % Timeslice of output used. %for G21_T2_GRL_upload.nc, this goes to 240.
site = [1655 -2282];%  [1572 -2111]; %U1358 [-9 -1564]; %U1521  %Set a core site location. In polar stereographic coordinates.
sed_mod = 2; %Predicts erosion rate depending on input type. Can: can use squared ice velocity (=0); read 'quarryrate' directly from ISM (=1); or use product of basal shear and ice velocity (=2).
hgrams = 0; % Option to make histogram plot at: point with most flowline endpoints (=1), specified site offshore from surface currents (=2). Only possible if seedmeth = 0.
surf_months = 2; % Number of surface current months used. 324 recommended.
bot_months = 1; % Number of bottom current months used. 24 recommended.
rep = 1; % Specify how many time bottom current/gravity flow methods are iterated. 2 recommended.
radius3 = 200; % Radius of exclusion zone around active volcanoes (km) for plotting and statistics.
%% Terrestrial component variables
seedmeth = 1; %Method use to seed flowlines. Either randomly generated (=0), or using an erosion rate threshold and weighted (=1).
s_points = 500; %Select the number of streamline seed locations used. Only used if seedmeth = 0. 32000 recommended.
seedscale = 1; %Can increase the seed locations for the terrestrial part by this proportion. Note that only the number specified in 's_points' are read into the surface current estimate. Only used if seedmeth = 0.
q_thresh = 3; %Erosion rate threshold for seed locations (mm/kyr). Only used if seedmeth = 1.
kval = 5.1E-10; % Value of quarrying coefficient. 0 uses spatially variable from Pollard and DeConto (2019), otherwise can specifiy a uniform value.
%% Surface current method variables
vweighting = 1; % Weighting of ice vs ocean velocities where no ocean velocity data. 1 recommended.
sed_dist = 1; % Type of distribution of basal sediment. 1 = linear, 2 = exponential, 3 = constant.
ocean_T = 0; % Uniform water temperature offset from the modern ocean reanalysis (degrees C).
basalsed = 4; % Thickness of basal sediment layer (m). 4 default.
decay = 5; % Decay constant used if exponential decay of debris. 5 default.
plataeu = 1e-5; % Minimum ('englacial') sediment yield, relative to base of ice column. 1e-5 default.
smth = 0; % Option to smooth output, specifies moving mean window size. Sometimes useful for a palaeo ice sheet. 0 to turn off.
radius4 = 40; %interpolation around surface velocity flowlines (km). 40 default.
surf_simp = 1; %Option to simplify surface current estimation by combining repeated streamlines.
mrg = 3; %Number of cells to blend modern and palaeo ocean currents over. 3 recommended.
mrg2 = 2; % Number of cells to blend ocean and ice velocities over. 2 recommended.
%% Bottom Current Variables
erosionv = 0.15; % Ocean velocity required to mobilise sediment (m/s). 0.15 default.
c0 = 1; % Starting sediment concentration (relative, so 1 is ok).
kappa = 0.41; % von Kármán constant. 0.41 default.
z0 = 0.0005; % Roughness Length (m). 0.0005 default.
vs = 3.3e-5; % Settling velocity (m/s). 0.000033 default.
susph = 15; % Suspended sediment thickness (m).
rho = 1025; % Seawater density (kg/m3). 1025 default.
t0 = 0.08; % Critical depositional stress (Pa). 0.08 default.
radius2 = 40; %Radius around bottom current flow paths (km). 40 default.
%% Gravity Flow variables
shelfmin = -1200; %Minimum depth of shelf break (m). -1200 default.
slopethresh = 1; %Slope threshold above which gravity flows occur (degrees).
radius5 = 300; %Distance of interpolation distance around gravity flows (km).
%% Variables used for analysis of results
grouping = 10; %depth interval grouping for result analyses (m). 10 default.
Ea = 135; % East Antarctic cutoff (degrees east).
We = -62; % West Antarctic cutoff (degrees west). Peninsula = ~-60.
radius = 200; % Maximum distance from coretop data in spatial analysis (km). 200 default.

%First check these inputs and options are in correct format.
parse_input(runname,plot_figs,W,nc1,model,palaeo,t,site,sed_mod, ...
    hgrams,surf_months,bot_months,rep,s_points,kval,vweighting,sed_dist,ocean_T,...
    basalsed,decay,plataeu,smth,erosionv,c0,kappa,z0,vs,susph,rho,t0,...
    radius2,shelfmin,grouping,Ea,We,radius,method_file,radius3,...
    slopethresh,radius4,invdistweight,undershelf,remshelf,forcebergs,...
    biasWAIS,shelfmelt,tuning,sed_dist_figs,sourcetrace_figs,radius5,...
    seedmeth,seedscale,q_thresh,surf_simp,mrg,mrg2)

%rng shuffle
%% Call Main Function
[eNd] = prov(runname,plot_figs,W,nc1,model,palaeo,t,site,sed_mod,hgrams,...
    surf_months,bot_months,rep,s_points,kval,vweighting,sed_dist,ocean_T,...
    basalsed,decay,plataeu,smth,erosionv,c0,kappa,z0,vs,susph,rho,t0,...
    radius2,shelfmin,grouping,Ea,We,radius,method_file,radius3,PC,...
    slopethresh,radius4,invdistweight,undershelf,remshelf,forcebergs,...
    biasWAIS,shelfmelt,tuning,sed_dist_figs,sourcetrace_figs,radius5,...
    seedmeth,seedscale,q_thresh,surf_simp,mrg,mrg2)

function parse_input(runname,plot_figs,W,nc1,model,palaeo,t,site,sed_mod, ...
    hgrams,surf_months,bot_months,rep,s_points,kval,vweighting,sed_dist,ocean_T,...
    basalsed,decay,plataeu,smth,erosionv,c0,kappa,z0,vs,susph,rho,t0,...
    radius2,shelfmin,grouping,Ea,We,radius,method_file,radius3,...
    slopethresh,radius4,invdistweight,undershelf,remshelf,forcebergs,...
    biasWAIS,shelfmelt,tuning,sed_dist_figs,sourcetrace_figs,radius5,...
    seedmeth,seedscale,q_thresh,surf_simp,mrg,mrg2)
%This subfunction checks all the inputs are suitable before executing the 
% main function.

assert(ischar(runname)==1,'Runname must be a string!')

assert(plot_figs==0 || plot_figs==1,'Plot_figs must be 0 or 1!')

assert(isnumeric(W)==1,'Number of parallel workers must be numeric!')

assert(ischar(nc1)==1,'ISM input must not be numeric!')
nc1_str = string(nc1);
assert(contains(nc1_str,".nc")==1,'ISM input must be a .nc file!')

assert(model==1 || model==2,'Model option must be 1 or 2!')

assert(palaeo==0 || palaeo==1,'Palaeo option must be 0 or 1!')

assert(isnumeric(t)==1,'t must be a number!')
assert(t>0 && round(t)==t,'t must be a positive integer!')

assert(isnumeric(site)==1,'Site location must be numeric!')

assert(size(site,2)==2 && size(site,1)==1,'Site location must be a 1 by 2 vector!')

assert(sed_mod==0 || sed_mod==1 || sed_mod==2,'Sed_mod option must be 0, 1 or 2!')

if seedmeth == 0
assert(hgrams==0 || hgrams == 1,'hgrams must be 0 or 1!')
end

assert(isnumeric(surf_months)==1 && round(surf_months)==surf_months && ...
    surf_months>0, 'surf_months must be a numeric, positive integer!')

assert(isnumeric(bot_months)==1 && round(bot_months)==bot_months && ...
    bot_months>0, 'bot_months must be a numeric, positive integer!')

assert(isnumeric(rep)==1 && round(rep)==rep && rep>0,...
    'rep must be a numeric, positive integer!')

if seedmeth == 0
assert(isnumeric(s_points)==1 && round(s_points)==s_points && ...
    s_points>0, 's_points must be a numeric, positive integer!')
end

assert(isnumeric(kval)==1 && kval>=0, 'kval must be numeric and positive!')

assert(isnumeric(vweighting)==1 && vweighting>=0, 'vweighting must be numeric and positive!')

assert(sed_dist==1 || sed_dist==2 || sed_dist==3,'Sed_dist option must be 1, 2 or 3!')

assert(isnumeric(ocean_T)==1,'ocean_T must be numeric!')

assert(isnumeric(basalsed)==1 && basalsed>0, 'basalsed must be numeric and positive!')

assert(isnumeric(decay),'decay must be numeric!')

assert(isnumeric(plataeu)==1 && plataeu>=0, 'plataeu must be numeric and positive!')

assert(isnumeric(smth)==1 && round(smth)==smth && ...
   smth>=0, 'Smoothing must be a numeric, positive integer!')

assert(isnumeric(erosionv)==1 && erosionv>0, 'erosionv must be numeric and positive!')

assert(isnumeric(c0)==1 && c0>=0, 'c0 must be numeric and positive!')

assert(isnumeric(kappa)==1 && kappa>=0, 'kappa must be numeric and positive!')

assert(isnumeric(z0)==1 && z0>=0, 'z0 must be numeric and positive!')

assert(isnumeric(vs)==1 && vs>=0, 'vs must be numeric and positive!')

assert(isnumeric(susph)==1 && susph>=0, ...
    'suspended sediment thickness must be numeric and positive!')

assert(isnumeric(rho)==1 && rho>=0,...
    'seawater density must be numeric and positive!')

assert(isnumeric(t0)==1 && t0>=0, ...
    'critical depsositional stress must be numeric and positive!')

assert(isnumeric(radius)==1 && radius>=0, 'radius must be numeric and positive!')
assert(isnumeric(radius2)==1 && radius2>=0, 'radius2 must be numeric and positive!')
assert(isnumeric(radius3)==1 && radius3>=0, 'radius3 must be numeric and positive!')
assert(isnumeric(radius4)==1 && radius4>=0, 'radius4 must be numeric and positive!')
assert(isnumeric(radius5)==1 && radius5>=0, 'radius5 must be numeric and positive!')

assert(isnumeric(shelfmin)==1 && shelfmin<=0, 'shelfmin must be numeric and negative!')

assert(isnumeric(grouping)==1 && round(grouping)==grouping && grouping>0,...
    'grouping must be a numeric, positive integer!')

assert(isnumeric(Ea)==1 && Ea>-180 && Ea<180,'Ea must be a numeric and between -180 and 180!')
assert(isnumeric(We)==1 && We>-180 && We<180,'We must be a numeric and between -180 and 180!')

if palaeo == 1
    assert(ischar(method_file)==1,'Method file used must not be numeric!')
    method_file_str = string(method_file);
    assert(contains(method_file_str,".mat")==1,'Method input must be a .mat file!')
    load(method_file)
    assert(exist('factr','var')==1,'Selected method file does not contain required variables')
    clear factr
    clear x0_out
    clear y0_out
end

assert(isnumeric(slopethresh)==1 && slopethresh>=0, 'slopethresh must be numeric and positive!')

assert(invdistweight==0 || invdistweight==1,'invdistweight option must be 0 or 1!')
assert(undershelf==0 || undershelf==1,'undershelf option must be 0 or 1!')
assert(remshelf==0 || remshelf==1,'remshelf option must be 0 or 1!')
assert(forcebergs==0 || forcebergs==1,'forcebergs option must be 0 or 1!')
assert(biasWAIS==0 || biasWAIS==1,'biasWAIS option must be 0 or 1!')
assert(shelfmelt==0 || shelfmelt==1,'invdistweight option must be 0 or 1!')
assert(tuning==0 || tuning==1,'tuning option must be 0 or 1!')
assert(sed_dist_figs==0 || sed_dist_figs==1,'sed_dist_figs option must be 0 or 1!')
assert(sourcetrace_figs==0 || sourcetrace_figs==1,'sourcetrace_figs option must be 0 or 1!')
assert(surf_simp==0 || surf_simp==1,'surf_simp option must be 0 or 1!')

assert(isnumeric(seedscale)==1 && round(seedscale)==seedscale && seedscale>0,...
    'seedscale must be a numeric, positive integer!')

assert(isnumeric(q_thresh)==1 && q_thresh>=0, 'q_thresh must be numeric and positive!')

assert(isnumeric(mrg)==1 && mrg>=0,...
    'merge interval must be numeric and positive!')
assert(isnumeric(mrg2)==1 && mrg2>=0,...
    'merge interval must be numeric and positive!')

