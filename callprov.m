function callprov
% This script uses the output from an ice sheet model to estimate
% subglacial erosion and predict the epsilon Nd values of debris at the ice
% sheet margin.
% This function contains the variables and options.
% The following functions are required:
% prov.m             - this reads in files and variables and calls the 
%                      other functions.
% terrestrial_part.m - this performs the terrestrial (i.e. sub-ice) parts
%                      of the calculations.
% surf_currents.m    - this estimates detritus transport in surface
%                      currents as an IRD proxy.
% bot_currents.m     - this estimates bottom current transport of detritus.
% grav_flows.m       - this uses bathymetry to route submarine gravity
%                      flows.
% modern_analysis.m  - this performs some analyses, makes plots and a file 
%                      storing the best fitting method to observations 
%                      across the domain file, if using a modern ice sheet.
% mappingplotter.m   - called with the above subfunctions. A short function
%                      to plot data nicely using the MATLAB mapping toolbox.
% The script uses functions from Antarctic Mapping Tools: 
% https://www.mathworks.com/matlabcentral/fileexchange/47638-antarctic-mapping-tools
% Greene, C.A. et al. (2017). Antarctic Mapping Tools for Matlab. 
% Computers & Geosciences, 104, 151–57, doi:10.1016/j.cageo.2016.08.003.

% Clear workspace/close figures
close all;
clear;
clc;

runname = 'newmethod'; %Name of this run. Used in figure file names.

%Select directory and add relevant paths.
PC = 0;% (0=personal, 1=BP, other=Imp.)

if PC == 0
    cd 'C:\Users\Jim\OneDrive - Imperial College London\MATLAB\Modelling';
    W = 2; %Workers in parfor loop
elseif PC == 1
    cd '/user/home/fe20161/MATLAB/provenance/all_prov/';
    addpath '/user/work/fe20161/Data_For_Matlab/';  %Contains input data files
    addpath '/user/home/fe20161/MATLAB/';
    addpath '/user/home/fe20161/MATLAB/provenance';
    addpath '/user/home/fe20161/MATLAB/provenance/all_prov';
    addpath '/user/home/fe20161/MATLAB/Ant_Mapping_Tools';
    W = 16; %Workers in parfor loop
else
    cd 'C:\Users\jwm17\OneDrive - Imperial College London\MATLAB\Modelling';
    addpath 'C:\Users\jwm17\OneDrive - Imperial College London\MATLAB';
    addpath 'C:\Users\jwm17\OneDrive - Imperial College London\MATLAB\Ant_Mapping_Tools';
    W = 4; %Workers in parfor loop
end

%% Input Parameters and Options
tuning = 1; %Option to 
shelfmelt = 1; %Option to set ocean velocities to ice flow velocities to simulate sub ice shelf melt.
biasWAIS = 0; %Option when collapsed WAIS to bias quarryrate here.
undershelf = 1; %Option to track provenance under ice shelves.
remshelf = 0; %Option to remove ice shelves from the provenacne estimate.
invdistweight = 0; %Option to use inverse distance weighting interpolation from ice sheet margin instead of full method.
plot_figs = 1; %Option to plot and save figures.
nc1 = 'Modern_control_deconto_fort.92.nc';% 'G21_T2_GRL_upload.nc'; % %Model file name. 
model = 1; %1 = P&D, 2 = PISM. For reading outputs from different ISMs.
palaeo = 0; %Is the model run modern (0) or not (1)?
forcebergs = 0; %Force Wilkes Subglacial Basin flow to mirror modern ice velocities.
method_file = 'SedTransMethCalibration_t1_11mar.mat'; % Method file used - only needed if paleo = 1, otherwise a new one is created.
t = 1; % Timeslice of output used. %for G21_T2_GRL_upload.nc, this goes to 240.
site = [1655 -2282];%  [1572 -2111]; %U1358 [-9 -1564]; %U1521  %Set a core site location.
sed_mod = 2; %Does nc1 contain 'quarryrate' output? 0 = no, use velocity^2; 1 = yes; 2 = no, use taub*velocity
hgrams = 0; % Option to make histogram plot at 1) point with most flowline endpoints, 2) specified site offshore from surface currents
surf_months = 12; % Number of surface current months used (1 to 324).
bot_months = 2; % Number of bottom current months used (1 to 24).
rep = 1; % Specify how many time bottom/gravity methods are iterated. 2 default.
radius3 = 200; % Exclusion zone around active volcanoes (km).
plotmonthly = 0; %Option to split the surface velocity data into 12 months to look at seasonal trends (1 = on, 0 = off).
%% Terrestrial component variables
s_points = 4000; %Select the number of streamline seed locations used.
kval = 5.1E-10; % Value of quarrying coefficient. 0 uses spatially variable from Pollard and DeConto (2019), otherwise can specifiy a uniform value.
%% Surface current method variables
vweighting = 1; % Weighting of ice vs ocean velocities where no ocean velocity data.
sed_dist = 1; % Type of distribution of basal sediment. 1 = linear, 2 = exponential, 3 = constant.
ocean_T = 0;% Mean water temperature (oC)
basalsed = 4; % Thickness of basal sediment layer (m). 8m in Christoffersen et al. (2010), but assume some will be lost in GZ, shelf (if present), so 6 default.
decay = 5; % Decay constant used if exponential decay of debris. 5 default.
plataeu = 1e-5; % Minimum ('englacial') sediment yield, relative to base of ice column. 0.00001 default.
smth = 0; % Option to smooth output, specifies moving mean window size. Sometimes useful for a palaeo ice sheet. 0 to turn off.
radius4 = 40; %km around surface velocity flowlines. 40 default.
%% Bottom Current Variables
erosionv = 0.15; % Ocean velocity to erode sediment (m/s). 0.15 defualt.
c0 = 1; % Starting sediment concentration (relative, so 1 is ok).
kappa = 0.41; % von Kármán constant. 0.41 default.
z0 = 0.0005; % Roughness Length (m). 0.0005 default.
vs = 3.3e-5; % Settling velocity (m/s). 0.000033 defualt.
susph = 15; % Suspended sediment thickness (m).
rho = 1025; % Seawater denisty (kg/m3). 1025 defualt.
t0 = 0.08; % Critical depositional stress (Pa). 0.08 default.
radius2 = 40; %Radius around bottom current flow paths (km). 40 default.
%% Gravity Flow variables
shelfmin = -1200; %Minimum depth of shelf break (m). -1200 default.
slopethresh = 1; %Slope threshold above which gravity flows occur (degrees).
%% Variables used for analysis of results
grouping = 10; %depth interval grouping for result analyses (m). 10 default.
Ea = 135; % East Antarctic cutoff (degrees east).
We = -62; % West Antarctic cutoff (degrees west). Peninsula = ~-60.
radius = 200; % Maximum distance from coretop data in spatial analysis (km). 200 default.
% eNd_filter = 100; % Maximum eNd discrepancy for best method match in analysis.

%rng shuffle
%% Call Main Function
[eNd] = prov(runname,plot_figs,W,nc1,model,palaeo,t,site,sed_mod,hgrams,surf_months,...
    bot_months,rep,s_points,kval,vweighting,sed_dist,ocean_T,basalsed,...
    decay,plataeu,smth,erosionv,c0,kappa,z0,vs,susph,rho,t0,radius2,shelfmin,...
    grouping,Ea,We,radius,method_file,radius3,PC,...
    slopethresh,plotmonthly,radius4,invdistweight,undershelf,remshelf,...
    forcebergs,biasWAIS,shelfmelt,tuning)

% surf_temps = 'global-reanalysis-phy-001-031-grepv2-monthly_1660130057544.nc';
% sst_orig = double(ncread(surf_temps,char('thetao_oras'),[1 1 1 1],[Inf Inf Inf 12]));
% lat = double(ncread(surf_temps,char('latitude')));
% long = double(ncread(surf_temps,char('longitude')));
% [lat,long] = meshgrid(lat,long);
%     x1 = double(ncread(nc1,'x1')*1e3);
%     y1 = double(ncread(nc1,'y1')*1e3);
%     [xs,ys] = meshgrid(x1,y1);
% 
% for m = 1:12
%     %Convert to polar stereo
%     [X,Y] = polarstereo_fwd(lat,long); 
%     size(X)
%     size(Y)
%     size(sst_orig)
%     % now interpolate to ice sheet model grid - uneven spacing, have to use
%     % ScatteredInterp. 
%     sst = sst_orig(:,:,m);
%     Fa = scatteredInterpolant(X(:),Y(:),sst(:)); 
%     sst = Fa(ys,xs); 
%     sst = inpaint_nans(sst,4);
%     sst = sst + ocean_T; %Add temperature offset, if required.
%     meltrate = real(60*60*24*2.08*10^-7.*(1.8+sst).^1.5);% Russell-Head (1980)
%     
%     subplot(3,4,m)
%     imagesc(meltrate)
%     colorbar
%     caxis([-.0001 .0001])
% end

