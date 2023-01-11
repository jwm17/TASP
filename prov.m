function [eNd] = prov(runname,plot_figs,W,nc1,model,paleo,t,site,sed_mod,hgrams,...
    surf_months,bot_months,rep,s_points,kval,vweighting,sed_dist,ocean_T,...
    basalsed,decay,plataeu,smth,erosionv,c0,kappa,z0,vs,susph,rho,t0,...
    radius2,shelfmin,grouping,Ea,We,radius,method_file,radius3,PC,...
    slopethresh,radius4,invdistweight,undershelf,remshelf,...
    forcewsb,biasWAIS,shelfmelt,tuning,sed_dist_figs,sourcetrace_figs,...
    radius5,seedmeth,seedscale,q_thresh,surf_simp,mrg,mrg2)
% This function uses an ice sheet model output (nc1) to calculate erosion
% rate over Antarctica. This is used to weight the generation of
% streamlines, which trace the expected epsilon Nd signature of sediments
% from their source to the ice sheet margin. The function then uses three
% methods to estimate the transport of sediment offshore and compares this
% to measurements from seafloor surface sediments. Variables are set in
% callprov.m.

%% Make figure filenames.
fig_no = 1; % Starting figure number
calib_file = ['SedTransMethCalibration_t' num2str(t) '_' runname '.mat'];
fig1 = ['Erosion_t' num2str(t) '_' runname '.pdf'];
fig2 = ['Ice_Streamlines_t' num2str(t) '_' runname '.png'];
fig2b = ['IDW_t' num2str(t) '_' runname '.pdf'];
fig4a = ['Surface_Interp_t' num2str(t) '_' runname '.pdf'];
fig4b = ['Sed_Delivery_Heatmap_t' num2str(t) '_' runname '.png'];
fig4c = ['Iceberg_Heatmap_t' num2str(t) '_' runname '.png'];
fig6 = ['Bottom_Interp_t' num2str(t) '_' runname '.pdf'];
fig7 = ['Grav_Flow_Paths_t' num2str(t) '_' runname '.pdf'];
fig8 = ['Grav_Flow_Interp_t' num2str(t) '_' runname '.pdf'];
fig9 = ['Coretop_Data_t' num2str(t) '_' runname '.pdf'];
fig10 = ['Calibration_factor_t' num2str(t) '_' runname '.pdf'];
fig11a = ['Best_eNd_match_t' num2str(t) '_' runname '.pdf'];
fig11b = ['Best_eNd_masked_t' num2str(t) '_' runname '.pdf'];
fig12 = ['Best_Method_t' num2str(t) '_' runname '.pdf'];
fig13 = ['Depth_method_dist_t' num2str(t) '_' runname '.pdf'];
fig13b = ['Waterdepth_method_t' num2str(t) '_' runname '.pdf'];
fig14 = ['Scatter_t' num2str(t) '_' runname '.pdf'];
fig15 = ['eNdmatch_histograms' num2str(t) '_' runname '.pdf'];
figx = ['sourcetrace_map' num2str(t) '_' runname '.png'];
figy = ['sourcetrace_histogram' num2str(t) '_' runname '.pdf'];

%% Define other filenames
plio_uvel = 'UVELcomb.mat';
plio_vvel = 'VVELcomb.mat';
eNd_map = 'merge7_krg.tif';
surf_sed_data = 'Coretop_Data_9.xls';
surf_interp = 'surf_krg2.tif';
volcanoes = 'Active_Volcanos.xlsx';
%Quarter degree, ocean surface, ORAS5 reanalysis. January 1993 to December 2019.
filename = 'global-reanalysis-phy-001-031-grepv2-monthly_1632926904430.nc';
%Quarter degree, all depths, ORAS5 reanalysis. 2019 and 2018 (all months).
filename2 = 'global-reanalysis-phy-001-031-grepv2-monthly_1629111297849.nc';
filename3 = 'bot2018_global-reanalysis-phy-001-031-grepv2-monthly_1636538864227.nc';
if model == 1
    modern_ism = 'Modern_control_deconto_fort.92.nc'; 
elseif model == 2
    modern_ism = 'G21_T2_GRL_upload.nc';
end

%Threshold in ice coverage (1E6 km^2) below which Pliocene ocean velocities are used.
extent_thresh = 12; 

%% Read in parameters from the model output NetCDF file
if model == 1
    ubot = ncread(nc1,'ubot',[1 1 t],[Inf Inf 1])'; % Get basal velocities at selected timestep
    vbot = ncread(nc1,'vbot',[1 1 t],[Inf Inf 1])';
    if paleo == 0
        % Remove unrealistic modern PIG flow directions.
        ubot(246:249,114:121) = 0;
        vbot(246:249,114:121) = 0;
    end
    maskwater =  ncread(nc1,'maskwater',[1 1 t],[Inf Inf 1])'; % Water mask
    x0_orig = ncread(nc1,char('x0')); % Polar stereographic coords from input file. 
    y0_orig = ncread(nc1,char('y0'));
    x1 = double(ncread(nc1,'x1')*1e3);
    y1 = double(ncread(nc1,'y1')*1e3); 
    alond = double(ncread(nc1,char('alond'))); % Longitude (degrees)
    alatd = double(ncread(nc1,char('alatd'))); % Latitude (degrees)
    h = double(ncread(nc1,'h',[1 1 t],[Inf Inf 1])); % Ice thickness data
    hb = ncread(nc1,char('hb'),[1 1 t],[Inf Inf 1]); % Bed elevation data
elseif model == 2
    vel = ncread(nc1,'velbase_mag',[1 1 t],[Inf Inf 1])';
    [ubot,vbot] = gradient(vel);
    ubot(isnan(ubot)) = 0;
    vbot(isnan(vbot)) = 0;
    %As approximating trajectories, need to scale velocities so match
    %magnitude.
    grad_mag = sqrt(ubot.^2+vbot.^2);
    ubot = ubot.*(vel./grad_mag);
    vbot = vbot.*(vel./grad_mag);    
%     ubot = ubot';
%     vbot = vbot';
    maskwater =  ncread(nc1,'mask',[1 1 t],[Inf Inf 1]); % Water mask
    maskwater(maskwater<=2) = 0; % Convert to 0=grounded 1=ungrounded
    maskwater(maskwater>2) = 1;
    maskwater = double(maskwater');

    x0_orig = ncread(nc1,char('x')); % Polar stereographic coords from input file. 
    y0_orig = ncread(nc1,char('y'));
    [x0_orig2,y0_orig2] = meshgrid(x0_orig,y0_orig);
    [alatd,alond] = polarstereo_inv(x0_orig2,y0_orig2);
    alond = alond';
    alatd = alatd';
    dx = -(x0_orig(1)-x0_orig(2))/2; % Half resolution
    x1 = x0_orig + dx; %Convert to centre points
    y1 = y0_orig + dx;
    x0_orig = (x0_orig/1000); %Convert to km
    y0_orig = (y0_orig/1000);
    %x is in metres, 20000, 0 ,20000 at centre.
    %x1 is metres, 5000 -5000 at centre. x0 is km, 10 0 -10 at centre (lhs
    h = double(ncread(nc1,'thk',[1 1 t],[Inf Inf 1])); % Ice thickness data
    hb = ncread(nc1,'topg',[1 1 t],[Inf Inf 1]); % Bed elevation data 
else
    disp('Error - Enter valid model type!')
    return
end

%Dummies for parfor loop
ubot2 = 0;
vbot2 = 0;

%Option to remove ice shelves from estimate. Sets ice thickness and velocities to 0.
if remshelf == 1 || shelfmelt == 1
    if shelfmelt == 1
        ubot2 = ubot;
        vbot2 = vbot;
    end
    ubot(maskwater == 1) = 0;
    vbot(maskwater == 1) = 0;
    if remshelf == 1
        h(imrotate(flipud(maskwater),270) == 1) = 0;
    end
end

[xs,ys] = meshgrid(x1,y1); %size of model domain
x0 = x0_orig;
y0 = y0_orig;
indi = h>0; % Indices with ice

%% Some other bits and pieces
[sy,sx] = size(vbot);
[x22,y22] = meshgrid(1:sx,1:sy); 

if model == 1
    vmag = sqrt(ubot.^2+vbot.^2); %ice basal velocity magnitude
else
    vmag = vel;
end

%Round ISM polar coordinates to nearest 10 km.
x0 = 10*round(x0./10); 
y0 = 10*round(y0./10);
% Find ISM resolution
resolution = round((abs(min(x0))+abs(max(x0)))./size(x0,1));
resolution = 10*round(resolution./10);

%Calculate ice extent
extent = (resolution^2*sum(indi,'all'))/1000000 %Million km2
%Read in velocity data if small ice sheet.
if extent <= extent_thresh
    data = load(plio_uvel);
    UVEL = data.UVEL;
    data = load(plio_vvel);
    VVEL = data.VVEL;
else %Need to make dummies to feed into bot_currents.m
    UVEL = 0;
    VVEL = 0;
end

%Create red-blue divergent colour map.
divmap = makeColorMap([1 0 0],[1 1 1],[0 0 1], 256);

%% Read in other external data
[pv,geodata] = readgeoraster(eNd_map); %Read input eNd map (pv).
% Read spreadsheet with core top locations
coresites = xlsread(surf_sed_data);
% Read in core top interpolation
[coretops,geodata2] = readgeoraster(surf_interp);
% Rotate so it matches everything else 
coretops = imrotate(coretops,270); 
% Read locations of active Antarctic volcanoes
avolcs = xlsread(volcanoes);
%Exclude Webber Nunatak as uncertain volcanic activity and diluted by high
%PIG sed flux.
avolcs = [avolcs(1:4,:); avolcs(6:end,:)];

% Get ocean velocity coordinate data.
lat = double(ncread(filename,char('latitude')));
long = double(ncread(filename,char('longitude')));
[lat,long] = meshgrid(lat,long);

%% Conversion of eNd input to model grid
xcoords = (-3332.75:10:3336.25)'; % Manually generated. Given in 'geodata'.
ycoords = (-3332.75:10:3336.25)';
xcoords = round(xcoords./10).*10;
ycoords = round(ycoords./10).*10;
scaling = resolution/20; %number of points either side to average over

if resolution == 10 % If resolutions match
    lower = find(xcoords==min(x0)); % find points
    upper = find(xcoords==max(x0));
    pv = pv(lower:upper,lower:upper); 
else % If resolutions don't match
    pvsmall = zeros(size(x0,1),size(y0,1)); %create blank matrix
    if model == 2
        pv = pv';
    end
    % Loops to convert 
    for i = 1:size(x0)
        for j = 1:size(y0)
            if x0(i)>min(xcoords) && x0(i)<max(xcoords) && y0(j)>min(ycoords) ...
                    && y0(j)<max(ycoords)
                pvsmall(i,j) = mean(mean(pv(find(xcoords==x0(i))-(scaling-1):find(xcoords==x0(i))+scaling,...
                    find(ycoords==y0(j))-(scaling-1):find(ycoords==y0(j))+scaling)));
            end
        end
    end
    pv = pvsmall;  
end
pv = imrotate(pv,270); %rotate so they match

% Get coordinates of selected site
sitex = site(1); 
sitey = site(2);
sitex = round(sitex./resolution).*resolution;
sitey = round(sitey./resolution).*resolution;
sitex_loc = find(sitex==x0);
sitey_loc = find(sitey==y0);

%% Conversion of core top interpolation to model grid
xcoords2 = (geodata2.XWorldLimits(1)/1000:10:geodata2.XWorldLimits(2)/1000)';
ycoords2 = (geodata2.YWorldLimits(1)/1000:10:geodata2.YWorldLimits(2)/1000)';
xcoords2 = round(xcoords2./10).*10;
ycoords2 = round(ycoords2./10).*10;

coretopssmall = zeros(size(x0,1),size(y0,1)); %create blank matrix
%Loops to convert to core top map to model grid
for i = 1:size(x0)
    for j = 1:size(y0)
        if x0(i)>min(xcoords2) && x0(i)<max(xcoords2) && y0(j)>min(ycoords2) ...
                && y0(j)<max(ycoords2)
            coretopssmall(i,j) = coretops(find(xcoords2==x0(i)),find(ycoords2==y0(j)));
        end
    end
end
coretops = coretopssmall;

%If a large ice sheet, set values beyond modern margin to modern marine sediment
%values.
if paleo == 1
    %Read in modern simulation to get ice extent.
    if model == 1
        mod_h = ncread(modern_ism,'h',[1 1 1],[Inf Inf 1])';
    else 
        mod_h = ncread(modern_ism,'thk',[1 1 1],[Inf Inf 1])';
    end
    % Find areas grounded for palaeo ice sheet and not covered by modern ice.
    newgrounded = zeros(sx,sy);
    newgrounded(maskwater == 0 & mod_h == 0) = 1;
    newgrounded = flipud(imrotate(newgrounded,90));
    pv(newgrounded == 1) = coretops(newgrounded == 1);
end
% Convert core sites and volcano locations to the model grid and create a 
% mask around them.
corelat = coresites(:,1);
corelon = coresites(:,2);
volclat = avolcs(:,1);
volclon = avolcs(:,2);
[coreX,coreY] = polarstereo_fwd(corelon,corelat);
[volcX,volcY] = polarstereo_fwd(volclon,volclat);

for n = 1:numel(coreX)
    corex = coreX(n)/1000;
    corey = coreY(n)/1000;
    corex = round(corex/resolution)*resolution;
    corey = round(corey/resolution)*resolution;
    if corex>min(x0) && corex<max(x0) && corey>min(y0) && corey<max(y0)
        corex_loc(n) = find(corex==x0);
        corey_loc(n) = find(corey==y0);
    end
end

for n = 1:numel(volcX)
    volcx = volcX(n)/1000; 
    volcy = volcY(n)/1000;
    volcx = round(volcx/resolution)*resolution;
    volcy = round(volcy/resolution)*resolution;
    if volcx>min(x0) && volcx<max(x0) && volcy>min(y0) && volcy<max(y0)
        volcx_loc(n) = find(volcx==x0);
        volcy_loc(n) = find(volcy==y0);
    end
end

%Select only points within 'radius' km of data
dist = zeros(size(corex_loc)); %vector containing distance to each coretop
volcs = zeros(size(volcx_loc)); %vector containing distance to each coretop
nrst_smpl = zeros(sx,sy); %Variable containing filter for points within 'radius' km
volc_msk = ones(sx,sy); %Variable containing volcano mask
for i = 1:sx
    for j = 1:sy
        for k = 1:size(corex_loc,2)
            if corex_loc(k)~=0 && corey_loc(k)~=0
                dist(k) = resolution*sqrt((xcoords(corex_loc(k))-xcoords(i)).^2+...
                    (ycoords(corey_loc(k))-ycoords(j)).^2);
            else
                dist(k) = 10^6; %set to arbitrary high value so definitely > radius
            end
        end
        if min(dist)<(resolution*radius)
            nrst_smpl(i,j) = 1; %set to 1 if close to a core top.
        end
        %Loop for volcano mask
        for k = 1:size(volcx_loc,2)
            if volcx_loc(k)~=0 && volcy_loc(k)~=0
                volcs(k) = resolution*sqrt((xcoords(volcx_loc(k))-xcoords(i)).^2+...
                    (ycoords(volcy_loc(k))-ycoords(j)).^2);
            else
                volcs(k) = 10^6; %set to arbitrary high value so definitely > radius3
            end
        end
        if min(volcs)<(resolution*radius3)
            volc_msk(i,j) = 0; %set to 1 if close to a volcano.
        end
    end
end
volc_msk(volc_msk==0) = NaN;

%Find ice shelves.
shelf = abs(indi-1)+flipud(imrotate(maskwater,90));
shelf(shelf==2) = 0;

%Call function to estimate provenance value under ice.
[average,fig_no,x2,y2,shelf_average,packages,quarry2] = terrestrial_part(sed_mod,t,vmag,...
    maskwater,nc1,sx,sy,kval,alond,alatd,biasWAIS,s_points,fig1,x22,y22,...
    ubot,vbot,pv,indi,undershelf,resolution,x0,y0,fig_no,plot_figs,...
    paleo,PC,fig2,hgrams,model,shelf,seedmeth,seedscale,q_thresh,shelfmelt);

%If just using a 'nearest-neighbour' type approach offshore...
if invdistweight == 1
    [plotx, ploty] = find(~isnan(average)); %Get x and y coordinates of endpoints
    e = 2; % Distance weight. 2 gives optimal results for modern ice sheet.
    ng = 10; % Number of neighbours used. 10 gives optimal results for modern ice sheet.
    invdisw = IDW(plotx,ploty,average(~isnan(average)),1:sx,1:sy,e,'ng',ng);
    invdisw = flipud(invdisw);
    plot_invdisw = invdisw;
    plot_invdisw(imrotate(indi,90)) = NaN; %Set ice to NaN
    plot_invdisw = imrotate(plot_invdisw,270); 
    
    subplot(2,1,1)
    mappingplotter(alatd,alond,plot_invdisw,-20,5,viridis);
    colormap(viridis)
    hold on
    freezeColors
    contourm(alatd,alond,double(flipud(imrotate(maskwater,90))),'LineColor','k')
    
    nrst_smpl_nans = nrst_smpl;
    nrst_smpl_nans(nrst_smpl_nans == 0) = NaN;
    plot_invdisw_dif = coretops - plot_invdisw;
    plot_invdisw_dif = plot_invdisw_dif.*nrst_smpl_nans;
    
    subplot(2,1,2)
    mappingplotter(alatd,alond,plot_invdisw_dif,-5,5,divmap);
    colormap(divmap)
    hold on
    freezeColors
    contourm(alatd,alond,double(flipud(imrotate(maskwater,90))),'LineColor','k')
   
    print(gcf,fig2b,'-dpdf','-bestfit');
else
%Otherwise, use a more realistic movement of marine detritus.
    %Call surface current estimate function.
    [oceantracks,fig_no,nrst_smpl_nans] = surf_currents(plot_figs,...
        fig_no,surf_months,sed_dist,sx,sy,PC,W,resolution,indi,extent_thresh,...
        lat,long,forcewsb,shelfmelt,ubot,vbot,ubot2,vbot2,maskwater,x22,y22,...
        x2,y2,basalsed,decay,plataeu,alatd,alond,runname,coretops,nrst_smpl,...
        ocean_T,filename,ys,xs,shelf,vweighting,average,sitex_loc,sitey_loc,...
        figx,figy,s_points,fig4b,fig4c,xcoords,ycoords,radius4,smth,...
        undershelf,divmap,fig4a,model,shelf_average,t,tuning,x0,y0,...
        sed_dist_figs,sourcetrace_figs,seedmeth,packages,surf_simp,...
        quarry2,paleo,nc1,mrg,mrg2);
    
    iter = 0; %Needs a dummy value even if not used.
    for repeat = 1:rep  % Iterate 'rep' times, using gravity flow output as bottom current input.
        %Call bottom current estimate function.
            [oceantracks2,fig_no] = bot_currents(sx,sy,bot_months,plot_figs,filename2,...
        filename3,ys,xs,extent,extent_thresh,erosionv,alatd,alond,kappa,hb...
        ,z0,rho,t0,c0,vs,susph,oceantracks,repeat,iter,runname,rep,xcoords,...
        ycoords,radius2,resolution,vbot,ubot,smth,nrst_smpl_nans,coretops,...
        divmap,maskwater,fig6,indi,x22,y22,model,UVEL,VVEL);
       %Call gravity flow estimate function.
            [oceantracks3,fig_no,iter] = grav_flows(x0_orig,y0_orig,sx,sy,hb,slopethresh,...
        maskwater,coretops,oceantracks,oceantracks2,shelfmin,nrst_smpl_nans,alatd,...
        alond,divmap,fig8,iter,repeat,rep,plot_figs,resolution,paleo,fig_no,fig7,...
        vbot,ubot,model,xcoords,ycoords,radius5);
    end
    
    %Switch site x and y values.
    sitex_loc_new = sy - sitey_loc;
    sitey_loc_new = sitex_loc;

end %IDW option

%% Analysis of Results  
if paleo == 0
    % If using the modern ice sheet
    if invdistweight == 0
        [fig_no,eNd,~] = modern_analysis(plot_figs,coretops,oceantracks,oceantracks2,...
            oceantracks3,sitex_loc,sitey_loc,x22,y22,maskwater,...
            fig_no,y0,fig9,fig10,fig11a,fig11b,fig12,fig13,fig13b,fig14,fig15,sx,sy,...
            x0,grouping,hb,Ea,We,nrst_smpl,coresites,calib_file,corex_loc,...
            corey_loc,alatd,alond,volc_msk,runname,PC,corelat,corelon,...
            tuning);
    else
         [fig_no,eNd,~] = modern_analysis(plot_figs,coretops,imrotate(invdisw,270),imrotate(invdisw,270),...
            invdisw,sitex_loc,sitey_loc,x22,y22,maskwater,...
            fig_no,y0,fig9,fig10,fig11a,fig11b,fig12,fig13,fig13b,fig14,fig15,sx,sy,...
            x0,grouping,hb,Ea,We,nrst_smpl,coresites,calib_file,corex_loc,...
            corey_loc,alatd,alond,volc_msk,runname,PC,corelat,corelon,...
            tuning);
    end
elseif paleo == 1 && invdistweight == 0
    % If using a non-modern ice sheet configuration
    % Load existing method calibration file
    load(method_file)

    if resolution == 10
        lower = find(x0_out==min(x0)); %if resolutions match, find points
        upper = find(x0_out==max(x0));
        factr = factr(lower:upper,lower:upper); 
    else %if resolutions don't match
        fctrsmall = zeros(size(x0,1),size(y0,1)); %create blank matrix
        %loops to convert 
        for i = 1:size(x0)
            for j = 1:size(y0)
                if x0(i)>min(x0_out) && x0(i)<max(x0_out) && y0(j)>min(y0_out) ...
                        && y0(j)<max(y0_out)
                    fctrsmall(i,j) = mean(mean(factr(find(x0_out==x0(i))-(scaling-1):find(x0_out==x0(i))+scaling,...
                        find(y0_out==y0(j))-(scaling-1):find(y0_out==y0(j))+scaling)));
                end
            end
        end
        factr = fctrsmall;  
    end
    
    if plot_figs == 1
        figure(fig_no);
        fig_no = fig_no + 1;
        factr = imrotate(factr,90);
        imagesc(factr)
        colorbar
    end
    
    %Merge output into single array and save it as .mat file.
    oceantracks(:,:,1) = imrotate(oceantracks,90);
    oceantracks(:,:,2) = imrotate(oceantracks2,90);
    oceantracks(:,:,3) = oceantracks3;
    
    %Get site of interest coordinates
    sitey_loc = sy-sitey_loc;
    oceantracks(sitey_loc,sitex_loc,:)

    wghtd = zeros(sx,sy); %weighted eNd output.
    for i = 1:sx
        for j = 1:sy
            if isnan(factr(i,j)) || factr(i,j)==0 % If inland 
                wghtd(i,j) = mean(oceantracks(i,j,~isnan(oceantracks(i,j,:))));
            else
                if round(factr(i,j)) == factr(i,j) % If an exact method
                   wghtd(i,j) = oceantracks(i,j,factr(i,j));
                elseif factr(i,j)>1 && factr(i,j)<2
                    lwr = oceantracks(i,j,1);
                    upr = oceantracks(i,j,2);
                    wghtd(i,j) = lwr-(lwr-upr)*(factr(i,j)-1);
                elseif factr(i,j)>2 && factr(i,j)<3
                    lwr = oceantracks(i,j,2);
                    upr = oceantracks(i,j,3);
                    wghtd(i,j) = lwr-(lwr-upr)*(factr(i,j)-2);  
                else
                    lwr = oceantracks(i,j,3);
                    upr = oceantracks(i,j,1);
                    wghtd(i,j) = lwr-(lwr-upr)*factr(i,j); 
                end
            end
        end
    end
    % Set NaNs to mean of other methods
    means = mean(oceantracks,3,'omitnan');
    wghtd(isnan(wghtd)) = means(isnan(wghtd));

    if smth > 0 %Option to apply smoothing to final result.
        wghtd = movmean_2d(wghtd,smth)
    end

    if plot_figs == 1
        figure(fig_no);
        fig_no = fig_no + 1;
        subplot(2,1,1)
        mappingplotter(alatd,alond,imrotate(wghtd,270),-20,5,viridis);
        colormap(viridis)
        hold on
        freezeColors
        contourm(alatd,alond,double(flipud(imrotate(maskwater,90))),'LineColor','k')
        
        subplot(2,1,2) %Plot difference from modern.
        mappingplotter(alatd,alond,nrst_smpl_nans.*(imrotate(wghtd,270)-coretops),-5,5,divmap);
        colormap(divmap)
        hold on
        freezeColors
        contourm(alatd,alond,double(flipud(imrotate(maskwater,90))),'LineColor','k') 
        
        if PC==1
            print(gcf,fig11a,'-dpdf','-fillpage')
        else
            set(gcf,'Position',[50 50 600 1200]);
            exportgraphics(gcf,fig11a,'Resolution',300);
        end
        close 
    end

    oceantracks(:,:,4) = wghtd; %Set best guess eNd.
    save([runname '_' num2str(t) '_output.mat'],'oceantracks'); %Save data.

    eNd = wghtd(sitey_loc,sitex_loc);
elseif paleo == 1 && invdistweight == 1
    oceantracks = plot_invdisw;
    save([runname '_' num2str(t) '_output_IDW.mat'],'oceantracks')
else
    disp('Enter valid palaeo option: 1 = Palaeo, 0 = Modern')
end