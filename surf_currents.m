function [oceantracks,fig_no,nrst_smpl_nans] = surf_currents(plot_figs,...
    fig_no,surf_months,sed_dist,sx,sy,PC,W,resolution,indi,extent_thresh,...
    lat,long,forcewsb,shelfmelt,ubot,vbot,ubot2,vbot2,maskwater,x22,y22,...
    x2,y2,basalsed,decay,plataeu,alatd,alond,runname,coretops,nrst_smpl,...
    ocean_T,filename,ys,xs,shelf,vweighting,average,sitex_loc,sitey_loc,...
    figx,figy,s_points,fig4b,fig4c,xcoords,ycoords,radius4,smth,...
    undershelf,divmap,fig4a,model,shelf_average,t,tuning,x0,y0,sed_dist_figs,...
    sourcetrace_figs,seedmeth,packages,surf_simp,quarry2,palaeo,nc1,mrg,mrg2)

% Estimate of provenance value using ocean surface currents.

disp(['Estimating surface current transport using modern ocean velocities. Using '...
    num2str(surf_months) ' months.'])

%Print a statement regarding selection of basal debris distribution.
if sed_dist == 1
    disp('Using linear decay of iceberg debris content.'); %days for the sediment to be lost
elseif sed_dist == 2
    disp('Using exponential decay of iceberg debris content.')
elseif sed_dist == 3
    disp('Using a constant basal sediment content')
else
    disp('Error - enter valid basal debris distribution option!')
end

oceantracks_out = zeros(sy,sx,surf_months);

%% Read in and process monthly data prior to parfor loop

[X,Y] = polarstereo_fwd(lat,long); %Convert to polar stereo
meltrate = zeros(sx,sy,surf_months); %variable to contain melt rate maps.
Uout_all = meltrate; %Contains ocean surface veolcity U component.
Vout_all = meltrate; %Contains ocean surface veolcity V component.

if shelfmelt == 1 && palaeo == 0
    %%% Read in melt rates under ice shelves
    shelfmelt_file = 'bb0448974g_3_1.h5'; % https://library.ucsd.edu/dc/object/bb0448974g
    shelfm = h5read(shelfmelt_file,'/w_b/'); %Read melt rates
    shelfm = shelfm./365.25; %Convert from m/yr to m/day
    shelfx = h5read(shelfmelt_file,'/x/'); %Read coordinates
    shelfy = h5read(shelfmelt_file,'/y/');
    shelfx = shelfx./1000; %convert m to km
    shelfy = shelfy./1000;
    %Convert to model grid 
    shelfm_model = zeros(sx,sy);
    for i = 1:sx
        for j = 1:sy
            if x0(i)>min(shelfx) && x0(i)<max(shelfx) && y0(j)>min(shelfy) ...
                    && y0(j)<max(shelfy)
                shelfm_model(i,j) = shelfm(find(shelfx==x0(i)),find(shelfy==y0(j)));
            end
        end
    end
    shelfm = shelfm_model;
elseif shelfmelt == 1 && palaeo == 1
    %Read shelf melt from NetCDF file
    shelfm = ncread(nc1,'oceanmelt',[1 1 t],[Inf Inf 1]); % Ocean melt rate data
    shelfm = shelfm./365.25; %Convert from m/yr to m/day
end

%Calculate extent.
extent = (resolution^2*sum(indi,'all'))/1000000; %Million km2

if extent <= extent_thresh  %If lots of retreat (i.e. collapsed WAIS)
    %Get latitude and longitude data
    ULAT = ncread('Pliomip2_ULATLONG_extracted.nc','ULAT');
    [~,y] = find(ULAT<=-60);
    topl = max(y);
    ULONG = ncread('Pliomip2_ULATLONG_extracted.nc','ULONG');
    %Get only correct
    ULONG = ULONG(:,1:topl);
    ULONG(ULONG>180) = ULONG(ULONG>180)-360;
    ULAT = ULAT(:,1:topl);
    %Convert to polar stereo
    [PLIOX,PLIOY] = polarstereo_fwd(ULAT,ULONG); 
end
 
disp('Reading in and processing monthly data...')
for month = 1:surf_months
    if tuning ~= 1
        disp(['Month ' num2str(month)])
    end
    %%% Read in ocean temperature data
    %Get correct month's data. Use most recent.
    surf_month_vec = 325-surf_months:324;
    startLoc = [1 1 1 surf_month_vec(month)];
    count = [Inf Inf Inf 1];
    surf_temps = 'global-reanalysis-phy-001-031-grepv2-monthly_1660130057544.nc';
    sst = double(ncread(surf_temps,char('thetao_oras'),startLoc,count));
    
    % now interpolate to ice sheet model grid - uneven spacing, have to use
    % ScatteredInterp. 
    Fa = scatteredInterpolant(X(:),Y(:),sst(:)); 
    sst = Fa(ys,xs); 
    sst = inpaint_nans(sst,4);
    sst = sst + ocean_T; %Add temperature offset, if required.
    % calculate melt rate in m/day Russell-Head (1980). Temporarily held in
    % variable 'meltrate2'.
    meltrate2 = real(60*60*24*2.08*10^-7.*(1.8+sst).^1.5);

    if shelfmelt == 1
        %Add ice shelf melt rate to ocean data. Some NaNs in shelfm, so set
        %these to the (interpolated) ocean data.
        meltrate_orig = meltrate2;
        meltrate2(imrotate(flipud(maskwater),270)==1 & indi==1) = ...
        shelfm(imrotate(flipud(maskwater),270)==1 & indi==1);
        meltrate2(isnan(meltrate2)) = meltrate_orig(isnan(meltrate2));
    end

    meltrate(:,:,month) = meltrate2;

    %%% Read in ocean velocity data
    V = double(ncread(filename,char('vo_oras'),startLoc,count));
    U = double(ncread(filename,char('uo_oras'),startLoc,count));  
 
    %convert velocity field to stereo
    [vx,vy] = uv2vxvy(lat,long,U,V);
       
    % now interpolate to ice sheet model grid - uneven spacing, have to use
    % ScatteredInterp. 
    Fa = scatteredInterpolant(X(:),Y(:),vx(:)); 
    Uout = Fa(ys,xs); 
    Fb = scatteredInterpolant(X(:),Y(:),vy(:)); 
    Vout = Fb(ys,xs);

    %%% Read in paleo ocean velocity data if small ice sheet
    if extent <= extent_thresh  %If lots of retreat (i.e. collapsed WAIS)
        data = load('UVELcomb.mat');
        UVEL = data.UVEL;
        data = load('VVELcomb.mat');
        VVEL = data.VVEL;
        VVEL = VVEL(:,:,1,surf_month_vec(month));
        UVEL = UVEL(:,:,1,surf_month_vec(month)); 
        %convert velocity field to stereo
        [PLIO_vx,PLIO_vy] = uv2vxvy(ULAT,ULONG,UVEL,VVEL);
        % now interpolate to ice sheet model grid - uneven spacing, have to use
        % ScatteredInterp. 
        Fa_PLIO = scatteredInterpolant(PLIOX(:),PLIOY(:),PLIO_vx(:)); 
        UPLIO = Fa_PLIO(ys,xs); 
        Fb_PLIO = scatteredInterpolant(PLIOX(:),PLIOY(:),PLIO_vy(:)); 
        VPLIO = Fb_PLIO(ys,xs); 
        
        if forcewsb > 0
            umod = ncread('Modern_control_deconto_fort.92.nc','ubot',[1 1 1],[Inf Inf 1])';
            vmod = ncread('Modern_control_deconto_fort.92.nc','vbot',[1 1 1],[Inf Inf 1])';
            vscaling = 1/5000; %Assume ice vel. in m/yr is ~5000x ocean vel in m/s
            if forcewsb == 1
                VPLIO(350:420,70:160) = vscaling.*vmod(70:160,350:420)';
                UPLIO(350:420,70:160) = vscaling.*umod(70:160,350:420)';
            elseif forcewsb == 2
                VPLIO = vscaling.*vmod';
                UPLIO = vscaling.*umod';
            end   
        end
        
        % Flag the margin between the modern/paleo datasets to blend later.
        moderndata = ~isnan(Vout);
        flag = zeros(sx,sy,surf_months);
        for i = mrg+1:sx-mrg % mrg = number of cells to blend over
            for j = mrg+1:sy-mrg
            subset = moderndata(i-mrg:i+mrg,j-mrg:j+mrg);
                if sum(subset,'all') ~= (1+2*mrg)^2 && sum(subset,'all') ~= 0
                    flag(i,j,month) = 1;
                end
            end
        end
        
        %Fill holes in modern ocean velocity data with Pliocene data
        Vout(isnan(Vout)) = VPLIO(isnan(Vout))/100;
        Uout(isnan(Uout)) = UPLIO(isnan(Uout))/100;
    end
    
    Uout(indi) = NaN;
    Vout(indi) = NaN;

    Uout_all(:,:,month) = Uout;
    Vout_all(:,:,month) = Vout;  
end

if surf_simp == 1
    [x2, y2] = find(packages>0);
    disp(['Reducing streamlines to ' num2str(numel(x2)) ' tracks.'])
end

disp('Starting parfor loops...')
if PC == 1
    %Set up parfor loop
    pc = parcluster('local');
    parpool(pc,W)
end

if extent>extent_thresh
    flag = zeros(1,1,surf_months); %Set flag variable for parfor (not used).
end

% for month = 1:surf_months
parfor (month = 1:surf_months,W)
    h22 = 0;
    if tuning ~= 1
        disp(['Month ' num2str(month)])
    end

    %Get month's U and V ocean surface velocity data
    Uout = Uout_all(:,:,month);
    Vout = Vout_all(:,:,month);
        
    %% Interpolate around margin to ensure all velocities non-zero
    %interpolate where ice velocities are 0
    if shelfmelt == 1
        ubotB = double(ubot2');
        vbotB = double(vbot2');
    else
        ubotB = double(ubot');
        vbotB = double(vbot');
    end
    ubotB(ubotB==0) = NaN;
    vbotB(vbotB==0) = NaN;
    vbotB = inpaint_nans(vbotB,4);
    ubotB = inpaint_nans(ubotB,4);
    
    if shelfmelt == 1
    % Set ocean velocities to shelf velocity where present
         Vout(shelf==1) = vbotB(shelf==1)./(365.25*24*60*60);
         Uout(shelf==1) = ubotB(shelf==1)./(365.25*24*60*60);
    end

    Uout_old = Uout;
    Vout_old = Vout;
    
    %Loop to blend together small band around the edge of the ocean velocity data
    %with ice velocities.
    for i = mrg2+1:(sy-mrg2)
        for j = mrg2+1:(sx-mrg2)
            Vsub = Vout_old((i-mrg2:i+mrg2),(j-mrg2:j+mrg2)); %Get subset
            naans = sum(sum(isnan(Vsub)));
            %If at least 1 data point but at least 1 NaN and the centre is a NaN
            if naans>0 && naans<(2*mrg2+1)^2 && isnan(Vout_old(i,j))
                Vsub_oc = inpaint_nans(Vsub,4); %fill gaps based on ocean velocities.
                Vsub_ice = vbotB((j-mrg2:j+mrg2),(i-mrg2:i+mrg2)); %interpolated ice velocities in subset.
                fact = abs(mean(mean(Vsub_oc)))/abs(mean(mean(Vsub_ice))); %scaling factor between ocean and ice velocities.
                Vsub_ice = fact.*Vsub_ice; %change ice velocities so similar magnitude to ocean.
                Vsub = (vweighting*Vsub_oc + Vsub_ice)./(vweighting+1); %Take mean of ocean and scaled ice velocities, weighting twice towards ocean velocities.
                Vout((i-mrg2:i+mrg2),(j-mrg2:j+mrg2)) = Vsub;
%                 nearmargin(i,j) = 1;
            end
            Usub = Uout_old((i-mrg2:i+mrg2),(j-mrg2:j+mrg2));
            naans = sum(sum(isnan(Usub)));
            if naans>0 && naans<(2*mrg2+1)^2 && isnan(Uout_old(i,j))
                Usub_oc = inpaint_nans(Usub,4);
                Usub_ice = ubotB((j-mrg2:j+mrg2),(i-mrg2:i+mrg2));
                fact = abs(mean(mean(Usub_oc)))/abs(mean(mean(Usub_ice)));
                Usub_ice = fact.*Usub_ice;
                Usub = (vweighting*Usub_oc + Usub_ice)./(vweighting+1);
                Uout((i-mrg2:i+mrg2),(j-mrg2:j+mrg2)) = Usub;
%                 nearmargin(i,j) = 1;
            end
        end
    end
    
    % Fill gaps using interpolation
    VoutC = Vout;
    UoutC = Uout;
    if extent <= extent_thresh
        flagn = flag(:,:,month);
        UoutC(flagn==1) = NaN; %Fill overlap by first setting to NaN.
        VoutC(flagn==1) = NaN;
    end
    VoutC = inpaint_nans(VoutC); 
    UoutC = inpaint_nans(UoutC); 

    % Set areas with ice to NaN to stop ocean flowlines inland.
    UoutC(indi) = Uout(indi);
    VoutC(indi) = Vout(indi); 

    Vout = VoutC;
    Uout = UoutC;

    %% Get ocean velocity streamlines
    % Using ice end points as ocean start points
    fig_no = 3; 
    if plot_figs == 1
        hold on
        contour(1:sx,1:sy,maskwater,1,'k')
        axis equal
        axis off
    end

if model == 1
    h22 = streamline(x22,y22,Uout',Vout',x2,y2,[.25,5000]);
elseif model == 2
    h22 = streamline(x22,y22,Uout,Vout,x2,y2,[.25,5000]);
else
    disp('Warning: Ice sheet model number (''model'') invalid!')
end

    ocean_eNd = zeros(size(Uout)); %records total eNd
    ocean_count = zeros(size(Uout)); %records relative sediment amount
    berg_no = ocean_count;
    sitetrace = ocean_count; %records origins of streamlines
    sitetraceB = ocean_count;
    savezone = ocean_count; %records map of area streamlines are saved over.
    ocean_log = zeros(1,11);

    Icemag = sqrt(Uout.^2+Vout.^2); %magnitude of iceberg velocity for meltout calc (m/s)
    Icemag = Icemag*60*60*24; %convert to m/day

    for k = 1:length(h22) % for each stream line
        xk = get(h22(k),'xdata'); % get its x data
        yk = get(h22(k),'ydata'); % get its y dat            
        xk = xk(~isnan(xk)); %remove NaNs if present
        yk = yk(~isnan(yk));

        icelost = 0; %reset cumulative amount of basal ice melted (m).
        icelost_max = 0; %reset maximum amount of basal ice melted (m).

        % Loops calculating relative sediment volume, eNd values and number
        % of iceberg tracks through each ocean cell.
        if xk(1)~= xk(end) % Occasionally a track just stops - ignore these.
            % Loop over each point in streamline
            for p = 2:(length(xk)-1)
                if ~isnan(xk(p)) && ~isnan(yk(p)) %if not NaN
                    if model == 1
                        xkc = round(xk(p));
                        ykc = round(yk(p));
                    else
                        ykc = round(xk(p));
                        xkc = round(yk(p));
                    end
                    %calculate distance from previous cell (m).
                    distance = resolution*1000*...
                          sqrt(abs(xk(p-1)-xk(p)).^2+abs(yk(p-1)-yk(p)).^2);
                    %Using the ocean surface velocity magnitude, calculate 
                    %the time elapsed after leaving the previous cell (days).
                    elapsed = distance/Icemag(xkc,ykc);
                    %calculate the amount of debris-rich ice lost for this
                    %'iceberg'(m). Can be negative!
                    icelost = icelost + meltrate(xkc,ykc,month)*elapsed;
                    
                    % Monitor maximum ice lost in case of refreezing.
                    if icelost > icelost_max
                       icelost_max = icelost;
                    end

                    %Get corresponding debris concentration
                    if icelost_max == icelost %If ice is not refrozen or melt rate 0
                        if (basalsed-icelost) <= plataeu %If all basal ice gone
                            D = plataeu;
                        else
                            if sed_dist == 1
                                D = 1 - (icelost/basalsed);
                            elseif sed_dist == 2
                                D = exp(-decay*(icelost/basalsed));
                            else
                                D = 1;
                            end
                            D = max(D,plataeu); % In case very little ice lost.
                        end
                    else %If refrozen ice (no debris). Should only be under shelves.
                       D = 0; 
                    end
                    % Calculate relative amount of debris deposited 
                    % (i.e. melted ice thickness* debris concentration).
                    if surf_simp == 1
                        %Include a terrestrial sediment volume weighting.
                        sed_prop = meltrate(xkc,ykc,month)*elapsed*D*...
                            packages(round(xk(1)),round(yk(1)));
                    else
                        if seedmeth == 0
                            %in this instance, the sediment volume
                            %weighting is provided by the repetition of
                            %streamlines.
                            sed_prop = meltrate(xkc,ykc,month)*elapsed*D;
                        else
                            %Otherwise, the inital terrestrial erosion rate
                            % must be used to weight sed volume.
                            sed_prop = meltrate(xkc,ykc,month)*elapsed*D*...
                                quarry2(k);
                        end
                    end
            
                    %add point to eNd map 
                    epsilonval = average(round(xk(1)),round(yk(1)));
                    ocean_eNd(xkc,ykc) = ocean_eNd(xkc,ykc) + epsilonval*sed_prop;
                    % add relative proportion of sediment
                    ocean_count(xkc,ykc) = ocean_count(xkc,ykc) + sed_prop; %Can be greater than 1 if lots of melt
                    berg_no(xkc,ykc) = berg_no(xkc,ykc) + 1; %Record number of 'iceberg' tracks for figure.
                    
                    if model == 2
                        ykc_temp = ykc;
                        ykc = xkc;
                        xkc = ykc_temp;
                    end
                    
                    % keep record of streamlines near site. Record eNd of
                    % start for histogram.
                    if alond(xkc,ykc)>143 && alond(xkc,ykc)<145 && alatd(xkc,ykc)<-63.9 && alatd(xkc,ykc)>-64.8 %(xkc == sitex_loc) && (ykc == sitey_loc)
                        savezone(xkc,ykc) = 1;
                        if epsilonval >=0
                           ocean_log(1) = ocean_log(1) + sed_prop;
                        elseif epsilonval <0 && epsilonval>-2
                            ocean_log(2) = ocean_log(2) + sed_prop;
                        elseif epsilonval <-2 && epsilonval>-4
                            ocean_log(3) = ocean_log(3) + sed_prop;
                        elseif epsilonval <-4 && epsilonval>-6
                            ocean_log(4) = ocean_log(4) + sed_prop;
                        elseif epsilonval <-6 && epsilonval>-8
                            ocean_log(5) = ocean_log(5) + sed_prop;
                        elseif epsilonval <-8 && epsilonval>-10
                            ocean_log(6) = ocean_log(6) + sed_prop;
                        elseif epsilonval <-10 && epsilonval>-12
                            ocean_log(7) = ocean_log(7) + sed_prop;
                        elseif epsilonval <-12 && epsilonval>-14
                            ocean_log(8) = ocean_log(8) + sed_prop;
                        elseif epsilonval <-14 && epsilonval>-16
                            ocean_log(9) = ocean_log(9) + sed_prop;
                        elseif epsilonval <-16 && epsilonval>-18
                            ocean_log(10) = ocean_log(10) + sed_prop;
                        elseif epsilonval <=-18
                            ocean_log(11) = ocean_log(11) + sed_prop;
                        end
                        %Variables to store streamlines crossing specified
                        %transect.
                        sitetrace(round(xk(1)),round(yk(1))) = sitetrace(round(xk(1)),round(yk(1)))+1;
                        sitetraceB(round(xk(1)),round(yk(1))) = sitetraceB(round(xk(1)),round(yk(1)))+sed_prop;
                    end
                end
            end
        end
    end
    
    oceantracks = ocean_eNd./ocean_count;
    oceantracks_out(:,:,month) = oceantracks; %Save data
    ocean_logout(:,month) = ocean_log;

    if sed_dist_figs == 1 %To make sedimentation distribution figure
        ocean_count_out(:,:,month) = ocean_count;
        berg_no_out(:,:,month) = berg_no;
    end
    if sourcetrace_figs == 1
    sitetrace_out(:,:,month) = sitetrace;
    sitetraceB_out(:,:,month) = sitetraceB;
    savezone_out(:,:,month) = savezone;
    end

    if plot_figs == 1
        fig3 = ['Surface_flowpaths_' runname num2str(month) '.png'];
        print(gcf,fig3,'-dpng','-r300');
        close
    else
        close all
    end    
end %end parfor

if plot_figs == 1
    if sourcetrace_figs == 1
        sitetrace = sum(sitetrace_out,3);
        sitetraceB = sum(sitetraceB_out,3);
        sitetrace(sitetrace == 0) = NaN;
        savezone = sum(savezone_out,3);
        maskwater_trace = maskwater;
        maskwater_traceB = maskwater;
        maskwater_trace(sitetrace>0) = 0;
        maskwater_traceB(sitetraceB>0) = 0;
        
        figure(fig_no)
        fig_no = fig_no + 1;
        subplot(2,2,1)
        
        mappingplotter(alatd,alond,savezone,0,1,viridis);
        axis equal
        axis off
        colorbar;
        hold on
        freezeColors
        contourm(alatd,alond,imrotate(flipud(double(maskwater_trace)),270),'LineColor','k')
        hold on
        [sitetracex, sitetracey] = find(sitetrace>0);
        
        subplot(2,2,2)
        contour(1:sx,1:sy,maskwater_trace,1,'k')
        axis equal
        axis off
        hold on
        scatter(sitetracex,sitetracey,4,sitetrace(sitetrace>0),'filled')
        colorbar;

        subplot(2,2,3)
        contour(1:sx,1:sy,maskwater_traceB,1,'k')
        hold on
        [sitetraceBx, sitetraceBy] = find(sitetraceB>0);
        scatter(sitetraceBx,sitetraceBy,4,sitetraceB(sitetraceB>0),'filled')
        axis equal
        axis off
        colorbar;
        hold on
        scatter(sitex_loc,sitey_loc,5,'x')

        print(gcf,figx,'-dpng','-r300')
        close
    end

    ocean_logout = sum(ocean_logout,2);
    % plot result
    figure(fig_no)
    fig_no = fig_no + 1;
    bar(ocean_logout)
    print(gcf,figy,'-dpdf','-bestfit')
    close
end

%Get  surface current data
oceantracks = oceantracks_out;
oceantracks_out(~isnan(oceantracks_out)) = 1; 
oceantracks_out(isnan(oceantracks_out)) = 0;
trackno = sum(oceantracks_out,3);
oceantracks(isnan(oceantracks)) = 0;
ocean_eNd = sum(oceantracks,3);
ocean_eNd(ocean_eNd==0) = NaN;
oceantracks = ocean_eNd./trackno;

if sed_dist_figs == 1
    % Get sediment density array
    oceancountno = sum(ocean_count_out,3);
    oceancountno(oceancountno==0) = NaN;

    if plot_figs == 1
        figure(fig_no)
        fig_no = fig_no + 1;
        %Some very high points can throw off colourbar, so set limit to improve figure.
        if surf_simp == 1
	        colmax = surf_months*1.5; %1.5 manually adjusted to make figure look nice.
        else
	        if	seedmeth == 0
	            colmax = surf_months*s_points/1000; %1000 manually adjusted to improve figure.
            else
                colmax = surf_months*numel(x2)/1000; %1000 manually adjusted to improve figure.
	        end 
        end
    	oceancountno(oceancountno>colmax) = colmax;
        %Rescale data to 0:255 for colourbar (size 256).
        colour_int = max(max(oceancountno))./255;
        oceancountno = ceil(oceancountno./colour_int);
        oceancountno(oceancountno<=1) = 1; %Make sure all paths with icebergs above 1.
        %Set NaNs to 0.
        oceancountno(isnan(oceancountno)) = 0;
        %Modify 1st colourbar value for white background
        hotmod = hot;
        hotmod(1,:) = [1 1 1];
        %Plot using subfunction and adjust figure
        mappingplotter(alatd,alond,oceancountno,0.96,255,hotmod);   
        set(gca,'ColorScale','log')
        colormap(hotmod)
        axis equal
        axis off
        colorbar;
        hold on
    %Plotting ice margin seems to cause trouble with plotting.
    %     freezeColors 
    %     contourm(alatd,alond,double(indi),'LineColor','k')
        title('Sediment delivery heatmap')
        print(gcf,fig4b,'-dpng')
        close
  
    % Iceberg path figure
        figure(fig_no)
        fig_no = fig_no + 1;
        % Get sediment density array
        bergcountno = sum(berg_no_out,3);
        bergcountno(bergcountno==0) = NaN;
        %Adjust maximum to improve plot
        colmax = max(max(bergcountno));
        bergcountno(bergcountno>(colmax/100)) = colmax/100; %100 manually adjusted to improve figure
        %Rescale data to 0:255 for colourbar (size 256).
        colour_int = max(max(bergcountno))./256;
        bergcountno = ceil(bergcountno./colour_int);
	
        %Set NaNs to 0.
        bergcountno(isnan(bergcountno)) = 0;
        %Modify 1st colourbar value for white background
        hotmod = hot;
        hotmod(1,:) = [1 1 1];
        %Plot using subfunction and adjust figure
        mappingplotter(alatd,alond,bergcountno,0.96,255,hotmod);   
        set(gca,'ColorScale','log')
        colormap(hotmod)
        axis equal
        axis off
        colorbar;
        hold on
        title('Iceberg path heatmap')
        
        print(gcf,fig4c,'-dpng')
        close

    end
end

%Select only points within 'radius2' km of data
[x,y] = find(~isnan(oceantracks));
surfcur_mask = zeros(sx,sy); %map of points within radius
% Loop over each grid cell, and for each one loop over each point
% filled by a streamline and get the distance from it.
for i = 1:sx
    for j = 1:sy
        for k = 1:size(x,1)
            if x(k)~=0 && y(k)~=0
                dist = resolution*sqrt((xcoords(x(k))-xcoords(i)).^2+...
                    (ycoords(y(k))-ycoords(j)).^2);
                if dist<(resolution*radius4)
                    surfcur_mask(i,j) = 1; % Set to 1 if close to a bottom current path.
                    break % Exit loop and move to next cell.
                end
            end
        end
    end
end
surfcur_mask(surfcur_mask==0) = NaN; %Set zeros to NaNs

if model == 2
    oceantracks = oceantracks';
end

%fill gaps between ocean streamlines, then set ice (grounded/no velocity) to NaN
oceantracks = inpaint_nans(oceantracks,4);
oceantracks = surfcur_mask.*oceantracks;
%Set ice to 0.
oceantracks(vbot'~=0) = NaN;
oceantracks(ubot'~=0) = NaN;
oceantracks(maskwater'==0) = NaN;
    
if smth > 0 %Option to apply a moving average smoothing.
    oceantracks = movmean_2d(oceantracks,smth);
end
oceantracks(vbot'~=0) = NaN;
oceantracks(ubot'~=0) = NaN;
oceantracks(maskwater'==0) = NaN;

% If making an estimate under ice shelves, add to surface current estimate.
if shelfmelt == 0 && undershelf == 1
    oceantracks(isnan(oceantracks)) = shelf_average(isnan(oceantracks));   
end

if plot_figs == 1
    % plot result
    figure(fig_no)
    fig_no = fig_no + 1;
    ax(1) = subplot(2,1,1);
    mappingplotter(alatd,alond,oceantracks,-20,5,viridis);
    colormap(ax(1),viridis)
    hold on
    freezeColors
    contourm(alatd,alond,double(flipud(imrotate(maskwater,90))),'LineColor','k')
    title('Surface Current \epsilon_{Nd} map')

    %Get model-data difference and limit to 'radius' of samples.
    surf_diff = coretops-oceantracks;
    nrst_smpl_nans = nrst_smpl;
    nrst_smpl_nans(nrst_smpl_nans == 0) = NaN;
    surf_diff = surf_diff.*nrst_smpl_nans;

    ax(2) = subplot(2,1,2);
    mappingplotter(alatd,alond,surf_diff,-5,5,viridis);
    colormap(ax(2),divmap)
    hold on
    freezeColors
    contourm(alatd,alond,double(flipud(imrotate(maskwater,90))),'LineColor','k') 
    set(gcf,'Position',[50 50 600 1200]);
    title('Seafloor surface sediment - Surface current estimate')
    
    print(gcf,fig4a,'-dpdf','-fillpage')
    close 
    
end

% If running on HPC cluster, run function to merge monthly surface output.
if PC == 1 && plot_figs == 1
    streamline_merger(runname,t,surf_months);
end
