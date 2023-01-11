function [average,fig_no,x2_out,y2_out,shelf_average,packages,quarry2]...
    = terrestrial_part(sed_mod,t,vmag,maskwater,nc1,sx,sy,kval,alond,alatd,biasWAIS,...
    s_points,fig1,x22,y22,ubot,vbot,pv,indi,undershelf,resolution,x0,y0,...
    fig_no,plot_figs,paleo,PC,fig2,hgrams,model,shelf,seedmeth,seedscale,...
    q_thresh,shelfmelt)

%Calculates erosion rate from ISM output and predicts provenance value at
%the ice sheet margin.

%% Get erosion rate
if sed_mod == 1
    quarryrate = ncread(nc1,'quarryrate'); %read in quarry rate
    quarryrate = quarryrate(:,:,t);
elseif sed_mod == 0  %approximation if there is no quarryrate from model.
    E = vmag.^-2;
    E(E>300) = 300;
    maskwater2 = abs(maskwater-1);
    maskwater2 = double(maskwater2);
    quarryrate = maskwater2.*E(:,:);
elseif sed_mod == 2
    % Uses heatb output (i.e. basal shear stress * basal velocity) to get
    % erosion rate.
    if model == 1
        heatb = ncread(nc1,'heatb',[1 1 t],[Inf Inf 1]);
    elseif model == 2
        taub = ncread('taud-correction-factor.nc','taub_mag');
        taub(isnan(taub)) = 0;
        taub(isinf(taub)) = 0; %At edge, sometimes hit inf.
        taub(taub>1) = 1;
        vmag(isnan(vmag)) = 0;
        heatb = taub.*vmag'*1E9; %wants to be Pa and m/yr
        heatb = inpaint_nans(heatb,4);
    end
    kmap = zeros(sx,sy);  
    if kval == 0
        % Use k (quarry coeff) values from D&P (2019). Sequentially: DML,
        % EL, MRL-PEL, QML-WL, G5l-OL, Ross, Bel-Amun, W Penins, Weddell.
        kv = 1E-10*[1.240558555 0.3053126918 0.4340391295 0.1377941166 ...
             0.2498541251 0.5382989050 0.1852198851 0.9329537193 0.3174665456];  
        for i = 1:sx
            for j = 1:sy
                if alond(i,j)>=0 && alond(i,j)<=30
                    kmap(i,j) = kv(1);
                elseif alond(i,j)>30 && alond(i,j)<=60
                    kmap(i,j) = kv(2);
                elseif alond(i,j)>60 && alond(i,j)<=94
                    kmap(i,j) = kv(3);
                elseif alond(i,j)>94 && alond(i,j)<=124
                    kmap(i,j) = kv(4);
                elseif alond(i,j)>124 && alond(i,j)<=165
                    kmap(i,j) = kv(5);
                elseif alond(i,j)>165 && alond(i,j)<=180
                    kmap(i,j) = kv(6);
                elseif alond(i,j)<-150
                    kmap(i,j) = kv(6);
                elseif alond(i,j)<-80 && alond(i,j)>=-150
                    kmap(i,j) = kv(7);
                elseif alond(i,j)<-50 && alond(i,j)>=-80
                    kmap(i,j) = kv(8);
                elseif alond(i,j)<0  && alond(i,j)>=-50
                    kmap(i,j) = kv(9);
                end
            end
        end  
        kmap = movmean_2d(kmap,9); %apply 9 point smoothing as in P&D 2019
    else
        kmap = kmap + kval; %Use specified uniform k value, kval.
    end
    quarryrate = 1000.*heatb.*kmap; %get erosion rate in mm/yr
    quarryrate = imrotate(flipud(quarryrate),270);
    % If WAIS collapse, MBL ice caps are slow flowing and not well
    % represented so add seeds disproportionately.
    if biasWAIS == 1
        for i = 1:sx
            for j = 1:sy
                if alond(i,j)<=-65
                    quarryrate(i,j) = quarryrate(i,j).*biasWAIS;
                end
            end
        end
    end
else
    disp('Error: Please enter valid sed_mod option')
end


%% Generate seed locations based on quarry rate
if seedmeth == 0
    % Randomly seed locations dependant on the erosion rate.

    totquarry = sum(sum(quarryrate)); %total quarrying rate over all cells
    fracquarry = quarryrate./totquarry; %fraction of quarrying coming from each cell
    sumquarry = zeros(sy,sx);%cumulative quarrying rate
    runtot = 0;
    
    %loop over each cell. make each one cumulatively higher.
    for i = 1:sy
        for j = 1:sx
            runtot = runtot + fracquarry(i,j);
            sumquarry(i,j) = runtot;
        end
    end

    %As this part is relatively quick, automatically use more seed locations 
    % to reduce randomness.
    s_points = s_points*seedscale;
    disp(['Randomly generating ' num2str(s_points) ' seed locations...'])
    startx = zeros(s_points,1);
    starty = zeros(s_points,1);
    pointmap = zeros(sx,sy);
    
    for n = 1:s_points
        dummy = zeros(sx,sy);
        number = rand(1,1); % generate random number
        for i = 1:sy %loop over each cell, assign 1 if less than random number
            for j = 1:sx
                if sumquarry(i,j)<number
                    dummy(i,j) = 1;
                else
                    dummy(i,j) = 0;
                end
            end
        end
        indmin = sum(sum(dummy)); %sum number of cells less than random number
        starty(n) = floor(indmin/sx); %y coordinate
        startx(n) = round(sx*((indmin/sx) - starty(n))); % x coordinate
        %Very occasionally a point is selected outside of limits - issue warning.
        if startx(n) <= 0 || startx(n) >= sx
            disp('Warning: X coordinate generated outside limits!')
            startx(n) = 1;
        end
        if starty(n) <= 0 || starty(n) >= sy
            disp('Warning: Y coordinate generated outside limits!')
            starty(n) = 1;
        end
        pointmap(startx(n),starty(n)) = pointmap(startx(n),starty(n)) + 1;
    end
    
    if plot_figs == 1
        %Plot seed locations alongside quarry rate
        figure(fig_no)
        fig_no = fig_no + 1;
        subplot(2,1,1)  
        if sed_mod == 1
            mappingplotter(alatd,alond,flipud(imrotate(pointmap,90)),0,max(s_points/1000,1),viridis);
        else
            mappingplotter(alatd,alond,pointmap,0,max(s_points/4000,1),viridis);
        end
        
        colormap(viridis)
        hold on
        freezeColors
        contourm(alatd,alond,double(flipud(imrotate(maskwater,90))),'LineColor','k')
    
        subplot(2,1,2)
        quarryrate_plot = quarryrate;
        quarryrate_plot(maskwater==1) = NaN;
        quarryrate_plot(quarryrate_plot<(8/255)) = 8/255;
        if paleo == 1
            mappingplotter(alatd,alond,quarryrate_plot,0,2,viridis);
        else
            mappingplotter(alatd,alond,flipud(imrotate(quarryrate_plot,90)),0,8,viridis);
        end
        viridis_mod = viridis;
        viridis_mod(1,:) = 1;
        colormap(viridis_mod)
        hold on
        freezeColors
        contourm(alatd,alond,double(flipud(imrotate(2*maskwater,90))),'LineColor','k')
        title('Erosion Rate (mm yr^-^1)')
        
    %     if PC == 1
    %         print(gcf,fig1,'-dpdf','-bestfit','-r300'); %This makes the figure on BP.
    %     else
    %         set(gcf,'Position',[50 50 600 1200]);
        exportgraphics(gcf,fig1,'Resolution',300); 
    %     end
        close
    end
elseif seedmeth == 1
    [starty,startx] = find(quarryrate>q_thresh);
    disp(['Using ' num2str(numel(startx)) ' seed locations...'])
else
    disp('Invalid seedmeth value specified!')
    return
end

if sed_mod == 1
    %switching x and y to make it line up - seems to be needed.
    starty2 = starty;
    startx2 = startx;
    starty = startx2;
    startx = starty2;
end

%% Generate Streamlines
disp('Calculating ice sheet flow lines...')

if plot_figs == 1
    figure(fig_no)
    fig_no = fig_no + 1;
    colormap(viridis)
    hold on
end

h22 = streamline(x22,y22,ubot,vbot,startx,starty); %gets streamlines
 
if plot_figs == 1
    set(h22,'linewidth',0.1,'color',[0.5 0.5 0.5])
    axis equal
    set(gca,'ydir','normal')
    hold on
end

%Create some variables to hold start and end locations
x1b = zeros(size(h22));
y1b = zeros(size(h22));
x2 = zeros(size(h22));
y2 = zeros(size(h22));
shelf_output = zeros(sx,sy); %Contains summed eNd values
shelf_packages = zeros(sx,sy); %Contains the number of ice flow lines to pass through the cell.
output = zeros(sx,sy);
packages = zeros(sx,sy);

for k = 1:length(h22) % for each stream line.
    xk = get(h22(k),'xdata'); % get its x data
    yk = get(h22(k),'ydata'); % get its y data
    x1b(k) = xk(1);  %get start and end points
    x2(k) = xk(end);
    y1b(k) = yk(1); 
    y2(k) = yk(end);
end

if model == 2
    pv = imrotate(pv,90);
end

idx = sub2ind(size(pv),x1b,y1b);
idx = idx(~isnan(idx)); 
x2 = x2(~isnan(x2));
y2 = y2(~isnan(y2));

if seedmeth == 0
    %Shrink number of output end locations back to 1/x th of seed locations.
    %Take first tenth as representative for ocean methods.
    x2_out = x2(1:s_points/seedscale);
    y2_out = y2(1:s_points/seedscale);
else
    x2_out = x2;
    y2_out = y2;
end
        
pv2 = pv(idx); % Set eNd values to start cell values
%Get erosion rate values.
if seedmeth == 1
    quarryrate2 = imrotate(flipud(quarryrate),270);
    quarry2 = quarryrate2(idx);
else
    quarry2 = zeros(size(pv2));
end

if shelfmelt == 0 && undershelf == 1
    for k = 1:length(h22) % for each stream line.
        xk = get(h22(k),'xdata'); % get its x data
        yk = get(h22(k),'ydata'); % get its y data
        xk_round = round(xk);
        yk_round = round(yk);
        for n = 1:length(xk)
            if shelf(xk_round(n),yk_round(n)) == 1
                shelf_packages(xk_round(n),yk_round(n)) = shelf_packages(xk_round(n),yk_round(n)) + 1;
                shelf_output(xk_round(n),yk_round(n)) = shelf_output(xk_round(n),yk_round(n)) + pv2(k);
            end
        end
    end
    shelf_average = shelf_output./shelf_packages; %eNd values under ice shelves.
else
    shelf_average = 0; %Dummy value
end

disp('Assigning end point provenance values...')

%round each flow line end point to nearest cell
x2b = round(x2);
y2b = round(y2);

%loop over each end point and determine how many in each cell and total
%value.
for m = 1:numel(x2b)
    if seedmeth == 0
        packages(x2b(m),y2b(m)) = packages(x2b(m),y2b(m)) + 1;
        output(x2b(m),y2b(m)) = output(x2b(m),y2b(m)) + pv2(m); 
    else
        packages(x2b(m),y2b(m)) = packages(x2b(m),y2b(m)) + quarry2(m);
        output(x2b(m),y2b(m)) = output(x2b(m),y2b(m)) + pv2(m)*quarry2(m); 
    end
end

endcells = sum(sum(packages>0)); %number of cells which have an endpoint in them.
average = output./packages; %mean value of each cell

%% Create histogram data
if hgrams == 1 && seedmeth == 0
    maxpac = ceil(max(max(packages)));
    histograms = 100.*ones(endcells,maxpac);
    histograms = [1000.*ones(endcells,1) histograms];
    
    coords = [x2b y2b];
    
    if model == 1
        pvsize = size(pv2,1);
    else
        pvsize = size(pv2,2);
    end
    
    loc = ones(pvsize,1); %contains end point location numbers
    p = 1;
    
    for m = 2:pvsize %loop over end points
        new = 1;
        for n = m-1:-1:1 %loop over previous end points
            if coords(m,1)==coords(n,1) && coords(m,2)==coords(n,2) %if repeated
                new = 0;
                loc(m) = loc(n); %set loc to first occurance of these coordiantes
            end
        end
        if new == 1 %if not repeated, put new number in loc
            p = p+1;
            loc(m) = p;
        end
    end
    
    for m = 1:size(pv2)
        added = 0;
        for n = 2:maxpac+1
            if histograms(loc(m),n) == 100 && histograms(loc(m),n-1)~= 100 && added == 0
                histograms(loc(m),n) = pv2(m);
                added = 1; %stops value filling to end of row
            end
        end
    end
end

%% Make other plots
if plot_figs == 1
    %Ice Streamline Figure
    caxis([-20 5])
    colormap(viridis)
    contour(1:sx,1:sy,maskwater,1,'k')
    title('Streamlines and end points')
    % Plot mean output for each cell at ice margin
    [plotx, ploty] = find(~isnan(average)); %Get x and y corrodinates of endpoints
    hold on
    scatter(plotx,ploty,20,average(~isnan(average)),'filled')
    colorbar
%     xlim([100 150]) %Just plot Amundsen Sea
%     ylim([200 300])
    axis equal
    axis off
   
    print(gcf,fig2,'-dpng','-r800'); 
    close
end

%% Histograms
if hgrams == 1 && seedmeth == 0
    %generate histogram at ice sheet cell with most end points (xind, yind)
%     [xind, yind] = find(packages == max(max(packages)));
    xind = 1572; %looking at specific point
    yind = -2111;
    xind = round(xind./resolution).*resolution;
    yind = round(yind./resolution).*resolution;
    xind = find(xind==x0);
    yind = find(yind==y0);
    
    for n = 1:size(pv2)
        if x2b(n) == xind && y2b(n) == yind
            most = loc(n);
        end
    end
    figure(fig_no)
    fig_no = fig_no + 1; %plot histogram
    scatter(xind,yind,10,'k','o','filled')

    figure(fig_no)
    fig_no = fig_no + 1; %plot histogram
    histogram(histograms(most,2:end),'BinWidth',1)%(min(find(histograms(most,:)==100))-1))) %selects which point - most
    print(gcf,'Histograms.pdf','-dpdf','-bestfit');
    close
end

disp('Terrestrial Component Finished.')