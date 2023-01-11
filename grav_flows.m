function [oceantracks3,fig_no,iter] = grav_flows(x0_orig,y0_orig,sx,sy,hb,slopethresh,...
    maskwater,coretops,oceantracks,oceantracks2,shelfmin,nrst_smpl_nans,alatd,...
    alond,divmap,fig8,iter,repeat,rep,plot_figs,resolution,paleo,fig_no,fig7,...
    vbot,ubot,model,xcoords,ycoords,radius5)
% Estimate debris movement due to gravity flows. This reads in several
% variables from prov.m and outputs an eNd value estimate, figure number and
% the iteration number if iterating with the bottom current method to
% represent interaction between these processes.

disp('Estimating offshore transport via gravity flows...')
    
[x3,y3] = meshgrid(x0_orig,y0_orig); %x and y coordinates 
% direction vectors (up, down, left, right, diagonals)
ish = [1 -1 0 0 -1 1 -1 1];
jsh = [0 0 -1 1 -1 -1 1 1];

%Create blank arrays.
slump_eNd = zeros(sx,sy);
slump_count = slump_eNd;
slopeedge = zeros(sx,sy);
slopeedgeNd = zeros(sx,sy);
bestoftwo = zeros(sx,sy);

% Get slope map (degrees). ArcGIS formula used:
% (https://desktop.arcgis.com/en/arcmap/10.3/tools/spatial-analyst-toolbox/how-slope-works.htm)
slope = zeros(sx,sy);
slopeedge = slope;
for i = 2:sx-1
    for j = 2:sy-1
        dzdx = ((hb(i+1,j-1)+2*hb(i+1,j)+...
            hb(i+1,j+1))-(hb(i-1,j-1)+2*hb(i-1,j)+hb(i-1,j+1)))/(8*(resolution*1000));
        dzdy = ((hb(i-1,j+1)+2*hb(i,j+1)+...
            hb(i+1,j+1))-(hb(i-1,j-1)+2*hb(i,j-1)+hb(i+1,j-1)))/(8*(resolution*1000));
        slope(i,j) = atan(sqrt(dzdx^2+dzdy^2))*57.29578; %Find slope and convert to degrees.
    end
end
slopeedge(slope>=slopethresh) = 1; %Get areas above slope threshold.
slopeedge(hb<-2000) = 0; %Exclude deep water areas as not relevant.
slopeedge(flipud(imrotate(maskwater,90)) == 0) = 0; %Exclude inland areas.

% Calculate values to use at start of gravity flow paths.
if paleo == 0
    % Get nearest match to core tops at shelf break
    difference = abs(oceantracks - coretops);
    difference2 = abs(oceantracks2 - coretops);
    method_match = min(difference,difference2);
    bestoftwo(method_match==difference) = oceantracks(method_match==difference);
    bestoftwo(method_match==difference2) = oceantracks2(method_match==difference2);
else
    % If paleo, assume an average of bottom/ surface currents.
    bestoftwo = (oceantracks + oceantracks2)./2;
    % Fill with one value if the other is NaN.
    bestoftwo(isnan(bestoftwo)) = oceantracks(isnan(bestoftwo));
    bestoftwo(isnan(bestoftwo)) = oceantracks2(isnan(bestoftwo));
end

[icoords,jcoords] = find(slopeedge == 1);

for q = 1:size(icoords,1) %Loop over cells on slope break
    j = jcoords(q);
    i = icoords(q);
    ivec = 0;
    jvec = 0;
    p = 0;
    %From starting point at ice margin, loop downslope unless edge
    %reached
    while j<sy && j>0 && i<sx && i>0
        p = p+1;
        %Set current point to value at shelf break
        slump_eNd(i,j) = slump_eNd(i,j) + bestoftwo(icoords(q),jcoords(q));
        slump_count(i,j) = slump_count(i,j) + 1;
        %Find steepest-descent neighbor
        zsmax = -1e20;
        for look = 1:8
            ii = max(1,min(sx, i + ish(look)));
            jj = max(1,min(sy, j + jsh(look)));
            if model == 1
                zx = sqrt((x3(i,j)-x3(ii,jj))^2+(y3(i,j)-y3(ii,jj))^2);
            else
                zx = sqrt((x3(j,i)-x3(jj,ii))^2+(y3(j,i)-y3(jj,ii))^2);
            end
            zs = (hb(i,j) - hb(ii,jj)) / max (1.e-20,zx);
            if zs>zsmax
                zsmax = zs; % slope, m/m
                iimax = ii; % x coord downslope
                jjmax = jj; % y coord downslope
            end
        end
        %record coordinate to check it's not repeated
        ivec(p) = i;
        jvec(p) = j;
        %move to next point
        i = iimax;
        j = jjmax;
        %Loop over previous points and if the same, exit loop.
        if p>=2
            for vec = 1:(size(ivec,2)-1)
                if ivec(vec) == i && jvec(vec) == j
                    i = 0;
                    j = 0;
                end
            end
        end
    end
end

oceantracks3 = slump_eNd./slump_count;

if plot_figs == 1
    figure(fig_no);
    fig_no = fig_no + 1;
    mappingplotter(alatd,alond,oceantracks3,-20,5,viridis);
    colormap(viridis)
    hold on
    freezeColors
    contourm(alatd,alond,double(flipud(imrotate(maskwater,90))),'LineColor','k') 
    title('Gravity flow path \epsilon_{Nd}')

    print(gcf,fig7,'-dpdf','-bestfit')
    close
end

%Select only points within 'radius5' km of data
[grav_x,grav_y] = find(~isnan(oceantracks3)); %Get indices where there's a flow path
grav_x = grav_x';
grav_y = grav_y';
grav_dist = zeros(size(grav_x)); %vector containing distance to each feature
grav_msk = zeros(sx,sy); %Variable containing mask
for i = 1:sx %Loop over each cell, deciding whether it's in the mask.
    for j = 1:sy
        for k = 1:size(grav_x,2)
            if grav_x(k)~=0 && grav_y(k)~=0
                grav_dist(k) = resolution*sqrt((xcoords(grav_x(k))-xcoords(i)).^2+...
                    (ycoords(grav_y(k))-ycoords(j)).^2);
            else
                grav_dist(k) = 10^6; %set to arbritary high value so definitely > radius
            end
        end
        if min(grav_dist)<(resolution*radius5)
            grav_msk(i,j) = 1; %set to 1 if close to a core top.
        end
    end
end

%fill gaps between ocean streamlines
oceantracks3 = flipud(imrotate(oceantracks3,90));
oceantracks3(oceantracks3==0) = NaN;
oceantracks3_orig = oceantracks3; %Save rotated uninterpolated version.
oceantracks3 = inpaint_nans(oceantracks3,4);

%Set shelf to NaN except where direct flow path.
oceantracks3((hb'>shelfmin)&(isnan(oceantracks3_orig))) = NaN;
oceantracks3 = flipud(imrotate(oceantracks3,90));

% Set ice (grounded/no velocity) to NaN
oceantracks3(vbot'~=0) = NaN;
oceantracks3(ubot'~=0) = NaN;
oceantracks3(maskwater'==0) = NaN;

% Apply maximum distance mask
oceantracks3 = oceantracks3.*grav_msk;
oceantracks3(grav_msk==0) = NaN;
%    % Exclude areas near active volcanoes
%     oceantracks3 = volc_msk.*oceantracks3;

oceantracks3 = imrotate(oceantracks3,90);

if plot_figs == 1
    figure(fig_no);
    fig_no = fig_no + 1;
    
    %Get model-data difference and limit to 'radius' of samples.
    grav_diff = coretops-imrotate(oceantracks3,270);
    grav_diff = grav_diff.*nrst_smpl_nans;

    ax(2) = subplot(2,1,2);
    mappingplotter(alatd,alond,grav_diff,-5,5,divmap);
    colormap(ax(2),divmap)
    hold on
    freezeColors
    contourm(alatd,alond,double(flipud(imrotate(maskwater,90))),'LineColor','k') 
    title(['Measured - Gravity Flow \epsilon_{Nd}'])

    ax(1) = subplot(2,1,1);
    mappingplotter(alatd,alond,imrotate(oceantracks3,270),-20,5,viridis);
    colormap(ax(1),viridis)
    hold on
    freezeColors
    contourm(alatd,alond,double(flipud(imrotate(maskwater,90))),'LineColor','k') 
    title('Gravity Flow \epsilon_{Nd}')
    set(gcf,'Position',[50 50 600 1200]);

    print(gcf,fig8,'-dpdf','-fillpage')
    close
end

% If iterating bottom currents/gravity flows, use gravity flow output
% as input where possible and bottom current elsewhere.
if repeat == 1 && rep > 1
    iter = imrotate(oceantracks3,270);
    iter(isnan(iter)) = oceantracks2(isnan(iter));
else
    iter = 0;
end