function [fig_no,eNd,rmse_surf] = modern_analysis(plot_figs,coretops,...
            oceantracks,oceantracks2,oceantracks3,sitex_loc,sitey_loc,...
            x22,y22,maskwater,fig_no,y0,fig9,fig10,fig11a,fig11b,fig12,...
            fig13,fig13b,fig14,fig15,sx,sy,x0,grouping,hb,Ea,...
            We,nrst_smpl,coresites,calib_file,corex_loc,corey_loc,alatd,...
            alond,volc_msk,runname,PC,corelat,corelon,tuning)
% If using the modern ice sheet, does some analysis of results, makes some
% plots and generates a new calibration file.

% Make NaNs in core top interpolation so it looks nicer.
coretopnans = coretops;
coretopnans(isnan(oceantracks)) = NaN;  

%Print eNd values at input core site
disp('At core site:')
Bottom_Currents = oceantracks2(sitex_loc,sitey_loc)
Iceberg_Rafting = oceantracks(sitex_loc,sitey_loc)
oceantracks3 = imrotate(oceantracks3,270);
Gravity_Flow = oceantracks3(sitex_loc,sitey_loc)

best = zeros(sx,sy); % Best eNd fit of three methods.
threemeth = best; % Closest matching method. 
eNdvals = [0 0 0];

for i = 1:sx
    for j = 1:sy
        eNdvals = [oceantracks(i,j) oceantracks2(i,j) oceantracks3(i,j)];
        act = coretops(i,j);
        nearest = find(min(abs(eNdvals-act))==abs(eNdvals-act));
        if size(nearest,2) == 0 %no mechanisms available
            threemeth(i,j) = NaN;
            best(i,j) = NaN;
        elseif size(nearest,2) == 2 %two mechanisms the same.
            threemeth(i,j) = mean(nearest);
            if threemeth(i,j) == 2 %If 1 and 3 the same, want to avoid getting 2.
                threemeth(i,j) = 0;
            end
            best(i,j) = eNdvals(nearest(1));
        elseif size(nearest,2) == 3 %all three mechanisms the same
            threemeth(i,j) = 4;
            best(i,j) = eNdvals(nearest(1));
        else %one mechanism matches
            threemeth(i,j) = nearest;
            best(i,j) = eNdvals(nearest);
        end
    end
end

if plot_figs == 1
    %Plot core top data
    figure(fig_no)
    fig_no = fig_no + 1;
    imagesc(imrotate(coretopnans,90))
    colormap(viridis);
    colorbar;
    caxis([-20 5])
    title('Interpolated Core Top (Actual) \epsilonNd')
    hold on
    scatter(corex_loc,numel(y0)-corey_loc,3,'filled','white');
    contour(x22,y22,flipud(maskwater),1,'linewidth',1,'color','k')
    axis equal
    axis off
    print(gcf,fig9,'-dpdf','-bestfit')
    close
end

%% Calculate best match
%Get the nearest method match, assuming a mixture of methods is ok. Also
%get the method matching value, factr (1 = IRD, 2 = Bottom, 3 = gravity).
best_new = zeros(sx,sy); % Nearest matching eNd values.
factr = zeros(sx,sy);

for i = 1:sx
    for j = 1:sy
        eNdvals = [oceantracks(i,j) oceantracks2(i,j) oceantracks3(i,j) coretops(i,j)]';
        [sorted, index] = sort(eNdvals);
        pos = find(index==4); %Find where coretop value is
        if sum(isnan(eNdvals)) == 3 %If all methods are NaNs, set to NaN.
            best_new(i,j) = NaN;
        else
            if isnan(oceantracks2(i,j)) && isnan(oceantracks3(i,j))
                best_new(i,j) = oceantracks(i,j); %If surface velocities are only non-NaNs, set this.
            else
                if pos == 4 %if coretop value is the biggest value
                    best_new(i,j) = sorted(3);
                    factr(i,j) = index(3);   
                elseif pos == 3 && (isnan(oceantracks2(i,j)) || isnan(oceantracks3(i,j))) %Or coretop is biggest but there's a NaN
                    best_new(i,j) = sorted(2);
                    factr(i,j) = index(2); 
                elseif pos == 1 %if smallest value
                    best_new(i,j) = sorted(2);
                    factr(i,j) = index(2);
                else
                    %get fraction of amount coretop value lies between adjacent
                    %methods
                    fr = (sorted(pos)-sorted(pos-1))/(sorted(pos+1)-sorted(pos-1)); %sorted(pos) should be coretops(i,j)
                    if abs(index(pos+1)-index(pos-1)) == 1 %if method reference number is conveniently 1 apart
                       factr(i,j) = index(pos-1)-(index(pos-1)-index(pos+1))*fr;
                    else %if method reference numbers are 1 and 3, set 3 to 0.
                        if index(pos-1) == 3
                            index(pos-1) = 0;
                        else
                            index(pos+1) = 0;
                        end
                        factr(i,j) = index(pos-1)-(index(pos-1)-index(pos+1))*fr;
                    end            
                    %eNd is somewhere between estimates and therefore exact.
                    best_new(i,j) = coretops(i,j);
                end    
            end
        end
    end
end

x0_out = x0;
y0_out = y0;
save(calib_file,'factr','x0_out','y0_out');

eNd = best_new(sitey_loc,sitex_loc);
save([runname '_output.mat'],'best_new','oceantracks','oceantracks2','oceantracks3'); %Save output

%Extrapolate 1 cell around data to capture locations near coast.
best_new_interp = inpaint_nans(best_new);
oceantracks_interp = inpaint_nans(oceantracks);
for n = 1:numel(corex_loc)
    if corex_loc(n)~=0 && corey_loc(n)~=0 && isnan(best_new(corex_loc(n),corey_loc(n)))
        best_new(corex_loc(n),corey_loc(n)) = best_new_interp(corex_loc(n),corey_loc(n));
        oceantracks(corex_loc(n),corey_loc(n)) = oceantracks_interp(corex_loc(n),corey_loc(n));
    end
end

%Filter to range of longitudes
selected = zeros(size(hb)); %map of points within the selected longitudes
for i = 1:size(hb,1)
    for j = 1:size(hb,2)
        if alond(i,j)>Ea || alond(i,j)<=We
            selected(i,j) = 1;
        end
    end
end

corelat62 = corelat(corelon>62);
corelon62 = corelat(corelon>62);
[coreX62,coreY62] = polarstereo_fwd(corelon62,corelat62);

if plot_figs == 1

    figure(fig_no);
    fig_no = fig_no + 1;
    ax(1) = subplot(2,1,1);
    mappingplotter(alatd,alond,best_new,-20,5,viridis);
    colormap(ax(1),viridis);
    hold on
    freezeColors
    contourm(alatd,alond,double(flipud(imrotate(maskwater,90))),'LineColor','k');
    title('Best match \epsilon_{Nd}')
    
    divmap = makeColorMap([1 0 0],[1 1 1],[0 0 1], 256); %Seems to not read in from main function, so remake this.,
    
    ax(2) = subplot(2,1,2);
    mappingplotter(alatd,alond,coretops-best_new,-5,5,divmap);
    colormap(ax(2),divmap);
    hold on
    freezeColors
    contourm(alatd,alond,double(flipud(imrotate(maskwater,90))),'LineColor','k');
    title('Difference from Measured')

    if PC ~= 1
        set(gcf,'Position',[50 50 600 1200]);
    end   

    print(gcf,fig11a,'-dpdf','-fillpage')
    close   
    
    %Apply mask to exclude points away from core tops or out of longitude
    %range.
    best_new_mask = selected.*best_new; %_mask
    coretops_mask = selected.*coretops; %_mask
    best_new_mask = nrst_smpl.*best_new_mask;
    coretops_mask = nrst_smpl.*coretops_mask;
    %set unwanted areas to NaN.
    best_new_mask(best_new_mask==0) = NaN;
    coretops_mask(coretops_mask==0) = NaN;
    
    figure(fig_no);
    fig_no = fig_no + 1;
    ax(1) = subplot(2,1,1); 
    mappingplotter(alatd,alond,best_new_mask,-20,5,viridis);
    colormap(ax(1),viridis);
    hold on
    freezeColors
    scatter(coreX62./600000,coreY62./600000,3,'filled','black')
    contourm(alatd,alond,double(flipud(imrotate(maskwater,90))),'LineColor','k') 
    title('Best match \epsilon_{Nd}')
    
    ax(2) = subplot(2,1,2);
    mappingplotter(alatd,alond,coretops_mask-best_new_mask,-5,5,divmap);
    colormap(ax(2),divmap);
    hold on
    freezeColors
    contourm(alatd,alond,double(flipud(imrotate(maskwater,90))),'LineColor','k');
    scatter(coreX62./60000,coreY62./60000,3,'filled','black')
    if PC ~= 1
        set(gcf,'Position',[50 50 600 1200]);
    end 
    title('Difference from Measured')

    print(gcf,fig11b,'-dpdf','-fillpage')
    close 
end

%% Do some analysis of depth vs transport mechanism.

dgroups = -5000:grouping:0; %Depth groups
threemeth_orig = threemeth;

%select only correct longitude range
threemeth = selected.*threemeth;
hb = selected.*hb;

%Apply filter for vicinity to core top data
threemeth = nrst_smpl.*threemeth;
hb = nrst_smpl.*hb;

%set unwanted areas to NaN.
hb(hb==0) = NaN; 
threemeth(threemeth==0) = NaN;

if plot_figs
    figure(fig_no)
    fig_no = fig_no + 1;
    subplot(2,1,1)
    mappingplotter(alatd,alond,threemeth_orig,0,3,viridis);
    colormap(viridis);
    hold on
    freezeColors
    contourm(alatd,alond,double(flipud(imrotate(maskwater,90))),'LineColor','k') 
    title('1=IRD, 2=Bottom Currents, 3=Gravity Flow')
    
    subplot(2,1,2)
    mappingplotter(alatd,alond,threemeth,0,3,viridis);
    colormap(viridis);
    hold on
    freezeColors
    contourm(alatd,alond,double(flipud(imrotate(maskwater,90))),'LineColor','k') 
    title('Nearest Method Match Subset')
    set(gcf,'Position',[50 50 600 1200]);
    
    print(gcf,fig12,'-dpdf','-fillpage')
    close
end

%Calculate scatter relationships to visualise model accuracy
coredata = coresites(:,3); %read in coretop eNd values
coredata(isnan(coredata)) = 0; %set NaN in header row to 0 (if present).

corex_loc_full = corex_loc;
corey_loc_full = corey_loc;
coredata_full = coredata;
corex_loc_v = NaN(size(corex_loc));
corey_loc_v = NaN(size(corex_loc));
coredata_v = NaN(size(corex_loc));

for n = 1:numel(corex_loc)
    if corex_loc(n)~=0 && corey_loc(n)~=0 && ...
            isnan(volc_msk(corex_loc(n),corey_loc(n)))
        corex_loc(n) = NaN;
        corey_loc(n) = NaN;
        coredata(n) = NaN;
        corex_loc_v(n) = corex_loc_full(n);
        corey_loc_v(n) = corey_loc_full(n);
        coredata_v(n) = coredata_full(n);
    end
end
corex_loc(isnan(corex_loc)) = [];
corey_loc(isnan(corey_loc)) = [];
coredata(isnan(coredata)) = [];
corex_loc_v(isnan(corex_loc_v)) = [];
corey_loc_v(isnan(corey_loc_v)) = [];
coredata_v(isnan(coredata_v)) = [];

grav_eNdvals = NaN(size(corex_loc));
surf_eNdvals = grav_eNdvals;
bot_eNdvals = grav_eNdvals;
best_eNdvals = grav_eNdvals;

%For areas away from radius3 km of active volcanoes.
for n = 1:numel(corex_loc)
    if corex_loc(n) ~= 0 %as values are 0 if site is outside of model domain.
        %Restrict to analysis longitudes
        if alond(corex_loc(n),corey_loc(n))>Ea || alond(corex_loc(n),corey_loc(n))<=We
            grav_eNdvals(n) = oceantracks3(corex_loc(n),corey_loc(n));
            surf_eNdvals(n) = oceantracks(corex_loc(n),corey_loc(n));
            bot_eNdvals(n) = oceantracks2(corex_loc(n),corey_loc(n));
            best_eNdvals(n) = best_new(corex_loc(n),corey_loc(n));
            best_meth(n) = threemeth(corex_loc(n),corey_loc(n));
            best_hb(n) = hb(corex_loc(n),corey_loc(n));
        end
    end
end
best_meth(isnan(best_meth)) = 4;
best_meth(isnan(best_eNdvals)) = NaN;

%Linear regression
size(coredata)
size(surf_eNdvals)
surf_mdl = fitlm(coredata,surf_eNdvals);
bot_mdl = fitlm(coredata,bot_eNdvals);
grav_mdl = fitlm(coredata,grav_eNdvals);
best_mdl = fitlm(coredata,best_eNdvals);
%Extract R^2
r2_surf = surf_mdl.Rsquared.Ordinary
r2_bot = bot_mdl.Rsquared.Ordinary
r2_grav = grav_mdl.Rsquared.Ordinary
r2_best = best_mdl.Rsquared.Ordinary

%Remove NaNs. Transpose model data.
coredata_surf = coredata(~isnan(surf_eNdvals));
coredata_bot = coredata(~isnan(bot_eNdvals));
coredata_grav = coredata(~isnan(grav_eNdvals));
if tuning ~= 1
    coredata_best = coredata(~isnan(best_eNdvals))
    best_eNdvals = best_eNdvals(~isnan(best_eNdvals))'
else
    coredata_best = coredata(~isnan(best_eNdvals));
    best_eNdvals = best_eNdvals(~isnan(best_eNdvals))';
end
surf_eNdvals = surf_eNdvals(~isnan(surf_eNdvals))';
bot_eNdvals = bot_eNdvals(~isnan(bot_eNdvals))';
grav_eNdvals = grav_eNdvals(~isnan(grav_eNdvals))';
best_hb = best_hb(~isnan(best_eNdvals))';
best_meth = best_meth(~isnan(best_meth))';

%Calculate root mean square error.
rmse_surf = sqrt(mean((surf_eNdvals(~isnan(surf_eNdvals))-...
    coredata_surf(~isnan(coredata_surf))).^2,'all'))
rmse_bot = sqrt(mean((bot_eNdvals(~isnan(bot_eNdvals))-...
    coredata_bot(~isnan(coredata_bot))).^2,'all'))
rmse_grav = sqrt(mean((grav_eNdvals(~isnan(grav_eNdvals))-...
    coredata_grav(~isnan(coredata_grav))).^2,'all'))
rmse_best = sqrt(mean((best_eNdvals(~isnan(best_eNdvals))-...
    coredata_best(~isnan(coredata_best))).^2,'all'))

grav_eNdvals_v = NaN(size(corex_loc_v));
surf_eNdvals_v = grav_eNdvals_v;
bot_eNdvals_v = grav_eNdvals_v;
best_eNdvals_v = grav_eNdvals_v;

%For areas within radius3 km of active volcanoes.
for n = 1:numel(corex_loc_v)
    if corex_loc_v(n) ~= 0 %as values are 0 if site is outside of model domain.
        %Restrict to analysis longitudes
        if alond(corex_loc_v(n),corey_loc_v(n))>Ea ||...
                alond(corex_loc_v(n),corey_loc_v(n))<=We
            grav_eNdvals_v(n) = oceantracks3(corex_loc_v(n),corey_loc_v(n));
            surf_eNdvals_v(n) = oceantracks(corex_loc_v(n),corey_loc_v(n));
            bot_eNdvals_v(n) = oceantracks2(corex_loc_v(n),corey_loc_v(n));
            best_eNdvals_v(n) = best_new(corex_loc_v(n),corey_loc_v(n));
            best_meth_v(n) = threemeth(corex_loc_v(n),corey_loc_v(n));
        end
    end
end
if exist('best_meth_v','var')
    best_meth_v(isnan(best_meth_v)) = 4;
    best_meth_v(isnan(best_eNdvals_v)) = NaN;

    best_meth_v(best_eNdvals_v==0) = NaN;
    
    best_meth_v = best_meth_v(~isnan(best_meth_v))';
end

surf_eNdvals_v(surf_eNdvals_v==0) = NaN;
bot_eNdvals_v(bot_eNdvals_v==0) = NaN;
grav_eNdvals_v(grav_eNdvals_v==0) = NaN;
best_eNdvals_v(best_eNdvals_v==0) = NaN;

%Remove NaNs. Transpose model data.
coredata_surf_v = coredata_v(~isnan(surf_eNdvals_v));
coredata_bot_v = coredata_v(~isnan(bot_eNdvals_v));
coredata_grav_v = coredata_v(~isnan(grav_eNdvals_v));
coredata_best_v = coredata_v(~isnan(best_eNdvals_v));
surf_eNdvals_v = surf_eNdvals_v(~isnan(surf_eNdvals_v))';
bot_eNdvals_v = bot_eNdvals_v(~isnan(bot_eNdvals_v))';
grav_eNdvals_v = grav_eNdvals_v(~isnan(grav_eNdvals_v))';
best_eNdvals_v = best_eNdvals_v(~isnan(best_eNdvals_v))';

%Difference vectors
dsurf = coredata_surf - surf_eNdvals;
dbot = coredata_bot - bot_eNdvals;
dgrav = coredata_grav - grav_eNdvals;
dbest = coredata_best - best_eNdvals;

% Put into column vectors.
depths = hb(:);
coef = threemeth(:);
% Corresponding depths of locations of points where each transport method is best match.
s1d = depths(coef==1);
s2d = depths(coef==2);
s3d = depths(coef==3);
% Round depths to intervals 'grouping'
s1d10 = grouping*round(s1d/grouping);
s2d10 = grouping*round(s2d/grouping);
s3d10 = grouping*round(s3d/grouping);
dpth10 = grouping*round(depths/grouping);
% Loop over depth intervals, counting number of points in each one.
n = 0;
for d = dgroups
    n = n+1;
    d1(n) = numel(find(s1d10==d));
    d2(n) = numel(find(s2d10==d));
    d3(n) = numel(find(s3d10==d));
    dpth(n) = numel(find(dpth10==d));
end

if plot_figs == 1
    %Plot scatter plots of actual vs modelled eNd
    figure(fig_no)
    fig_no = fig_no + 1;
    subplot(2,2,1)
    scatter(coredata_surf,surf_eNdvals,'x','k')
    hold on
    scatter(coredata_surf_v,surf_eNdvals_v,'x','r')
    hold on
    plot([-30 30],[-30 30])
    xlim([-21 6])
    ylim([-21 6])
    title(['R^2 = ' num2str(r2_surf)])
    xlabel('Core top values')
    ylabel('Modelled surface current values')
    subplot(2,2,2)
    scatter(coredata_bot,bot_eNdvals,'x','k')
    hold on
    scatter(coredata_bot_v,bot_eNdvals_v,'x','r')
    hold on
    plot([-30 30],[-30 30])
    xlim([-21 6])
    ylim([-21 6])
    title(['R^2 = ' num2str(r2_bot)])
    xlabel('Core top values')
    ylabel('Modelled bottom current values')
    subplot(2,2,3)
    scatter(coredata_grav,grav_eNdvals,'x','k')
    hold on
    scatter(coredata_grav_v,grav_eNdvals_v,'x','r')
    hold on
    plot([-30 30],[-30 30])
    xlim([-21 6])
    ylim([-21 6])
    title(['R^2 = ' num2str(r2_grav)])
    xlabel('Core top values')
    ylabel('Modelled gravity flow values')
    subplot(2,2,4)
    scatter(coredata_best,best_eNdvals,'x','k')
    hold on
    scatter(coredata_best_v,best_eNdvals_v,'x','r')
    colormap(viridis);
    hold on
    plot([-30 30],[-30 30])
    xlim([-21 6])
    ylim([-21 6])
    title(['R^2 = ' num2str(r2_best)])
    xlabel('Core top values')
    ylabel('Closest Mathching Method')

    print(gcf,fig14,'-dpdf','-bestfit')
    close
    
    %Make offset bar plot
    edges = -10:10;
    
    figure(fig_no)
    fig_no = fig_no + 1;
    subplot(2,2,1)
    histogram(dsurf,edges,'FaceColor','k')
    ylim([0 35])
    xlabel('Actual - Predicted \epsilon_{Nd}')
    ylabel('Number of sites')
    subplot(2,2,2)
    histogram(dbot,edges,'FaceColor','k')
    ylim([0 35])
    xlabel('Actual - Predicted \epsilon_{Nd}')
    ylabel('Number of sites')
    subplot(2,2,3)
    histogram(dgrav,edges,'FaceColor','k')
    ylim([0 35])
    xlabel('Actual - Predicted \epsilon_{Nd}')
    ylabel('Number of sites')
    subplot(2,2,4)
    histogram(dbest,edges,'FaceColor','k')
    ylim([0 35])
    xlabel('Actual - Predicted \epsilon_{Nd}')
    ylabel('Number of sites')
    
    print(gcf,fig15,'-dpdf','-bestfit')
    close
 
    %Get % of values within certain distance from core tops.
    disp('Percentages with certain eNd thresholds')
    within_1 = round(100*sum(dbest<1&dbest>-1)/numel(dbest))
    within_2 = round(100*sum(dbest<2&dbest>-2)/numel(dbest))
    within_3 = round(100*sum(dbest<3&dbest>-3)/numel(dbest))
    within_5 = round(100*sum(dbest<5&dbest>-5)/numel(dbest))
    
    figure(fig_no)
    fig_no = fig_no + 1;
    scatter(dbest,best_hb,'k','filled')
    xlabel('\epsilon_{Nd} Deviation')
    ylabel('Water depth (m)')
    print(gcf,fig13,'-dpdf','-bestfit')
    close
    
    % Plot depth analyses.
    figure(fig_no)
    fig_no = fig_no + 1;
    % Ocean depth profile for comparison
    subplot(3,2,1)
    [curve,~,~] = fit(dgroups',100*dpth','smoothingspline','SmoothingParam',0.005);
    plot(curve,'k',dgroups,100*dpth,'w','LineSpec2','w')
    legend('off')
    title('Ocean depth distribution')
    ylabel('Area (km^3)')
    xlabel('Water Depth (m)')
    ylim([0 50000]);
    xlim([-5000 0]);
    % Tranport method vs depth
    subplot(3,2,2)
    hold on
    [curve,~,~] = fit(dgroups',100.*d1','smoothingspline','SmoothingParam',0.005);
    plot(curve,'b',dgroups,100.*d1,'w','LineSpec2','w')
    [curve,~,~] = fit(dgroups',100.*d2','smoothingspline','SmoothingParam',0.005);
    plot(curve,'g',dgroups,100.*d2,'w','LineSpec2','w')
    [curve,~,~] = fit(dgroups',100.*d3','smoothingspline','SmoothingParam',0.005);
    plot(curve,'m',dgroups,100.*d3,'w','LineSpec2','w')
    legend({'Iceberg Rafting','Bottom Currents','Gravity Flows'},'Location','northwest');
    ylim([0 30000]);
    xlim([-5000 0]);
    title('Best transport method match compared to depth')
    ylabel('Area (km^3)')
    xlabel('Water Depth (m)')
    % Normalised i.e. fraction of points for each method at each depth interval
    subplot(3,2,3)
    hold on
    [curve,~,~] = fit(dgroups',(d1./sum(d1))','smoothingspline','SmoothingParam',0.005);
    plot(curve,'b',dgroups,d1./sum(d1),'w','LineSpec2','w')
    [curve,~,~] = fit(dgroups',(d2./sum(d2))','smoothingspline','SmoothingParam',0.005);
    plot(curve,'g',dgroups,d2./sum(d2),'w','LineSpec2','w')
    [curve,~,~] = fit(dgroups',(d3./sum(d3))','smoothingspline','SmoothingParam',0.005);
    plot(curve,'m',dgroups,d3./sum(d3),'w','LineSpec2','w')
    legend({'Iceberg Rafting','Bottom Currents','Gravity Flows'},'Location','northwest');
    title('Normalised')
    ylabel('Fraction')
    xlabel('Water Depth (m)')
    ylim([0 0.03]);
    xlim([-5000 0]);
    %Repeat for just <1000 m
    subplot(3,2,4)
    [curve,~,~] = fit(dgroups',100*dpth','smoothingspline','SmoothingParam',0.005);
    plot(curve,'k',dgroups,100*dpth,'w','LineSpec2','w')
    legend('off')
    title('Ocean depth distribution')
    ylabel('Area (km^3)')
    xlabel('Water Depth (m)')
    ylim([0 50000]);
    xlim([-1000 0]);
    % Tranport method vs depth
    subplot(3,2,5)
    hold on
    [curve,~,~] = fit(dgroups',100.*d1','smoothingspline','SmoothingParam',0.005);
    plot(curve,'b',dgroups,100.*d1,'w','LineSpec2','w')
    [curve,~,~] = fit(dgroups',100.*d2','smoothingspline','SmoothingParam',0.005);
    plot(curve,'g',dgroups,100.*d2,'w','LineSpec2','w')
    [curve,~,~] = fit(dgroups',100.*d3','smoothingspline','SmoothingParam',0.005);
    plot(curve,'m',dgroups,100.*d3,'w','LineSpec2','w')
    legend({'Iceberg Rafting','Bottom Currents','Gravity Flows'},'Location','northwest');
    ylim([0 30000]);
    xlim([-1000 0]);
    title('Best transport method match compared to depth')
    ylabel('Area (km^3)')
    xlabel('Water Depth (m)')
    % Normalised i.e. fraction of points for each method at each depth interval
    subplot(3,2,6)
    hold on
    [curve,~,~] = fit(dgroups',(d1./sum(d1))','smoothingspline','SmoothingParam',0.005);
    plot(curve,'b',dgroups,d1./sum(d1),'w','LineSpec2','w')
    [curve,~,~] = fit(dgroups',(d2./sum(d2))','smoothingspline','SmoothingParam',0.005);
    plot(curve,'g',dgroups,d2./sum(d2),'w','LineSpec2','w')
    [curve,~,~] = fit(dgroups',(d3./sum(d3))','smoothingspline','SmoothingParam',0.005);
    plot(curve,'m',dgroups,d3./sum(d3),'w','LineSpec2','w')
    legend({'Iceberg Rafting','Bottom Currents','Gravity Flows'},'Location','northwest');
    title('Normalised')
    ylabel('Fraction')
    xlabel('Water Depth (m)')
    ylim([0 0.03]);
    xlim([-1000 0]);

    print(gcf,fig13b,'-dpdf','-bestfit')
    close
end   
