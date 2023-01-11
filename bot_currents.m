function [oceantracks2,fig_no] = bot_currents(sx,sy,bot_months,plot_figs,filename2,...
    filename3,ys,xs,extent,extent_thresh,erosionv,alatd,alond,kappa,hb...
    ,z0,rho,t0,c0,vs,susph,oceantracks,repeat,iter,runname,rep,xcoords,...
    ycoords,radius2,resolution,vbot,ubot,smth,nrst_smpl_nans,coretops,...
    divmap,maskwater,fig6,indi,x22,y22,model,UVEL,VVEL)
% Estimate debris movement due to bottom current transport where velocity 
% exceeds a threshold. This function reads in several variables from prov.m
% and outputs an eNd value estimate and figure number.

disp('Estimating transport using bottom currents...')

oceantracks2 = zeros(sx,sy); %Create blank array and set to NaN.
oceantracks2(oceantracks2==0) = NaN;

%Calculate some bits before entering loop.
lat2 = double(ncread(filename2,char('latitude')));
long2 = double(ncread(filename2,char('longitude')));
[lat2,long2] = meshgrid(lat2,long2); %Gets lat lon in grid.
[X2,Y2] = polarstereo_fwd(lat2,long2); % Convert ocean reanalysis coordinates to polar stereo

for month = 1:bot_months
    disp(['Month ' num2str(month)])
    if plot_figs == 1
        fig_no = 5;
    end

    count = [Inf Inf Inf 1];
    if bot_months>12 && month<=12 %Read 2018 data for first 12 months if more than 12 in total
        startLoc = [1 1 1 month];
        V2 = double(ncread(filename3,char('vo_oras'),startLoc,count));
        U2 = double(ncread(filename3,char('uo_oras'),startLoc,count));
    elseif bot_months>12 && month>12 %Read 2019 data after 12 months
        startLoc = [1 1 1 month-12];
        V2 = double(ncread(filename2,char('vo_oras'),startLoc,count));    
        U2 = double(ncread(filename2,char('uo_oras'),startLoc,count));
    else %Read 2019 data only if less than 12 months in total.
        startLoc = [1 1 1 month];
        V2 = double(ncread(filename2,char('vo_oras'),startLoc,count));
        U2 = double(ncread(filename2,char('uo_oras'),startLoc,count));
    end

    %Create blank velocity arrays
    [sX,sY,~] = size(V2);
    U_bot = zeros(sX,sY);
    V_bot = U_bot;

    for i = 1:sX
        for j = 1:sY
            if isempty(max(find(~isnan(U2(i,j,:)))))
                U_bot(i,j) = NaN;
                V_bot(i,j) = NaN;
            else
                U_bot(i,j) = U2(i,j,max(find(~isnan(U2(i,j,:)))));
                V_bot(i,j) = V2(i,j,max(find(~isnan(V2(i,j,:)))));
            end
        end
    end    

    [vt,ut] = uv2vxvy(lat2,long2,U_bot,V_bot); %Convert velocities to polar stereo 

    Fa2 = scatteredInterpolant(X2(:),Y2(:),vt(:)); 
    Fb2 = scatteredInterpolant(X2(:),Y2(:),ut(:));  
    botU = Fa2(ys,xs);
    botV = Fb2(ys,xs);

    if extent <= extent_thresh  %If lots of retreat (i.e. collapsed WAIS)

        ULAT = ncread('Pliomip2_ULATLONG_extracted.nc','ULAT');
        [~,y] = find(ULAT<=-60);
        topl = max(y);
        ULONG = ncread('Pliomip2_ULATLONG_extracted.nc','ULONG');
        %Get only correct latitude data
        ULONG = ULONG(:,1:topl);
        ULONG(ULONG>180) = ULONG(ULONG>180)-360;
        ULAT = ULAT(:,1:topl);

        UVEL_temp = UVEL;
        VVEL_temp = VVEL;
        VVEL_temp = VVEL_temp(:,:,:,month);
        UVEL_temp = UVEL_temp(:,:,:,month);
        %Create blank velocity arrays
        [sX_PLIO,sY_PLIO,~] = size(VVEL_temp);
        U_botPLIO = zeros(sX_PLIO,sY_PLIO);
        V_botPLIO = U_botPLIO;

        for i = 1:sX_PLIO
            for j = 1:sY_PLIO
                if isempty(max(find(~isnan(UVEL_temp(i,j,:)))))
                    U_botPLIO(i,j) = NaN;
                    V_botPLIO(i,j) = NaN;
                else
                    U_botPLIO(i,j) = UVEL_temp(i,j,max(find(~isnan(UVEL_temp(i,j,:)))));
                    V_botPLIO(i,j) = VVEL_temp(i,j,max(find(~isnan(VVEL_temp(i,j,:)))));
                end
            end
        end    
        %Convert to polar stereo
        [PLIOX,PLIOY] = polarstereo_fwd(ULAT,ULONG); 
        %Convert velocities to polar stereo 
        [vtPLIO,utPLIO] = uv2vxvy(PLIOX,PLIOY,U_botPLIO,V_botPLIO);

        Fa2PLIO = scatteredInterpolant(PLIOX(:),PLIOY(:),vtPLIO(:)); 
        Fb2PLIO = scatteredInterpolant(PLIOX(:),PLIOY(:),utPLIO(:));  
        botUPLIO = Fa2PLIO(ys,xs);
        botVPLIO = Fb2PLIO(ys,xs); 

        % Flag the margin between the datasets to blend later.
        moderndata = ~isnan(botV);
        mrg = 3;
        flag = zeros(sx,sy);
        for i = mrg+1:sx-mrg
            for j = mrg+1:sy-mrg
            subset = moderndata(i-mrg:i+mrg,j-mrg:j+mrg);
                if sum(subset,'all') ~= (1+2*mrg)^2 && sum(subset,'all') ~= 0
                    flag(i,j) = 1; 
                end
            end
        end

        %Fill holes in modern data
        botV(isnan(botV)) = botVPLIO(isnan(botV))/100;
        botU(isnan(botU)) = botUPLIO(isnan(botU))/100;
    end

    % Interpolate gaps
    botUB = botU;
    botVB = botV;
    if extent <= extent_thresh
        botUB(flag==1) = NaN; %Fill overlap 
        botVB(flag==1) = NaN;
    end
    botVB = inpaint_nans(botVB); 
    botUB = inpaint_nans(botUB); 
    botUB(indi) = NaN;
    botVB(indi) = NaN;
    botV = botVB;
    botU = botUB;

    % Get start locations where velocity above erosion threshold.
    botmag = sqrt(botV.^2+botU.^2);
    er_dep = zeros(sx,sy);
    er_dep(botmag>(erosionv/2.5)) = 1;
    er_dep(alatd>-60) = 0; % Exclude areas out of reanalysis (above 60oS).
    [row,col] = find(er_dep == 1);

    ocean_eNd2 = zeros(size(botU)); % Records total eNd.
    ocean_count2 = zeros(size(botU)); % Records number of tracks through cell.

    % Using ice end points as ocean start points
    if plot_figs == 1
        figure(fig_no)
        fig_no = fig_no + 1;
        hold on
        contour(1:sx,1:sy,maskwater,1,'k')
        axis equal
        axis off
    end
    
    if model == 1
          h22 = streamline(x22,y22,botU',botV',row,col,[.25,5000]);
    elseif model == 2
          h22 = streamline(x22,y22,botU,botV,row,col,[.25,5000]);
    else
        disp('Ice sheet model number (''model'') invalid!')
        return
    end

    for k = 1:length(h22) % for each stream line.
        xk = get(h22(k),'xdata'); % get its x data
        yk = get(h22(k),'ydata'); % get its y dat            
        xk = xk(~isnan(xk)); %remove NaNs
        yk = yk(~isnan(yk));
        if model == 1
            xk = round(xk); %round to exact cells
            yk = round(yk);
        else
            xk_temp = xk;
            xk = round(yk); %round to exact cells
            yk = round(xk_temp);
        end
        vmags = zeros(length(xk),1);
        for n = 1:length(xk)
            vmags(n) = botmag(xk(n),yk(n));
        end
        %loops adding up eNd values and number of tracks through each ocean
        %cell.
        if xk(1)~= xk(end) %if not a line that just stops
            sed_conc = 1; % Concentration of sediment
            distance = 0; %distance sediment has travelled.
            current_age = 0; %cumulative age of the current flowpath               
            p = 1;
            %loop over each point in streamline.
            while p<length(xk) && sed_conc>0
                p = p+1;
                if ~isnan(xk(p)) && ~isnan(yk(p)) %if not NaN
                    if model == 2 %swap back over
                        Cd = real(kappa^2/(log(-hb(yk(p),xk(p))/z0)-1)^2); %Drag coefficient   
                    else
                        Cd = real(kappa^2/(log(-hb(xk(p),yk(p))/z0)-1)^2); %Drag coefficient          
                    end
                    tw = rho*Cd*botmag(xk(p),yk(p))^2; % Shear stress (N/m^2) 
                    prob = max(1-tw/t0,0); %Probability of a particle sticking to the bed
                    %get new distance travelled
                    distance = distance + resolution*1000*...
                        sqrt(abs(xk(p-1)-xk(p)).^2+abs(yk(p-1)-yk(p)).^2);
                    elapsed2 = distance./botmag(xk(p),yk(p)); %time elapsed at this grid cell (days)
                    current_age = elapsed2 + current_age; %cumulative age of flowline (seconds)
                    sed_conc_old = sed_conc;
                    % Calculate new sediment concentration. Velocity 
                    %changes can give negative deposition, so set to 
                    %previous if larger.
                    sed_conc = min(c0*exp(-prob*current_age*vs/susph),sed_conc_old);
                    sed_dropped = sed_conc_old - sed_conc;
                    if month == 1 && repeat == 1
                        % If first bottom current estimate, use surface
                        % current estimate.
                        if model == 2
                        ocean_eNd2(xk(p),yk(p)) = ocean_eNd2(xk(p),yk(p))...
                            + oceantracks(yk(1),xk(1))*sed_dropped;
                        else
                        ocean_eNd2(xk(p),yk(p)) = ocean_eNd2(xk(p),yk(p))...
                            + oceantracks(xk(1),yk(1))*sed_dropped;
                        end
                    elseif month == 1 && repeat > 1
                        % If first bottom current estimate after at
                        % least 1 gravity flow iteration, use gravity 
                        % method beyond shelf.
                        ocean_eNd2(xk(p),yk(p)) = ocean_eNd2(xk(p),yk(p))...
                            + iter(xk(1),yk(1))*sed_dropped;
                    else
                        % Otherwise, use previous month output filled
                        % with surface current estimate where NaN.
                        ocean_eNd2(xk(p),yk(p)) = ocean_eNd2(xk(p),yk(p))...
                            + oceantracks2_filled(xk(1),yk(1)).*sed_dropped;
                    end
                    ocean_count2(xk(p),yk(p)) = ocean_count2(xk(p),yk(p)) + sed_dropped; 
                end
            end
        end
    end  

    oceantracks2_old = oceantracks2; % Save old output
    oceantracks2 = ocean_eNd2./ocean_count2; % Calculate new map
    % Where NaNs, set to previous month.
    oceantracks2(isnan(oceantracks2)) = oceantracks2_old(isnan(oceantracks2));
    % When no previous bottom current data, use surface current estimate.
    oceantracks2_filled = oceantracks2;
    oceantracks2_filled(isnan(oceantracks2_filled)) = ...
        oceantracks(isnan(oceantracks2_filled));

    if plot_figs == 1 && (month/6 == round(month/6))
        fig5 = ['Bottom_flowpaths_' runname num2str(month) '.png'];
        print(gcf,fig5,'-dpng','-r800')
        close 
    end  
end  

if repeat == rep %Only interpolate bottom flow paths if last bottom/gravity iteration.
    %Select only points within 'radius2' km of data
    [x, y] = find(~isnan(oceantracks2));
    dist = 10^6; % Vector containing distance to each coretop
    botcur_mask = zeros(sx,sy); %map of points within radius
    % Loop over each grid cell, and for each one loop over each point
    % filled by a streamline and get the distance from it.
    for i = 1:sx
        for j = 1:sy
            for k = 1:size(x,1)
                if x(k)~=0 && y(k)~=0
                    dist = resolution*sqrt((xcoords(x(k))-xcoords(i)).^2+(ycoords(y(k))-ycoords(j)).^2);
                    if dist<(resolution*radius2)
                        botcur_mask(i,j) = 1; % Set to 1 if close to a bottom current path.
                        break % Exit loop and move to next cell.
                    end
                end
            end
        end
    end
    botcur_mask(botcur_mask==0) = NaN; %Set zeros to NaNs
end

if model == 2
    oceantracks2 = oceantracks2';
end

%fill gaps between ocean streamlines, then set ice (grounded/no velocity) to NaN
if repeat == rep
    oceantracks2 = inpaint_nans(oceantracks2,4);
    oceantracks2 = botcur_mask.*oceantracks2;
end
oceantracks2(vbot'~=0) = NaN;
oceantracks2(ubot'~=0) = NaN;
oceantracks2(maskwater'==0) = NaN;

if (repeat == rep) && (smth > 0) %Option to apply a moving average smoothing.
    oceantracks2 = movmean_2d(oceantracks2,smth);
end

if plot_figs == 1
    figure(fig_no)
    fig_no = fig_no + 1;
    ax(1) = subplot(2,1,1);
    mappingplotter(alatd,alond,oceantracks2,-20,5,viridis);
    colormap(ax(1),viridis)
    hold on
    freezeColors
    contourm(alatd,alond,double(flipud(imrotate(maskwater,90))),'LineColor','k') 
    title('Bottom Current \epsilon_{Nd} map')

    %Get model-data difference and limit to 'radius' of samples.
    bot_diff = coretops-oceantracks2;
    bot_diff = bot_diff.*nrst_smpl_nans;

    ax(2) = subplot(2,1,2);
    mappingplotter(alatd,alond,bot_diff,-5,5,divmap);
    colorbar
    colormap(ax(2),divmap)
    hold on
    freezeColors
    contourm(alatd,alond,double(flipud(imrotate(maskwater,90))),'LineColor','k') 
    title('Actual - Bottom Current \epsilon_{Nd}')
    set(gcf,'Position',[50 50 600 1200]);

    print(gcf,fig6,'-dpdf','-fillpage')
    close
end
