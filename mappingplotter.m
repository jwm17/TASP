function [k] = mappingplotter(lat,lon,in,cmin,cmax,colormap)

axesm('MapProjection','stereo','grid','on','frame','on','glinewidth'...
    ,0.25,'glinestyle','-','origin',[-90 0],'MapLatLimit',[-90 -62],...
    'MeridianLabel','off','MLabelParallel','north')
k = pcolorm(lat,lon,in);
set(k,'edgecolor','none');
% colormap(colormap)
colorbar; 
 if exist('cmin','var') == 1 && exist('cmax','var') == 1
     caxis([cmin cmax])
 end