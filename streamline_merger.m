function streamline_merger(filein,t,m)

skipped = 0;
[sx,sy] = size(imread(['Surface_flowpaths_' filein '1.png']));
sy = sy/3; %For some reason sy is tripled.
heatmap = zeros(sx,sy);

for month = 1:m
  %  month
    if isfile(['Surface_flowpaths_' filein num2str(month) '.png'])
        data = imread(['Surface_flowpaths_' filein num2str(month) '.png']);
        if month == 1
            total = data;
        else
            %Need if to ignore all black images sometimes produced.
            if max(data(:)) == 0
                    disp('Skipped')
                    skipped = skipped + 1;
            else
                for n = 1:3
                        total(:,:,n) = min(total(:,:,n),data(:,:,n));
                       heatmap = heatmap + im2double(data(:,:,n));
                end
            end
        end
    end
end

heatmap = (((m-1)*3)- heatmap)/2;
heatmap(heatmap==0) = NaN;
heatmapcp = heatmap;

%Loop over cells near grounding line and set to black.
for i = 2:sx-1
    for j = 2:sy-1
        subset = heatmapcp(i-1:i+1,j-1:j+1);
        if max(max(subset)) == 1.5*(m-1)
            heatmap(i-1:i+1,j-1:j+1) = 1.5*(m-1)*ones(3,3);
        end
    end
end

xmin = min(find(~isnan(max(heatmap,[],1))));
xmax = max(find(~isnan(max(heatmap,[],1))));
ymax = max(find(~isnan(max(heatmap,[],2))));
ymin = min(find(~isnan(max(heatmap,[],2))));

heatmap = heatmap - 1.5*skipped;
heatmap = heatmap*(255/m); %Adjust so 1:256 scale for colourbar.
maximum = max(max(heatmap));
heatmap(heatmap<maximum/256) = ceil(maximum/256);
heatmap(heatmap == 0) = NaN;

heatmapsmall = heatmap(ymin:ymax,xmin:xmax);

figure(2)
imagesc(total)
imwrite(total,['Surface_flowpaths_merged_' filein '.png']);
disp(['Finished. ' num2str(skipped) ' months skipped.'])
close
