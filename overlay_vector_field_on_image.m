function overlay_vector_field_on_image(vector_sheet,hair_image,bw_map)

%vector_sheet: path to the csv or xlsx sheet with the output from
%OrientationJ vector field.
%ex: 'Vector field rst P4 18_05_18 M1072 n1.xlsx'

%hair_image: path to the image of the hairs upon which you want the vectors
%overlayed
%ex: '18_05_18 M1072 1 P4 rst rst invert sharpen scaled to control.tif'

%bw_map: black and white image where white represents regions to flip
%vectors and black represents regions to not flip
%ex: 'test_map.tif'

%ex of running: 
%overlay_vector_field_on_image('Vector field rst P4 18_05_18 M1072 n1.xlsx','18_05_18 M1072 1 P4 rst rst invert sharpen scaled to control.tif','test_map.tif')

t = readtable(vector_sheet);
I = imread(hair_image);

if ndims(I) == 2
    I = ind2rgb(imadjust(I),gray(single(intmax('uint16'))));
end

map = imread(bw_map);

x_offset = t.X(1) - .5;
y_offset = t.Y(1) - .5;

map_i = zeros(size(map),'int8');
map_i(map>0) = -1;
map_i(map==0) = 1;

indic = @(x,y,dx) dx*double(map_i(y,x));
rev_dx = rowfun(indic,table(round(t.X-x_offset),round(t.Y-y_offset),t.DX),'OutputFormat','uniform');
rev_dy = rowfun(indic,table(round(t.X-x_offset),round(t.Y-y_offset),t.DY),'OutputFormat','uniform');

%{
figure; 
imshow(I,[])
hold on;
quiver(t.X,t.Y,rev_dx/2,rev_dy/2,.25,'y','LineWidth',2);
quiver(t.X,t.Y,-rev_dx/2,-rev_dy/2,.25,'y','LineWidth',2,'ShowArrowHead','off');
hold off;
axis off;

[~,fname,~] = fileparts(hair_image);
exportgraphics(gcf,[fname '_vec_overlay.tif']);
saveas(gcf,[fname '_vec_overlay.fig']);
close all;
%}

n = 720;        % two colors defined per angle for an even distribution
% offsets the hsv colormap to center 0° at cyan
cmap = hsv2rgb(horzcat(((0:n)./n)', repelem(1, n+1)', repelem(1, n+1)'));  % 0° = cyan

[~,fname,~] = fileparts(hair_image);
[~,vs_fname,~] = fileparts(vector_sheet);
updated_orientation = -(rad2deg(atan2(rev_dy,rev_dx))-90)+90; %0 = west, 90 = north, 180 = east, 270 = south
updated_orientation(updated_orientation<0) = updated_orientation(updated_orientation<0) + 360;
updated_t = table(t.X,t.Y,t.Slice,rev_dx,rev_dy,updated_orientation,t.Coherency,t.Energy,'VariableNames',{'X','Y','Slice','Updated_DX','Updated_DY','Updated_Orientation','Coherencey','Energy'});
writetable(updated_t,[vs_fname '_updated.xlsx']);

X = reshape(t.X,[length(unique(t.X)) length(unique(t.Y))]);
Y = reshape(t.Y,[length(unique(t.X)) length(unique(t.Y))]);
C = reshape(updated_t.Updated_Orientation,[length(unique(t.X)) length(unique(t.Y))]);

Y(end+1,:) = Y(end,:);
Y(:,end+1) = Y(:,end)+12;
X(end+1,:) = X(end,:) + 12;
X(:,end+1) = X(:,end);
C(end+1,:) = 0;
C(:,end+1) = 0;

%{
f = figure;
ax1 = axes;
imshow(I,[],'Parent',ax1); hold on;
axis off;
ax2 = copyobj(ax1,f);
axes(ax2);
s = pcolor(X'-x_offset,Y'-y_offset,C');
colormap(cmap)
s.FaceColor = 'flat';
s.EdgeColor = 'none';
s.FaceAlpha = .3;
quiver(t.X,t.Y,rev_dx/2,rev_dy/2,.25,'k','LineWidth',2);
quiver(t.X,t.Y,-rev_dx/2,-rev_dy/2,.25,'k','ShowArrowHead','off','LineWidth',2);
hlink = linkprop([ax1,ax2],{'XLim','YLim','ZLim','Position'});
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
% Give each one its colormap
colormap(ax1);
colormap(ax2,cmap);
caxis(ax2, [0 360]);
axis off;
%colorbar
%colorbar('Ticks', [0,90,180,270,360], 'TickLabels', [-180,-90,0,90,180]);

exportgraphics(gcf,[fname '_colored_vec_overlay.tif']);
saveas(gcf,[fname '_colored_vec_overlay.fig']);
close all;
%}

f = figure;
ax1 = axes;
imshow(I,[],'Parent',ax1); hold on;
axis off;
ax2 = copyobj(ax1,f);
axes(ax2);
s = pcolor(X-x_offset,Y-y_offset,C);
colormap(cmap)
s.FaceColor = 'flat';
s.EdgeColor = 'none';
s.FaceAlpha = .3;
hlink = linkprop([ax1,ax2],{'XLim','YLim','ZLim','Position'});
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
% Give each one its colormap
colormap(ax1);
colormap(ax2,cmap);
caxis(ax2, [0 360]);
axis off;
%colorbar
%colorbar('Ticks', [0,90,180,270,360], 'TickLabels', [-180,-90,0,90,180]);

exportgraphics(gcf,[fname '_colored_overlay.tif']);
saveas(gcf,[fname '_colored_overlay.fig']);
close all;

end