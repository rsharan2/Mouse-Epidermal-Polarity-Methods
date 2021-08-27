function hair_local_order_analysis(vector_sheet,hair_image,radius)

%vector_sheet: path to the xlsx sheet with the output from
%overlay_vector_field_on_image
%ex: 'Vector field rst P4 18_05_18 M1072 n1_updated.xlsx'

%hair_image: path to the image of the hairs upon which you want the vectors
%overlayed
%ex: '18_05_18 M1072 1 P4 rst rst invert sharpen scaled to control.tif'

%radius: radius to consider neighborhood for local order calculation
%ex: 60

%ex of running: 
%hair_local_order_analysis('Vector field rst P4 18_05_18 M1072 n1_updated.xlsx','18_05_18 M1072 1 P4 rst rst invert sharpen scaled to control.tif',60)

t = readtable(vector_sheet);
I = imread(hair_image);

x_offset = t.X(1) - .5;
y_offset = t.Y(1) - .5;

NN_all = rangesearch([t.X t.Y],[t.X t.Y],radius,'SortIndices',true);

%order_calc = @(x) sum((cosd(t(x,:).Updated_Orientation(2:end) - t(x,:).Updated_Orientation(1))+1)./2)/(length(x)-1);
%order_param = cellfun(order_calc,NN_all);

%perform pairwise subtraction for each follicle in neighborhood for each
%neighborhood
order_prep = @(x) t(x,:).Updated_Orientation - t(x,:).Updated_Orientation';
order_prep_out = cellfun(order_prep,NN_all,'UniformOutput',false);

%perform rest of calcuations on the upper trianglular matrix excluding diagonal
order_calc = @(A) mean((cosd(A(triu(true(size(A)),1)))+1)./2);
order_param = cellfun(order_calc,order_prep_out);

t.Local_Order = order_param;

[~,vs_fname,~] = fileparts(vector_sheet);
writetable(t,[vs_fname '_local_order_radius_' num2str(radius) '.xlsx']);

X = reshape(t.X,[length(unique(t.X)) length(unique(t.Y))]);
Y = reshape(t.Y,[length(unique(t.X)) length(unique(t.Y))]);
C = reshape(t.Local_Order,[length(unique(t.X)) length(unique(t.Y))]);

Y(end+1,:) = Y(end,:);
Y(:,end+1) = Y(:,end)+12;
X(end+1,:) = X(end,:) + 12;
X(:,end+1) = X(:,end);
C(end+1,:) = 0;
C(:,end+1) = 0;

%{
f = figure;
ax1 = axes; 
imshow(I,[]);
hold on;
s = pcolor(X'-x_offset,Y'-y_offset,C');
s.FaceColor = 'flat';
s.EdgeColor = 'none';
quiver(t.X,t.Y,t.Updated_DX/2,t.Updated_DY/2,.5,'k');
quiver(t.X,t.Y,-t.Updated_DX/2,-t.Updated_DY/2,.5,'k','ShowArrowHead','off');
ax1.Visible = 'off';
ax1.XTick = [];
ax1.YTick = [];
% Give each one its colormap
colormap(hot);
caxis([0 1]);
axis off;
%colorbar

[~,fname,~] = fileparts(hair_image);
exportgraphics(gcf,[fname '_local_order_radius_' num2str(radius) '_vec.tif'],'Resolution',500);
saveas(gcf,[fname '_local_order_radius_' num2str(radius) '_vec.fig']);
close all;
%}

f = figure;
ax1 = axes;
imshow(I,[]);
hold on;
s = pcolor(X'-x_offset,Y'-y_offsetC');
colormap(hot)
s.FaceColor = 'flat';
s.EdgeColor = 'none';
ax1.Visible = 'off';
ax1.XTick = [];
ax1.YTick = [];
% Give each one its colormap
colormap(hot);
caxis([0 1]);
axis off;
%colorbar

[~,fname,~] = fileparts(hair_image);
exportgraphics(gcf,[fname '_local_order_radius_' num2str(radius) '.tif'],'Resolution',500);
saveas(gcf,[fname '_local_order_radius_' num2str(radius) '.fig']);
close all;

end
