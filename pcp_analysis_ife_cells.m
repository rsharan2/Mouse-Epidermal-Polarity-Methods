%% *Analysis of polarity vectors*
%% *Loading images and data*
%%
%load images and data
top_dir = ['ceph' filesep '16bit celsr' filesep '16bit celsr-new'];

if ispc
    top_dir = '.';
end

imgs = dir([top_dir filesep '*.tif']);
imgs_names = {imgs.name};
disp(imgs_names)
area_and_perimeter = table;
pca_mag_and_orientation = table;

warn_id = 'MATLAB:table:ModifiedAndSavedVarnames';
warning('off',warn_id);

for i=1:length(imgs_names)
    img = imgs_names{i};
    [~,img_name,~] = fileparts(img);
    
    genotype = "";
    if contains(img,'fzd6','IgnoreCase',true)
        genotype = "fzd6KO";
    elseif contains(img,'rst','IgnoreCase',true)
        genotype = "rst";
    elseif contains(img,'b6','IgnoreCase',true)
        genotype = "wt";
    else
        disp([ 'Unknown genotype: "' img '"' newline char(9) 'Does not contain b6, fzd6, or rst'])       
    end
    results = dir([top_dir filesep img_name filesep 'Result_*']);
    ap_tbl = readtable([top_dir filesep img_name filesep 'Result_03' filesep 'Cell_Area_Perimeter.csv']);
    Genotype = array2table(repelem(genotype,height(ap_tbl),1), 'VariableNames', {'Genotype'});
    Image_name = array2table(repelem(convertCharsToStrings(img_name),height(ap_tbl),1),'VariableNames', {'Image_name'});
    ap_tbl = [ap_tbl Genotype Image_name];
    area_and_perimeter = [area_and_perimeter; ap_tbl];
    
    mo_tbl = readtable([top_dir filesep img_name filesep 'Result_03' filesep 'PCA_Cell-by-Cell_Polarity.csv']);
    if height(ap_tbl) ~= height(mo_tbl)
        disp(['Different number of rows between Area_Perimeter and Polarity CSVs: ' img])
    end
    Genotype = array2table(repelem(genotype,height(mo_tbl),1), 'VariableNames', {'Genotype'});
    Image_name = array2table(repelem(convertCharsToStrings(img_name),height(mo_tbl),1),'VariableNames', {'Image_name'});
    mo_tbl = [mo_tbl Genotype Image_name];
    
    pca_mag_and_orientation = [pca_mag_and_orientation; mo_tbl];
end
%save('stupid_tables.mat','area_and_perimeter','pca_mag_and_orientation');
combined_tbl = join(area_and_perimeter,pca_mag_and_orientation,'Keys',[1 5],'KeepOneCopy','Genotype');
areas = table2array(area_and_perimeter(:,2));
%area_threshold = mean(areas) + 2*std(areas);
area_threshold = median(areas) + 6*mad(areas);
%%
wt_tbl = combined_tbl(combined_tbl.Genotype == "wt",:);
fzko_tbl = combined_tbl(combined_tbl.Genotype == "fzd6KO",:);
rst_tbl = combined_tbl(combined_tbl.Genotype == "rst",:);
warning('on',warn_id);

%%
neigh_p_vec_all = [];
neigh_p_ang_all = [];
neigh_local_order_all = [];
Centroid_X_all = [];
Centroid_Y_all = [];
n = 720;        % two colors defined per angle for an even distribution
% offsets the hsv colormap to center 0° at cyan
cmap = hsv2rgb(horzcat([(0:round(n/2)-1)./n (round(n/2):n)./n ]', repelem(1, n+1)', repelem(1, n+1)'));  % 0° = cyan
% adds zeros to correspond to negative numbers (ie NaNs, to color black)
cmap = vertcat([0 0 0; 0 0 0], cmap);

neighborhood_radius = 40;

outputfolder = [top_dir filesep 'output' filesep 'r_' num2str(neighborhood_radius)];
%%
if ~exist(outputfolder, 'dir')
    mkdir(outputfolder)
end

for i=1:length(imgs_names)
    img = imgs_names{i};
    [~,img_name,~] = fileparts(img);
    results = dir('Result_*');
    id_im = imread([top_dir filesep img_name filesep 'Result_03' filesep 'Cell_Identity.tif']);
    labeled_by_pca = zeros(size(id_im));
    labeled_by_pca(id_im==0) = NaN; 
    opacity_by_pca = zeros(size(id_im));
    
    labeled_by_local_order = zeros(size(id_im));
    labeled_by_local_order(id_im==0) = NaN; 

    centroids = regionprops('table',id_im,'Centroid');
    c = centroids.Centroid;
    
    neighborhoods = rangesearch(c,c,neighborhood_radius);
    
    curr_tbl = combined_tbl(combined_tbl.Image_name == img_name,:);
    neigh_p_vec = zeros(1,height(curr_tbl));
    neigh_p_ang = zeros(1,height(curr_tbl));
    neigh_local_order = zeros(1,height(curr_tbl));

    parfor n=1:length(neighborhoods)
        %tic;
        curr_cells = curr_tbl(find(ismember(curr_tbl.CellIdentity,neighborhoods{n})),:);
        %filter cells too big
        curr_id = neighborhoods{n}(1);
        temp_ang = zeros(size(id_im));
        temp_vec = zeros(size(id_im));
        temp_local_order = zeros(size(id_im));
        p_vec = NaN;
        p_ang = NaN;
        n_local_order = NaN;
        if curr_cells(curr_cells.CellIdentity == curr_id,:).CellArea_pixel_2_ < area_threshold && any(curr_cells.CellArea_pixel_2_ < area_threshold) && any(curr_cells.PCAMagnitude > .02)
            %filter out cells that are unpolarized and too big
            curr_cells = curr_cells(curr_cells.CellArea_pixel_2_ < area_threshold & curr_cells.PCAMagnitude > .02,:);
            curr_mags = curr_cells.PCAMagnitude;
            curr_angles = curr_cells.PCAAngle___;
            
            %calculate all pairwise comparisons
            A = curr_cells.PCAAngle___-curr_cells.PCAAngle___'; 
            
            %max angle difference cannot be greater than 90, if so, take the supplementary angle of the difference
            %cos is symmetric around 0, so no need to worry about signs
            A(abs(A) > 90) = (180 - abs(A(abs(A) > 90)));
            
            n_local_order = mean(cosd(A(triu(true(size(A)),1))).^2);
            p_vec = mean(vecnorm([(curr_mags.*(cosd(2*curr_angles)))'; (curr_mags.*(sind(2*curr_angles)))']));
            p_ang = .5*atan2d(mean((curr_mags.*(sind(2*curr_angles)))),mean((curr_mags.*(cosd(2*curr_angles)))));
        end
        
        neigh_p_vec(n) = p_vec;
        neigh_p_ang(n) = p_ang;
        neigh_local_order(n) = n_local_order;
        temp_ang(id_im==curr_id)= p_ang;
        temp_vec(id_im==curr_id)= p_vec;
        temp_local_order(id_im == curr_id) = n_local_order;
        labeled_by_pca = labeled_by_pca + temp_ang;
        opacity_by_pca = opacity_by_pca + temp_vec;
        labeled_by_local_order = labeled_by_local_order + temp_local_order;
        %toc
    end
    
    neigh_p_ang_all = [neigh_p_ang_all neigh_p_ang];
    neigh_p_vec_all = [neigh_p_vec_all neigh_p_vec];
    neigh_local_order_all = [neigh_local_order_all neigh_local_order];
    Centroid_X_all = [Centroid_X_all c(:,1)];
    Centroid_Y_all = [Centroid_Y_all c(:,2)];
    
    save([outputfolder filesep img_name '_r_' num2str(neighborhood_radius) '_workspace.mat'],'curr_tbl','labeled_by_pca','opacity_by_pca','labeled_by_local_order');
end
disp('done with calculations, now making tables')
save([outputfolder filesep 'r_' num2str(neighborhood_radius) '_workspace.mat'],'combined_tbl','neigh_p_vec_all','neigh_p_ang_all','neigh_local_order');
return_tbl = table;
return_tbl.CellIdentity = combined_tbl.CellIdentity;
return_tbl.Image_name = combined_tbl.Image_name;
return_tbl.Centroid_x = Centroid_X_all;
return_tbl.Centroid_y = Centroid_Y_all;
return_tbl.Neighborhood_mag = neigh_p_vec_all';
return_tbl.Neighborhood_ang = neigh_p_ang_all';
return_tbl.Neighborhood_local_order = neigh_local_order_all';
%%
combined_neigh_tbl = join(combined_tbl,return_tbl,'LeftKeys',[1 5], 'RightKeys', [1 2]);
wt_neigh_tbl = combined_neigh_tbl(combined_neigh_tbl.Genotype == "wt",:);
fzko_neigh_tbl = combined_neigh_tbl(combined_neigh_tbl.Genotype == "fzd6KO",:);
rst_neigh_tbl = combined_neigh_tbl(combined_neigh_tbl.Genotype == "rst",:);
save([outputfolder filesep 'r_' num2str(neighborhood_radius) '_combined_tbl.mat'],'combined_neigh_tbl');
%% normalize magnitudes by all data
radius = ['r_' num2str(neighborhood_radius)];
load([outputfolder filesep radius '_combined_tbl.mat']);
wmats = dir([outputfolder filesep '*_r_' num2str(neighborhood_radius) '_workspace.mat']);
map = colorcet('C7','shift',.75);
map = vertcat([0 0 0; 0 0 0], map);
for m=1:length(wmats)
    load([wmats(m).folder filesep wmats(m).name]);
    [~,img_name,~] = fileparts(wmats(m).name);
    im_name = img_name(1:strfind(wmats(m).name,['_' radius])-1);
    opacity_by_pca(opacity_by_pca==0) = NaN;
    opacity_by_pca(isnan(opacity_by_pca)) = 1;
    %clip values less than -88, as below that in cmap is black for "Zero" and "NaN" value
    labeled_by_pca(labeled_by_pca <= -88) = -88;
    h=figure;set(h,'Visible','on');
    s = imshow(labeled_by_pca,[]);
    colormap(map);
    caxis([-90 90])
    colorbar;
    alpha(s,imadjust(opacity_by_pca,stretchlim(rmmissing(combined_neigh_tbl.Neighborhood_mag)),[0 1]));
    saveas(h,[outputfolder filesep 'r_' num2str(neighborhood_radius) '_' im_name '_' radius '_pca_mag_norm_by_all.tif'])
    %{
    h=figure; set(h,'Visible','on');
    s = imshow(labeled_by_local_order,[]);
    colormap('hot');
    colorbar;
    caxis([0 1]);
    saveas(h,[outputfolder filesep 'r_' num2str(neighborhood_radius) '_' im_name '_' radius '_local_order.tif'])
    %}
end

%% Plots
%{

%% Compare distributions of wt and rst
figure; hold on; 
histogram(wt_neigh_tbl.PCAMagnitude); 
histogram(fzko_neigh_tbl.PCAMagnitude); 
histogram(rst_neigh_tbl.PCAMagnitude);
legend({'WT','Fzd6KO','rst'})
saveas(gcf,'pca_mag_distributions.fig');
close all;

%% Rose plot of original data
figure
ax1 = subplot(1,3,1,polaraxes);
polarhistogram(ax1,[deg2rad(wt_neigh_tbl.PCAAngle___); deg2rad(wt_neigh_tbl.PCAAngle___)+pi])
rlim([0 4000])
title('wt PCA angle ')

ax2 = subplot(1,3,2,polaraxes);
polarhistogram(ax2,[deg2rad(fzko_neigh_tbl.PCAAngle___); deg2rad(fzko_neigh_tbl.PCAAngle___)+pi])
rlim([0 4000])
title('fz6KO PCA angle')

ax3 = subplot(1,3,3,polaraxes);
polarhistogram(ax3,[deg2rad(rst_neigh_tbl.PCAAngle___); deg2rad(rst_neigh_tbl.PCAAngle___)+pi])
rlim([0 4000])
title('rst PCA angle')
sgtitle('PCA angle')
saveas(gcf,'pca_angle_3plot.fig');

figure; polaraxes; 
hold on; 
polarhistogram([deg2rad(wt_neigh_tbl.PCAAngle___); deg2rad(wt_neigh_tbl.PCAAngle___)+pi])
polarhistogram([deg2rad(fzko_neigh_tbl.PCAAngle___); deg2rad(fzko_neigh_tbl.PCAAngle___)+pi])
polarhistogram([deg2rad(rst_neigh_tbl.PCAAngle___); deg2rad(rst_neigh_tbl.PCAAngle___)+pi])
legend({'WT','Fzd6KO','rst'})
saveas(gcf,'pca_angle_overlay_plot.fig');
close all;

%}

%% *Combined Plots*
%%
%{
%Determine area outliers
areas = table2array(area_and_perimeter(:,2));
histogram(areas)
std(areas)
mean(areas)
histogram((areas-mean(areas)))
outliers = (areas-mean(areas))>3*std(areas);
num_outliers = sum(outliers)
area_threshold = mean(areas) + 3*std(areas)
%%
pca_mag = table2array(pca_mag_and_orientation(:,2));
histogram(pca_mag)

pca_angle = table2array(pca_mag_and_orientation(:,3));
polarhistogram(deg2rad(pca_angle))
%% *WT Plots*
%%
areas = table2array(wt_tbl(:,2));
%histogram(areas)
std(areas)
mean(areas)
%histogram((areas-mean(areas)))
outliers = (areas-mean(areas))>3*std(areas);
num_outliers = sum(outliers)

pca_mag = table2array(wt_tbl(:,6));
%histogram(pca_mag)
pca_angle = table2array(wt_tbl(:,7));
%polarhistogram(deg2rad(pca_angle))
%% *Fz6KO Plots*
%%
areas = table2array(fzko_tbl(:,2));
%histogram(areas)
std(areas)
mean(areas)
%histogram((areas-mean(areas)))
outliers = (areas-mean(areas))>3*std(areas);
num_outliers = sum(outliers)

pca_mag = table2array(fzko_tbl(:,6));
%histogram(pca_mag)

pca_angle = table2array(fzko_tbl(:,7));
%polarhistogram(deg2rad(pca_angle))
%% *rst Plots*
%%

areas = table2array(rst_tbl(:,2));
%histogram(areas)
std(areas)
mean(areas)
%histogram((areas-mean(areas)))
outliers = (areas-mean(areas))>3*std(areas);
num_outliers = sum(outliers)
pca_mag = table2array(rst_tbl(:,6));
%histogram(pca_mag)
pca_angle = table2array(rst_tbl(:,7));
%polarhistogram(deg2rad(pca_angle))
%% *Comparison*
%%
figure
ax1 = subplot(1,3,1,polaraxes);
polarhistogram(ax1,deg2rad(table2array(wt_tbl(:,7))))
rlim([0 2000])
title('WT PCA angle')

ax2 = subplot(1,3,2,polaraxes);
polarhistogram(ax2,deg2rad(table2array(fzko_tbl(:,7))))
rlim([0 2000])
title('Fz6KO PCA angle')

ax3 = subplot(1,3,3,polaraxes);
polarhistogram(ax3,deg2rad(table2array(rst_tbl(:,7))))
rlim([0 2000])
title('rst PCA angle')
%}
%%
%{
figure;
histogram(neigh_p_vec_all)
figure;
polarhistogram(deg2rad(neigh_p_ang_all))
%}
%%
%{
figure
ax1 = subplot(1,3,1,polaraxes);
polarhistogram(ax1,deg2rad(table2array(wt_neigh_tbl(:,9))))
rlim([0 2000])
title('WT PCA weighted angle ')

ax2 = subplot(1,3,2,polaraxes);
polarhistogram(ax2,deg2rad(table2array(fzko_neigh_tbl(:,9))))
rlim([0 2000])
title('Fz6KO PCA weighted angle')

ax3 = subplot(1,3,3,polaraxes);
polarhistogram(ax3,deg2rad(table2array(rst_neigh_tbl(:,9))))
rlim([0 2000])
title('rst PCA weighted angle')
sgtitle('PCA angle neighboorhood weighted by magnitude')
%}
%%
