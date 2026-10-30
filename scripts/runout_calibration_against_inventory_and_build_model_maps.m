% Validation for Runout Model
% Description: This code imports a local landslide inventory and compares
% runout model outputs, generated using Grfin Tools, to inventoried
% landslide extents.
%
% Author: Nicolas Mathews
% Date: 11/5/2025

clear variables
close all
clc
homedir=cd;

%% Add Paths
cd('..');
functions_path=addpath(genpath('functions'));
geospatial_path=addpath(genpath('geospatial'));
cd(homedir);

%% Options

% Name of Runout Folder/File
study_area_name='POW2026_runout';
nongrowth_raster_name='POW_a1_011_a2_08_BCx100_growthoff.tif';

% Runout Parameters
slope=[15,20,25]; % entrainment slope threshold
growth_rate=[5,10,15,20,25,30];

% Name of Landslide Inventory File
landslide_inventory_filename='POW2026_TNF_areas_rasterized.tif';

% Name of Watershed Raster
watershed_name='POW2026_huc14.tif';

% DEM Name
DEM_name='POW2026_5m.tif';

%% Simulated Map Creation Options

% Toggle On/Off (simulations must exist for this to work)
create_simulated_maps=1; % 1: create simulated maps; 2: do not create simulated maps

% Options
study_area_name_sim='POW2026_runout';
nongrowth_raster_name_sim='POW_model_a1_011_a2_08_divide0070_growthoff.tif';

%% Import Runout Rasters  (Inventory)

% Navigate to directory
cd('..');
cd('grfin_setup');
cd(study_area_name);

% Load Non-Growth Runout Raster
cd(nongrowth_raster_name(1:end-4))
[runout_nogrowth,R]=readgeoraster(nongrowth_raster_name);
cd('..');

% Load Growth Runout Rasters
runout=cell(length(slope),length(growth_rate));

% Initiate Counter
total=numel(slope)*numel(growth_rate);
counter=0;
msg=sprintf(['Importing Inventory-Based Runout Rasters: %0' num2str(length(num2str(total))) 'd/%d'], counter, total);
fprintf('%s',msg);

for i=1:length(slope)
    for j=1:length(growth_rate)

        % Counter
        counter=counter+1;
        fprintf(repmat('\b', 1, numel(msg)));
        msg=sprintf(['Importing Inventory-Based Runout Rasters: %0' num2str(length(num2str(total))) 'd/%d'], counter, total);
        fprintf('%s',msg);

        cd(['POW_a1_011_a2_08_BCx0070_slope' num2str(slope(i)) '_growth' sprintf('%02d',growth_rate(j))]);

        [runout_raster,~]=readgeoraster(['POW_a1_011_a2_08_BCx0070_slope' num2str(slope(i)) '_growth' sprintf('%02d',growth_rate(j)) '.tif']);
        runout_raster(runout_raster<0)=0;

        % Add Non-Growth Runout
        runout_raster(runout_nogrowth>0)=20;

        runout{i,j}=zeros(size(runout_raster));
        runout{i,j}(runout_raster>0)=2;

        cd('..');
    end
end

fprintf('... Complete.\n\n\n');

% Return to Home Directory
cd(homedir);

%% Import Runout Rasters (Simulated)

if create_simulated_maps==1

    % Navigate to directory
    cd('..');
    cd('grfin_setup');
    cd(study_area_name);

    % Load Non-Growth Runout Raster
    cd(nongrowth_raster_name_sim(1:end-4))
    [runout_nogrowth_sim,R]=readgeoraster(nongrowth_raster_name_sim);
    cd('..');

    % Load Growth Runout Rasters
    runout_sim=cell(length(slope),length(growth_rate));

    % Initiate Counter
    total=numel(slope)*numel(growth_rate);
    counter=0;
    msg=sprintf(['Importing Simulated Runout Rasters: %0' num2str(length(num2str(total))) 'd/%d'], counter, total);
    fprintf('%s',msg);

    for i=1:length(slope)
        for j=1:length(growth_rate)

            % Counter
            counter=counter+1;
            fprintf(repmat('\b', 1, numel(msg)));
            msg=sprintf(['Importing Simulated Runout Rasters: %0' num2str(length(num2str(total))) 'd/%d'], counter, total);
            fprintf('%s',msg);

            cd(['POW_model_a1_011_a2_08_slope' num2str(slope(i)) '_growth' sprintf('%02d',growth_rate(j))]);

            [runout_raster_sim,~]=readgeoraster(['POW_model_a1_011_a2_08_slope' num2str(slope(i)) '_growth' sprintf('%02d',growth_rate(j)) '.tif']);
            runout_raster_sim(runout_raster_sim<0)=0;

            % Add Non-Growth Runout
            runout_raster_sim(runout_nogrowth_sim>0)=20;

            runout_sim{i,j}=zeros(size(runout_raster_sim));
            runout_sim{i,j}(runout_raster_sim>0)=2;

            cd('..');
        end
    end

    fprintf('... Complete.\n\n\n');

    % Return to Home Directory
    cd(homedir);
end

%% Import Landslide Inventory Raster

% Navigate to Directory
cd('..')
cd('geospatial')

% Import Landslide Inventory
[inventory,~]=readgeoraster(landslide_inventory_filename);

% Turn NaNs into Zeros
inventory(isnan(inventory))=0;

% Return to Home Directory
cd(homedir)

%% Import Watershed Raster

% Navigate to Directory
cd('..')
cd('geospatial')

% Import Watershed Raster
[watershed,~]=readgeoraster(watershed_name);
watershed=single(watershed);
watershed(watershed<1 | watershed>120)=NaN;

% Return to Home Directory
cd(homedir)

%% Calculate Similarity Index For Inventory vs. Modeled Runout (new)

% Create cells
watershed_id=unique(watershed(~isnan(watershed)));
FN=cell(1,length(watershed_id));
FP=cell(1,length(watershed_id));
TP=cell(1,length(watershed_id));
similarity_index=cell(1,length(watershed_id));

% Initiate Counter
total=length(watershed_id)*length(slope)*length(growth_rate);
counter=0;
msg=sprintf(['Calculating Similarity: %0' num2str(length(num2str(total))) 'd/%d'], counter, total);
fprintf('%s',msg);

for k=1:length(watershed_id)

    % Create Watershed Mask
    watershed_mask=watershed==(watershed_id(k));

    % Create Landslide Map in Watershed Mask
    inventory_watershed=zeros(size(watershed_mask));
    inventory_watershed(watershed_mask)=inventory(watershed_mask);

    for i=1:length(slope)
        for j=1:length(growth_rate)

            % Counter
            counter=counter+1;
            fprintf(repmat('\b', 1, numel(msg)));
            msg=sprintf(['Calculating Similarity: %0' num2str(length(num2str(total))) 'd/%d'], counter, total);
            fprintf('%s',msg);

            % Initiate Variables
            runout_watershed=zeros(size(watershed_mask));
            runout_watershed(watershed_mask)=runout{i,j}(watershed_mask);

            % Add Runout and Inventory for Similarity Calculation
            similarity_raster=inventory_watershed+runout_watershed;

            % Count Runout, Inventory, and Overlap Cells (1: inventory only; 2: runout only; 3: overlap/agreement)
            FN{k}(i,j)=sum(similarity_raster==1,'all'); % false negative: inventoried slide, but no modeled runout
            FP{k}(i,j)=sum(similarity_raster==2,'all'); % false positive: modeled runout, but no inventoried slide
            TP{k}(i,j)=sum(similarity_raster==3,'all'); % true positive: model and inventory overlap

            % Calculate Similarity Index (Gorr et al., 2022; Heiser et al., 2017)
            alpha_t{k}(i,j)=TP{k}(i,j)/(TP{k}(i,j)+FP{k}(i,j)+FN{k}(i,j));
            beta_t{k}(i,j)=FN{k}(i,j)/(TP{k}(i,j)+FP{k}(i,j)+FN{k}(i,j));
            gamma_t{k}(i,j)=FP{k}(i,j)/(TP{k}(i,j)+FP{k}(i,j)+FN{k}(i,j));

            similarity_index{k}(i,j)=alpha_t{k}(i,j)-beta_t{k}(i,j)-gamma_t{k}(i,j);
        end
    end
end

fprintf('... Complete.\n\n\n');

%% Find Best Parameters for Each Watershed

for k=1:length(watershed_id)

[similarity_max_watershed(k),similarity_max_watershed_idx(k)]=max(similarity_index{k}(:),[],'omitnan');
[similarity_i(k), similarity_j(k)] = ind2sub(size(similarity_index{k}),similarity_max_watershed_idx(k));

alpha_t_opt(k)=alpha_t{k}(similarity_max_watershed_idx(k)); % true positive
beta_t_opt(k)=beta_t{k}(similarity_max_watershed_idx(k)); % false negative
gamma_t_opt(k)=gamma_t{k}(similarity_max_watershed_idx(k)); % false positive

slope_watershed(k)=slope(similarity_i(k));
growth_rate_watershed(k)=growth_rate(similarity_j(k));
end

% Make Arrays Vertical
similarity_max_watershed=similarity_max_watershed';
slope_watershed=slope_watershed';
growth_rate_watershed=growth_rate_watershed';

% Insert NaN Where No Landslides Exist
slope_watershed(isnan(similarity_max_watershed))=NaN;
growth_rate_watershed(isnan(similarity_max_watershed))=NaN;

%% Build Optimized Runout Raster (Inventoried Slides)

% Initialize Runout Raster
runout_opt=zeros(size(runout_raster));

% Make NaN (Watersheds With No Landslides) Default to Highest Slope and Lowest Growth
similarity_i_raster_opt=similarity_i;
similarity_j_raster_opt=similarity_j;

similarity_i_raster_opt(isnan(similarity_max_watershed))=find(slope==max(slope));
similarity_j_raster_opt(isnan(similarity_max_watershed))=find(growth_rate==min(growth_rate));

for k=1:length(watershed_id)

    % Create Watershed Mask
    watershed_mask=watershed==(watershed_id(k));

    % Assign Runout to Watershed
    runout_opt(watershed_mask)=runout{similarity_i_raster_opt(k),similarity_j_raster_opt(k)}(watershed_mask);
end

%% Build Optimized Raster with Range of Growth (Will Not Work for 10:growth, 20:nogrowth Rasters)

% Initialize Runout Raster
runout_opt_growth=zeros(size(runout_raster));

% Make NaN (Watersheds With No Landslides) Default to Highest Slope and Lowest Growth
similarity_i_raster_opt=similarity_i;
similarity_j_raster_opt=similarity_j;

similarity_i_raster_opt(isnan(similarity_max_watershed))=find(slope==max(slope));
similarity_j_raster_opt(isnan(similarity_max_watershed))=find(growth_rate==min(growth_rate));

for k=1:length(watershed_id)
k
    % Create Watershed Mask
    watershed_mask=watershed==(watershed_id(k));

    for j=fliplr(1:length(growth_rate))

        % Assign Runout to Watershed
        runout_mask_current=zeros(size(watershed_mask));
        runout_mask_current(watershed_mask)=runout{similarity_i_raster_opt(k),j}(watershed_mask);
        runout_mask_current=logical(runout_mask_current);
        
        runout_opt_growth(runout_mask_current)=growth_rate(j);
    end
end




%% Build Optimized Runout Raster (Simulated Landslides)

if create_simulated_maps==1

    % Initialize Runout Raster
    runout_opt_sim=zeros(size(runout_raster_sim));

    % Make NaN (Watersheds With No Landslides) Default to Highest Slope and Lowest Growth
    similarity_i_raster_opt=similarity_i;
    similarity_j_raster_opt=similarity_j;

    similarity_i_raster_opt(isnan(similarity_max_watershed))=find(slope==max(slope));
    similarity_j_raster_opt(isnan(similarity_max_watershed))=find(growth_rate==min(growth_rate));

    for k=1:length(watershed_id)
        k
        % Create Watershed Mask
        watershed_mask=watershed==(watershed_id(k));

        % Assign Runout to Watershed
        runout_opt_sim(watershed_mask)=runout_sim{similarity_i_raster_opt(k),similarity_j_raster_opt(k)}(watershed_mask);
    end
end

%% Build Optimized Raster with Range of Growth (Simulated Landslides; Will Not Work for 10:growth, 20:nogrowth Rasters)

% Initialize Runout Raster
runout_opt_growth_sim=zeros(size(runout_raster_sim));

% Make NaN (Watersheds With No Landslides) Default to Highest Slope and Lowest Growth
similarity_i_raster_opt=similarity_i;
similarity_j_raster_opt=similarity_j;

similarity_i_raster_opt(isnan(similarity_max_watershed))=find(slope==max(slope));
similarity_j_raster_opt(isnan(similarity_max_watershed))=find(growth_rate==min(growth_rate));

for k=1:length(watershed_id)
    k
    % Create Watershed Mask
    watershed_mask=watershed==(watershed_id(k));

    for j=fliplr(1:length(growth_rate))

        % Assign Runout to Watershed
        runout_mask_current=zeros(size(watershed_mask));
        runout_mask_current(watershed_mask)=runout_sim{similarity_i_raster_opt(k),j}(watershed_mask);
        runout_mask_current=logical(runout_mask_current);
        
        runout_opt_growth_sim(runout_mask_current)=growth_rate(j);
    end
end

%% Save Rasters

% Navigate to Directory
cd('..')
cd('outputs')
cd('calibrated_runout')

% Export Runout Maps for Inventoried Slides
saveraster(DEM_name,runout_opt,'POW2026_optimized_runout',R);
saveraster(DEM_name,runout_opt_growth,'POW2026_optimized_runout_growth',R);

% Export Runout Maps for Modeled/Simulated Slides
if create_simulated_maps==1
    saveraster(DEM_name,runout_opt_sim,'POW2026_model_optimized_runout',R);
    saveraster(DEM_name,runout_opt_growth_sim,'POW2026_model_optimized_runout_growth',R);
end