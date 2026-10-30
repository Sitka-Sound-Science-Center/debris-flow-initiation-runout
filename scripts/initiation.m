% Landslide Initiation Model for Southeast Alaska - Tuning and Validation
% Desctription: This code is just to develop validation and tuning method
% for adapting landslide initiation analysis. It does not iterate through
% parameters to find the best combination.
%
% Author: Nicolas Mathews
% Date: 10/13/2025

clear variables
close all
clc
homedir=cd;

%% Add Paths
cd('..');
functions_path=addpath(genpath('functions'));
data_path=addpath(genpath('data'));
geospatial_path=addpath(genpath('geospatial'));
outputs_path=addpath(genpath('outputs'));
config_grfin_path=addpath(genpath('config_grfin'));
cd(homedir);

%% Options

% Soil Parameterization Option
soil_param_option=2; % 1: use parameters below; 2: use soil lookup table

% Soil Parameters (Soil Parameterization Option #1)
c_uniform=9.5; % effective cohesion (kPa)
c_roots_uniform=4.9; % root cohesion (kPa) from Wu et al (1979)
phi_uniform=38.2; % effective internal friction angle (degrees)
gamma_uniform=10.9; % unit weight of soil (kN/m³)
gamma_w=9.81; % unit weight of water (kN/m³)

% Maximum Soil Depth
depth_max=1.5; % maximum soil depth (m)

% ROC Validation On?
ROC_AUC_on=1; % 1: on; 0: off

% Select DEM Name
dem_path='POW2026_5m.tif'; % DEM file name

% TNF Landslide Points Shapefile Name (Used for Validation)
initiation_pts_name='TNF_initiation_POW2026.shp';

% Soil ID Map Name
soil_id_name='POW2026_soil_id.tif';

%% Data Export Options
export_source_areas=1; % 1: export as GeoTIFF | 0: do not export
export_source_areas_suffix='_soil_series_strength_4.9kparoots'; % suffix for raster name

% Initiation Table Output (Input for Grfin Tools)
export_initiation_table=1; % 0: no | 1:yes, export as CSV
init_table_name=[dem_path(1:end-4) '_soil_series_strength_4.9kparoots' '.csv']; % suffix for file name

% Factor of Safety Map
export_fs_map=1; % 1: export as GeoTIFF | 0: do not export

%% Figures Settings
fs_contour_fig_on=0; % 1: on; 0: off

%% Import DEM

fprintf('Loading DEM and Calculating Derivatives...');
   
% Navigate to DEM
cd('..')
cd('geospatial')

% Import DEM
[z,R]=readgeoraster(dem_path);

% Elevations Less Than Zero are NaN
z(z<0)=NaN;

% Return to Home Directory
cd(homedir)

%% Calculate DEM Derivatives

% Extract DEM Characteristics
x_cellsize=R.CellExtentInWorldX; % size of each raster cell, x-dir
x_ext=x_cellsize:x_cellsize:R.RasterExtentInWorldX; % full extent of raster, x-dir
y_cellsize=R.CellExtentInWorldY; % size of each raster cell, y-dir
y_ext=y_cellsize:y_cellsize:R.RasterExtentInWorldY; % full extent of raster, x-dir

% World Coordinates
x_world_ext=R.XWorldLimits;
y_world_ext=R.YWorldLimits;
x_world=single(x_world_ext(1)+x_cellsize:x_cellsize:x_world_ext(2))';
y_world=flipud(single(y_world_ext(1)+y_cellsize:y_cellsize:y_world_ext(2))');

% Calculate Slope and Aspect
[slope,aspect,dx,dy]=gradient_prince(z,x_cellsize,y_cellsize);

fprintf(' Complete.\n\n');

%% Import Soil Depth Model (Create if Nonexistent)

% Navigate to Soil Depth Directory
cd('..')
cd('geospatial')

% Check for existing soil depth raster
soil_depth_file_check=isfile([dem_path(1:end-4) '_soil_depth.mat']);

% Create or Load Soil Depth Layer
if soil_depth_file_check==1
    fprintf('Loading Soil Depth Raster...');
    load([dem_path(1:end-4) '_soil_depth.mat']); % load soil depths from hillslope evolution model
    fprintf(' Complete.\n\n');
elseif soil_depth_file_check==0
    depth=soil_depth(z,x_cellsize,5000,dem_path,0); % run non-linear hillslope evolution model (Roering, 2008)

    % Save soil depth raster
    fprintf('Saving Soil Depth Raster...');
    save([dem_path(1:end-4) '_soil_depth.mat'],'depth');
    fprintf(' Complete.\n\n');
end

% Correct Strange Soil Depth Artifacts
depth(depth>depth_max)=depth_max;

% Return to Home Directory
cd(homedir)

%% Load TWI-Groundwater Relationship (from Swiss Data)

% Read TWI-Groundwater Table
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames');
twi_table=readtable("twi_rinderer_2014.csv");

% Assign Data to Variables
twi_data=twi_table.TWI_ln_m__;
gw_depth_data=twi_table.MedianGroundwaterLevelRelativeToSoilDepth;

% Limit GW Depth to Realistic Values
gw_depth_data(gw_depth_data<0)=0;
gw_depth_data(gw_depth_data>1)=1;

%% Load TWI Raster

% Navigate to Hydrology Directory
cd('..')
cd('geospatial')

% Import Topographic Wetness Index (TWI) Raster
twi=readgeoraster([dem_path(1:end-4) '_topographic_index.tif']);

% Return to Home Directory
cd(homedir)

%% Load Inventory Shapefile (Initiation Points)

% Navigate to Inventory Directory
cd('..')
cd('geospatial')

% IMPORT LANDSLIDE POINTS HERE % TNF_initiation_mitkofkuplind
initiation_pts=shaperead(initiation_pts_name);

% Extract Landslide Coordinates
x_init=[initiation_pts.X];
y_init=[initiation_pts.Y];

% Create Intrinsic Coordinates
[col_init,row_init] = worldToIntrinsic(R, x_init, y_init);

% Return to Home Directory
cd(homedir);

%% Import Soil ID Map, Import Soil Strength Lookup Table, And Build Raster

if soil_param_option==1 % Uniform Soil Parameters

    fprintf('Assigning Uniform Soil Strength...');    

    % Build Soil Parameter Rasters
    c=ones(size(z))*(c_uniform+c_roots_uniform);
    phi=ones(size(z))*phi_uniform;
    gamma=ones(size(z))*gamma_uniform;

    % Make NaNs Match DEM
    c(isnan(z))=NaN;
    phi(isnan(z))=NaN;
    gamma(isnan(z))=NaN;

fprintf(' Complete.\n\n');

elseif soil_param_option==2

    fprintf('Assigning Soil Strength Using Soil Series...');
    

    % Navigate to Soil ID Directory
    cd('..')
    cd('geospatial')

    % Import Soil ID Raster
    soil_id=readgeoraster(soil_id_name);
    soil_id(isnan(soil_id))=0;

    % Return to Home Directory
    cd(homedir);

    % Navigate to Soil Strength Lookup Table
    cd('..')
    cd('data')

    % Import Lookup Table
    soil_lookup=readmatrix('soil_properties_lookup_table.csv');

    % Return to Home Directory
    cd(homedir);

    % Sample Soil Series from Inventory Points
    soil_series_init = interp2(soil_id, col_init, row_init, 'linear');

    % Find Most Prevalent Landslide Soil Series and Assign its Parameters to Undefined Soils
    soil_series_mode=mode(soil_series_init(~isnan(soil_series_init) & soil_series_init~=0));

    phi_undefined=soil_lookup(soil_series_mode,4);
    c_undefined=soil_lookup(soil_series_mode,5);
    gamma_undefined=soil_lookup(soil_series_mode,10);

    % Build Soil Parameter Rasters
    phi=ones(size(z))*phi_undefined;
    c=ones(size(z))*(c_undefined+c_roots_uniform);
    gamma=ones(size(z))*gamma_undefined;

    for i=1:height(soil_lookup) % loop through soil series
        phi(soil_id==i)=soil_lookup(i,4);
        c(soil_id==i)=soil_lookup(i,5)+c_roots_uniform;
        gamma(soil_id==i)=soil_lookup(i,10); % mean bulk unit weight (kN/m3)
    end

    % Make NaNs Match DEM
    c(isnan(z))=NaN;
    phi(isnan(z))=NaN;
    gamma(isnan(z))=NaN;
    
    fprintf(' Complete.\n\n');
end

%% Fit TWI - Saturation to Logistic Equation (Currently for Swiss Data)

% Function Definition
model=@(p, x) p(1)./(1+exp(-p(2)*(x-p(3))));

% Initial parameter guess: [a, b, c]
p0 = [0, 0, 0]; % [L, k, x0]

% Evaluate Best Fit
lsq_options = optimoptions('lsqcurvefit', 'Display', 'off');
[p_fit,r_fit_log]=lsqcurvefit(model,p0,twi_data,gw_depth_data,[],[],lsq_options);

% Generate Fitted Values
twi_fit=0:0.2:max(max(twi));
gw_depth_fit=model(p_fit,twi_fit);

% Compute R^2
gw_depth_fit_rsq=model(p_fit,twi_data);
gw_depth_data_mean=mean(gw_depth_data);  % Mean of observed data
ss_total=sum((gw_depth_data-gw_depth_data_mean).^2);  % Total sum of squares
ss_residual=sum((gw_depth_data-gw_depth_fit_rsq).^2);  % Residual sum of squares
r_squared=1-(ss_residual/ss_total);

%% Compute Saturation Ratio Using Topographic Index (Nonlinear or Linear Fit)

mw=interp1(twi_fit,gw_depth_fit,twi);
mw(mw>1)=1; % limit saturation ratio to 1

%% Compute Stability Using Infinite Slope
fs=(c./(gamma.*depth.*sind(slope)))+((cosd(slope).*tand(phi))./sind(slope)).*(1-mw.*(gamma_w./gamma));

% Make Sure that NaNs are Reflected in FS Raster
fs(isnan(z))=nan;

%% Create Binary Landslide Map (Source Areas)
lsio=false(size(fs));
lsio=fs<=1;

%% Group Unstable Clusters Into Discrete Landslides
ls_cc=bwconncomp(lsio);
ls_cluster_idx=ls_cc.PixelIdxList;

%% Find Volumes of Each Landslide
ls_vol=single(zeros(length(ls_cluster_idx),1)); % landslide volume

% Calculate Volume of Each Landslide Cluster
for i=1:length(ls_vol)
    ls_vol(i)=sum(depth(ls_cluster_idx{i})*x_cellsize*y_cellsize); % landslide volume (m^3)
end

%% Statistics of Landslide Volumes
ls_vol_5th=prctile(ls_vol,5);
ls_vol_95th=prctile(ls_vol,95);
ls_vol_mean=mean(ls_vol);

%% -- VALIDATION: INITIATION -- %%
if ROC_AUC_on==1

% INFILL THIS SECTION TO CALL ROC FUNCTION

else
end

%% -- EXPORT DATA -- %

%% Find Location of Maximum Elevation of Each Landslide (Source Location for Debris Flow Initiation)
ls_max_elev=zeros(length(ls_cluster_idx),1); % max. landslide elevation
ls_max_elev_idx=zeros(length(ls_cluster_idx),1); % index of max. landslide elevation
ls_max_elev_x=zeros(length(ls_cluster_idx),1); % x-location of max. landslide elevation
ls_max_elev_y=zeros(length(ls_cluster_idx),1); % y-location of max. landslide elevation

% Get Index of Pixel with Highest Elevation for Each Landslide
for i=1:length(ls_max_elev)
    [ls_max_elev(i),ls_max_elev_idx_loc(i)]=max(z(ls_cluster_idx{i}));
    ls_max_elev_idx(i)=ls_cluster_idx{i}(ls_max_elev_idx_loc(i));
end

% Convert Index into Subscripts
[y_elev_idx,x_elev_idx]=ind2sub(size(z),ls_max_elev_idx);

% Create Arrays of Initiation Coordinates
ls_init_x=x_world(x_elev_idx);
ls_init_y=y_world(y_elev_idx);

%% Create and Export Table of Debris Flow Initiation Points as CSV (x, y, volume, designator [optional])
init_table=zeros(length(ls_cluster_idx),3);

init_table(:,1)=ls_init_x;
init_table(:,2)=ls_init_y;
init_table(:,3)=ls_vol;

% Export Table
if export_initiation_table==1

    % Navigate to Source Areas Directory
    cd('..');
    cd('grfin_setup');
    cd('initiation_points');

    writematrix(init_table,init_table_name);

    % Navigate to Home Directory
    cd(homedir)
else
end

%% Export Source Area Raster
if export_source_areas==1

    % Navigate to Source Areas Directory
    cd('..');
    cd('outputs');
    cd('source_areas');

    % Save Source Area Raster
    saveraster(dem_path,lsio,[dem_path(1:end-4) '_source_areas' export_source_areas_suffix '.tif'],R);

    % Navigate to Home Directory
    cd(homedir);
end

%% Export FS Raster
if export_fs_map==1

    % Navigate to Source Areas Directory
    cd('..');
    cd('outputs');
    cd('factor_of_safety');

    % Save Factor of Safety Raster
    saveraster(dem_path, fs, [dem_path(1:end-4) '_fs' export_source_areas_suffix '.tif'], R);

    % Navigate to Home Directory
    cd(homedir);
end

%% -- FIGURES -- %%

%% Plot TWI vs. Groundwater Depth Data
figure; hold on;

plot(twi_data,gw_depth_data,'k*');
plot(twi_fit,gw_depth_fit,'-k');

grid on;
xlabel('Topographic Wetness Index [ln(m)]')
ylabel('Median Saturation Ratio, m_w [ ]')

ylim([0 1.05]);
xlim([2 10]);

%% Soil Series Distribution Histogram Plot

% Flatten raster and remove NaNs
soil_id_hist=soil_id;
soil_id_hist(soil_id_hist==0)=NaN;
vals = soil_id_hist(:);
vals = vals(isfinite(vals));

% Get unique soil IDs and counts
[soil_ids, ~, idx] = unique(vals);
counts = accumarray(idx, 1);

% Convert to probability
prob = counts / sum(counts);

% Sort by probability (descending)
[prob_sorted, order] = sort(prob, 'descend');
soil_ids_sorted = soil_ids(order);

% Map soil IDs to names
soil_names = strings(size(soil_ids_sorted));

for i = 1:length(soil_ids_sorted)
    switch soil_ids_sorted(i)
        case 11
            soil_names(i) = "Tolstoi";
        case 1
            soil_names(i) = "Karta";
        case 6
            soil_names(i) = "Wadleigh";
        case 5
            soil_names(i) = "Ulloa";
        case 4
            soil_names(i) = "Tokeen";
        case 9
            soil_names(i) = "Kupreanof";
        case 10
            soil_names(i) = "Traitors";
        otherwise
            soil_names(i) = "Unknown";
    end
end

% Plot
figure
bar(prob_sorted)
xticks(1:length(soil_ids_sorted))
xticklabels(soil_names)
xtickangle(45)
ylabel('Probability')
xlabel('Landslide Source Soil Series')
grid on;

%% Contour Plot (Depth vs. mw with FS shading, slope 35–45)
if fs_contour_fig_on==1

% Flatten rasters into vectors and convert to double
depth_vec = double(depth(:));
mw_vec    = double(mw(:));
fs_vec    = double(fs(:));
slope_vec = double(slope(:));

% Remove NaNs and apply slope filter (35–40°)
valid = ~isnan(depth_vec) & ...
        ~isnan(mw_vec) & ...
        ~isnan(fs_vec) & ...
        ~isnan(slope_vec) & ...
        slope_vec >= 35 & slope_vec <= 45;

depth_vec = depth_vec(valid);
mw_vec    = mw_vec(valid);
fs_vec    = fs_vec(valid);

% Clip FS to [0, 5]
fs_vec(fs_vec < 0) = 0;
fs_vec(fs_vec > 5) = 5;


% Interpolate FS onto a regular grid

% Define grid resolution
depth_grid = linspace(min(depth_vec), max(depth_vec), 150);
mw_grid    = linspace(min(mw_vec), max(mw_vec), 150);

[Depth, MW] = meshgrid(depth_grid, mw_grid);

% Interpolate FS values onto grid
FS_grid = griddata(depth_vec, mw_vec, fs_vec, ...
                   Depth, MW, 'natural');


% Create contour plot

figure

% Filled contours
contourf(Depth, MW, FS_grid, 30, ...
         'LineColor', 'none');

% Color settings
ylorrd=[255 255 204;
        255 237 160;
        254 217 118;
        254 178  76;
        253 141  60;
        252  78  42;
        227  26  28;
        189   0  38;
        128   0  38]/255;
cmap_fs=interp1(linspace(0,1,size(ylorrd,1)),ylorrd,linspace(0,1,256));
cmap_fs=flipud(cmap_fs); % reverse so FS = 1 is dark red

colormap(cmap_fs)
caxis([1 5])

% Colorbar
% cb = colorbar;
% cb.Label.String = 'Factor of Safety (FS)';

% Axis labels
xlabel('Colluvium Depth (m)')
ylabel('Saturation Ratio (m_w)')

% Add FS = 1 failure threshold
hold on

[C,h] = contour(Depth, MW, FS_grid, [1 1], ...
                'k-', 'LineWidth', 2);

clabel(C,h, ...
       'FontSize', 10, ...
       'Color', 'k');

hold off

% Formatting
grid on
set(gca, 'Color', 'w')
box on

else
end