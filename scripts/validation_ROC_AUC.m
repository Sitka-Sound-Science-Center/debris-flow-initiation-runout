%% -- VALIDATION: INITIATION -- %%
% Description: Compares 12 versions of the infinite slope model using the
% same inventory and background points. Background points are only drawn
% from steep terrain (slope >= slope_thresh).

homedir=cd;

%% Add Paths
cd('..');
functions_path=addpath(genpath('functions'));
cd(homedir);

%% Validation Options

% Steep Terrain Threshold
slope_thresh=25; % background points only drawn from slopes >= this value (degrees)

% Filter Inventory Points by Slope?
filter_inventory=0; % 1: only use inventory points on slopes >= slope_thresh | 0: use all inventory points

% Buffer Options
r_m=10; % buffer radius around inventory points (m)
buffer_bg=0; % 1: apply same buffer to background points | 0: use single pixel for background points

% Background Point Options
excl_m=20; % minimum distance between background and inventory points (m)
num_rand_pts=5000; % number of background points per draw
n_draws=20; % number of background point draws

% Replace Infinite FS (Zero Colluvium Depth)
fs_cap=100;

% Model Options
depth_opts={depth,0.5,1.0,1.5}; % colluvium depth options (first is modeled colluvium)
depth_lbl={'Variable','0.5 m','1.0 m','1.5 m'};
sat_opts={0.5,mw,1}; % saturation ratio options
sat_lbl={'Uniform 0.5','TWI (median)','Saturated'};

% Reference Model
ref_s=2; % reference saturation (TWI-based)
ref_d=1; % reference depth (modeled colluvium)

% Set Random Seed for Reproducibility
rng(1);

%% Convert Distances to Pixels

cell_m=mean([x_cellsize,y_cellsize]); % cell size (m)
r_pix=round(r_m/cell_m); % buffer radius (pixels)
r_pix_bg=r_pix*buffer_bg; % background buffer radius (pixels), 0 = single pixel
ex_pix=round(excl_m/cell_m); % exclusion distance (pixels)

%% Create Steep Terrain Mask

steep=(slope>=slope_thresh) & ~isnan(z) & ~isnan(fs);

%% Process Inventory Points

% Convert Inventory Coordinates to Rows and Columns
[col_i,row_i]=worldToIntrinsic(R,x_init,y_init);
row_i=round(row_i(:));
col_i=round(col_i(:));

% Remove Points Outside of DEM
in_dem=row_i>=1 & row_i<=size(fs,1) & col_i>=1 & col_i<=size(fs,2);
row_i=row_i(in_dem);
col_i=col_i(in_dem);
lin_i_all=sub2ind(size(fs),row_i,col_i);

% Check Number of Inventory Points on Steep Terrain
on_steep=steep(lin_i_all);
fprintf('Initiation points on slopes >= %d deg: %d of %d (%.0f%%)\n', ...
    slope_thresh,sum(on_steep),numel(on_steep),100*mean(on_steep));

% Select Inventory Points
if filter_inventory==1
    keep_i=on_steep;
    fprintf('Using only inventory points on steep terrain\n\n');
elseif filter_inventory==0
    keep_i=~isnan(fs(lin_i_all)); % all points with valid FS
    fprintf('Using all %d inventory points with valid FS\n\n',sum(keep_i));
end

row_i=row_i(keep_i);
col_i=col_i(keep_i);
N_init=numel(row_i);

%% Create Background Points (Steep Cells Away from Inventory Points)

fprintf('Creating Background Points...');

% Exclude Cells Near Any Inventory Point
near_ls=false(size(fs));
near_ls(lin_i_all)=true;
near_ls=imdilate(near_ls,strel('disk',ex_pix,0));

% Pool of Valid Background Cells
bg_pool=find(steep & ~near_ls);

% Draw All Background Points Before Evaluating Models (Same Points for Every Model)
bg_idx=zeros(num_rand_pts,n_draws);
for i=1:n_draws
    bg_idx(:,i)=randsample(bg_pool,num_rand_pts);
end

% Convert Indices to Rows and Columns
[row_bg_all,col_bg_all]=ind2sub(size(fs),bg_idx(:));

fprintf(' Complete. Background pool: %d steep cells\n\n',numel(bg_pool));

clear near_ls bg_pool

%% Establish Model Evaluation Order (Reference Model First)

n_depth=numel(depth_opts); % number of depth options
n_sat=numel(sat_opts); % number of saturation options
n_models=n_depth*n_sat; % total number of models
ref_k=(ref_s-1)*n_depth+ref_d; % reference model index

model_order=zeros(n_models,2); % [saturation option, depth option]
model_order(1,:)=[ref_s ref_d];
count=1;

for s=1:n_sat
    for d=1:n_depth
        if ~(s==ref_s && d==ref_d)
            count=count+1;
            model_order(count,:)=[s d];
        end
    end
end

%% Evaluate Each Model

% Preallocate
fs_min_init=nan(N_init,n_models); % minimum FS within buffer, inventory points
fs_min_bg=nan(num_rand_pts*n_draws,n_models); % minimum FS within buffer, background points
mdl_depth=strings(n_models,1); % depth label for each model
mdl_sat=strings(n_models,1); % saturation label for each model
hit_rate=nan(n_models,1); % fraction of inventory points with FS <= 1
area_flag=nan(n_models,1); % fraction of steep terrain with FS <= 1
area_match=nan(n_models,1); % fraction of steep terrain flagged at reference hit rate
t_match=nan(n_models,1); % FS threshold matching reference hit rate

for i=1:n_models

    % Model Options
    s=model_order(i,1);
    d=model_order(i,2);
    k=(s-1)*n_depth+d; % model index
    mdl_depth(k)=depth_lbl{d};
    mdl_sat(k)=sat_lbl{s};

    fprintf('Evaluating Model %d of %d (Depth: %s, Saturation: %s)...', ...
        i,n_models,depth_lbl{d},sat_lbl{s});

    % Compute Factor of Safety
    fs_model=inf_slope_fs(c,gamma,phi,slope,depth_opts{d},sat_opts{s},gamma_w,z,fs_cap);

    % Check that Reference Model Reproduces FS Raster
    if k==ref_k
        fs_check=fs;
        fs_check(isinf(fs_check))=fs_cap;
        fprintf(' Max. difference from FS raster: %.2e...', ...
            max(abs(fs_model(steep)-fs_check(steep))));
        clear fs_check
    end

    % Extract Minimum FS for Inventory and Background Points
    fs_min_init(:,k)=buffer_min(fs_model,row_i,col_i,r_pix);
    fs_min_bg(:,k)=buffer_min(fs_model,row_bg_all,col_bg_all,r_pix_bg);

    % Hit Rate and Flagged Area on Steep Terrain
    fs_steep=fs_model(steep);
    hit_rate(k)=mean(fs_min_init(:,k)<=1);
    area_flag(k)=mean(fs_steep<=1);

    % Area Required to Match Reference Model Hit Rate
    t_match(k)=quantile(fs_min_init(:,k),hit_rate(ref_k));
    area_match(k)=mean(fs_steep<=t_match(k));

    clear fs_model fs_steep

    fprintf(' Complete.\n');
end

fprintf('\n');

%% ROC/AUC Performance Evaluation (Repeated Draws)

% Labels (1: landslide, 0: background)
labels=[ones(N_init,1); zeros(num_rand_pts,1)];

% Preallocate
AUC_all=nan(n_draws,n_models); % AUC for each draw and model
TPR_fs1=nan(n_draws,n_models); % true positive rate at FS = 1
FPR_fs1=nan(n_draws,n_models); % false positive rate at FS = 1
roc1=cell(n_models,1); % ROC curves from first draw (for plotting)

% Valid Inventory Points (Valid in All Models)
ok_init=all(~isnan(fs_min_init),2);

for i=1:n_draws

    % Background Points for This Draw
    rows_draw=(i-1)*num_rand_pts+(1:num_rand_pts);

    % Points Valid in All Models
    ok=[ok_init; all(~isnan(fs_min_bg(rows_draw,:)),2)];
    lab=labels(ok);

    for k=1:n_models
        scores=[fs_min_init(:,k); fs_min_bg(rows_draw,k)];
        scores=scores(ok);

        % Lower FS = more unstable, so negate for perfcurve
        [fpr,tpr,thr,AUC_all(i,k)]=perfcurve(lab,-scores,1);

        % TPR and FPR at FS = 1
        TPR_fs1(i,k)=mean(scores(lab==1)<=1);
        FPR_fs1(i,k)=mean(scores(lab==0)<=1);

        % Save Curves from First Draw
        if i==1
            roc1{k}=struct('FPR',fpr,'TPR',tpr,'T',thr, ...
                'FPR1',FPR_fs1(i,k),'TPR1',TPR_fs1(i,k));
        end
    end
end

%% Reference Model Outputs (Used Later in Script)

FPR=roc1{ref_k}.FPR;
TPR=roc1{ref_k}.TPR;
T=roc1{ref_k}.T;
AUC=mean(AUC_all(:,ref_k));

% Background Point Coordinates (First Draw)
[x_bg,y_bg]=intrinsicToWorld(R,col_bg_all(1:num_rand_pts),row_bg_all(1:num_rand_pts));

% Compute Best Threshold (Youden's J)
J=TPR-FPR;
[~,idx_best]=max(J);
best_threshold=T(idx_best);
best_fs_cutoff=-best_threshold; % convert back to original FS scale

%% Summary Table

roc_summary=table(mdl_depth,mdl_sat, ...
    round(mean(AUC_all)',3),round(min(AUC_all)',3),round(max(AUC_all)',3), ...
    round(100*mean(TPR_fs1)',1),round(100*mean(FPR_fs1)',1), ...
    round(100*min(FPR_fs1)',1),round(100*max(FPR_fs1)',1), ...
    round(100*area_flag,1),round(100*area_match,1), ...
    'VariableNames',{'Depth','Saturation','AUC_mean','AUC_min','AUC_max', ...
    'TPR_FS1_pct','FPR_FS1_pct','FPR_FS1_min','FPR_FS1_max', ...
    'SteepAreaFS1_pct','AreaAtRefHitRate_pct'});

% Print Summary
fprintf('--- Steep Terrain Validation (Slope >= %d deg, %d Inventory Points, %d Draws x %d Background Points) ---\n', ...
    slope_thresh,N_init,n_draws,num_rand_pts);

if buffer_bg==1
    fprintf('Buffer: %d m around inventory and background points\n',r_m);
elseif buffer_bg==0
    fprintf('Buffer: %d m around inventory points only; background points use single pixel\n',r_m);
end

fprintf('TPR_FS1 = %% of inventory points with FS <= 1\n');
fprintf('FPR_FS1 = %% of background points with FS <= 1 (mean, min, max over draws)\n');
fprintf('SteepAreaFS1 = %% of all steep terrain with FS <= 1\n');
fprintf('AreaAtRefHitRate = %% of steep terrain flagged to match reference hit rate (%.1f%%)\n\n', ...
    100*hit_rate(ref_k));

disp(roc_summary)

% Export Summary Table
% writetable(roc_summary,'roc_summary_steep.csv');

%% Print Values for Figure Caption

fprintf('\n--- Figure Caption Values (Slope >= %d deg, %d Draws) ---\n',slope_thresh,n_draws);

depth_leg={'Modeled colluvium','0.5 m uniform','1.0 m uniform','1.5 m uniform'};

for s=1:n_sat
    for d=1:n_depth
        k=(s-1)*n_depth+d;
        fprintf('%-18s %-13s AUC = %.3f (%.3f-%.3f)\n',depth_leg{d},sat_lbl{s}, ...
            mean(AUC_all(:,k)),min(AUC_all(:,k)),max(AUC_all(:,k)));
    end
end

fprintf('Inventory points with no modeled colluvium within buffer: %.1f%%\n\n', ...
    100*mean(fs_min_init(:,ref_k)>=fs_cap));

%% -- FIGURES AND PLOTS -- %%

%% Plot ROC Curves for All Models

% Figure Size (inches)
fig_w=6.5; % width (3.5: single column; 7: full page)
fig_h=4.0; % height

% Line Colors (Colluvium Depth)
depth_col=[0.165 0.471 0.839;  % modeled colluvium (blue)  #2a78d6
           0.949 0.604 0.431;  % 0.5 m (light orange)      #f29a6e
           0.831 0.337 0.110;  % 1.0 m (orange)            #d4561c
           0.541 0.208 0.063]; % 1.5 m (dark orange)       #8a3510

% Line Styles (Saturation)
sat_ls={':','-','--'}; % dotted: mw = 0.5 | solid: TWI-based | dashed: saturated
sat_leg={'m_w = 0.5','TWI (median)','saturated'};

% Open Figure
figure;
set(gcf,'Units','inches');
set(gcf,'Position',[1 1 fig_w fig_h]);
hold on;

% Plot ROC Curve for Each Model
h_roc=gobjects(n_models,1);
leg_txt=cell(n_models,1);

for s=1:n_sat
    for d=1:n_depth
        k=(s-1)*n_depth+d;
        h_roc(k)=plot(roc1{k}.FPR,roc1{k}.TPR,sat_ls{s}, ...
            'Color',depth_col(d,:),'LineWidth',1.2+0.8*(d==1));
        leg_txt{k}=sprintf('%s, %s (AUC = %.3f)',depth_leg{d},sat_leg{s},mean(AUC_all(:,k)));
    end
end

% Plot Random Guessing Line
plot([0 1],[0 1],':','Color',[0.5 0.5 0.5],'LineWidth',1,'HandleVisibility','off');

% Axis Settings
grid on;
box on;
xlim([0 1]);
ylim([0 1]);
xlabel('False Positive Rate');
ylabel('True Positive Rate');
legend(h_roc,leg_txt,'Location','southeast','FontSize',8);
set(gca,'FontSize',9);

hold off

%% Plot AUC Grid (Rows: Colluvium Depth, Columns: Saturation)

% Figure Size (inches)
fig_w_b=4.5; % width
fig_h_b=3.2; % height

% Text Options (labels can be added in PowerPoint instead)
show_axis_labels=0; % 1: show row/column names and axis titles | 0: hide
show_cell_values=1; % 1: show AUC in each cell | 0: hide

% AUC Grid (Depth x Saturation)
AUC_grid=reshape(mean(AUC_all,1),n_depth,n_sat);

% Color Limits (Data Range Rounded to Nearest 0.01)
clim_auc=[floor(100*min(AUC_grid(:)))/100, ceil(100*max(AUC_grid(:)))/100];
% clim_auc=[0.5 1]; % anchor at random guessing (reduces contrast)

% Row and Column Labels
row_lbl={'Modeled','Uniform, 0.5 m','Uniform, 1.0 m','Uniform, 1.5 m'};
col_lbl={'Uniform ({\itm_w} = 0.5)','TWI-based (median {\itm_w})','Fully saturated ({\itm_w} = 1)'};

% Green Color Ramp (Light to Dark = Low to High AUC)
color_stops=[217 240 211;  % #d9f0d3
             166 219 160;  % #a6dba0
              90 174  97;  % #5aae61
              27 120  55;  % #1b7837
               0  68  27]/255; % #00441b
cmap_auc=interp1(linspace(0,1,size(color_stops,1)),color_stops,linspace(0,1,256));

% Open Figure
figure;
set(gcf,'Units','inches');
set(gcf,'Position',[1 1 fig_w_b fig_h_b]);

% Plot AUC Grid
imagesc(AUC_grid);
colormap(gca,cmap_auc);
caxis(clim_auc);
axis tight;
hold on;

% White Gridlines Between Cells
for i=0.5:1:n_sat+0.5
    plot([i i],[0.5 n_depth+0.5],'w-','LineWidth',2);
end
for i=0.5:1:n_depth+0.5
    plot([0.5 n_sat+0.5],[i i],'w-','LineWidth',2);
end

% AUC Value in Each Cell (White Text on Dark Cells, Dark Text on Light Cells)
if show_cell_values==1
    for d=1:n_depth
        for s=1:n_sat
            auc_val=AUC_grid(d,s);
            auc_frac=(auc_val-clim_auc(1))/diff(clim_auc);
            if auc_frac>0.55
                txt_col=[1 1 1];
            else
                txt_col=[0.1 0.1 0.1];
            end
            text(s,d,sprintf('%.3f',auc_val),'HorizontalAlignment','center', ...
                'VerticalAlignment','middle','Color',txt_col,'FontSize',9,'FontWeight','bold');
        end
    end
end

% Axis Settings
set(gca,'XTick',1:n_sat,'YTick',1:n_depth,'TickLength',[0 0], ...
    'FontSize',8,'XAxisLocation','top','Box','off');
colbar=colorbar;
colbar.FontSize=8;

if show_axis_labels==1
    set(gca,'XTickLabel',col_lbl,'YTickLabel',row_lbl);
    xtickangle(0);
    xlabel('Saturation');
    ylabel('Colluvium Depth');
    colbar.Label.String='AUC';
elseif show_axis_labels==0
    set(gca,'XTickLabel',[],'YTickLabel',[]);
end

hold off