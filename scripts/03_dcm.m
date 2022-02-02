%% 03_dcm

% author: Casey Paquola
% email: casey.paquola@gmail.com
% feel free to get in touch if you have any questions

for init_project = 1
   
    % set directories
    GH          = 'C:/Users/cpaq3/Desktop/';    % github repo location
    myDir       = 'C:/Users/cpaq3/OneDrive/';   % personal directory for housing extra things
    homeDir     = [GH '/DMN'];                  % home directory
    bbwDir      = [GH '/BigBrainWarp'];         % BigBrainWarp location
    
    % require fsaverage5 subject directory from Freesurfer
    fsAv        = [myDir '/useful_scripts/fsaverage5/'];  % fsaverage5 directory
    
    % addpaths - requires micaopen (github), BigBrainWarp (github), 
    % freesurfer matlab toolbox, and gifti toolbox for matlab
    addpath(genpath([GH '/micaopen']));
    addpath([bbwDir '/scripts/']);
    addpath(genpath([homeDir '/scripts/']));
    addpath(genpath([myDir '/useful_scripts/gifti-master']));
    addpath(genpath([myDir 'useful_scripts/freesurfer_matlab']));
    addpath(genpath([GH '/BrainSpace/matlab/']))
    
    % load colourmaps - using colorbrewer and scientific colours 
    load([myDir '/useful_scripts/colourmaps/colorbrewer.mat'])
    scicol = [myDir '/useful_scripts/colourmaps/ScientificColourMaps7'];
    scicol_names = ["devon", "lajolla", "lapaz", "batlow", "roma", "vik", "berlin", "bam", "cork"];
    for ii = 1:length(scicol_names)
        load(strcat(scicol, '/', scicol_names(ii), '/', scicol_names(ii), '.mat'));    
    end
    cmap_types = [127 140 172; 139 167 176; 171 186 162; 218 198 153; 253 211 200; 252 229 252]/255;
    [~, ~, atlas7_ctb] = read_annotation([fsAv '/lh.Yeo2011_7Networks_N1000.annot']);
    cmap_atlas7 = atlas7_ctb.table(:,1:3)/255;
    cmap_dmn = [0.8 0.8 0.8; cmap_atlas7(end,:)];
    [~, ~, atlas17_ctb] = read_annotation([fsAv '/lh.Yeo2011_17Networks_N1000.annot']);
    cmap_atlas17 = atlas17_ctb.table(:,1:3)/255;
    
    % load surfaces
    FS = SurfStatAvSurf({[fsAv '/surf/lh.pial'], [fsAv '/surf/rh.pial']});
    FSinf = SurfStatAvSurf({[fsAv '/surf/lh.inflated'], [fsAv '/surf/rh.inflated']});
    
end

for load_parcellations_and_cytoarchitectonic_map = 1
      
    parc_name='Schaefer2018_400Parcels_7Networks_order';  % what parcellation the connectomes are in
    
    % load atlases on fsaverage5
    atlas7_fs = annot2classes([fsAv '/label/lh.Yeo2011_7Networks_N1000.annot'], ...
        [fsAv '/label/rh.Yeo2011_7Networks_N1000.annot'], 0); % from CBIG github (https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation)
    econ_fs = annot2classes([fsAv '/label/lh.economo.annot'],[fsAv '/label/rh.economo.annot'],0); % from Dutch Connectome Lab (http://www.dutchconnectomelab.nl/)
    [~,~,econ_ctb] = read_annotation([fsAv '/label/lh.economo.annot']); 
    econ_ctb.types = [0;0;2;3;4;3;3;3;2;2;3;3;3;4;5;6;6;6;5;4;6;6;4;4;6;6;6;2;1;1;2;1;2;3;2;3;4;3;3;2;1;1;2;4;5]; % hard coded based on table data in Garcia-Cabezas (2021)
    econ_types_fs = BoSurfStatMakeParcelData(econ_ctb.types([1 3:45]),FS,econ_fs);
    econ_types_fs(econ_types_fs==0) = nan;
    [sch400, ~] = annot2classes([fsAv '/lh.Schaefer2018_400Parcels_7Networks_order.annot'],...
        [fsAv '/rh.Schaefer2018_400Parcels_7Networks_order.annot'],1);
    sch400_cort = [2:201 203:402];
    
    % calculate modal atlas7 and econ type assignment for each parcel
    usch400 = unique(sch400);
    for ii = 1:length(usch400)
        sch400_atlas7(ii) = mode(atlas7_fs(sch400==usch400(ii)));
        sch400_econ_types(ii) = mode(econ_types_fs(sch400==usch400(ii)));
    end
    sch400_atlas7_cort = sch400_atlas7(sch400_cort);
    sch400_econ_types_cort = sch400_econ_types(sch400_cort);
    
    % load data-driven cytoarchitectural axis
    %(requires same parcellations on bigbrain)
    load([homeDir '/output/downsample_bigbrain_100k.mat'], 'nn_bb', 'BB100')
    load([homeDir '/output/bigbrain_embedding_100k_thresh.mat'], 'embedding', 'results')
    BB = SurfStatAvSurf({[bbwDir '/spaces/tpl-bigbrain/tpl-bigbrain_hemi-L_desc-pial_sym.obj'], ...
     [bbwDir '/spaces/tpl-bigbrain/tpl-bigbrain_hemi-R_desc-pial_sym.obj']}); % BigBrain pial surface
    parc_name='Yeo2011_7Networks_N1000';  
    lh = gifti([bbwDir '/spaces/tpl-bigbrain/tpl-bigbrain_hemi-L_desc-' parc_name '.label.gii']);
    rh = gifti([bbwDir '/spaces/tpl-bigbrain/tpl-bigbrain_hemi-R_desc-' parc_name '.label.gii']);
    atlas7_bb = [lh.cdata; rh.cdata]; 
    toMap = zeros(1,length(unique(nn_bb)));
    toMap(ismember(unique(nn_bb),unique(nn_bb(atlas7_bb==7)))) = embedding(:,1) * -1;
    E1 = BoSurfStatMakeParcelData(toMap, BB, nn_bb); % upsample to whole cortex
    
    % parcellate data-driven cytoarchitectural axis
    sch400_bb = annot2classes([homeDir '/utilities/tpl-bigbrain_hemi-L_desc-Schaefer2018_400Parcels_7Networks_order.annot'], ...
        [homeDir '/utilities/tpl-bigbrain_hemi-R_desc-Schaefer2018_400Parcels_7Networks_order.annot'], 1);
    E1_parc = grpstats(E1, sch400_bb);
    E1_parc_cort = E1_parc(sch400_cort);
      
    % spin the atlas
    run_spins = 0;
    if run_spins == 1;
        sch400_centroids_onsphere = grpstats(FSinf.coord', sch400);
        perm_id = rotate_parcellation(sch400_centroids_onsphere(1:end/2,:), ...
            sch400_centroids_onsphere((end/2)+1:end,:), 10000);
        save([homeDir '/utilities/sch400_spin.mat'], 'perm_id')
    else
       load([homeDir '/utilities/sch400_spin.mat'], 'perm_id')
    end
    E1_spun = E1(perm_id);
    sch400_atlas7_spun = sch400_atlas7(perm_id);
    
end

for load_data = 1

    % for more details on construction of the rDCM model see the 
    % the nested preprocessing scripts [homeDir '/scripts/rDCM/'] or 
    % TAPAS (https://github.com/translationalneuromodeling/tapas)
    dataset = 'mics';
    load([homeDir '/data/' dataset '_rDCM_sch400.mat'], 'results');
    mat = results.mean_Amatrix_allSubjects;
    mat = abs(mat) .* ~eye(size(mat));
    
    % compile afferent (aff) and efferent (eff) connectivity of the DMN
    mat_aff_node_type = grpstats(mat(sch400_atlas7_cort==8,sch400_atlas7_cort~=8)', ...
        sch400_econ_types_cort(sch400_atlas7_cort~=8));
    mat_eff_node_type = grpstats(mat(sch400_atlas7_cort~=8,sch400_atlas7_cort==8), ...
        sch400_econ_types_cort(sch400_atlas7_cort~=8));
    
    aff_dmn_node = mean(mat(sch400_atlas7_cort==8,sch400_atlas7_cort~=8)'); 
    std(aff_dmn_node)/mean(aff_dmn_node)
    eff_dmn_node = mean(mat(sch400_atlas7_cort~=8,sch400_atlas7_cort==8));
    std(eff_dmn_node)/mean(eff_dmn_node)

    
end

for comparison_with_cyto = 1

    % Axis and interaction effects
    x = E1_parc_cort(sch400_atlas7(sch400_cort)==8); % E1 within DMN
    aff_r = corr(x, mean(mat(sch400_atlas7(sch400_cort)==8,sch400_atlas7(sch400_cort)~=8),2));
    eff_r = corr(x, mean(mat(~sch400_atlas7(sch400_cort)~=8,sch400_atlas7(sch400_cort)==8)',2));
    aff_types_r = corr(x, mat_aff_node_type');
    eff_types_r = corr(x, mat_eff_node_type');
    
    % recalculate with spun atlas
    mat_exp = nan(length(unique(sch400)),length(unique(sch400)));  % expand matrix to correspond with spun indicies
    mat_exp(sch400_cort,sch400_cort) = mat; 
    aff_perm = zeros(length(perm_id), length(aff_r));
    eff_perm = zeros(length(perm_id), length(aff_r));
    for kk = 1:length(perm_id)
        keep_idx = sch400_atlas7_spun(:,kk)==8;
        % average 
        aff_perm(kk) = corr(E1_spun(keep_idx), ...
            nanmean(mat_exp(keep_idx,~keep_idx'),2),'rows', 'pairwise');
        eff_perm(kk) = corr(E1_spun(keep_idx), ...
            nanmean(mat_exp(keep_idx' & sch400_econ_types>0, keep_idx))','rows', 'pairwise');
        % type-specific connectivity
        aff_types_perm(kk,:) = corr(E1_spun(keep_idx), ...
            grpstats(mat_exp(keep_idx,~keep_idx')',sch400_econ_types(~keep_idx'))',...
            'rows', 'pairwise');
        eff_types_perm(kk,:) = corr(E1_spun(keep_idx), ...
            grpstats(mat_exp(~keep_idx',keep_idx),sch400_econ_types(~keep_idx'))',...
            'rows', 'pairwise');
    end
    [aff_r aff_types_r; ... 
        1 - (sum(aff_r < aff_perm)/length(perm_id)) 1 - (sum(aff_types_r < aff_types_perm)/length(perm_id)); ...
        eff_r eff_types_r; ... 
        1 - (sum(eff_r < eff_perm)/length(perm_id)) 1 - (sum(eff_types_r < eff_types_perm)/length(perm_id))]
    
    % comparison of correlation coefficients between afferent and efferent
    % connectivity
    [hyp,p,z] = mengz(eff_r, aff_r, ...
        corr(mean(mat(sch400_atlas7(sch400_cort)==8,sch400_atlas7(sch400_cort)~=8),2), ...
            mean(mat(sch400_atlas7(sch400_cort)~=8,sch400_atlas7(sch400_cort)==8)',2)), ...
            length(x))
        
    % k-fold validation of stability of correlations
    c = cvpartition(find(sch400_atlas7_cort==8),'KFold',10,'Stratify',false);
    r_cv = [];
    dmn_idx = find(sch400_atlas7(sch400_cort)==8);
    for ii = 1:10
       r_cv(ii,:,1) = corr(x(training(c,ii)==1), ...
           mat_aff_node_type(:,training(c,ii)==1)');
       r_cv(ii,:,2) = corr(x(training(c,ii)==1), ...
           mat_eff_node_type(:,training(c,ii)==1)');
    end
   
    for figure_time = 1
        
        f = figure('units','centimeters','outerposition',[0 0 20 20]);
        
        cmap_r = flipud(interp_colormap(colorbrewer.div.RdBu{1,9}/255,23));
              
        a(2) = axes('position', [0.33 0.8 0.15 0.15])
        x = E1_parc(sch400_atlas7==8);
        mdl = fitlm(x, mean(mat_aff_node_type));
        Xnew = linspace(min(x), max(x), 100)';
        [ypred,yci] = predict(mdl,Xnew);
        fill([Xnew; flipud(Xnew)]', [yci(:,1); flipud(yci(:,2))], [0.7 0.7 0.7], 'LineStyle','none')
        alpha 0.5
        hold on
        plot(Xnew, ypred,'LineWidth',1,'Color','k')  
        scatter(E1_parc(sch400_atlas7==8), mean(mat_aff_node_type), 5, 'k', ...
            'filled')
        xlim(a(2), [min(Xnew)-0.1 max(Xnew)+0.2])
        ylim(a(2), [min(mean(mat_aff_node_type))-0.0001 max(mean(mat_aff_node_type))+0.00015])

        a(3) = axes('position', [0.53 0.8 0.15 0.15])
        for ii = 1:6
            bar(ii, median(r_cv(:,ii,1)), 'FaceColor', batlow(ii*40,:))
            hold on
        end
        alpha 0.5
        for ii = 1:6
            plot([ii ii], [min(r_cv(:,ii,1)) max(r_cv(:,ii,1))], 'k')
        end
        xlim([0.5 6.5])
        ylim([-0.7 0.2])

        a(3) = axes('position', [0.33 0.6 0.15 0.15])
        mdl = fitlm(E1_parc(sch400_atlas7==8), mean(mat_eff_node_type));
        Xnew = linspace(min(E1_parc), max(E1_parc), 100)';
        [ypred,yci] = predict(mdl,Xnew);
        fill([Xnew; flipud(Xnew)]', [yci(:,1); flipud(yci(:,2))], [0.7 0.7 0.7], 'LineStyle','none')
        alpha 0.5
        hold on
        plot(Xnew, ypred,'LineWidth',1,'Color','k')  
        scatter(E1_parc(sch400_atlas7==8), mean(mat_eff_node_type), 5, 'k', ...
            'filled')
        xlim(a(3), [min(Xnew)-0.1 max(Xnew)+0.2])
        ylim(a(3), [min(mean(mat_eff_node_type))-0.0001 max(mean(mat_eff_node_type))+0.00015])
        
        a(3) = axes('position', [0.53 0.6 0.15 0.15])
        for ii = 1:6
            bar(ii, median(r_cv(:,ii,2)), 'FaceColor', batlow(ii*40,:))
            hold on
        end
        alpha 0.5
        for ii = 1:6
            plot([ii ii], [min(r_cv(:,ii,2)) max(r_cv(:,ii,2))], 'k')
        end
        xlim([0.5 6.5])
        ylim([-0.7 0.2])
        
        exportfigbo(f, [homeDir '/figures/dcm_global_glm_' lower(dataset) '.png'], 'png', 5)
    end
    
end

for imbalance_between_networks_and_types = 1
   
    % compile afferent (aff) and efferent (eff) connectivity for each
    % functional network
    mat_aff_net = []; mat_eff_net = [];
    for ii = 1:7
        net_idx = sch400_atlas7(sch400_cort)==ii+1;
        mat_aff_net(:,net_idx) = grpstats(mat(net_idx,~net_idx & sch400_econ_types_cort>0)', ...
            sch400_econ_types_cort(~net_idx & sch400_econ_types_cort>0));
        mat_eff_net(:,net_idx) = grpstats(mat(~net_idx & sch400_econ_types_cort>0,net_idx), ...
            sch400_econ_types_cort(~net_idx & sch400_econ_types_cort>0));
    end
    
    % divergence from balance across types
    mat_aff_net_mean = grpstats(mat_aff_net', sch400_atlas7(sch400_cort));
    mat_eff_net_mean = grpstats(mat_eff_net', sch400_atlas7(sch400_cort));
    for ii = 1:7
        pVecNull = ones(1,6)/6;
        pVecNet = mat_aff_net_mean(ii,:)/(sum(mat_aff_net_mean(ii,:)));
        kl(ii,1) = nansum(pVecNull .* (log2(pVecNull)-log2(pVecNet)));
        pVecNet = mat_eff_net_mean(ii,:)/(sum(mat_eff_net_mean(ii,:)));
        kl(ii,2) = nansum(pVecNull .* (log2(pVecNull)-log2(pVecNet)));
    end
    % permutation testing - spin testing
    kl_perm = [];
    for kk = 1:size(perm_id,2)
        perm_val_aff = grpstats(mat_aff_net', sch400_atlas7(perm_id(sch400_cort,kk)));
        perm_val_eff = grpstats(mat_eff_net', sch400_atlas7(perm_id(sch400_cort,kk)));
        for ii = 1:7
            pVecNull = ones(1,6)/6;
            pVecNet = perm_val_aff(ii,:)/(sum(perm_val_aff(ii,:)));
            kl_perm(ii,kk,1) = nansum(pVecNull .* (log2(pVecNull)-log2(pVecNet)));
            pVecNet = perm_val_eff(ii,:)/(sum(perm_val_eff(ii,:)));
            kl_perm(ii,kk,2) = nansum(pVecNull .* (log2(pVecNull)-log2(pVecNet)));
        end
    end
    kl_perm(isinf(kl_perm)) = nan;
    
    % key stats
    [kl(:,1) sum(kl(:,1)>kl_perm(:,:,1),2)/10000 ...
        kl(:,2) sum(kl(:,2)>kl_perm(:,:,2),2)/10000]'
    
end

for extended_model = 1
   
    load([homeDir '/data/' dataset '_rDCM_sch400_extended.mat'], 'results');
    intra_dmn = [zeros(1,14) sch400_atlas7(sch400_cort)==8]==1;
    sch400_econ_types_sc = [7:13 7:13 sch400_econ_types_cort];
    mat = results.mean_Amatrix_allSubjects;
    mat = abs(mat) .* ~eye(size(mat));
    
    % dmn only
    mat_aff_node_type = grpstats(mat(intra_dmn,~intra_dmn)', sch400_econ_types_sc(~intra_dmn));
    mat_eff_node_type = grpstats(mat(~intra_dmn,intra_dmn), sch400_econ_types_sc(~intra_dmn));
    
    % Axis and interaction effects
    x = E1_parc_cort(sch400_atlas7(sch400_cort)==8); % E1 within DMN
    aff_r = corr(x, mean(mat(intra_dmn,~intra_dmn),2));
    eff_r = corr(x, mean(mat(~intra_dmn,intra_dmn)',2));
    aff_types_r = corr(x, mat_aff_node_type');
    eff_types_r = corr(x, mat_eff_node_type');
    [aff_r aff_types_r(1:6); eff_r eff_types_r(1:6)]
    aff_perm = zeros(length(perm_id), length(aff_r));
    eff_perm = zeros(length(perm_id), length(aff_r));
    % recalculate with spun atlas
    mat_exp = nan(length(unique(sch400)),length(unique(sch400)));
    mat_exp(sch400_cort,sch400_cort) = mat(15:end,15:end);
    for kk = 1:length(perm_id)
        keep_idx = sch400_atlas7_spun(:,kk)==8;
        % average 
        aff_perm(kk) = corr(E1_spun(keep_idx), ...
            nanmean(mat_exp(keep_idx,~keep_idx'),2),'rows', 'pairwise');
        eff_perm(kk) = corr(E1_spun(keep_idx), ...
            nanmean(mat_exp(keep_idx' & sch400_econ_types>0, keep_idx))','rows', 'pairwise');
        % type-specific connectivity
        aff_types_perm(kk,:) = corr(E1_spun(keep_idx), ...
            grpstats(mat_exp(keep_idx,~keep_idx')',sch400_econ_types(~keep_idx'))',...
            'rows', 'pairwise');
        eff_types_perm(kk,:) = corr(E1_spun(keep_idx), ...
            grpstats(mat_exp(~keep_idx',keep_idx),sch400_econ_types(~keep_idx'))',...
            'rows', 'pairwise');
    end
    1 - (sum(aff_r < aff_perm)/length(perm_id))
    1 - (sum(eff_r < eff_perm)/length(perm_id))
    1 - (sum(aff_types_r(1:6) < aff_types_perm)/length(perm_id))
    1 - (sum(eff_types_r(1:6) < eff_types_perm)/length(perm_id))
    
    % divergence from balance across types
    mat_aff_net = []; 
    mat_eff_net = []; 
    for ii = 1:7
        net_idx = [zeros(1,14) (sch400_atlas7(sch400_cort)==ii+1)]==1;
        mat_aff_net(:,net_idx) = grpstats(mat(net_idx,~net_idx & (sch400_econ_types_sc<7))', ...
            sch400_econ_types_sc(~net_idx & (sch400_econ_types_sc<7)));
        mat_eff_net(:,net_idx) = grpstats(mat(~net_idx & sch400_econ_types_sc<7,net_idx), ...
            sch400_econ_types_sc(~net_idx & (sch400_econ_types_sc<7)));
    end
    mat_aff_net_mean = grpstats(mat_aff_net', [zeros(1,14) sch400_atlas7(sch400_cort)]);
    mat_eff_net_mean = grpstats(mat_eff_net', [zeros(1,14) sch400_atlas7(sch400_cort)]);
    for ii = 1:7
        pVecNull = ones(1,6)/6;
        pVecNet = mat_aff_net_mean(ii+1,:)/(sum(mat_aff_net_mean(ii+1,:)));
        kl(ii,1) = nansum(pVecNull .* (log2(pVecNull)-log2(pVecNet)));
        pVecNet = mat_eff_net_mean(ii+1,:)/(sum(mat_eff_net_mean(ii+1,:)));
        kl(ii,2) = nansum(pVecNull .* (log2(pVecNull)-log2(pVecNet)));
    end
    % permutation testing - spin testing
    kl_perm = [];
    for kk = 1:size(perm_id,2)
        atlas_spun = [zeros(1,14) sch400_atlas7(perm_id(sch400_cort,kk))];
        atlas_spun(atlas_spun<2) = nan; % ignore non-isocortex
        perm_val_aff = grpstats(mat_aff_net', atlas_spun);
        perm_val_eff = grpstats(mat_eff_net', atlas_spun);
        for ii = 1:7
            pVecNull = ones(1,6)/6;
            pVecNet = perm_val_aff(ii,:)/(sum(perm_val_aff(ii,:)));
            kl_perm(ii,kk,1) = nansum(pVecNull .* (log2(pVecNull)-log2(pVecNet)));
            pVecNet = perm_val_eff(ii,:)/(sum(perm_val_eff(ii,:)));
            kl_perm(ii,kk,2) = nansum(pVecNull .* (log2(pVecNull)-log2(pVecNet)));
        end
    end
    kl_perm(isinf(kl_perm)) = nan;
    
    % key stats
    [kl(:,1) sum(kl(:,1)>kl_perm(:,:,1),2)/10000 ...
        kl(:,2) sum(kl(:,2)>kl_perm(:,:,2),2)/10000]'
    
end