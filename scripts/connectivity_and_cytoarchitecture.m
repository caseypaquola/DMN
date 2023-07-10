%% connectivity_and_cytoarchitecture

for init_project = 1
   
    % set directories
    GH          = '/Users/cpaquola/Documents/GitHub/';          % github repo location
    homeDir     = [GH '/DMN'];                               % home directory
    bbwDir      = [GH '/BigBrainWarp'];                      % BigBrainWarp location
    
    % preset fsaverage5 subject directory from Freesurfer
    fsDir       = '/Applications/freesurfer/'; 
    fsAv        = [fsDir '/subjects/fsaverage5/'];           % fsaverage5 directory
    
    % addpaths of github repos - micaopen, BigBrainWarp, gifti toolbox,
    % munesoft, BrainSpace, toolbox_fast_marching
    addpath([homeDir '/scripts/']);                         % from project-specific GitHub                     
    addpath(genpath([GH '/micaopen/surfstat']));            % surfstat toolbox, housed in micaopen (from MICA Lab)         
    addpath([bbwDir '/scripts/']);                          % BigBrainWarp code                     
    addpath(genpath([GH '/gifti']));                        % gifti toolbox
    addpath(genpath([GH '/BrainSpace/matlab/']));           % matlab component of the BrainSpace toolbox (from MICA Lab)
    
    % load colourmaps - using colorbrewer and scientific colours 
    cmap_types = [127 140 172; 139 167 176; 171 186 162; 218 198 153; 253 211 200; 252 229 252]/255;
    [~, ~, atlas7_ctb] = read_annotation([homeDir '/utilities/tpl-fsaverage5/tpl-fsaverage5_hemi-L_desc-Yeo_7Networks.annot']);
    cmap_atlas7 = atlas7_ctb.table(:,1:3)/255;
    cmap_dmn = [0.8 0.8 0.8; cmap_atlas7(end,:)];
    
    % load surfaces
    FS = SurfStatAvSurf({[fsAv '/surf/lh.pial'], [fsAv '/surf/rh.pial']});    % fsaverage5 pial surface
    FSinf = SurfStatAvSurf({[fsAv '/surf/lh.inflated'], [fsAv '/surf/rh.inflated']});
    
end

for load_parcellations_and_cytoarchitectonic_map = 1
      
    parc_name='Schaefer2018_400Parcels_7Networks_order';  % what parcellation the connectomes are in
    
    % load atlases on fsaverage5
    atlas7_fs = annot2classes([homeDir '/utilities/tpl-fsaverage5/tpl-fsaverage5_hemi-L_desc-Yeo_7Networks.annot'], ...
        [homeDir '/utilities/tpl-fsaverage5/tpl-fsaverage5_hemi-R_desc-Yeo_7Networks.annot'], 0); % orig. from CBIG github (https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation)
    atlas17_fs = annot2classes([homeDir '/utilities/tpl-fsaverage5/tpl-fsaverage5_hemi-L_desc-Yeo_17Networks.annot'], ...
        [homeDir '/utilities/tpl-fsaverage5/tpl-fsaverage5_hemi-R_desc-Yeo_17Networks.annot'], 0); % orig. from CBIG github
    econ_fs = annot2classes([homeDir '/utilities/tpl-fsaverage5/tpl-fsaverage5_hemi-L_desc-economo.annot'], ...
        [homeDir '/utilities/tpl-fsaverage5/tpl-fsaverage5_hemi-R_desc-economo.annot'],0); % from Dutch Connectome Lab (http://www.dutchconnectomelab.nl/)
    [~,~,econ_ctb] = read_annotation([homeDir '/utilities/tpl-fsaverage5/tpl-fsaverage5_hemi-L_desc-economo.annot']); 
    econ_ctb.types = [0;0;2;3;4;3;3;3;2;2;3;3;3;4;5;6;6;6;5;4;6;6;4;4;6;6;6;2;1;1;2;1;2;3;2;3;4;3;3;2;1;1;2;4;5]; % hard coded based on table data in Garcia-Cabezas (2021)
    econ_types_fs = BoSurfStatMakeParcelData(econ_ctb.types([1 3:45]),FS,econ_fs);
    dmn_fs = atlas7_fs==8;
    [sch400, ~] = annot2classes([fsAv '/label/lh.' parc_name '.annot'],...
        [fsAv '/label/rh.' parc_name '.annot'],1);
    sch400_cort = [2:201 203:402];
    hemi_id = [2:201; 203:402]; % excludes cortical wall
    hemi_id_nocort = [1:200; 201:400]; % when cortical wall is already excluded
   
    % calculate modal atlas7 and econ type assignment for each parcel
    usch400 = unique(sch400);
    for ii = 1:length(usch400)
        sch400_atlas7(ii) = mode(atlas7_fs(sch400==usch400(ii)));
        sch400_econ_types(ii) = mode(econ_types_fs(sch400==usch400(ii)));
    end
    sch400_atlas7_dmn = sch400_atlas7==8;
    sch400_atlas7_cort = sch400_atlas7(sch400_cort);
    sch400_atlas7_cort_dmn = sch400_atlas7_cort==8;
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
    sch400_bb = annot2classes([homeDir '/utilities/tpl-bigbrain/tpl-bigbrain_hemi-L_desc-Schaefer2018_400Parcels_7Networks_order.annot'], ...
        [homeDir '/utilities/tpl-bigbrain/tpl-bigbrain_hemi-R_desc-Schaefer2018_400Parcels_7Networks_order.annot'], 1);
    E1_parc = grpstats(E1, sch400_bb);
    E1_parc_cort = E1_parc(sch400_cort);
    E1_parc_dmn = E1_parc(sch400_atlas7_dmn);
      
    % spin the atlas
    run_spins = 0;
    if run_spins == 1
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
    
    % load data collated in "compute_navigation_efficiency" (tractography-based) or
    % "compute_functional_flow" (fMRI-based DCM)
    data_type = 'enav'; % can be enav, afferent or efferent
    dataset = 'mics'; % can be mics or hcp
    load([homeDir '/data/' data_type '_' dataset '.mat'], 'conn_mat')

end


for compare_conn_and_cyto = 1

    % AXIS EFFECT
    % relationship between E1 and average connectivity acros parcels
    r = corr(E1_parc_dmn, conn_extrinsic(sch400_atlas7_cort_dmn)')
    E1_parc_spun = E1_parc(perm_id); % going to spin E1, as well as the DMN (used as a mask)
    r_perm = zeros(size(perm_id,2), 1);
    for p = 1:size(sch400_atlas7_spun,2)
        r_perm(p,:) = corr(E1_parc_spun(sch400_atlas7_spun(:,p)==8,p), conn_extrinsic(sch400_atlas7_cort_dmn)');
    end
    spun_p = sum(r > squeeze(r_perm)) / size(perm_id,2)

    
    % INTERACTION WITH TYPE
    % correlation of E1 with type-specific connectivity 
    type_conn_extrinsic_dmn = type_conn_extrinsic(:,sch400_atlas7_cort_dmn)';
    r = corr(E1_parc_dmn, type_conn_extrinsic_dmn)
    % k-fold validation of stability of type-specific correlations
    c = cvpartition(find(sch400_atlas7_cort_dmn==1),'KFold',10,'Stratify',false);
    r_cv = [];
    for ii = 1:10
       r_cv(ii,:) = corr(E1_parc_dmn(training(c,ii)==1), type_conn_extrinsic_dmn(training(c,ii)==1,:));
    end
    % spin test
    % evaluates whether the correlation of type-specific Enav with E1 is
    % independent of the size/spatial autocorrelation of the DMN
    E1_parc_spun = E1_parc(perm_id); % going to spin E1, as well as the DMN (used as a mask)
    r_perm = zeros(size(perm_id,2), 6);
    for p = 1:size(sch400_atlas7_spun,2)
        r_perm(p,:) = corr(E1_parc_spun(sch400_atlas7_spun(:,p)==8,p), type_conn_extrinsic_dmn);
    end
    spun_p = sum(r > squeeze(r_perm)) / size(perm_id,2)

end

for balance_across_types = 1
   
    % export for ridgeplots in R
    tbl = table(type_conn_extrinsic(1,:)',type_conn_extrinsic(2,:)',...
        type_conn_extrinsic(3,:)',type_conn_extrinsic(4,:)',...
        type_conn_extrinsic(5,:)',type_conn_extrinsic(6,:)',...
        sch400_atlas7_cort', ...
        'VariableNames', {'type1', 'type2', 'type3', 'type4', 'type5', 'type6', 'net'});
    writetable(tbl, [homeDir '/output/' data_type '_all_networks_by_type.csv'])

    % divergence from balance across types
    for ii = 1:7
        pVecNull = ones(1,6)/6;
        pVecNet = type_net_avg(2:end,ii)'/(sum(type_net_avg(2:end,ii)));
        kl(ii,1) = nansum(pVecNull .* (log2(pVecNull)-log2(pVecNet)));
    end
        
    % permutation testing - extending matrix
    conn_whole_ext = nan(402,402);
    conn_whole_ext(sch400_cort,sch400_cort) = conn_mat;
    % spin testing
    kl_perm = [];
    for kk = 1:size(perm_id,2)
        disp(kk)
        perm_idx = perm_id(:,kk);
        perm_conn_net_avg = grpstats(conn_whole_ext,sch400_atlas7(perm_idx))';
        perm_type_net_avg = grpstats(perm_conn_net_avg, sch400_econ_types);
        for ii = 1:7
            pVecNull = ones(1,6)/6;
            pVecNet = perm_type_net_avg(2:end,ii)'/(sum(perm_type_net_avg(2:end,ii)));
            kl_perm(ii,kk,1) = nansum(pVecNull .* (log2(pVecNull)-log2(pVecNet)));
        end
    end
    writematrix([kl(:,1)'; kl_perm(:,:,1)'], [homeDir '/output/' data_type '_kl_perm.csv']); % for pretty plotting in R

    % key stats (divergence and spun p-value)
    [kl(:,1) sum(kl(:,1)>kl_perm(:,:,1),2)/10000]'
    
end
