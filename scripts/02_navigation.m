%% 02_navigation

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
    sch400_atlas7_spin = sch400_atlas7(perm_id);
    
end

for load_data = 1
    
    dataset = 'mics';
    
    % mics tractography
    load([homeDir '/data/' dataset '_tractography_sch400.mat'], 'sc_group')
   
end

for navigation_calculation = 1
   
    % create distance matrix
    sch400_centroids = grpstats(FS.coord', sch400);
    D = squareform(pdist(sch400_centroids));
    
    % organise parcels and vertices by hemisphere
    hemi_id = [2:201; 203:402]; % excludes cortical wall
    hemi_id_nocort = [1:200; 201:400]; % when cortical wall is already excluded
    hemi_id_vert = [1:length(FS.coord)/2; (length(FS.coord)/2)+1:length(FS.coord)];
    Enav_whole = nan(400,400);
    Enav_net_grp = [];
    for h = 1:2 % for each hemisphere
        % navigation paths
        sc_intra = sc_group(hemi_id_nocort(h,:),hemi_id_nocort(h,:)); % connectivity within hemisphere
        L = -log10(sc_intra./(max(sc_intra(:)) + min(sc_intra(sc_intra>0)))); % log transform weights (see Seguin et al., 2019)
        L(isinf(L)) = 0;
        nav = navigate(L,D(hemi_id(h,:),hemi_id(h,:))); 
        
        % calculate navigation effiency
        Enav = 1./(nav.pl_MS);
        Enav(isinf(Enav)) = 0;
        net_avg = [];     
        
        % whole enav matrix
        Enav_whole(hemi_id_nocort(h,:),hemi_id_nocort(h,:)) = Enav;
        
        % network-average Enav for each node
        Enav_net_grp(:,hemi_id_nocort(h,:)) = grpstats(Enav,sch400_atlas7(hemi_id(h,:)));
        
        for jj = 1:7 % for each network
            % define the network
            net_idx = sch400_atlas7(hemi_id(h,:))==jj+1;
            % calculate average extrinsic Enav
            Enav_all(:,sch400_atlas7==jj+1 & ismember(1:402,hemi_id(h,:))) = mean(Enav(~net_idx,net_idx));
            % calculate type-average extrinsic Enav
            Enav_type_net(:,sch400_atlas7==jj+1 & ismember(1:402,hemi_id(h,:))) = ...
                grpstats(Enav(net_idx,~net_idx & sch400_econ_types(hemi_id(h,:))>0)',...
                sch400_econ_types(sch400_atlas7~=jj+1 & sch400_econ_types>0 & ismember(1:402,hemi_id(h,:))), 'mean');
            % calculate type-average Enav (intrinsic and extrinsic)
            Enav_type_net_all(:,sch400_atlas7==jj+1 & ismember(1:402,hemi_id(h,:))) = ...
                grpstats(Enav(net_idx,sch400_econ_types(hemi_id(h,:))>0)',...
                sch400_econ_types(sch400_econ_types>0 & ismember(1:402,hemi_id(h,:))), 'mean');
            % calculate type-average intrinsic Enav
            for ii = 1:6 % loops through in a different way because not all types will be represented within each network
                if nnz(net_idx & sch400_econ_types(hemi_id(h,:))==ii)
                    Enav_type_net_int(ii,sch400_atlas7==jj+1 & ismember(1:402,hemi_id(h,:))) = ...
                        mean(Enav(net_idx,net_idx & sch400_econ_types(hemi_id(h,:))==ii),2);
                end
            end
        end
    end
    save([homeDir '/output/enav_' dataset '.mat'], 'Enav_whole', 'Enav_net_avg', 'Enav_all', ...
        'Enav_type_net', 'Enav_type_net_all')

    % axis effect: relationship between average Enav and E1
    tbl = table(E1_parc(sch400_atlas7==8), ...
        Enav_type_net(:,sch400_atlas7==8)', ...
        Enav_type_net_all(:,sch400_atlas7==8)', ...
        Enav_type_net_int(:,sch400_atlas7==8)', ...
        'VariableNames', {'E1', 'econ_avg_nodmn', 'econ_avg', 'econ_avg_dmn'});
    [r,p] = corr(tbl.E1, Enav_all(sch400_atlas7==8)')
    
    % interaction: correlation of type-specific Enav with E1
    r = corr(tbl.E1, tbl.econ_avg_nodmn);
    
    % k-fold validation of stability of correlations
    c = cvpartition(find(sch400_atlas7==8),'KFold',10,'Stratify',false);
    r_cv = [];
    for ii = 1:10
       r_cv(ii,:) = corr(tbl.E1(training(c,ii)==1), tbl.econ_avg_nodmn(training(c,ii)==1,:));
    end
    
    % spin test
    % evaluates whether the correlation of type-specific Enav with E1 is
    % independent of the size/spatial autocorrelation of the DMN
    E1_parc_spun = E1_parc(perm_id); % going to spin E1, as well as the DMN (used as a mask)
    r_perm = zeros(size(perm_id,2), 6);
    for p = 1:size(sch400_atlas7_spin,2)
        r_perm(p,:) = corr(E1_parc_spun(sch400_atlas7_spin(:,p)==8,p), tbl.econ_avg_nodmn);
    end
    p = sum(r > squeeze(r_perm)) / size(perm_id,2);
    
end

for figure_Enav = 1
            
    f = figure('units','centimeters','outerposition',[0 0 20 20]);
    
    cmap_r = flipud(interp_colormap(colorbrewer.div.RdBu{1,9}/255,23));
    
    toMap = zeros(1,length(unique(sch400)));
    clim = [min(Enav_all(sch400_atlas7==8)) max(Enav_all(sch400_atlas7==8))];
    toMap(sch400_atlas7==8) = Enav_all(sch400_atlas7==8);
    OnSurf = BoSurfStatMakeParcelData(toMap, FS, sch400);
    BoSurfStat_calibrate4Views(OnSurf, FSinf, ...
        [0.05 0.86 0.1 0.1; 0.05 0.78 0.1 0.1; 0.15 0.78 0.1 0.1; 0.15 0.86 0.1 0.1], 1:4, ...
        clim, [1,1,1; lajolla])
    
    a(2) = axes('position', [0.33 0.8 0.15 0.15])
    mdl = fitlm(tbl.E1, Enav_all(sch400_atlas7==8)');
    Xnew = linspace(min(tbl.E1), max(tbl.E1), 100)';
    [ypred,yci] = predict(mdl,Xnew);
    fill([Xnew; flipud(Xnew)]', [yci(:,1); flipud(yci(:,2))], [0.7 0.7 0.7], 'LineStyle','none')
    alpha 0.5
    hold on
    plot(Xnew, ypred,'LineWidth',1,'Color','k')  
    scatter(tbl.E1, Enav_all(sch400_atlas7==8)', 5, 'k', ...
        'filled')
    colormap(a(2), lajolla)
    caxis(a(2), clim)
    xlim(a(2), [min(Xnew)-0.1 max(Xnew)+0.2])
    
    a(3) = axes('position', [0.53 0.8 0.15 0.15])
    for ii = 1:6
        bar(ii, median(r_cv(:,ii)), 'FaceColor', batlow(ii*40,:))
        hold on
    end
    alpha 0.5
    for ii = 1:6
        plot([ii ii], [median(r_cv(:,ii))-std(r_cv(:,ii)) median(r_cv(:,ii))+std(r_cv(:,ii))], ...
            'k')
    end
    xlim([0.5 6.5])
    ylim([-0.7 0.25])
        
    exportfigbo(f, [homeDir '/figures/enav_glm_' dataset '.png'], 'png', 5)
  
    
end

for navigation_calculation_all_networks = 1
   
    % calculate strong connectivity
    Enav_net_grp = grpstats(Enav_whole,sch400_atlas7(sch400_cort))';
    Enav_nt_grp = grpstats(Enav_net_grp, sch400_econ_types(sch400_cort));
    
    % divergence from balance across types - using network-average per
    % extrinsic node
    for ii = 1:7
        pVecNull = ones(1,6)/6;
        pVecNet = Enav_nt_grp(:,ii)'/(sum(Enav_nt_grp(:,ii)));
        kl(ii,1) = nansum(pVecNull .* (log2(pVecNull)-log2(pVecNet)));
    end
        
    % permutation testing - extending matrix
    Enav_whole_ext = nan(402,402);
    Enav_whole_ext(sch400_cort,sch400_cort) = Enav_whole;
    % spin testing
    kl_perm = [];
    for kk = 1:size(perm_id,2)
        disp(kk)
        perm_idx = perm_id(:,kk);
        perm_Enav_net_avg = grpstats(Enav_whole_ext,sch400_atlas7(perm_idx))';
        perm_Enav_nt_avg = grpstats(perm_Enav_net_avg, sch400_econ_types);
        for ii = 1:7
            pVecNull = ones(1,6)/6;
            pVecNet = perm_Enav_nt_avg(:,ii)'/(sum(perm_Enav_nt_avg(:,ii)));
            kl_perm(ii,kk,1) = nansum(pVecNull .* (log2(pVecNull)-log2(pVecNet)));
        end
    end

    % key stats
    [kl(:,1) sum(kl(:,1)>kl_perm(:,:,1),2)/10000]'
    
end
