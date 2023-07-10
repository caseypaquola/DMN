% compute_navigation_efficiency

for init_project = 1
   
    % set directories
    GH          = '/Users/cpaquola/Documents/GitHub/';       % github repo location
    homeDir     = [GH '/DMN'];                               % home directory
    bbwDir      = [GH '/BigBrainWarp'];                      % BigBrainWarp location
    
    % preset fsaverage5 subject directory from Freesurfer
    fsDir       = '/Applications/freesurfer/'; 
    fsAv        = [fsDir '/subjects/fsaverage5/'];           % fsaverage5 directory
    
    % addpaths of github repos - micaopen, BigBrainWarp, gifti toolbox,
    % munesoft, BrainSpace, toolbox_fast_marching
    addpath([homeDir '/scripts/']);                         % from project-specific GitHub                     
    addpath(genpath([GH '/micaopen/surfstat']));            % surfstat toolbox, housed in micaopen (from MICA Lab)         
    addpath(genpath([GH '/gifti']));                        % gifti toolbox

    % load surfaces
    FS = SurfStatAvSurf({[fsAv '/surf/lh.pial'], [fsAv '/surf/rh.pial']});    % fsaverage5 pial surface
    FSinf = SurfStatAvSurf({[fsAv '/surf/lh.inflated'], [fsAv '/surf/rh.inflated']});
    
end

for load_parcellation_and_atlases = 1

    % parcellation
    parc_name='Schaefer2018_400Parcels_7Networks_order';  % what parcellation the connectomes are in
    [sch400, ~] = annot2classes(...
        [homeDir '/utilities/tpl-fsaverage5/tpl-fsaverage5_hemi-L_' parc_name '.annot'],...
        [homeDir '/utilities/tpl-fsaverage5/tpl-fsaverage5_hemi-R_' parc_name '.annot'],1);
    sch400_cort = [2:201 203:402]; % exclude non-cortical parcels

    % load atlases
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

    % calculate modal atlas7 and econ type assignment for each parcel
    usch400 = unique(sch400);
    for ii = 1:length(usch400)
        sch400_atlas7(ii) = mode(atlas7_fs(sch400==usch400(ii)));
        sch400_econ_types(ii) = mode(econ_types_fs(sch400==usch400(ii)));
    end
    sch400_atlas7_cort = sch400_atlas7(sch400_cort);
    sch400_econ_types_cort = sch400_econ_types(sch400_cort);


end

for load_data = 1
    
    % input: matrix of tract-strength between each pair of parcels, derived from diffusion-weighted
    % imaging
    dataset = 'hcp';
    load([homeDir '/data/' dataset '_tractography_sch400.mat'], 'sc_group');

end

for compute = 1

    % create distance matrix
    sch400_centroids = grpstats(FS.coord', sch400);
    D = squareform(pdist(sch400_centroids));
    
    % organise parcels and vertices by hemisphere
    hemi_id = [2:201; 203:402]; % excludes cortical wall
    hemi_id_nocort = [1:200; 201:400]; % when cortical wall is already excluded
    hemi_id_vert = [1:length(FS.coord)/2; (length(FS.coord)/2)+1:length(FS.coord)];
    conn_mat = nan(400,400);
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
        conn_mat(hemi_id_nocort(h,:),hemi_id_nocort(h,:)) = Enav;
        
        for jj = 1:7 % for each network
            % define the network
            net_idx = sch400_atlas7(hemi_id(h,:))==jj+1;
            % calculate average extrinsic Enav
            conn_extrinsic(:,sch400_atlas7==jj+1 & ismember(1:402,hemi_id(h,:))) = mean(Enav(~net_idx,net_idx));
            % calculate type-average extrinsic Enav
            type_conn_extrinsic(:,sch400_atlas7==jj+1 & ismember(1:402,hemi_id(h,:))) = ...
                grpstats(Enav(net_idx,~net_idx & sch400_econ_types(hemi_id(h,:))>0)',...
                sch400_econ_types(sch400_atlas7~=jj+1 & sch400_econ_types>0 & ismember(1:402,hemi_id(h,:))), 'mean');
        end
    end
    conn_extrinsic = conn_extrinsic(sch400_atlas7_cort);
    type_conn_extrinsic = type_conn_extrinsic(:,sch400_atlas7_cort);

    % calculate average connectivity to each type for each network
    conn_net_avg = grpstats(conn_mat,sch400_atlas7_cort)';
    conn_type_net_avg = grpstats(conn_net_avg, sch400_econ_types(sch400_cort)); % note: first row is type=0 and should be ignored

    save([homeDir '/output/enav_' dataset '.mat'], 'conn_extrinsic', 'type_conn_extrinsic', 'conn_type_net_avg')

end
