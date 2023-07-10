% compute_functional_flow

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

    % for more details on construction of the rDCM model see the 
    % the nested preprocessing scripts [homeDir '/scripts/rDCM/'] or 
    % TAPAS (https://github.com/translationalneuromodeling/tapas)
    dataset = 'hcp';
    load([homeDir '/data/' dataset '_rDCM_sch400.mat'], 'results');
    mat = results.mean_Amatrix_allSubjects;
    mat = abs(mat) .* ~eye(size(mat));
    
    % compile afferent
    conn_extrinsic = mean(mat(sch400_atlas7_cort==8,sch400_atlas7_cort~=8),2);
    type_conn_extrinsic = grpstats(mat(sch400_atlas7_cort==8,sch400_atlas7_cort~=8)', ...
        sch400_econ_types_cort(sch400_atlas7_cort~=8));
    conn_type_avg = zeros(7,400);
    for ii = 1:7
        net_idx = sch400_atlas7(sch400_cort)==ii+1;
        utypes = unique(sch400_econ_types_cort(~net_idx));
        conn_type_avg(utypes+1,net_idx) = grpstats(mat(net_idx,~net_idx)', ...
            sch400_econ_types_cort(~net_idx));
    end
    conn_type_net_avg = grpstats(conn_type_avg', sch400_atlas7_cort)';
    conn_type_net_avg = conn_type_net_avg(2:end,:); % ignore any type=0 entries
    save([homeDir '/output/afferent_' dataset '.mat'], 'conn_extrinsic', 'type_conn_extrinsic', 'conn_type_net_avg')

    % compile efferent
    conn_extrinsic = mean(mat(sch400_atlas7_cort~=8,sch400_atlas7_cort==8))';
    type_conn_extrinsic = grpstats(mat(sch400_atlas7_cort~=8,sch400_atlas7_cort==8), ...
        sch400_econ_types_cort(sch400_atlas7_cort~=8));
    conn_type_avg = zeros(7,400);
    for ii = 1:7
        net_idx = sch400_atlas7(sch400_cort)==ii+1;
        utypes = unique(sch400_econ_types_cort(~net_idx));
        conn_type_avg(utypes+1,net_idx) = grpstats(mat(~net_idx,net_idx), ...
            sch400_econ_types_cort(~net_idx));
    end
    conn_type_net_avg = grpstats(conn_type_avg', sch400_atlas7_cort)';
    conn_type_net_avg = conn_type_net_avg(2:end,:); % ignore any type=0 entries
    save([homeDir '/output/efferent_' dataset '.mat'], 'conn_extrinsic', 'type_conn_extrinsic', 'conn_type_net_avg')
    
end


