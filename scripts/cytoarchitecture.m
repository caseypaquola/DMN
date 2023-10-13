%% cytoarchitecture

% author: Casey Paquola
% email: casey.paquola@gmail.com
% feel free to get in touch if you have any questions

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
    addpath([fsDir '/matlab/'])
    addpath(genpath([GH '/micaopen/surfstat']));            % surfstat toolbox, housed in micaopen (from MICA Lab)         
    addpath(genpath([GH '/BrainStat/brainstat_matlab/']));
    addpath(genpath([GH '/BrainSpace/matlab/']));
    addpath([bbwDir '/scripts/']);                          % BigBrainWarp code                     
    addpath(genpath([GH '/gifti']));                        % gifti toolbox
    addpath(genpath([GH '/toolbox_fast_marching']));        % fast marching toolbox (caution: compiling may not work on Apple Silicon chips)
    % also requires fs_LR_32k github repo
    
    % load colourmaps - using colorbrewer and scientific colours 
    load([GH '/munesoft/colourmaps/colorbrewer.mat'])
    scicol = [GH '/munesoft/colourmaps/ScientificColourMaps7'];
    scicol_names = ["devon", "lajolla", "lapaz", "batlow", "roma", "vik", "berlin", "bam", "cork"];
    for ii = 1:length(scicol_names)
        load(strcat(scicol, '/', scicol_names(ii), '.mat'));    
    end
    cmap_types = [127 140 172; 139 167 176; 171 186 162; 218 198 153; 253 211 200; 252 229 252]/255;
    [~, ~, atlas7_ctb] = read_annotation([homeDir '/utilities/tpl-fsaverage5/tpl-fsaverage5_hemi-L_desc-Yeo_7Networks.annot']);
    cmap_atlas7 = atlas7_ctb.table(:,1:3)/255;
    cmap_dmn = [0.8 0.8 0.8; cmap_atlas7(end,:)];
    
    % load surfaces
    BB = SurfStatAvSurf({[bbwDir '/spaces/tpl-bigbrain/tpl-bigbrain_hemi-L_desc-pial_sym.obj'], ...
     [bbwDir '/spaces/tpl-bigbrain/tpl-bigbrain_hemi-R_desc-pial_sym.obj']}); % BigBrain pial surface
    BBinf = SurfStatAvSurf({[homeDir '/utilities//tpl-bigbrain/tpl-bigbrain_hemi-L.inflated'], ...
     [homeDir '/utilities//tpl-bigbrain/tpl-bigbrain_hemi-R.inflated']});                   % inflated BigBrain surface
    BBinf.coord(1,1:end/2) = BBinf.coord(1,1:end/2) - 100;                    % need to inflated BigBrain surface shift for visualisation because mris_inflate (a freesurfer tool used to create this) changes the origin
    load([homeDir '/utilities//tpl-bigbrain/tpl-bigbrain_desc-wall.mat'], 'BB_wall');       % predefined cortical wall for BigBrain
    FS = SurfStatAvSurf({[fsAv '/surf/lh.pial'], [fsAv '/surf/rh.pial']});    % fsaverage5 pial surface
    FSinf = SurfStatAvSurf({[fsAv '/surf/lh.inflated'], [fsAv '/surf/rh.inflated']});
    
end

for load_parcellations = 1
      
    % load atlases on fsaverage5
    atlas7_fs = annot2classes([homeDir '/utilities/tpl-fsaverage5/tpl-fsaverage5_hemi-L_desc-Yeo_7Networks.annot'], ...
        [homeDir '/utilities/tpl-fsaverage5/tpl-fsaverage5_hemi-R_desc-Yeo_7Networks.annot'], 0); % orig. from CBIG github (https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation)
    econ_fs = annot2classes([homeDir '/utilities/tpl-fsaverage5/tpl-fsaverage5_hemi-L_desc-economo.annot'], ...
        [homeDir '/utilities/tpl-fsaverage5/tpl-fsaverage5_hemi-R_desc-economo.annot'],0); % from Dutch Connectome Lab (http://www.dutchconnectomelab.nl/)
    [~,~,econ_ctb] = read_annotation([homeDir '/utilities/tpl-fsaverage5/tpl-fsaverage5_hemi-L_desc-economo.annot']); 
    econ_ctb.types = [0;0;2;3;4;3;3;3;2;2;3;3;3;4;5;6;6;6;5;4;6;6;4;4;6;6;6;2;1;1;2;1;2;3;2;3;4;3;3;2;1;1;2;4;5]; % hard coded based on table data in Garcia-Cabezas (2021)
    econ_types_fs = BoSurfStatMakeParcelData(econ_ctb.types([1 3:45]),FS,econ_fs);
    dmn_fs = atlas7_fs==8;
    
    % load atlases on bigbrain
    parc_name='Yeo2011_17Networks_N1000';  
    lh = gifti([bbwDir '/spaces/tpl-bigbrain/tpl-bigbrain_hemi-L_desc-' parc_name '.label.gii']);
    rh = gifti([bbwDir '/spaces/tpl-bigbrain/tpl-bigbrain_hemi-R_desc-' parc_name '.label.gii']);
    atlas17_bb = [lh.cdata; rh.cdata];
    parc_name='Yeo2011_7Networks_N1000';  
    lh = gifti([bbwDir '/spaces/tpl-bigbrain/tpl-bigbrain_hemi-L_desc-' parc_name '.label.gii']);
    rh = gifti([bbwDir '/spaces/tpl-bigbrain/tpl-bigbrain_hemi-R_desc-' parc_name '.label.gii']);
    atlas7_bb = [lh.cdata; rh.cdata]; 
    econ_bb = annot2classes([homeDir '/utilities/tpl-bigbrain/tpl-bigbrain_hemi-L_desc-economo.annot'], ...
        [homeDir '/utilities/tpl-bigbrain/tpl-bigbrain_hemi-R_desc-economo.annot'], 0);
    econ_types_bb = BoSurfStatMakeParcelData(econ_ctb.types([1 3:45]),BB,econ_bb);
    dmn_bb = atlas7_bb==7;
    
end

for overlap_of_types_and_networks = 1
    
    % spin the fsaverage5 sphere
    % creates nperm random rotations of vertex assignements
    % takes a few minutes to run
    nperm = 10000;
    spheres{1} = SurfStatReadSurf([fsAv '/surf/lh.sphere']);
    spheres{2} = SurfStatReadSurf([fsAv '/surf/rh.sphere']);
    datas{1} = [1:length(spheres{1}.coord)]'; datas{2} = [1:length(spheres{2}.coord)]';
    perm_id = spin_permutations(datas,spheres,nperm); % takes a few minutes
    perm_id = [squeeze(perm_id{1}); squeeze(perm_id{2})+length(spheres{1}.coord)];
    
    % proportion in each type (all networks)
    econ_prop = crosstab(econ_types_fs,atlas7_fs); % chi squared test for equivalence of proportions
    supp_table_1 = (econ_prop(2:end,2:end)./repmat(sum(econ_prop(2:end,2:end)),6,1))'
    for ii = 1:7
        [~,econ_chi(ii)] = crosstab(econ_types_fs,atlas7_fs==ii+1); % chi squared test for equivalence of proportions
        [~,p(ii),ks(ii)] = kstest2(econ_types_fs(atlas7_fs==ii+1), econ_types_fs(dmn_fs)); % do distributions differ between dmn and other networks?
    end
    
    % proportion in each type (within dmn)
    econ_prop_p = [];
    for ii = 1:7
        parfor p = 1:nperm  % takes a a minute or two
            [econ_prop_perm(:,:,p),econ_chi_perm(p),~] = crosstab(econ_types_fs,atlas7_fs(perm_id(:,p))==ii+1);
        end
        [econ_prop,econ_chi] = crosstab(econ_types_fs,atlas7_fs==ii+1); % chi squared test for equivalence of proportions
        prop_zeros(:,ii) = mean(squeeze(econ_prop_perm(:,2,:))==0,2);
        econ_prop_p(:,ii) = 1 - (sum(econ_prop(:,2)<squeeze(econ_prop_perm(:,2,:)),2)/nperm);
    end
        
    for figure_types = 1
        
        f = figure('units','centimeters','outerposition',[0 0 20 20]);
        
        % DMN on bigbrain
        BoSurfStat_calibrate4Views(double(dmn_fs), FS, ...
            [0.05 0.77 0.15 0.15; 0.05 0.65 0.15 0.15; ...
            0.2 0.65 0.15 0.15; 0.2 0.77 0.15 0.15], ...
            1:4, [0 1], cmap_dmn)
        
        % cortical types
        OnSurf = econ_types_fs;
        OnSurf(isnan(OnSurf)) = 0;
        BoSurfStat_calibrate2Views(OnSurf, FS, ...
            [0.35 0.65 0.15 0.15], [0.35 0.77 0.15 0.15], ...
            1, 2, [0 max(econ_types_fs)], [0.8 0.8 0.8; cmap_types])
        a(11) = axes('position', [0.365 0.64 0.12 0.01]);
        imagesc(1:6); axis off
        colormap(a(11), cmap_types)
        a(1) = axes('position', [0.54 0.73 0.13 0.17]);
        freq_dmn = tabulate(econ_types_fs(dmn_fs));
        hb = bar(1:6, freq_dmn(2:end,3), 'FaceColor',  cmap_atlas7(end,:), 'BarWidth', 0.9);
        hold on
        xlim([0.3 6.7])
        ylim([0 max(freq_dmn(:,3))+2])
      
        exportfigbo(f, [homeDir '/figures/figure1a.png'], 'png', 5)
        
    end
        
end

for data_driven_axis = 1
    
    run_downsample = 0;
    run_embedding = 0;
    run_alternative_dim_red = 1;
    
    if run_downsample == 1
        % downsample BB using mesh decimation
        numFaces    = 100000;
        patchNormal = patch('Faces', BB.tri, 'Vertices', BB.coord.','Visible','off');
        Sds         = reducepatch(patchNormal,numFaces);
        [~,bb_downsample]  = intersect(patchNormal.Vertices,Sds.vertices,'rows');
        BB100.tri     = double(Sds.faces);
        BB100.coord   = Sds.vertices';

        % For each vertex on BB, find nearest neighbour on BB100, via mesh neighbours
        nn_bb = zeros(1,length(BB.coord));
        edg = SurfStatEdg(BB);
        parfor ii = 1:length(BB.coord)
            nei = unique(edg(sum(edg==ii,2)==1,:));
            if isempty(nei) && ismember(ii,bb_downsample)
                nn_bb(ii) = ii;
            else
                while sum(ismember(nei, bb_downsample))==0
                    nei = [unique(edg(sum(ismember(edg,nei),2)==1,:)); nei];
                end
                matched_vert = nei(ismember(nei, bb_downsample));
                if length(matched_vert)>1  % choose the mesh neighbour that is closest in Euclidean space
                    n1 = length(matched_vert);
                    d = sqrt(sum((repmat(BB.coord(1:3,ii),1,n1) - BB.coord(:,matched_vert)).^2));
                    [~, idx] = min(d);
                    nn_bb(ii) = matched_vert(idx);
                else
                    nn_bb(ii) = matched_vert;
                end
            end
        end
        save([homeDir '/output/downsample_bigbrain_100k.mat'], 'nn_bb', 'BB100')
    else
        load([homeDir '/output/downsample_bigbrain_100k.mat'], 'nn_bb', 'BB100')
    end

    % load microstructure profiles of BigBrain, select DMN and perform
    BB_MP_vert = reshape(dlmread([bbwDir '/spaces/tpl-bigbrain/tpl-bigbrain_desc-profiles.txt']),[],50)';
    BB_MP_vert = (BB_MP_vert*-1) + max(BB_MP_vert(:));  % Invert values so high values ~ more staining
    [BB_MPC, BB_MP] = build_mpc(BB_MP_vert(:,dmn_bb), nn_bb(dmn_bb)); % downsample to 100k mesh simultaneously

    if run_embedding == 1
        % diffusion map embedding
        BB_MPC = (BB_MPC + BB_MPC.')/2; % symmetricise
        BB_MPC(BB_MPC<0) = 0;
        [embedding, results] = mica_diffusionEmbedding(BB_MPC,'nComponents',10); % expect to run for 15minutes
        save([homeDir '/output/bigbrain_embedding_100k_thresh.mat'], 'embedding', 'results')
    else
        load([homeDir '/output/bigbrain_embedding_100k_thresh.mat'], 'embedding', 'results')
    end

    if run_alternative_dim_red == 1
        approaches = {'pca', 'le', 'dm'};
        thresh = [50:10:90];
        for ii = 1:numel(approaches)
            embedding_alt = GradientMaps('kernel','na','approach',approaches{ii});
            for jj = 1:numel(thres)
                embedding_alt = embedding_alt.fit(BB_MPC, 'sparsity', thresh(jj));
                alt_r(ii,jj) = corr(embedding(:,1), embedding_alt.gradients{1}(:,1));
                alt_lambda(ii,jj) = embedding_alt.gradients.
            end
        end
    end

    % project back to the cortical surface of BB
    toMap = zeros(1,length(unique(nn_bb)));
    toMap(ismember(unique(nn_bb),unique(nn_bb(dmn_bb)))) = embedding(:,1);
    E1 = BoSurfStatMakeParcelData(toMap, BB, nn_bb) .* BB_wall;
        
    for figure_datadriven_axis = 1
        
        % diffusion embedding based
        f = figure('units','centimeters','outerposition',[0 0 20 20]);
        cmap_grad = roma;
        cmap_grad((length(roma)/2)-2:(length(roma)/2)+2,:) = repmat([0.7 0.7 0.7],5,1);
        
        % first gradient
        clim = [-prctile(abs(E1),99) prctile(abs(E1),99)];
        BoSurfStat_calibrate4ViewsNoShade(E1, BBinf, ...
            [0.05 0.77 0.15 0.15; 0.05 0.65 0.15 0.15; ...
            0.2 0.65 0.15 0.15; 0.2 0.77 0.15 0.15], ...
            1:4, clim, cmap_grad)
        a(1) = axes('position', [0.15 0.65 0.1 0.01]);
        imagesc(1:100); axis off
        colormap(a(1), flipud(roma))
        
        % profiles
        a(4) = axes('position', [0.55 0.65 0.16 0.25]);
        vert_to_cmap = ceil(rescale(embedding(:,1), 1, length(flipud(roma))));
        nbins = 100;
        grad_bins = discretize(embedding(:,1),nbins);
        for ii = 1:nbins
            if sum(grad_bins==ii) > 1
                this_colour = roma(round(mean(vert_to_cmap(grad_bins==ii))),:);
            else
                this_colour = roma(vert_to_cmap(grad_bins==ii),:);
            end
            plot(mean(BB_MP(:,grad_bins==ii),2), flip(1:50), 'Color', this_colour, 'LineWidth', 0.2)
            hold on
        end
        xlim([5000 25000])
        ylim([1 51])
        
        exportfigbo(f, [homeDir '/figures/bigbrain_E1.png'], 'png', 5)
        
        f = figure('units','centimeters','outerposition',[0 0 20 20]);
        SurfStatViewData(E1,BBinf)
        colormap(cmap_grad)
        SurfStatColLim(clim)
        exportfigbo(f, [homeDir '/figures/bigbrain_E1_upclose.png'], 'png', 5)
        
    end
    
    % write out for bigbrainwarping
    writematrix(E1(1:end/2)', [homeDir '/output/tpl-bigbrain_hemi-L_desc-DMN.txt'])
    writematrix(E1((end/2)+1:end)', [homeDir '/output/tpl-bigbrain_hemi-R_desc-DMN.txt'])
    
end

for data_driven_axis_by_types = 1
    
    % downsample types map to BB100
    unn_bb = unique(nn_bb(dmn_bb));
    for ii = 1:length(unn_bb)
        econ_types_bb_dmn(ii) = mode(econ_types_bb(nn_bb==unn_bb(ii)));
    end
    
    % difference between types and eigenvector
    % scaling types to E1, taking into account proporitions of each type
    vert_steps = double(range(embedding(econ_types_bb_dmn'>0,1))) / sum(econ_types_bb_dmn>0);
    type_step = [];
    for ii = 1:6
        type_step(ii) = (sum(econ_types_bb_dmn==ii)*vert_steps);
    end
    type_val = cumsum(type_step) + min(embedding(econ_types_bb_dmn'>0,1));
    % reassign values
    type_rescaled = econ_types_bb_dmn;
    for ii = 1:6
        type_rescaled(type_rescaled==ii) = type_val(ii);
    end
    
    f = figure('units','centimeters','outerposition',[0 0 20 20]);
    
    % show rescaled types on surface
    toMap = zeros(1,length(unique(nn_bb)));
    toMap(ismember(unique(nn_bb),unique(nn_bb(dmn_bb)))) = type_rescaled;
    TOnSurf = BoSurfStatMakeParcelData(toMap, BB, nn_bb) .* BB_wall;
    BoSurfStat_calibrate4ViewsNoShade(TOnSurf, BBinf, ...
        [0.55 0.77 0.15 0.15; 0.55 0.65 0.15 0.15; ...
        0.7 0.65 0.15 0.15; 0.7 0.77 0.15 0.15], ...
        1:4, clim, flip(roma))
    a(1) = axes('position', [0.625 0.63 0.1 0.01])
    imagesc(unique(type_rescaled)); axis off
    colormap(a(1), flip(roma))
    caxis(a(1), clim)
    
    % subtract from eigenvector
    cmap_dev = flipud(bam);
    BoSurfStat_calibrate4ViewsNoShade(E1-TOnSurf, BBinf, ...
        [0.05 0.47 0.15 0.15; 0.05 0.35 0.15 0.15; ...
        0.2 0.35 0.15 0.15; 0.2 0.47 0.15 0.15], ...
        1:4, [max(abs(E1-TOnSurf))*-1 max(abs(E1-E1))], cmap_dev)
    
    exportfigbo(f, [homeDir '/figures/bigbrain_gradient_types.png'], 'png', 5)
       
    % writing out to R for raincloud plot
    writetable(table(embedding(econ_types_bb_dmn'>0,1), ...
        econ_types_bb_dmn(econ_types_bb_dmn'>0)', ...
        'VariableNames', {'E1', 'type'}), ...
        [homeDir '/output/bigbrain_embedding_type.csv'])
    
end

for subregion_analysis = 1
   
    run_subregioning = 0;
    if run_subregioning == 1
        % get spatially contiguous subregion clusters
        rand_data = rand(10, length(BB.coord));
        slm = SurfStatLinMod(rand_data,1,BB);
        slm = SurfStatT(slm, ones(10,1));
        slm.t = double((dmn_bb)');
        [peak, clus, clusid] = SurfStatPeakClus(slm, ones(1,length(dmn_bb)), 0.5);
        figure; SurfStatViewData(clusid, BBinf); colormap([0.7 0.7 0.7; parula]); caxis([0 12]) % inspect the ordering of subregions
        clusid(clusid>11) = 0; % remove small clusters
        % split dpfc and vpfc on left hemisphere
        clusid(clusid==1 & BBinf.coord(1,:)<-95 & BBinf.coord(2,:)<145 & BBinf.coord(3,:)<20) = 12;
        figure; SurfStatViewData(clusid, BBinf); colormap([0.7 0.7 0.7; parula]); caxis([0 12]) % inspect the ordering of subregions
        % note: SurfStat orders clusters by size so the order will be maintained when rerun
        clus_names = {'left dpfc', 'right dpfc', 'left ltg', 'right ltg', 'left pmc',...
            'left ipl', 'right pmc', 'right ipl', 'right vlpfc', 'left mtl', 'right mtl', 'left vlpfc'};
        sr_names = {'mtl','pmc',  'ipl', 'ltg', 'vlpfc', 'dpfc'};
        clusid_groups = zeros(size(clusid));
        for ii = 1:length(sr_names)
            sr_groups(contains(clus_names, sr_names{ii})) = ii;
            clusid_groups(ismember(clusid,find(sr_groups==ii))) = ii;
        end
        lh_reg = find(contains(clus_names, 'left'));
        tbl = table(E1(clusid_groups~=0)', clusid_groups(clusid_groups~=0)', ...
            'VariableNames', {'E1', 'reg'});
        writetable(tbl, [homeDir '/output/bigbrain_embedding_regions.csv']);
        save([homeDir '/output/bigbrain_subregions.mat'], 'clusid', 'clus_names', 'sr_names', 'clusid_groups', 'sr_groups')
    
        % extract surface mesh of subregion and calculate robustness parametesr
        % also make flatmap of each subregion for visualisation
        edg = SurfStatEdg(BB);
        pk_thresh = prctile(E1,95);
        vl_thresh = prctile(E1,5);
        sm_kern = [2 4 6 8];
        E1_contrast = zeros(size(E1s));
        for sr = 1:12
            disp('next subregion')
            % extract surface mesh of subregion
            % convert region to 2D image
            sr_list = find(clusid==sr); % list vertices in the subregion
            sr_list = unique(BB.tri(sum(ismember(BB.tri,sr_list),2)==3,:)); % remove vertices that aren't fully triangulated with the subregion
            % create a graph from the mesh of SFG
            sr_faces = BB.tri(sum(ismember(BB.tri,sr_list),2)==3,:);
            sr_surf = [];
            sr_surf.tri = zeros(size(sr_faces));
            for ii = 1:length(sr_list)
               sr_surf.coord(:,ii) = BBinf.coord(:,sr_list(ii));
               sr_surf.coord_inf(:,ii) = BBinf.coord(:,sr_list(ii));
               sr_surf.tri(sr_faces==sr_list(ii)) = ii;
            end
            all_sr_list{sr} = sr_list;
            all_sr_surf{sr} = sr_surf;
    
            % calculate roughness parameters
            param{sr} = roughness_parameters_surface(sr_surf, E1(sr_list), pk_thresh, vl_thresh);

            % isomap flattening
            % get geodesic distance between points on subregion mesh
            % (http://www.numerical-tours.com/matlab/meshdeform_3_flattening/)
            n = length(sr_surf.coord);
            D = zeros(n);
            % REQUIRES https://github.com/gpeyre/matlab-toolboxes/tree/master/toolbox_fast_marching
            % must be mex compiled on local system (caution: unlikely to work on
            % Apple Silicon chips)
            parfor ii=1:n
                if rem(ii,100)==0
                    disp(ii)
                end
                D(:,ii) = perform_fast_marching_mesh(sr_surf.coord,sr_surf.tri,ii);
            end
            D(isinf(D)) = max(max(D(~isinf(D))))+10; % replace inf with a very high number
            D = (D+D')/2; % symmetrise
            GD{sr} = D;
            % centred matrix
            J = eye(n) - ones(n)/n;
            W = -J*(D.^2)*J; 
            % diagonalise
            [U,S] = eig(W);
            S = diag(S);
            [S,I] = sort(S,'descend');
            U = U(:,I);
            %isomap 
            vertexF{sr} = U(:,1:2)' .* repmat(sqrt(S(1:2)), [1 n]);
        end
        save([homeDir '/output/subregion_parameters.mat'], 'vertexF', 'all_sr_list', 'all_sr_surf', 'param')
    else
        load([homeDir '/output/bigbrain_subregions.mat'], 'clusid', 'clus_names', 'sr_names', 'clusid_groups', 'sr_groups')
        load([homeDir '/output/subregion_parameters.mat'], 'vertexF', 'all_sr_list', 'all_sr_surf', 'param')
    end
    
    % smoothness: fit of E1 to spatial axes
    for sr = 1:12
        sr_list = all_sr_list{sr}; %  vertices in the subregion
        for ii = 1:4 % tested with different basis function complexities (linear, quadratic, cubic, quartic)
            mdl = polyfitn(vertexF{sr}', E1(sr_list)',ii);
            R2(sr,ii) = mdl.AdjustedR2;
            RMSE(sr,ii) = mdl.RMSE;
        end
    end
    tbl = table(R2(:,2), R2(:,3), R2(:,4), sr_groups', ...
        'VariableNames', {'ord2', 'ord3', 'ord4', 'sr'});
    writetable(tbl, [homeDir '/output/subregion_polyfit.csv']); % export for pretty plotting in R
    % anova test for differences between subregions in smoothness
    for ii = 1:4 % iterating through model complexities
        [p(ii),tbl,stats] = anova1(R2(:,ii), sr_groups', 'off'); 
        F(ii) = tbl{2,5};
    end
   
    % waviness: roughness parameter
    for ii = 1:length(clus_names)
        Wf(ii) = param{ii}.Wf; % note: waviness factor is the last roughness parameter
    end
    tbl = table(Wf', sr_groups', ...
        'VariableNames', {'param', 'sr'}); 
    writetable(tbl, [homeDir '/output/subregion_Wf.csv']); % export for pretty plotting in R
    [p,tbl,stats] = anova1(mean(param_all(end,:)',2), sr_groups', 'off');  % anova test for differences between subregions in waviness
        
    % visualise 3D landscape for both hemispheres
    f = figure('units','centimeters','outerposition',[0 0 20 20]);
    hemi_reg = [find(contains(clus_names, 'left'));find(contains(clus_names, 'right'))];
    for s = 1:length(hemi_reg)
        for ii = 1:2
            sr = hemi_reg(ii,s);
            disp(sr)
            sr_pos = sr_groups(sr);
            
            % get subregion mesh details
            sr_list = all_sr_list{sr}; %  vertices in the subregion
            sr_surf = all_sr_surf{sr};
            reg_embed = E1(sr_list) * -1;
          
            % plot flat map
            a(sr) = axes('position', [(sr_pos*0.166)-0.166 0.6-(ii*0.2) 0.15 0.15]);
            trisurf(sr_surf.tri, vertexF{sr}(1,:), vertexF{sr}(2,:), reg_embed, ...
                reg_embed, 'EdgeColor', 'none');
            view(-90, 90)
            xlim([min(min(vertexF{sr}')) max(max(vertexF{sr}'))])
            ylim([min(min(vertexF{sr}')) max(max(vertexF{sr}'))])
            colormap(a(sr), flipud(roma))
            caxis(a(sr), clim)
            axis off
            
            % plot 3D landscape
            a(1) = axes('position', [(sr_pos*0.166)-0.15 1-(ii*0.2) 0.145 0.145]);
            [X,Y] = meshgrid([min(vertexF{sr}(1,:)) max(vertexF{sr}(1,:))],...
                [min(vertexF{sr}(2,:)) max(vertexF{sr}(2,:))]);
            Z = zeros(size(X));
            azSurf = surf(X,Y,Z);
            set(azSurf, 'FaceAlpha', 0.5);
            set(azSurf, 'FaceColor', [0.9 0.9 0.9]);
            hold on
            tSurf = trisurf(sr_surf.tri, vertexF{sr}(1,:), vertexF{sr}(2,:), reg_embed, ...
                reg_embed, 'EdgeColor', 'none');
            view(168, 20)
            colormap(a(1), flipud(roma))
            caxis(a(1), clim)
            xlim([min(vertexF{sr}(1,:)) max(vertexF{sr}(1,:))]);
            ylim([min(vertexF{sr}(2,:)) max(vertexF{sr}(2,:))]);  
        end
    end
    exportfigbo(f, [homeDir '/figures/subregion_landscapes.png'], 'png', 4)
    
end

for functional_decoding_types = 1


    % neurosynth decoding of subregions of each type
    clear top_features
    for ii = 1:6
        rand_data = rand(10, length(FS.coord));
        slm = SurfStatLinMod(rand_data,1,FS);
        slm = SurfStatT(slm, ones(10,1));
        slm.t = double(econ_types_fs==ii);
        [peak, clus, clusid] = SurfStatPeakClus(slm, ones(1,length(slm.t)), 0.5);
        lh_clus = unique(clusid(1:end/2));
        lh_clus(lh_clus==0) = [];
      
        for jj = 1:length(lh_clus)
            [correlation, feature] = meta_analytic_decoder(double(clusid==lh_clus(jj)), ...
                'template', 'fsaverage5');
            top_features{ii,jj,1} = feature(correlation>2)';
            top_features{ii,jj,2} = correlation(correlation>2);
        end
    end

    % manually create select_features that excludes non-psych terms
    save([homeDir '/output/neurosynth_select_features.mat'], 'select_features')

    % 
    for ii = 1:6
        for jj = 1:size(select_features,2)
            if jj == 1
                top_features_type{ii,1} = top_features{ii,jj,1}(select_features{ii,jj});
                top_features_type{ii,2} = top_features{ii,jj,2}(select_features{ii,jj});
                top_features_type{ii,3} = [ones(length(select_features{ii,jj}),1)*jj];
            else
                top_features_type{ii} = [top_features_type{ii,1}; top_features{ii,jj,1}(select_features{ii,jj})];
                top_features_type{ii,2} = [top_features_type{ii,2}; top_features{ii,jj,2}(select_features{ii,jj})];
                top_features_type{ii,3} = [top_features_type{ii,3}; [ones(length(select_features{ii,jj}),1)*jj]];
            end
        end
    end


    f = figure('units','centimeters','outerposition',[0 0 20 20]);

    for ii = 1:6

        BoSurfStat_calibrate4ViewsNoShade(double(econ_types_fs==ii), FSinf, ...
        [(ii*0.167)-0.167 0.87 0.15 0.15; (ii*0.167)-0.167 0.76 0.15 0.15; ...
        0.05 0.8 0 0; 0.05 0.8 0 0], ...
        1:4, [0 1], [1 1 1; cmap_types(ii,:)])

        if ii > 3
            a(ii) = axes('position', [(ii*0.167)-0.167 0.6 0.15 0.15]);
            wc1 = wordcloud(top_features_type{ii,1},top_features_type{ii,2});
            wc1.HighlightColor = cmap_types(ii,:);
            wc1.Color = cmap_types(ii,:);
        else
            for jj = 1:length(top_features_type{ii})
                a(ii) = axes('position', [(ii*0.167)-0.167 0.75-(jj*0.15) 0.15 0.15]);
                wc1 = wordcloud(top_features_type{ii,1}(top_features_type{ii,3}==jj),...
                    top_features_type{ii,2}(top_features_type{ii,3}==jj));
                wc1.HighlightColor = cmap_types(ii,:);
                wc1.Color = cmap_types(ii,:);
            end
        end
    end
    exportfigbo(f, [homeDir '/figures/neurosynth_type_clouds.png'], 'png', 10)


end