%% 01_cytoarchitecture

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
    % freesurfer matlab toolbox, gifti toolbox for matlab and BrainSpace
    % (github)
    addpath(genpath([GH '/micaopen']));
    addpath([bbwDir '/scripts/']);
    addpath([homeDir '/scripts/']);
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
    BB = SurfStatAvSurf({[bbwDir '/spaces/tpl-bigbrain/tpl-bigbrain_hemi-L_desc-pial_sym.obj'], ...
     [bbwDir '/spaces/tpl-bigbrain/tpl-bigbrain_hemi-R_desc-pial_sym.obj']}); % BigBrain pial surface
    BBinf = SurfStatAvSurf({[homeDir '/utilities/tpl-bigbrain_hemi-L.inflated'], ...
     [homeDir '/utilities/tpl-bigbrain_hemi-R.inflated']});                   % inflated BigBrain surface
    BBinf.coord(1,1:end/2) = BBinf.coord(1,1:end/2) - 100;                    % need to shift for visualisation because mris_inflate (a freesurfer tool used to create this) changes the origin
    load([homeDir '/utilities/tpl-bigbrain_desc-wall.mat'], 'BB_wall');                      % predefined cortical wall for BigBrain
    FS = SurfStatAvSurf({[fsAv '/surf/lh.pial'], [fsAv '/surf/rh.pial']});
    
end

for load_parcellations = 1
      
    % load atlases on fsaverage5
    atlas7_fs = annot2classes([fsAv '/label/lh.Yeo2011_7Networks_N1000.annot'], ...
        [fsAv '/label/rh.Yeo2011_7Networks_N1000.annot'], 0); % from CBIG github (https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation)
    econ_fs = annot2classes([fsAv '/label/lh.economo.annot'],[fsAv '/label/rh.economo.annot'],0); % from Dutch Connectome Lab (http://www.dutchconnectomelab.nl/)
    [~,~,econ_ctb] = read_annotation([fsAv '/label/lh.economo.annot']); 
    econ_ctb.types = [0;0;2;3;4;3;3;3;2;2;3;3;3;4;5;6;6;6;5;4;6;6;4;4;6;6;6;2;1;1;2;1;2;3;2;3;4;3;3;2;1;1;2;4;5]; % hard coded based on table data in Garcia-Cabezas (2021)
    econ_types_fs = BoSurfStatMakeParcelData(econ_ctb.types([1 3:45]),FS,econ_fs);
    
    % load atlases on bigbrain
    parc_name='Yeo2011_17Networks_N1000';  
    lh = gifti([bbwDir '/spaces/tpl-bigbrain/tpl-bigbrain_hemi-L_desc-' parc_name '.label.gii']);
    rh = gifti([bbwDir '/spaces/tpl-bigbrain/tpl-bigbrain_hemi-R_desc-' parc_name '.label.gii']);
    atlas17_bb = [lh.cdata; rh.cdata];
    parc_name='Yeo2011_7Networks_N1000';  
    lh = gifti([bbwDir '/spaces/tpl-bigbrain/tpl-bigbrain_hemi-L_desc-' parc_name '.label.gii']);
    rh = gifti([bbwDir '/spaces/tpl-bigbrain/tpl-bigbrain_hemi-R_desc-' parc_name '.label.gii']);
    atlas7_bb = [lh.cdata; rh.cdata]; 
    econ_bb = annot2classes([homeDir '/utilities/tpl-bigbrain_hemi-L_desc-economo.annot'], ...
        [homeDir '/utilities/tpl-bigbrain_hemi-R_desc-economo.annot'], 0);
    econ_types_bb = BoSurfStatMakeParcelData(econ_ctb.types([1 3:45]),BB,econ_bb);
    
end

for overlap_analysis_on_fsaverage5 = 1
    
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
        [~,p(ii),ks(ii)] = kstest2(econ_types_fs(atlas7_fs==ii+1), econ_types_fs(atlas7_fs==8));
    end
    
    % proportion in each type (within dmn)
    econ_chi_perm = [];
    [econ_prop,econ_chi,p] = crosstab(econ_types_fs,atlas7_fs==8); % chi squared test for equivalence of proportions
    parfor p = 1:nperm  % takes a a minute or two
        [econ_prop_perm(:,:,p),econ_chi_perm(p),~] = crosstab(econ_types_fs,atlas7_fs(perm_id(:,p))==8);
    end
    econ_prop_p(:,1) = 1 - (sum(econ_prop(:,2)'<squeeze(econ_prop_perm(:,2,:))')/nperm);
    econ_prop_p(:,2) = 1 - (sum(econ_prop(:,2)'>squeeze(econ_prop_perm(:,2,:))')/nperm);
    
    for figure_types = 1
        
        f = figure('units','centimeters','outerposition',[0 0 20 20]);
        
        % DMN on bigbrain
        BoSurfStat_calibrate4Views(double(atlas7_fs==8), FS, ...
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
        freq_dmn = tabulate(econ_types_fs(atlas7_fs==8));
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
    [BB_MPC, BB_MP] = build_mpc(BB_MP_vert(:,atlas7_bb==7), nn_bb(atlas7_bb==7)); % downsample to 100k mesh simultaneously

    if run_embedding == 1
        % diffusion map embedding
        BB_MPC = (BB_MPC + BB_MPC.')/2;
        BB_MPC(BB_MPC<0) = 0;
        [embedding, results] = mica_diffusionEmbedding(BB_MPC,'nComponents',10); % expect to run for 15minutes
        save([homeDir '/output/bigbrain_embedding_100k_thresh.mat'], 'embedding', 'results')
    else
        load([homeDir '/output/bigbrain_embedding_100k_thresh.mat'], 'embedding', 'results')
    end
    
    for figure_datadriven_axis = 1
        
        % diffusion embedding based
        ii = 1;
        f = figure('units','centimeters','outerposition',[0 0 20 20]);
        cmap_grad = roma;
        cmap_grad((length(roma)/2)-2:(length(roma)/2)+2,:) = repmat([0.7 0.7 0.7],5,1);
        
        % first gradient
        toMap = zeros(1,length(unique(nn_bb)));
        toMap(ismember(unique(nn_bb),unique(nn_bb(atlas7_bb==7)))) = embedding(:,ii);
        E1 = BoSurfStatMakeParcelData(toMap, BB, nn_bb) .* BB_wall;
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
        vert_to_cmap = ceil(rescale(embedding(:,ii), 1, length(flipud(roma))));
        nbins = 100;
        grad_bins = discretize(embedding(:,ii),nbins);
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
        exportfigbo(f, [homeDir '/figures/bigbrain_gradient_' num2str(g) '_upclose.png'], 'png', 5)
        
    end
    
    % write out for bigbrainwarping
    writematrix(E1(1:end/2)', [homeDir '/output/tpl-bigbrain_hemi-L_desc-DMN.txt'])
    writematrix(E1((end/2)+1:end)', [homeDir '/output/tpl-bigbrain_hemi-R_desc-DMN.txt'])
    
end

for data_driven_axis_by_types = 1
    
    % downsample types map to BB100
    unn_bb = unique(nn_bb(atlas7_bb==7));
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
    toMap(ismember(unique(nn_bb),unique(nn_bb(atlas7_bb==7)))) = type_rescaled;
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
   
    % get subregion clusters
    rand_data = rand(10, length(BB.coord));
    slm = SurfStatLinMod(rand_data,1,BB);
    slm = SurfStatT(slm, ones(10,1));
    slm.t = double((atlas7_bb==7)');
    [peak, clus, clusid] = SurfStatPeakClus(slm, ones(1,length(atlas7_bb==7)), 0.5);
    figure; SurfStatViewData(clusid, BBinf); colormap([0.7 0.7 0.7; parula]); caxis([0 12]) % inspect the ordering of subregions
    clusid(clusid>11) = 0; % remove small clusters
    % split dpfc and vpfc on left hemisphere
    clusid(clusid==1 & BBinf.coord(1,:)<-95 & BBinf.coord(2,:)<145 & BBinf.coord(3,:)<20) = 12;
    figure; SurfStatViewData(clusid, BBinf); colormap([0.7 0.7 0.7; parula]); caxis([0 12]) % inspect the ordering of subregions
    % note: SurfStat orders clusters by size so the order will be maintained when rerun
    clus_names = {'left dpfc', 'right dpfc', 'left ltg', 'right ltg', 'left pmc',...
        'left ipl', 'right pmc', 'right ipl', 'left mtl', 'right vlpfc', 'right mtl', 'left vlpfc'};
    sr_names = {'mtl', 'ltg', 'ipl', 'pmc', 'vlpfc', 'dpfc'};
    clusid_groups = zeros(size(clusid));
    for ii = 1:length(sr_names)
        sr_groups(contains(clus_names, sr_names{ii})) = ii;
        clusid_groups(ismember(clusid,find(sr_groups==ii))) = ii;
    end
    lh_reg = find(contains(clus_names, 'left'));
    tbl = table(E1(clusid_groups~=0)', clusid_groups(clusid_groups~=0)', ...
        'VariableNames', {'E1', 'reg'});
    writetable(tbl, [homeDir '/output/bigbrain_embedding_regions.csv']);

    % extract surface mesh of subregion and calculate robustness parametesr
    % also make flatmap of each subregion for visualisation
    edg = SurfStatEdg(BB);
    for sr = 2:12
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
        all_sr_surf{sr} = sr_surf;
        
        % calculate roughness parameters
        param{sr} = roughness_parameters_surface(sr_surf, E1(sr_list)');
        
        % isomap flattening
        % get geodesic distance between points on subregion mesh
        % (http://www.numerical-tours.com/matlab/meshdeform_3_flattening/)
        n = length(sr_surf.coord);
        D = zeros(n);
        % REQUIRES https://github.com/gpeyre/matlab-toolboxes/tree/master/toolbox_fast_marching
        % must be mex compiled on local system
        parfor i=1:n
            disp(i)
            D(:,i) = perform_fast_marching_mesh(sr_surf.coord,sr_surf.tri,i);
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
    save([homeDir '/output/subregion_parameters.mat'], 'vertexF', 'param')
    
    % visualise 3D landscape for left hemisphere regions
    f = figure('units','centimeters','outerposition',[0 0 20 20]);
    for s = 1:length(lh_reg)
        sr = lh_reg(s);
        disp(sr)
        
        % get subregion mesh details
        sr_list = find(clusid==sr); % list vertices in the subregion
        sr_list = unique(BB.tri(sum(ismember(BB.tri,sr_list),2)==3,:)); % remove vertices that aren't fully triangulated with the subregion
        sr_surf = all_sr_surf{sr};
        
        % extract embedding for the subregion
        dmn_parc = [];
        for ii = 1:length(sr_list)
            dmn_parc(ii) = find(unique(nn_bb(atlas7_bb==7))==nn_bb(sr_list(ii)));
        end
        reg_embed = embedding(dmn_parc,1);
        
        % plot flat map
        a(sr) = axes('position', [0.05 1-(s*0.16) 0.15 0.15]);
        trisurf(sr_surf.tri, vertexF{sr}(1,:), vertexF{sr}(2,:), reg_embed, ...
            reg_embed, 'EdgeColor', 'none');
        view(-90, 90)
        xlim([min(min(vertexF{sr}')) max(max(vertexF{sr}'))])
        ylim([min(min(vertexF{sr}')) max(max(vertexF{sr}'))])
        colormap(a(sr), flip(roma))
        caxis(a(sr), clim)
        axis off
        
        % plot 3D landscape
        a(1) = axes('position', [0.25 1-(s*0.166) 0.15 0.15]);
        [X,Y] = meshgrid([min(vertexF{sr}(1,:)) max(vertexF{sr}(1,:))],...
            [min(vertexF{sr}(2,:)) max(vertexF{sr}(2,:))]);
        Z = zeros(size(X));
        azSurf = surf(X,Y,Z);
        set(azSurf, 'FaceAlpha', 0.5);
        set(azSurf, 'FaceColor', [0.9 0.9 0.9]);
        hold on
        tSurf = trisurf(sr_surf.tri, vertexF{sr}(1,:), vertexF{sr}(2,:), reg_embed, ...
            reg_embed, 'EdgeColor', 'none');
        view(-18, 31)
        colormap(a(1), flip(roma))
        caxis(a(1), clim)
        xlim([min(vertexF{sr}(1,:)) max(vertexF{sr}(1,:))]);
        ylim([min(vertexF{sr}(2,:)) max(vertexF{sr}(2,:))]);  
    end
    exportfigbo(f, [homeDir '/figures/subregion_landscapes.png'], 'png', 4)
    
    % collate roughness parameters across subregions
    param_all = [];    
    todo = ones(1,length(clus_names));
    %for ii = 1:length(clus_names) % normal run
    for ii = [9 11 6 8 5 7 3 4 10 12 1 2] % alternate: ordered according to figure 1    
        param_all(:,ii) = table2array(param{ii})';
    end
    param{1}.Properties.VariableNames % parameter name
    cov = std(param_all,[],2)./mean(param_all,2); % coefficient of variation
    
end
