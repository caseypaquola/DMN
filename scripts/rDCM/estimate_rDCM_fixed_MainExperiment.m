function estimate_rDCM_fixed_MainExperiment(dataset,subject,model,scale)
% Infer whole-brain effective connectivity for the resting-state dataset of
% the HCP and MICS dataset using regression DCM.
% 
% This function prepares the respective DCM model. Futher, it applies rDCM 
% to the model, which means assuming a fully (all-to-all) connected model,
% which is then pruned during model inversion to yield a sparse
% representation of network connectivity.
% 
% Input:
%   dataset             -- dataset: (1) HCP, (2) MICS
%   subject             -- subject to analyze
%   model               -- model to analyze
%   scale               -- scale BOLD signal time series
%   
% Output:
% 

% ----------------------------------------------------------------------
% 
% stefan_fraessle@gmx.de
%
% Author: Stefan Fraessle, TNU, UZH & ETHZ - 2021
% Copyright 2021 by Stefan Fraessle <stefan_fraessle@gmx.de>
%
% Licensed under GNU General Public License 3.0 or later.
% Some rights reserved. See COPYING, AUTHORS.
% 
% ----------------------------------------------------------------------


% get the path
m_path = mfilename('fullpath');
m_path = m_path(1:find(m_path=='/',1,'last'));

try
    load(fullfile(m_path,'ConfigFile.mat'))
catch err
    disp('Need to specify FilenameInfo containing paths!')
    rethrow(err)
end


% add the rDCM toolbox
addpath(genpath(fullfile(FilenameInfo.TAPAS_Path,'rDCM')))
addpath(FilenameInfo.SPM_Path)

% initialize spm
spm_jobman('initcfg')

% scaling or no scaling of data
scale_name = {'_noScale',''};

% names of the different datasets
dataset_name = {'HCP','MICS'};


% load the time series data
disp('Loading data...')
filename_temp = dir(fullfile(FilenameInfo.LongTermStorage_Path, dataset_name{dataset},'*rsfmri_ts_subcortical.mat'));
temp          = load(fullfile(FilenameInfo.LongTermStorage_Path, dataset_name{dataset}, filename_temp(1).name));
ts            = temp.ts;

% number of subjects
NrSub = size(ts,3);

% asign subcortical region names
labels{1,1} = '7Networks_LH_SubCort_Thalamus';
labels{2,1} = '7Networks_LH_SubCort_Caudate';
labels{3,1} = '7Networks_LH_SubCort_Putamen';
labels{4,1} = '7Networks_LH_SubCort_Pallidum';
labels{5,1} = '7Networks_LH_SubCort_Hippocampus';
labels{6,1} = '7Networks_LH_SubCort_Amygdala';
labels{7,1} = '7Networks_LH_SubCort_Accumbens';
labels{8,1} = '7Networks_RH_SubCort_Thalamus';
labels{9,1} = '7Networks_RH_SubCort_Caudate';
labels{10,1} = '7Networks_RH_SubCort_Putamen';
labels{11,1} = '7Networks_RH_SubCort_Pallidum';
labels{12,1} = '7Networks_RH_SubCort_Hippocampus';
labels{13,1} = '7Networks_RH_SubCort_Amygdala';
labels{14,1} = '7Networks_RH_SubCort_Accumbens';

% get the Schaefer region names
[~,labelsSchaefer,~,~,~,~] = textread(fullfile(FilenameInfo.Parcellation_path, 'Schaefer2018', 'Schaefer_labels.txt'),'%s%s%s%s%s%s');

% asign the Schaefer regions as well
for NrR = 1:length(labelsSchaefer)
    labels{NrR+14} = labelsSchaefer{NrR};
end


% asign the foldernames
pre_foldername = FilenameInfo.Data_Path;

% define the different names
foldername = fullfile(pre_foldername,dataset_name{dataset});


% define the filename
filename = ['DCM_RestingState_model' num2str(model)];



% choose which subjects to analyze
if ( isempty(subject) )
    subject_analyze = 1:NrSub;
else
    subject_analyze = subject;
end


% iterate overs subjects
for subject = subject_analyze
    
    % asign subject name
    numid = ['00' num2str(subject)];
    Subject = ['sub_' numid(end-2:end)];
    fprintf(['\nProcessing: ' Subject ' (' dataset_name{dataset} ') \n'])
    
    
    % clear relevant files
    clear DCM 
    clear output 
    clear options
    
    
    % define the standard option settings
    options.scale      = scale;
    options.estimateVL = 0;

    % asign the scale
    scale = options.scale;
    
    
    % display progress
    disp(['Found number of regions: ' num2str(size(ts(:,:,subject),2))])
    
    % number of regions
    DCM.n = size(ts(:,:,subject),2);

    % number of time points
    DCM.v = size(ts(:,:,subject),2);

    % specify the TR
    if ( dataset ==  1 )
        DCM.Y.dt  = 0.72;
    elseif ( dataset == 2 )
        DCM.Y.dt  = 0.6;
    end
    
    % specify the Y component of the DCM file
    DCM.Y.X0 = ones(DCM.v,1);

    % asign the data to the Y structure
    for i = 1:DCM.n
        DCM.Y.y(:,i)  = ts(:,i,subject);
        DCM.Y.name{i} = labels{i};
    end

    % define the covariance matrix
    DCM.Y.Q = spm_Ce(ones(1,DCM.n)*DCM.v);
    
    
    % specify empty driving input
    DCM.U.u     = zeros(size(DCM.Y.y,1)*16, 1);
    DCM.U.name  = {'null'};
    DCM.U.dt    = DCM.Y.dt/16;
    
    
    % DCM parameters
    DCM.delays = repmat(DCM.Y.dt/2,DCM.n,1);

    % asign echo time
    if ( dataset ==  1 )
        DCM.TE = 0.033;
    elseif ( dataset == 2 )
        DCM.TE = 0.03;
    end
    
    
    % DCM options
    DCM.options.nonlinear  = 0;
    DCM.options.two_state  = 0;
    DCM.options.stochastic = 0;
    DCM.options.nograph    = 0;


    % connectivity structure of respective model:
    %   Model 1: No driving inputs
    if ( model == 1 )
        DCM.a = ones(DCM.n);
        DCM.b = zeros(DCM.n,DCM.n,size(DCM.U.u,2));
        DCM.c = zeros(DCM.n,size(DCM.U.u,2));
        DCM.d = zeros(DCM.n,DCM.n,0);
    end


    % check whether the folder exists
    if ( ~exist(fullfile(foldername, Subject, 'firstlevel_dcm'),'dir') )
        mkdir(fullfile(foldername, Subject, 'firstlevel_dcm'))
    end


    % detrend and scale the data
    if ( options.scale )

        % detrend data
        DCM.Y.y = spm_detrend(DCM.Y.y);

        % scale data
        scale_factor   = max(max((DCM.Y.y))) - min(min((DCM.Y.y)));
        scale_factor   = 4/max(scale_factor,4);
        DCM.Y.y        = DCM.Y.y*scale_factor;

    end


    % set the prior type of the endogenous matrix
    DCM.options.wide_priors = 0;


    % display the progress
    disp(['Evaluating subject ' num2str(subject) ' - model ' num2str(model)])
    disp(' ')


    % shift the input regressor slightly
    options.u_shift = 0;
    options.filter_str = 4;
    
    % do not compute the signal for time reasons
    options.compute_signal = 0;


    % empirical analysis
    type = 'r';
    
    
    % get time
    currentTimer = tic;

    
    % run the rDCM analysis
    [output, options] = tapas_rdcm_estimate(DCM, type, options, 1);

    % output elapsed time
    time_rDCM = toc(currentTimer);
    disp(['Elapsed time is ' num2str(time_rDCM) ' seconds.'])


    % store the estimation time
    output.time_rDCM = time_rDCM;


    % create the directory
    if ( ~exist(fullfile(foldername, Subject, 'firstlevel_dcm', ['regressionDCM' scale_name{scale+1} '_connectome']),'dir') )
        mkdir(fullfile(foldername, Subject, 'firstlevel_dcm', ['regressionDCM' scale_name{scale+1} '_connectome']))
    end

    % save the estimated results
    if ( ~isempty(output) )
        save(fullfile(foldername, Subject, 'firstlevel_dcm', ['regressionDCM' scale_name{scale+1} '_connectome'], [filename '_rDCM.mat']),'output','-v7.3')
    end

    % save the options file
    if ( ~isempty(options) )
        save(fullfile(foldername, Subject, 'firstlevel_dcm', ['regressionDCM' scale_name{scale+1} '_connectome'], [filename '_options.mat']),'options')
    end

end

end
