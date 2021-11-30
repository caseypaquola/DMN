function get_rDCM_parameter_estimates_fixed_MainExperiment(dataset,model,scale)
% Get the individual parameter estimates from the whole-brain Dynamic Causal 
% Models (DCMs) for the resting state dataset of the HCP and MICS dataset.
% Parameter estimates have been computed using regression DCM (rDCM).
% 
% This function reads the DCM files that have been estimated using
% regression DCM. The function stores the individual parameter estimates
% for the sparse effective connectivity patterns.
% 
% Input:
%   dataset             -- dataset: (1) HCP, (2) MICS
%   model               -- model to analyze
%   scale               -- scaling of the data: (0) no scale, (1) scale
%   
% Output:
% 

% ----------------------------------------------------------------------
% 
% stefan_fraessle@gmx.de
%
% Author: Stefan Fraessle, TNU, UZH & ETHZ - 2016
% Copyright 2018 by Stefan Fraessle <stefan_fraessle@gmx.de>
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


% define the different names
scale_name      = {'_noScale',''};

% names of the different datasets
dataset_name = {'HCP','MICS'};


% asign the foldernames
pre_foldername = FilenameInfo.Data_Path;

% define the different names
foldername = fullfile(pre_foldername,dataset_name{dataset});

% get all subjects
Subject_List = dir(fullfile(foldername,'sub*'));


% define the filename
filename = ['DCM_RestingState_model' num2str(model) '_rDCM.mat'];

    
% iterate over all subjects
for subject = 1:length(Subject_List)

    % subject name
    Subject = Subject_List(subject).name;

    % set the directory
    foldername_rDCM = fullfile(foldername, Subject, 'firstlevel_dcm', ['regressionDCM' scale_name{scale+1} '_connectome']);

    % get the data
    try

        % load the DCM
        temp           = load(fullfile(foldername_rDCM,filename));
        output         = temp.output;

        % display subject name
        disp(['Subject: ' Subject ' (' dataset_name{dataset} ' - ' num2str(subject) ') - found'])

    catch

        % load a dummy result to get network size
        foldername_dummy = fullfile(foldername, 'sub_001', 'firstlevel_dcm', ['regressionDCM' scale_name{scale+1} '_connectome']);
        file_dummy = dir(fullfile(foldername_dummy, filename));
        temp = load(fullfile(foldername_dummy,file_dummy(1).name));

        % define dummy connectivity and driving input matrices
        output.Ep.A = NaN(size(temp.output.Ep.A,1),size(temp.output.Ep.A,2));
        output.Ep.C = NaN(size(temp.output.Ep.C,1),size(temp.output.Ep.C,2));
        
        % define a dummy array
        output.logF = NaN;

        % define a dummy time
        output.time_rDCM = NaN;

        % display subject name
        disp(['Subject: ' Subject ' (' dataset_name{dataset} ' - ' num2str(subject) ') - missing (wrong)'])

    end

    % clear the DCM file
    clear temp
    

    % define cells
    if ( subject == 1 )
        A_Matrix_AllSubjects    = cell(size(output.Ep.A,1),size(output.Ep.A,2));
        C_Matrix_AllSubjects    = cell(size(output.Ep.C,1),size(output.Ep.C,2));

        F_AllSubjects         	= cell(1,1);

        Time_rDCM_AllSubjects 	= cell(1,1);
    end

    % asign the a values for the endogenous parameters in each subject
    for int = 1:size(output.Ep.A,1)
        for int2 = 1:size(output.Ep.A,2)
            A_Matrix_AllSubjects{int,int2} = [A_Matrix_AllSubjects{int,int2}; output.Ep.A(int,int2)];
        end 
    end

    % asign the a values for the driving parameters in each subject
    for int = 1:size(output.Ep.C,1)
        for int2 = 1:size(output.Ep.C,2)
            C_Matrix_AllSubjects{int,int2} = [C_Matrix_AllSubjects{int,int2}; output.Ep.C(int,int2)];
        end 
    end


    % asign the negative free energy
    F_AllSubjects{1,1}           = [F_AllSubjects{1,1}; output.logF];


    % asign the run-time
    Time_rDCM_AllSubjects{1,1} = [Time_rDCM_AllSubjects{1,1}; output.time_rDCM];


    % asign the subject name
    if ( subject == 1 )
        results.AllSubjects{subject}        = Subject;
    else
        results.AllSubjects{end+1}          = Subject;
    end


    % asign the region names
    if ( subject == 1 )
        
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
       	
        % reorient the labels
        results.AllRegions = labels';
    end
end
    
% asign the results
results.A_Matrix_AllSubjects               = A_Matrix_AllSubjects;
results.C_Matrix_AllSubjects               = C_Matrix_AllSubjects;
results.F_AllSubjects                      = F_AllSubjects;
results.Time_rDCM_AllSubjects              = Time_rDCM_AllSubjects;


% display whether subject have to be re-run
A_Matrix_IncludedSubjects = A_Matrix_AllSubjects{int,int2};
disp(['Number of missing subjects: ' num2str(sum(isnan(A_Matrix_IncludedSubjects)))])

if ( sum(isnan(A_Matrix_IncludedSubjects)) > 0 )

    % find which subjects are missing
    temp = [];
    for int = 1:length(vector)
        if (vector(int)==1 && ~isfinite(results.A_Matrix_AllSubjects{1,1}(int)))
            temp = [temp, int];
        end
    end

    % output subject numbers that are missing
    fprintf(['Number of missing subjects: ' repmat('%d ',1,length(temp)) '\n'],temp)
    
end

% create the results folder
if ( ~exist(fullfile(FilenameInfo.Data_Path,dataset_name{dataset},'grouplevel_dcm',['regressionDCM' scale_name{scale+1} '_connectome']),'dir') )
    mkdir(fullfile(FilenameInfo.Data_Path,dataset_name{dataset},'grouplevel_dcm',['regressionDCM' scale_name{scale+1} '_connectome']))
end

% save the estimated result
if ( ~isempty(results) )
    save(fullfile(FilenameInfo.Data_Path,dataset_name{dataset},'grouplevel_dcm',['regressionDCM' scale_name{scale+1} '_connectome'], ['PosteriorParameterEstimates_model' num2str(model) '.mat']), 'results', '-v7.3')
end

end