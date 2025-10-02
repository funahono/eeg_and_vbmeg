clear
close all

% path
vbmeg_toolbox_dir = '/home/honoka/vbmeg3_0_0_a_2';
addpath_dir = {'func_show', 'func_tf_analysis', 'func_vbmeg', 'func_else', 'func_post_vbmeg', 'func_fmri'};
% check if file exist
vbmeg_init_dir_filename = fullfile(vbmeg_toolbox_dir, 'vbmeg.m');
if exist(vbmeg_init_dir_filename, 'file')~=2
    warning(['cannot find  :', vbmeg_init_dir_filename])
end

% change dir
cd(vbmeg_toolbox_dir)
vbmeg;


%% step 1 : set parameter
% make structure
p = struct;
p_log = struct;

% dir
p.read.eeg_dir        = '/media/honoka/HDD2/MATLAB/EEGMat/practice';
p.read.freesurfer_dir = '/media/honoka/HDD1/Funatsuki/mni_icbm152_t1_tal_nlin_asym_09c_fs';
p.read.mri_dir        = fullfile(p.read.freesurfer_dir, 'mni_icbm152_nlin_asym_09c');

p.save.dir            = '/media/honoka/HDD2/MATLAB/vbmeg_analysis';

% file comment
p.read.eeg_file_comment = '__1_4_7';
p.save.dir_comment = 'ear_ref_car_standard_brain';

% parameter
day_sub = '20250313_B93';
which_finger = 'Rindex';

% eeg
p.divide_num  = 4;
p.resamp_flag = true;
p.car_flag    = false; % common average reference

% fMRI, MRI, Neurosynth
p.read.mri_filename = 'mni_icbm152_t1_tal_nlin_asym_09c';
p.fmri_flag         = false; % fMRI data setting
p.fmri_meta_flag    = false; % neurosynth setting

% estimate source current
p.Tperiod         = 100;         % time windows for calculating source current [ms]
p.Next            = 50;          % time step for calculating source current [ms]
p.time_noise      = [-3 -1.5];   % Period for baseline [s]
p.prior_weight    = 0.0001;      % Relative influence of prior information (0-1)
p.variance_reduce = [1];

p.tf_trial_average     = false;
p.curr_trial_save_flag = true;

%% step 2 : decide project name

% for read eeg
p.read.eeg_filename = [day_sub, '_', which_finger, p.read.eeg_file_comment, '.eeg.mat'];
p.read.eeg_dir_filename = fullfile(p.read.eeg_dir, day_sub, p.read.eeg_filename);
disp('-- read dir & file name (eeg) --------------------------------------');
disp(p.read.eeg_dir_filename);

% for read mri data
p.read.mri_dir_filename = fullfile(p.read.mri_dir, [p.read.mri_filename, '.nii']);
disp('-- read dir & file name (mri) ---------------------------------------');
disp(p.read.mri_dir_filename);

% for save current
% get today time
today = datetime;
today.Format = 'yyyyMMdd';
today = char(today);

% save dir of all
p.save.dirname = [day_sub,'_', which_finger, '_', today, '_', p.save.dir_comment];
p.save.dir = fullfile(p.save.dir, p.save.dirname);
disp('-- saving dir -------------------------------------------------------');
disp(p.save.dir);

%eeg dir & filename
p.save.eeg_dir                = fullfile(p.save.dir, 'eeg');
p.save.eeg_dir_filename       = fullfile(p.save.eeg_dir, p.read.eeg_filename);

% brain dir & filename
p.save.brain_dir_filename     = fullfile(p.save.dir, 'brain', [p.read.mri_filename, '.brain.mat']);
p.save.area_dir_filename      = fullfile(p.save.dir, 'brain', [p.read.mri_filename, '.area.mat']);
p.save.act_dir_filename       = fullfile(p.save.dir, 'brain', [p.read.mri_filename, '.act.mat']);

% other dir
p.save.neurosynth_dir         = fullfile(p.save.dir, 'neurosynth');
p.save.fig_dir                = fullfile(p.save.dir, 'figure');
p.save.bayes_dir              = fullfile(p.save.dir, 'bayes');
p.save.current_dir            = fullfile(p.save.dir, 'current');

% other dir & filename
p.save.head_dir_filename      = fullfile(p.save.dir, 'head', [p.read.mri_filename, '.head.mat']);
p.save.basis_dir_filename     = fullfile(p.save.dir, 'basis', [p.save.dirname, '.basis.mat']);


%% step 3 : add path & check file exist

% make save dir
if isdir(p.save.dir)
    cd(p.save.dir)
    error(['this "p.save.dir" is already exist : ', p.save.dir])
else
    mkdir(p.save.dir)
end

% make eeg dir
mkdir(p.save.eeg_dir)

% move read eeg file
if ~isfile(p.read.eeg_dir_filename) % 2 — name is a file with extension .m, .mlx, or .mlapp, or name is the name of a file with a non-registered file extension (.mat, .fig, .txt).
    error(['cannot find this "p.read.eeg_dir_filename" : ', p.read.eeg_dir_filename])
end

copyfile(p.read.eeg_dir_filename, [p.save.eeg_dir,'/']);

% check mri data already exist or not
if ~isfile(p.read.mri_dir_filename)
    error(['cannot find  :  ', p.read.mri_dir_filename])
end

% check fmri data already exist or not
if p.fmri_meta_flag && p.fmri_flag % when both flag is true
    error('only one of fmri data must be choosen')

elseif p.fmri_flag % when only fmri_flag is true
    error('please setting fmri analysis')

elseif p.fmri_meta_flag % when only fmri_meta_flag is true
    error('please setting fmri analysis')
end


%% step 4 : construct brain model (.brain.mat, .act.mat, .area.mat)
% set parameter
brain_param = vb_set_brain_parm(); % Set parameters for leadfield calculation

brain_param.FS_left_file            = fullfile(p.read.freesurfer_dir, 'bem', 'lh.smoothwm.asc');
brain_param.FS_right_file           = fullfile(p.read.freesurfer_dir, 'bem', 'rh.smoothwm.asc');

brain_param.FS_left_infl_file       = fullfile(p.read.freesurfer_dir, 'bem', 'lh.inflated.asc');
brain_param.FS_right_infl_file      = fullfile(p.read.freesurfer_dir, 'bem', 'rh.inflated.asc');

brain_param.FS_left_curv_file       = fullfile(p.read.freesurfer_dir, 'bem', 'lh.curv.asc');
brain_param.FS_right_curv_file      = fullfile(p.read.freesurfer_dir, 'bem', 'rh.curv.asc');

brain_param.FS_sphere_key = 'sphere.reg';

brain_param.FS_left_sphere_file     = fullfile(p.read.freesurfer_dir, 'bem', ['lh.sphere.reg.asc']);
brain_param.FS_right_sphere_file    = fullfile(p.read.freesurfer_dir, 'bem', ['rh.sphere.reg.asc']);

brain_param.FS_left_label_file      = fullfile(p.read.freesurfer_dir, 'label', 'lh.cortex.label');
brain_param.FS_right_label_file     = fullfile(p.read.freesurfer_dir, 'label', 'rh.cortex.label');

brain_param.registration_mode       = 'FS'; % 'FS' : Freesurder sphere file os used

brain_param.analyze_file            = p.read.mri_dir_filename;
brain_param.brain_file              = p.save.brain_dir_filename;
brain_param.area_file               = p.save.area_dir_filename;
brain_param.act_file                = p.save.act_dir_filename;

vb_job_brain(brain_param); % Convert FreeSurfer data to VBMEG format.

%{
function [V,F,xx,BV_index,Vinfo] =vb_job_brain(varargin)

Standard brain information.
V         : Cortical vertex point cordinate (SPM_Right_m)
F         : Patch index structure
xx        : Normal vector to cortical surface
xxA       : area asigned for each vertex ([m^2])

BV_index  : original vertex index in BV/FS corresponding to brain
Vinfo     : Vertex dimension structure
%}

%% step 5 : save parameter
p.save.param_dir_filename = fullfile(p.save.dir, 'param.mat'); %%ここがおかしそう
%'/media/honoka/HDD2/MATLAB/vbmeg_analysis/20250313_B93_Rindex_20250814_ear_ref_car_standard_brain/media/honoka/HDD2/MATLAB/vbmeg_analysis' 
save(p.save.param_dir_filename, 'p')

%% step 6 : import fmri
% いつかかきます

%% step 7 : eeg preprocessing(.eeg.mat)
if p.resamp_flag || p.car_flag % p.resamoleflag or p.car_flag = true
    disp('preprocessing eeg')

    load(p.read.eeg_dir_filename, 'eeg_data', 'EEGinfo', 'Measurement')

    if p.resamp_flag
        fs_old = EEGinfo.SampleFrequency;
        [n_chan, n_tsample, n_trial] = size(eeg_data);
        eeg_data_dsamp = NaN(n_chan, round(n_tsample/p.divide_num)+1, n_trial);

        for now_chan_num = 1:n_chan % each channnel
            for now_trial_num = 1:n_trial % each trial
                eeg_data_dsamp(now_chan_num, :, now_trial_num) = ...
                decimate(eeg_data(now_chan_num, :, now_trial_num),p.divide_num); % y = decimate(idata,r,varargin) %resample
            end
        end

        EEGinfo.SampleFrequency = fs_old/p.divide_num;
        EEGinfo.Pretrigger      = round(EEGinfo.Pretrigger/p.divide_num)+1;
        EEGinfo.Nsample         = size(eeg_data_dsamp,2);
        eeg_data = eeg_data_dsamp;

    end

    if p.car_flag
        n_chan = size(eeg_data,1);
        channel_average = mean(eeg_data,1);
        eeg_data = eeg_data - repmat(channel_average,n_chan,1,1);
    end

    save(p.save.eeg_dir_filename, 'eeg_data', "EEGinfo", 'Measurement')

end

%% step 8 : head model (.head.mat) & leadfield (.basis.mat)
% construct 3-shell surface model
% detail : prepare_leadfield_eeg.m

% Head model
% set input files
head_param.brain_file       = p.save.brain_dir_filename;
% Absolute path
head_param.analyze_file     = p.read.mri_dir_filename;
head_param.fs_scalp_file    = fullfile(p.read.freesurfer_dir, 'bem', 'outer_skin_surface.asc');
head_param.fs_scalp_file    = fullfile(p.read.freesurfer_dir, 'bem', 'inner_skull_surface.asc');
head_param.gray_file        = fullfile(p.read.mri_dir, ['c1', p.read.mri_filename, '.nii']);
% Set output file ( absolute path )
head_param.head_file        = p.save.head_dir_filename;
head_param.Nsurf            = 3; %3-shell head model in the EEG data
% Construct 3-shell surface model
vb_job_head_3shell(head_param); % purpose : Make 1 shell or 3 shell surface model for MEG/EEG

% Leadfield
% set input files (relative path from proj_root)
basis_param.brain_file      = head_param.brain_file;
basis_param.area_file       = '';
basis_param.head_file       = head_param.head_file;
basis_param.device          = 'BRAINAMP'; % EEG device name
basis_param.bem_mode        = 4; % 3-shell BEM
basis_param.Basis_mode      = 1; % single dipole perpendicular to the cortical surface
basis_param.meg_file      = p.save.eeg_dir_filename; % EEGinfoからsensor座標を取得するために必要
% set output file (absolute path)
basis_param.basis_file = p.save.basis_dir_filename;
% make leadfield matrix
vb_job_leadfield(basis_param)
% take common average
load(basis_param.basis_file,'basis')
basis = basis - repmat(mean(basis,2),1,size(basis,2)); % channel num回，行列（頂点数分のデータの）平均の作業を繰り返す
vb_save(basis_param.basis_file, 'basis')

%% step 9 : estimate variance & current (.bayes.mat, .basis.mat)
for now_reduce = p.variance_reduce
    reduce_str = replace(num2str(now_reduce),'.','_'); % 0.5とかを0_5に変えてくれる
    p.save.bayes_dir_filename = fullfile(p.save.bayes_dir,['reduce_', reduce_str],[p.save.dirname,'.bayes.mat']);
    p.save.current_dir_filename = fullfile(p.save.current_dir, ['reduce_', reduce_str], [p.save.dirname,'.curr.mat']);
    p.reduce = now_reduce;

    %% estimate current variance and source currents from EEG data
    % set input files and parameter
    time_noise = p.time_noise;
    bayes_parm = vb_set_bayes_default_parameters; %  Set default value for bayes prior and estimation parameters % std or subj
    bayes_parm.brainfile = p.save.brain_dir_filename;
    bayes_parm.areafile = p.save.area_dir_filename;
    bayes_parm.actfile = p.save.act_dir_filename;
    bayes_parm.megfile{1,1} = p.save.eeg_dir_filename; % tkit ver: {1,1}, API ver:
    bayes_parm.basisfile = p.save.basis_dir_filename;
    bayes_parm.megfile_baseline = bayes_parm.megfile;

    bayes_parm.spatial_smoother ='subj'; % 'std': Standard brain , 'subj': Individual brain

    bayes_parm.area_key = 'Cortex'; % current area
    bayes_parm.act_key{1,1} = 'Uniform'; % 'Prior for current variance' %'Uniform' or p.fmri_contrase_name or p.neurosynth_term
    bayes_parm.prior_weight = p.prior_weight; % Relative influence of prior information (0-1)

    % Load time information
    [~, ~, time_info] = vb_load_meg_data(bayes_parm.megfile{1});
    time = time_info.time;
    % time_info
    %   .time : [-4, -3.9984, .... 3.9984, 4]
    %   .sample_frequency 625

    % Set time window for stimulus-evoked activity
    bayes_parm.twin_meg = [1, length(time_info.time)]; % time window for analysis (all time)
    bayes_parm.Tperiod = p.Tperiod; % 100[ms] % time windows for calculating source current
    bayes_parm.Tnext = p.Next; % 50[ms] % time step for calculating source current

    % Set time window for baseline
    [~, from] = min(abs(time - time_noise(1)));
    [~, to] = min(abs(time - time_noise(2)));
    bayes_parm.twin_noise = [from, to];
    bayes_parm.twin_baseline = [from, to];

    bayes_parm.patch_norm = ON; % patch area normalization is applied
    bayes_parm.reduce = p.reduce;

    % Set output file
    bayes_parm.bayesfile = p.save.bayes_dir_filename;

    % Set parameter for noise convarience estimation
    bayes_parm = vb_set_noise_estimation_model(bayes_parm);

    % Show figure and parameter
    bayes_show_title = 'estimate current variance';
    bayes_show_param = {
        'prior weight', bayes_parm.prior_weight;
        'act key', bayes_parm.act_key;
        'analyzing time ( twin meg )', bayes_parm.twin_meg;
        'length of time window (Tperiod)', bayes_parm.Tperiod;
        'step of time window (Tnext)', bayes_parm.Tnext;
        'beta (twin_noise, twin_baseline)', bayes_parm.twin_noise;
        'patch_norm', bayes_parm.patch_norm;
    };

    figure();
    ax = axes('Position',[0,0,1,1]);
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';

    num_param = size(bayes_show_param,1);
    line_between = 1/(num_param +2);

    text(0.01,1-line_between,['\bf ' bayes_show_title]);

    for now_num = 1:num_param
        now_y = 1-line_between*(now_num+1);
        if isnumeric(bayes_show_param{now_num,2})
            text(0.1,now_y, append(bayes_show_param{now_num,1}, ':', num2str(bayes_show_param{now_num,2})), 'interpreter', 'none');
        else
            text(0.1,now_y, [bayes_show_param{now_num,1}, ':', bayes_show_param{now_num,2}], 'interpreter', 'none');
        end
    end

    vb_job_vb(bayes_parm) %  Estimate parameters of hierarchical Bayse model using variational Bayes
    %↑ ここでEEG時系列データ使われてるかも

    % Estimate source current
    current_parm.megfile = bayes_parm.megfile;
    current_parm.bayesfile = bayes_parm.bayesfile;
    current_parm.trial_average = OFF;

    data = vb_load_meg_data(p.save.eeg_dir_filename);
    current_parm.twin_meg = [1 size(data, 2)];
    current_parm.Tperiod = size(data, 2);

    current_parm.currfile = p.save.current_dir_filename;
    if p.curr_trial_save_flag
        current_parm.jactdir = 'EachTrial';
    end

    vb_job_current(current_parm);

end

%% step 10 : save each trial vertex current
area_filename = [p.read.mri_filename, '_HCP_MMP1.area.mat'];
% p.save.area_dir_filename      = fullfile(p.save.dir, 'brain', [p.read.mri_filename, '.area.mat']);
choose_key = {'R_3b_ROI', 'L_3b_ROI', 'R_4_ROI', 'L_4_ROI'};

    %% area file & act file
    disp('decide ROI vertex')
    ROI_vertex_num = [];
    area_file_path = fullfile(p.save.dir, 'brain', area_filename);
    % area_dir_filename = fullfile(p.save.dir, 'brain', area_filename);
    key_list = vb_get_keyset_area(p.save.area_dir_filename); % This function is used to get the set of cortical area IDs. 

    if ~exist(area_file_path, 'file')
        warning('HCP_MMP1 area file is not found!!!!!!!')
    end

    for now_num = 1:length(choose_key)
        now_key = choose_key{1,now_num};
        area = vb_get_area(area_file_path, now_key); % This function is used to get cortical area data. 
        ROI_vertex_num = [ROI_vertex_num; area.Iextract];
    end

    %% load current
    curr_type = 1;
    ave_mode = OFF;
    reduce = p.variance_reduce;
    reduce_str = replace(num2str(reduce), '.', '_');
    curr_dir_filename = fullfile(p.save.current_dir, ['reduce_', reduce_str],[p.save.dirname, '.curr.mat']);

    disp('load current now ...')
    % [Jinfo, Zacts,ix_act] = vb_load_current(curr_dir_filename, curr_type, average_mode, [], ROI_vertex_num); %load estimated current
    [Jinfo, Zacts,ix_act] = vb_load_current(curr_dir_filename, curr_type, ave_mode, [], ROI_vertex_num);

    if ~isequal(ix_act, ROI_vertex_num) || (length(ROI_vertex_num)~=size(Zacts,1))
        error('check');
    end

    disp('save current now ...')
    save_dir = fullfile(p.save.current_dir, ['reduce_', reduce_str], 'EachVertex');

    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end

    save_dir_filename = fullfile(save_dir, 'Z_current_%d.mat');

    for now_num = 1:length(ROI_vertex_num)
        vertex_num = ROI_vertex_num(now_num);
        Zact = Zacts(now_num,:,:);
        now_save_dir_filename = sprintf(save_dir_filename, vertex_num);
        save(now_save_dir_filename, "Zact", 'Jinfo', 'vertex_num');
    end
