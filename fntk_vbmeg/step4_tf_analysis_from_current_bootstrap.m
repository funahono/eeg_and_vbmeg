clear
close all

% path
addpath(genpath('/home/honoka/vbmeg3_0_0_a_2'));
addpath(genpath('/home/honoka/programdir/kansuu'));


param.vbmeg_analysis_dir = '/media/honoka/HDD2/MATLAB/vbmeg_analysis';
proj_name = {'20250313_B93_Rindex_20250909_ear_ref_car_standard_brain'};
% /media/honoka/HDD2/MATLAB/vbmeg_analysis/20250313_B93_Rindex_20250909_ear_ref_car_standard_brain/tf_map_n_cycle15

param.mri_filename = 'mni_icbm152_t1_tal_nlin_asym_09c'; % omit 'nii'
param.tf_map_dir_comment = 'n_cycle_15';

% decide ROI vertex
param.brain_atlas = 'HCP_MMP1';
param.ROI_area_key = {'L_4_ROI','R_4_ROI','L_3b_ROI','R_3b_ROI'};
% Primary Motor Cortex : 'L_4_ROI','R_4_ROI'
% Primary Sensory Cortex : 'L_3b_ROI','R_3b_ROI'

param.how_many_vertex = 400;

param.alpha = [0.01, 0.05]; %信頼しない区間

for now_proj_num = 1:length(proj_name)
    now_proj_name = proj_name{1,now_proj_num};
    tf_map_dir = fullfile(param.vbmeg_analysis_dir, now_proj_name, ['tf_map_', param.tf_map_dir_comment]);

    if exist(tf_map_dir, 'dir') ~= 7
        disp(tf_map_dir)
        error('this tf map_dir cannot find ><')
    end
end

clear now_proj_num now_proj_name tf_map_dir

for now_proj_num = 1:length(proj_name)
    param.proj_name = proj_name{1,now_proj_num};
    func_tf_analysis_from_currents(param)
end