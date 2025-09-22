import os
import scipy.io
import numpy as np
import sys
sys.path.append('/home/honoka/programdir/kansuu')
from func_tf_analysis_from_current import fntk_current_tf_analysis

import mne
import matplotlib.pyplot as plt

class parameter:
    pass

param = parameter()

param.vbmeg_analysis_dir = '/media/honoka/HDD2/MATLAB/vbmeg_analysis'
proj_name = ['20250313_B93_Rindex_20250909_ear_ref_car_standard_brain']

param.brain_atlas         = 'HCP_MMP1'
param.ROI_area_name       = ['R_4_ROI','L_4_ROI','R_3b_ROI','L_3b_ROI']
param.curr_reduce_ratio   = 1

param.tf_method           = 'multitaper'
param.f_border            = [6., 40.]
param.f_step              = 1
param.ref                 = [-3., -1.]
param.ref_method          = 'all' # 'trial' or 'all'

tfa_n_cycle               = [15]

for now_proj_name in proj_name:
    param.proj_name = now_proj_name

    for now_n_cycle in tfa_n_cycle:
        param.n_cycles  = now_n_cycle
        param.tf_map_dir_comment = 'n_cycle'+str(param.n_cycles)
        fntk_current_tf_analysis(param)
