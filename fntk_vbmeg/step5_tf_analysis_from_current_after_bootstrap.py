import os
import scipy.io
import numpy as np
import mne
import sys
import matplotlib.pyplot as plt

sys.path.append('/home/honoka/programdir/kansuu')
from func_tf_analysis_from_current import fntk_save_tf_map

class parameter:
    pass

param = parameter()

param.vbmeg_analysis_dir = '/media/honoka/HDD2/MATLAB/vbmeg_analysis'
proj_names = ['20250313_B93_Rindex_20250909_ear_ref_car_standard_brain']

param.brain_atlas         = 'HCP_MMP1'
param.ROI_area_name       = ['R_4_ROI','L_4_ROI','R_3b_ROI','L_3b_ROI']
# param.curr_reduce_ratio   = 1

# param.tf_method           = 'multitaper'
# param.f_border            = [6., 40.]
# param.f_step              = 1
# param.ref                 = [-3., -1.]
# param.ref_method          = 'all' # 'trial' or 'all'

# tfa_n_cycle               = [15]

param.tf_map_dir_comment = 'n_cycle_15'

param.cmap = 'seismic'
param.vmin = -1
param.vmax = +1
param.tmin = -3
param.tmax = +7
param.savefig_format = 'png'

param.t_task = [0, 3]

for now_proj_name in proj_names:
    param.proj_name = now_proj_name
    fntk_save_tf_map(param)