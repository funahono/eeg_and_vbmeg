# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 00:01:18 2024

@author: Hono
"""

import mne
from mne.preprocessing import read_ica
import scipy
import os

##############実験毎に入力するところ#########################################

# 作業ディレクトリを指定
data_dir = 'C:/Users/Hono/Desktop/datadir'
day_sub = '20241223_B90'

#除外コンポーネントの指定
exclude_list = [[0,1,2,3,4,5,6,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24],
                [0,1,2,3,4,5,6,8,9,10,11,13,14,15,16,17,18,19,20,21,22,23,24],
                [0,1,2,3,4,5,6,8,9,10,12,13,14,15,16,17,18,19,20,21,22,23,24],
                [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24],
                [0,1,4,5,6,7,8,9,10,11,12,13,15,16,17,18,19,20,21,22,23,24],
                [0,2,3,4,5,6,7,8,11,12,13,14,15,16,17,18,19,20,21,22,23,24],
                [0,1,2,3,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24],
                [0,1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19,20,21,22,23,24],
                [0,1,2,3,4,5,7,9,10,12,13,14,15,16,17,18,19,20,21,22,23,24]]


# exclude_list = [[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24],
#                 [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24],
#                 [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24],
#                 [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24],
#                 [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24],
#                 [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24],
#                 [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24],
#                 [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24],
#                 [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]]#0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24

#paradim info
preparation_time = 3 #準備時間
task_time = 3 #課題時間
rest_time = 4 #休息時間
event_id = 8 #刺激提示タイミングのトリガ
trial_num = 9 #試行回数
########################################################################

#フォルダの生成
save_folder =f'{day_sub}_analysis_ref_EAR_CAR'
if not os.path.exists(save_folder):
    os.makedirs(save_folder)
    
for d in range(trial_num):
    
    
    # save exclude component data (.txt)
    with open(f'{save_folder}/{day_sub}_exclude_component_memomeo.txt', 'w') as f:
        f.write('exclude component :\n')
        components = exclude_list[d]
        f.write(f'{d+1:04}:[components]\n')
        f.close()
        
    # read fif data
    raw = mne.io.read_raw_fif(f'{save_folder}/{day_sub}_000{d+1}_ch_reject_raw.fif', preload=True)
    
    # exclude & apply ica
    ica = read_ica(f'{save_folder}/{day_sub}_000{d+1}_ica.fif')
    ica.exclude = exclude_list[d]
    data = ica.apply(raw)
    data.plot()
    
    # save ica　(.fif)
    data.save(f'{save_folder}/{day_sub}_000{d+1}_reconst_raw.fif', overwrite = True)
    
    #trigの取得，インベントリの列指定
    trig = mne.events_from_annotations(data)
    trig = trig[0]
    
    #epoching & plot
    tmin = (-1)*preparation_time - 2 #フィルタの影響を考慮して-2
    tmax = task_time + rest_time + 2 #フィルタの影響を考慮して+2
    epochs = mne.Epochs(data, trig, event_id, tmin, tmax) 
    epochs.plot()
    
    # epoch data dict設定
    epoch_data = epochs.get_data()    
    epoch_dict = {'epoch_data':epoch_data, 
                  'sfreq':epochs.info['sfreq'], 
                  'times':epochs.times, 
                  'event_trigger':epochs.events.T,
                  'ch_list':epochs._channel_type_idx,
                  'ch_names':epochs.ch_names, 
                  'drop_epoch_log':epochs.drop_log}
    
    #save epoch data (.mat),(.fif)
    scipy.io.savemat(f'{save_folder}/{day_sub}_000{d+1}_epo.mat', epoch_dict, format='5', long_field_names=True)
    epochs.save(f'{save_folder}/{day_sub}_000{d+1}_epo.fif', overwrite=True)