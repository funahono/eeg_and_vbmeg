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
data_dir = '/media/honoka/HDD2/Experiment'
day_sub = '20250313_B93'

#残すコンポーネントの指定
nokosu_list = [[3,4],#1
               [13,21],#2
               [8],#3
               [10,12],#4
               [5,8],#5
               [4,10,11],#6
               [3,6,7,8],#7
               [4,12],#8
               [4,6,7,9,11],#9
               []]#10

#paradim info
preparation_time = 3 #準備時間
task_time = 3 #課題時間
rest_time = 4 #休息時間
event_id = 8 #刺激提示タイミングのトリガ
trial_list = [1, 2, 3, 4, 5, 6, 7, 8, 9]  # 処理する試行番号を直接リストに
########################################################################

# make exclude list
trial_num =len(trial_list)
ch_list = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24] for d in range(trial_num+1)]

exclude_list =[[ch for ch in ch_list[i] 
                    if ch not in nokosu_list[i]] 
                    if nokosu_list[i]  
                    else ch_list[i]
                    for i in range(trial_num+1) ]

#print(len(ch_list))
#print(len(exclude_list))

#フォルダの生成
save_folder =f'/media/honoka/HDD2/Experiment/{day_sub}/{day_sub}_analysis_ref_EAR_CAR'
if not os.path.exists(save_folder):
    os.makedirs(save_folder)
    
for d in trial_list:
    
    
    # save exclude component data (.txt)
    with open(f'{save_folder}/{day_sub}_exclude_component_memomeo.txt', 'w') as f:
        f.write('exclude component :\n')
        components = exclude_list[d-1]
        f.write(f'{d:04}:[components]\n')
        f.close()
        
    # read fif data
    raw = mne.io.read_raw_fif(f'{save_folder}/{day_sub}_00{d:02}_ch_reject_raw.fif', preload=True)
    
    # exclude & apply ica
    ica = read_ica(f'{save_folder}/{day_sub}_00{d:02}_ica.fif')
    ica.exclude = exclude_list[d-1]
    data = ica.apply(raw)
    #data.plot()
    
    # save ica　(.fif)
    data.save(f'{save_folder}/{day_sub}_00{d:02}_reconst_raw.fif', overwrite = True)
    
    #trigの取得，インベントリの列指定
    trig = mne.events_from_annotations(data)
    trig = trig[0]
    
    #epoching & plot
    tmin = (-1)*preparation_time - 2 #フィルタの影響を考慮して-2
    tmax = task_time + rest_time + 2 #フィルタの影響を考慮して+2
    epochs = mne.Epochs(data, trig, event_id, tmin, tmax) 
    #epochs.plot()
    
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
    scipy.io.savemat(f'{save_folder}/{day_sub}_00{d:02}_epo.mat', epoch_dict, format='5', long_field_names=True)
    epochs.save(f'{save_folder}/{day_sub}_00{d:02}_epo.fif', overwrite=True)