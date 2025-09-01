# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 00:23:40 2024

@author: Hono
"""

##mneパッケージのインポート
import mne
import os

##############実験毎に入力するところ#########################################
# 作業ディレクトリを指定
data_dir = '/media/honoka/HDD2/Experiment'
day_sub = '20250313_B93'

#実験パラダイム
trial_num = 9 #試行回数

##ノイズの多いチャンネルの指定
reject_chan = ['TP9','T7'] #頻出チャンネル：'TP9','TP8','FT10'
########################################################################

#フォルダの生成
save_folder = f'{day_sub}_analysis_ref_EAR_CAR'
if not os.path.exists(save_folder):
    os.makedirs(save_folder)

  
# save reject channnel data (.txt)
with open(f'{save_folder}/{day_sub}_reject_ch_memomeo.txt', 'w') as f:
    f.write('reject channel : '+ ' '.join(map(str,reject_chan)))
    f.close()

for d in range(trial_num):
    ##入力データ
    vhdr_fname = f'{data_dir}/{day_sub}/{day_sub}_000{d+1}.vhdr'
    
    
    ##BrainVisionのEEG読み込む
    raw = mne.io.read_raw_brainvision(vhdr_fname,preload=True)
    
    
    ##チャンネルの除去,チャネルがなければ警告
    raw.drop_channels('EAR', on_missing='warn')
    raw.drop_channels(reject_chan, on_missing='warn')
    #raw.plot()
    
    # 特定のチャンネルタイプを設定する
    raw.set_channel_types({'EMG_right': 'emg',
                           'EMG_left': 'emg',
                           'EOG_V': 'eog',
                           'EOG_H': 'eog' })
    
    ##CAR空間フィルタ
    raw.set_eeg_reference(ref_channels = 'average', projection = False)
    
    #raw.plot()
    
    #montageの設定
    montage = mne.channels.make_standard_montage('standard_1020')
    raw.set_montage(montage)
    
    raw.save(f'{save_folder}/{day_sub}_000{d+1}_ch_reject_raw.fif', overwrite = True)
    
    #ICAのコンポーネントの表示
    filt_data = raw.copy().filter(l_freq=0.1, h_freq=None)
    before_ica_data = filt_data.copy().drop_channels(['EMG_right','EMG_left','EOG_V','EOG_H'])
    
    ica = mne.preprocessing.ICA(n_components=25,random_state=90,max_iter='auto')
    ica.fit(before_ica_data)
    ica.plot_components()
    ica.plot_sources(raw)
    
    ica.save(f'{save_folder}/{day_sub}_000{d+1}_ica.fif', overwrite=True)
    
    
    # いらない独立成分をメモ