# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 00:23:40 2024

@author: Hono
"""

##mneパッケージのインポート
import mne
import os
import matplotlib.pyplot as plt

##############実験毎に入力するところ#########################################
# 作業ディレクトリを指定
data_dir = '/media/honoka/HDD2/Experiment'
day_sub = '20250313_B93'

# 実験パラダイム
trial_list = [1, 2, 3, 4, 5, 6, 7, 8, 9]  # 処理する試行番号を直接リストに

# ノイズの多いチャンネルの指定
reject_chan = ['T7','TP9']  # 頻出チャンネル：'TP9','TP8','FT10'
########################################################################

# フォルダの生成
save_folder = f'/media/honoka/HDD2/Experiment/{day_sub}/{day_sub}_analysis_ref_EAR_CAR'
if not os.path.exists(save_folder):
    os.makedirs(save_folder)

# save reject channel data (.txt)
with open(f'{save_folder}/{day_sub}_reject_ch_memomeo.txt', 'w') as f:
    f.write('reject channel : ' + ' '.join(map(str, reject_chan)))
    f.close()

for d in trial_list:
    
    print(f'-------------------------- {d} / {len(trial_list)} ----------------------------------')
    # 試行番号のログ
    trial_num = f"Trial {d:02}"
    
    # 入力データ
    vhdr_fname = f'{data_dir}/{day_sub}/{day_sub}_00{d:02}.vhdr'

    ## BrainVisionのEEG読み込む
    raw = mne.io.read_raw_brainvision(vhdr_fname, preload=True)

    ## チャンネルの除去, チャネルがなければ警告
    raw.drop_channels('EAR', on_missing='warn')
    raw.drop_channels(reject_chan, on_missing='warn')
    # raw.plot()

    # 特定のチャンネルタイプを設定する
    raw.set_channel_types({'EMG_right': 'emg',
                           'EMG_left': 'emg',
                           'EOG_V': 'eog',
                           'EOG_H': 'eog'})

    ## CAR空間フィルタ
    raw.set_eeg_reference(ref_channels='average', projection=False)

    # raw.plot()

    # montageの設定
    montage = mne.channels.make_standard_montage('standard_1020')
    raw.set_montage(montage)

    raw.save(f'{save_folder}/{day_sub}_00{d:02}_ch_reject_raw.fif', overwrite=True)

    # ICAのコンポーネントの表示
    filt_data = raw.copy().filter(l_freq=0.1, h_freq=None)
    before_ica_data = filt_data.copy().drop_channels(['EMG_right', 'EMG_left', 'EOG_V', 'EOG_H'])

    ica = mne.preprocessing.ICA(n_components=25, random_state=90, max_iter='auto')
    ica.fit(before_ica_data)
    fig = ica.plot_sources(raw, title=f'ICA Sources - Trial {d:02} - wave')
    # ica.plot_components(title=f'ICA Sources - Trial {d:02}- topograph')

    plt.show(block=True)
    
    input("Press Enter to continue to the next components...")
    plt.close(fig)

    ica.save(f'{save_folder}/{day_sub}_00{d:02}_ica.fif', overwrite=True)
