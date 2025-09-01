# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 22:42:33 2024

@author: Hono
"""
# mneパッケージのインポート
import mne

##############実験毎に入力するところ#########################################
# 作業ディレクトリの指定
data_dir = 'C:/Users/Hono/Desktop/datadir'
day_sub = '20250221_B94'

# 実験パラダイム
trial_list = [10,11,12]  # 処理する試行番号を直接リストに
########################################################################

for d in trial_list:
    # 入力データ
    vhdr_fname = f'{data_dir}/{day_sub}/{day_sub}_00{d:02}.vhdr' 
    
    # BrainVisionのEEGを読み込む
    raw = mne.io.read_raw_brainvision(vhdr_fname)
    
    # データをプロットする
    raw.plot()

# 怖い電極をメモ
