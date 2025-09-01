
########## import library #########
import mne
import numpy as np
import os
import scipy.io

########## set parameter ###########
#dir
read_dir = '/media/honoka/HDD2/practice'
save_dir = '/media/honoka/HDD2/MATLAB/EEGMat/practice'

#read file
preprocessing_comment = 'analysis_ref_EAR_CAR'
channel_label_filename = 'label64ch.txt'

#parameter
day_sub = '20250313_B93'
pos_filename = '20250313_B93.pos.mat'
set_nums = [1,4,7]
which_finger = 'Rindex'

########## load data ##########
# eeg epoch data
epoch_list = []

for d in set_nums:
    epoch_fif_filename = f'{day_sub}_{d:04}_epo.fif'
    epoch_fif_dir_filename = os.path.join(read_dir, day_sub, f'{day_sub}_{preprocessing_comment}',epoch_fif_filename)
    now_epoch_list = mne.read_epochs(epoch_fif_dir_filename, preload = True)
    epoch_list.append(now_epoch_list)

epoch_list = mne.concatenate_epochs(epoch_list) # 各試行を1つのオブジェクトに連結してくれる
epoch_eeg_only = epoch_list.copy()
epoch_eeg_only = epoch_eeg_only.pick_types(eeg=True)

# .pos.mat file
pos_dir_filename = os.path.join(read_dir, day_sub, f'{day_sub}_pos', pos_filename)
pos_dot_mat = scipy.io.loadmat(pos_dir_filename)

# EEG channel label file
channel_label_dir_filename = os.path.join(read_dir,channel_label_filename)
with open(channel_label_dir_filename) as label:
    channel_label = label.read()
channel_label = np.array(channel_label.split())

########## adjust data ########
# reject bad channel from pos
epoch_channel = epoch_eeg_only.info['ch_names']
chan_idx = []
for now_epoch_channel in epoch_channel:
    now_chan_idx = np.where(channel_label == now_epoch_channel)[0][0]
    chan_idx.append(now_chan_idx)
channel_coord = pos_dot_mat['pos'][chan_idx,:]
if not np.allclose(chan_idx, np.sort(chan_idx)) : raise ValueError('please check souce code')

# print reject channel names
not_selected_chan = np.delete(channel_label,chan_idx)
print('====== channel names not selected ======')
print(not_selected_chan)
print('========================================')

########## for vbmeg ###########

# make class
class EEGinfo_struct:
    def __init__(self):
        self.SampleFrequency = 0
        self.Device = ''
        self.Nrepeat = 0
        self.Measurement = ''
        self.Pretrigger = 0
        self.Nsample = 0
        self.Coord = []
        self.Nchannel = 0

# get eeg data
eeg_get_data = epoch_eeg_only.get_data(units='V') # Nrepeat*Nchannel*Nsample
eeg_data = eeg_get_data.transpose(1,2,0) #Nchannel*Nsample*Nrepeat

# info for vbmeg
EEGinfo = EEGinfo_struct()
EEGinfo.SampleFrequency = float(epoch_eeg_only.info['sfreq'])
EEGinfo.Device = 'BASIC'
EEGinfo.Nchannel = float(eeg_data.shape[0])
EEGinfo.Nsample = float(eeg_data.shape[1])
EEGinfo.Nrepeat = float(eeg_data.shape[2])
EEGinfo.Measurement = 'EEG'
EEGinfo.Pretrigger = -1*float(epoch_eeg_only.times[0])*EEGinfo.SampleFrequency
EEGinfo.Coord = channel_coord

if eeg_data.shape != (EEGinfo.Nchannel, EEGinfo.Nsample, EEGinfo.Nrepeat) : raise ValueError('check')

# save
save_trial_name =''
for set_num  in set_nums:
    save_trial_name = save_trial_name +'_'+ '_'.join(str(set_num))

save_eeg_mat_dir_filename = os.path.join(save_dir,day_sub,f'{day_sub}_{which_finger}_{save_trial_name}.eeg.mat')
save_info_mat_dir_filename = os.path.join(save_dir,day_sub,f'{day_sub}_{which_finger}_{save_trial_name}_info.mat')

if os.path.exists(save_eeg_mat_dir_filename) == False:
    scipy.io.savemat(save_eeg_mat_dir_filename,mdict={'eeg_data':eeg_data, 'EEGinfo':EEGinfo,'Measurement':EEGinfo.Measurement})
    print('------------------------- save *.eeg.mat -------------------------')
else:
    print(f'This file have already exist : {save_eeg_mat_dir_filename}')

if os.path.exists(save_info_mat_dir_filename) == False:
    scipy.io.savemat(save_info_mat_dir_filename, mdict = {'ch_names':epoch_channel})
    print('------------------------- save *_info.mat -------------------------')
else:
    print(f'This file have already exist : {save_info_mat_dir_filename}')