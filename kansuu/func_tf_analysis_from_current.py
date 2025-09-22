import os
import numpy as np
import glob
import scipy.io
import mne

def fntk_make_ROI_vertices(param):
    # process of act.mat
    act_dir_filename = glob.glob(os.path.join(param.proj_dir, 'brain', '*'+param.brain_atlas+'.act.mat'))

    if not act_dir_filename:
        raise FileNotFoundError('cannot find act.mat')
    
    act = scipy.io.loadmat(act_dir_filename[0])['Act'][0,0]
    act_label_name = np.concatenate([act['label_name'][0,0][now_num,0] for now_num in range(len(act['label_name'][0,0]))])
    
    # process of area.mat
    area_dir_filename = np.char.replace(act_dir_filename,'.act.mat','.area.mat')[0]
    area = scipy.io.loadmat(area_dir_filename)['Area']
    ROI_vertices = np.empty((0,1), dtype = int)

    for now_ROI_area_name in param.ROI_area_name:
        area_idx = np.where(act_label_name == now_ROI_area_name)[0]
        if area[area_idx,0][0]['key'][0,0] != now_ROI_area_name:
            raise ValueError('area name is not correct')
        ROI_vertices = np.concatenate([ROI_vertices,area[area_idx,0][0]['Iextract'][0,0]], axis = 0)
    
    return ROI_vertices

def fntk_power2erds(power_instance, ref):

    print('\ncalculate ERD/S now ... \n')

    # Extract the power from the ref interval
    # (n_epochs, n_channels, n_freqs, n_times)
    power_data = power_instance.data
    n_epoch, n_channels, n_freq, n_times = power_data.shape
    power_times = power_instance.times
    ref_idx = [np.argsort(np.abs(power_times - now_num))[0] for now_num in ref ]
        # numpy.argsort(a, axis=-1, kind='quicksort', order=None) : Returns the indices that would sort an array.
    ref_power = power_data[:, :, :, ref_idx[0]:ref_idx[1]]

    # Combination of all trials & Time-axis direction average
    # will be (n_channels, n_freqs)
    ref_power = np.concatenate([ref_power[now_num, :, :, :] for now_num in range(ref_power.shape[0])], axis = 2 )
    ref_power = np.mean(ref_power, axis = 2)

    # Extend (n_channels, n_freqs) to (n_epochs, n_channels, n_freqs, n_times)
    ref_power = np.stack([ref_power for _ in range(n_epoch)], axis=0)
    ref_power = np.stack([ref_power for _ in range(n_times)], axis=3)
        # numpy.stack(arrays, axis=0, out=None, *, dtype=None, casting='same_kind') : Join a sequence of arrays along a new axis.

    # calculate ERD/S
    erds_data = power_data / ref_power -1

    info = power_instance.info
    erds_instance = mne.time_frequency.EpochsTFR(info = info, data = erds_data, times = power_times, freqs = power_instance.freqs)

    return erds_instance

def fntk_current_tf_analysis(param):
    param.proj_dir = os.path.join(param.vbmeg_analysis_dir,param.proj_name)

    # act file & area file
    ROI_vertices = fntk_make_ROI_vertices(param)
    curr_dir_filename = os.path.join(param.proj_dir, 'current','reduce_'+str(param.curr_reduce_ratio).replace('.','_'), 'EachVertex', 'Z_current_{}.mat')

    for now_loop_num in range(len(ROI_vertices)):
        print('\n[ ' , str(now_loop_num), ' / ' , str(len(ROI_vertices)), ' ]\n')
        now_vertex_num = ROI_vertices[now_loop_num][0]

        # load source current ( Njact, Nsample, Ntrial )
        curr_dot_mat = scipy.io.loadmat(curr_dir_filename.format(str(now_vertex_num)))

        # change to (n_epoch, n_channels, n_times)
        zact_data = curr_dot_mat['Zact'].transpose(2,0,1)

        # make epoch
        sfreq = curr_dot_mat['Jinfo'][0,0]['SampleFreq'][0,0].astype(np.float64)
        tmin = np.round((-1)*curr_dot_mat['Jinfo'][0,0]['Pretrigger'][0,0].astype(np.float64)/sfreq, decimals = 2)
        info = mne.create_info(ch_names = ['Zact_'+ str(now_vertex_num)],
                               ch_types = ['csd'], # 'csd' : Current source density (scaled by 1000 to plot in mV/m²)
                               sfreq = sfreq)
        epoch = mne.EpochsArray(zact_data, info, tmin = tmin)

        tfa_freqs = np.arange(param.f_border[0], param.f_border[1]+1., param.f_step)

        # time frequenvy
        if param.tf_method == 'multitaper':
            power = mne.time_frequency.tfr_multitaper(epoch, tfa_freqs, param.n_cycles, average = False, return_itc = False, use_fft = True)
            # mne.time_frequency.tfr_multitaper(inst, freqs, n_cycles, time_bandwidth=4.0, use_fft=True, return_itc=True, decim=1, n_jobs=None, picks=None, average=True, *, verbose=None)
        else:
            raise ValueError('this tf_method is not implemented yet >< : ' + param.tf_method)

        # calculate erds
        if param.ref_method == 'trial':
            erds_instance = power.apply_baseline(baseline = tuple(param.ref), mode = 'percent')
            # apply_baseline(self, baseline, mode="mean", verbose=None):
            # mode : 'mean' | 'ratio' | 'logratio' | 'percent' | 'zscore' | 'zlogratio'
            # - 'percent' : subtracting the mean of baseline values followed by dividing by the mean of baseline values

        elif param.ref_method == 'all':
            erds_instance = fntk_power2erds(power_instance = power, ref = param.ref)
            # calculate ERD/S using power imstance

        else:
            raise ValueError('this ref_method is not implemented yet >< : ' + param.ref_method)


        # plot figure with bootstrap

        # save erds data
        save_dir = os.path.join(param.proj_dir, 'tf_map_' + param.tf_map_dir_comment, str(now_vertex_num))
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)
            print('make save_dir now : ', save_dir)
        
        save_dir_filename = os.path.join(save_dir, 'erds_instance_'+ str(now_vertex_num)+'-tfr.h5')

        erds_data = erds_instance.data
        erds_instance.save(save_dir_filename, overwrite = False)

        save_tf_data_dir_filename = os.path.join(save_dir, 'tf_analysis_data_' + str(now_vertex_num) + '.mat')
        save_tf_info_dir_filename = os.path.join(save_dir, 'tfa_info_'+ str(now_vertex_num) + '.mat')
        
        scipy.io.savemat(save_tf_data_dir_filename, mdict = {'erds_data' : erds_data, 'sfreq' : erds_instance.info['sfreq']})
        scipy.io.savemat(save_tf_info_dir_filename, mdict = {'param' : param})
