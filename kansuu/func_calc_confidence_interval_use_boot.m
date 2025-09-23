function [sig_cl, sig_cu] = func_calc_confidence_interval_use_boot(erds_data_dir_filename, alpha)

    % disp(erds_data_dir_filename)
    m = matfile(erds_data_dir_filename);
    erds_data = m.erds_data;

    [n_epochs, n_channels, n_freqs, n_times] = size(erds_data);

    sig_cl = NaN(n_channels,n_freqs,n_times,length(alpha));
    sig_cu = NaN(n_channels,n_freqs,n_times,length(alpha));

    for now_num_chan = 1:n_channels
        for now_num_freqs = 1:n_freqs
            disp(append(' [', num2str(now_num_freqs), ' / ',  num2str(n_freqs), ' ]'))
            [~, now_cl, now_cu] = bootts(squeeze(erds_data(:,now_num_chan, now_num_freqs, :)), 300, alpha);
            sig_cl(now_num_chan,now_num_freqs,:,:) = now_cl';
            sig_cu(now_num_chan,now_num_freqs,:,:) = now_cu';
        end
    end
end