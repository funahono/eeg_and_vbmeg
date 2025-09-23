function func_tf_analysis_from_currents(param)
    disp(param);
    proj_dir = fullfile(param.vbmeg_analysis_dir, param.proj_name);

    ROI_vertex =[];
    area_dir_filename = fullfile(proj_dir, 'brain', append(param.mri_filename, '_', param.brain_atlas, '.area.mat'));

    for now_num = 1:length(param.ROI_area_key)
        now_key = param.ROI_area_key{now_num};
        area = vb_get_area(area_dir_filename, now_key);
        ROI_vertex = [ROI_vertex; area.Iextract];
    end

    clear now_key area now_num area_dir_filename

    chosen_vertex = [];
    memo_dir_filename = fullfile(proj_dir, ['tf_map', '_', param.tf_map_dir_comment], '%d', 'chosen.txt');
    save_dir_filename = fullfile(proj_dir, ['tf_map', '_', param.tf_map_dir_comment], '%d', 'tf_analysis_sig_%d.mat');

    now_loop_num = 1;
    loop_flag = true;

    while loop_flag
        now_memo_dir_filename = sprintf(memo_dir_filename, ROI_vertex(now_loop_num, 1));
        now_save_dir_filename = sprintf(save_dir_filename, ROI_vertex(now_loop_num, 1), ROI_vertex(now_loop_num,1));

        if (exist(now_memo_dir_filename, 'file') ~= 2) && (exist(now_save_dir_filename, 'file') ~= 2)
            chosen_vertex = [chosen_vertex; ROI_vertex(now_loop_num,1)];
            fileID = fopen(now_memo_dir_filename, 'w');
            fclose(fileID);
        else
            disp(ROI_vertex(now_loop_num,1))
        end

        if (now_loop_num == length(ROI_vertex)) || length(chosen_vertex) == param.how_many_vertex
            loop_flag = false;
        end

        now_loop_num = now_loop_num + 1;

    end

    clear now_loop_num loop_flag now_memo_dir_filename now_save_dir_filename

    erds_data_dir_filename = fullfile(proj_dir, ['tf_map_',param.tf_map_dir_comment], '%d', 'tf_analysis_data_%d.mat');
    disp(erds_data_dir_filename)
    now_loop_num = 0;
    
    for now_chosen_vertex = chosen_vertex'
        now_loop_num = now_loop_num + 1;
        disp([' -----------------------  [ ', num2str(now_loop_num), ' / ', num2str(length(chosen_vertex)), ' ]  -----------------------'])

        now_erds_dir_filename = sprintf(erds_data_dir_filename, now_chosen_vertex, now_chosen_vertex);

        disp(now_erds_dir_filename)
        disp('--------------------------------------------------------')
        [sig_cl, sig_cu] = func_calc_confidence_interval_use_boot(now_erds_dir_filename, param.alpha);
        
        alpha = param.alpha;

        if exist(sprintf(save_dir_filename, now_chosen_vertex, now_chosen_vertex), 'file') == 2
            error('this sig file already exist ><')
        end

        save(sprintf(save_dir_filename, now_chosen_vertex, now_chosen_vertex), 'sig_cu', 'sig_cl', 'alpha');
        delete(sprintf(memo_dir_filename, now_chosen_vertex))
        delete(now_erds_dir_filename)

    end
end