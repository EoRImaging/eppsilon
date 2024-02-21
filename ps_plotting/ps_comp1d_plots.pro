pro ps_comp1d_plots, folder_names, obs_info, ps_foldernames = ps_foldernames, $
    cube_types, pols, uvf_options = uvf_options, freq_options = freq_options, $
    ps_options = ps_options, plot_options = plot_options, $
    binning_1d_options = binning_1d_options, names = names, $
    all_type_pol = all_type_pol, $
    plot_filebase = plot_filebase, save_path = save_path, savefilebase = savefilebase, $
    quiet = quiet, window_num = window_num

    plot_2d_options = create_plot_2d_options()

    compare_plot_prep, folder_names, obs_info,  cube_types, pols, 'comp_1d', compare_files, $
        ps_foldernames = ps_foldernames, $
        uvf_options = uvf_options, freq_options = freq_options, ps_options = ps_options, $
        plot_options = plot_options, plot_2d_options = plot_2d_options, $
        binning_1d_options = binning_1d_options, $
        plot_filebase = plot_filebase, $
        save_path = save_path, savefilebase = savefilebase, $
        full_compare = all_type_pol

    if plot_options.pub then font = 1 else font = -1

    for i=0, n_elements(compare_files.wedge_amp) do begin

        if plot_options.pub then plotfiles_use = compare_files.plotfiles_1d[i]
        window_num = 1+i

        for j=0, compare_files.n_cubes-1 do begin

            if compare_files.n_cubes gt 1 then begin
                if j eq 0 then begin
                    if compare_files.n_cubes eq 6 then begin
                    ncol = 3
                    nrow = 2
                    endif else begin
                    nrow = 2
                    ncol = ceil(compare_files.n_cubes/nrow)
                    endelse
                    start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}
                    undefine, positions, pos_use
                endif else begin
                    pos_use = positions[*,j]
                endelse
            endif

            if j eq compare_files.n_cubes-1 and tag_exist(plot_options, 'note') then begin
                note_use = plot_options.note
            endif else begin
                note_use = ''
            endelse

            input_files = [compare_files.input_savefile1[j,i], compare_files.input_savefile2[j,i]]
            if n_elements(names) gt 0 then begin
                if n_elements(names) ne 2 then message, 'names must have two elements'
                names_use = names
            endif else begin
                if n_elements(folder_names) eq 2 or n_elements(ps_foldernames) eq 2 then begin
                    names_use = strsplit(obs_info.diff_note, '-', /extract)
                endif else begin
                    names_use = strsplit(plot_options.note, 'and', /regex, /extract)
                endelse
            endelse
            if i lt n_elements(compare_files.wedge_amp) then begin
                wedge_name_use = binning_1d_options.wedge_names[i]
            endif else begin
                wedge_name_use = "no wedge cut"
            endelse

            kpower_1d_plots, input_files, window_num=window_num, $
                start_multi_params = start_multi_params, multi_pos = pos_use, $
                names=names_use, plot_options = plot_options, note = note_use, $
                plotfile = plotfiles_use, title = compare_files.titles[j] + ' ' + wedge_name_use

            if j eq 0 and compare_files.n_cubes gt 1 then begin
            positions = pos_use
            undefine, start_multi_params
            endif

        endfor

    endfor

    if keyword_set(plot_options.pub) and compare_files.n_cubes gt 1 then begin
        if plot_options.png then begin
            if plot_options.pdf then delete_ps_use = 0 else delete_ps_use = plot_options.delete_ps
            cgps_close, /png, delete_ps = delete_ps_use, density = 600
        endif
        if plot_options.pdf then begin
            if not plot_options.png then cgps_close
            cgps2pdf, plotfiles_use, delete_ps=plot_options.delete_ps
        endif
    endif


end
