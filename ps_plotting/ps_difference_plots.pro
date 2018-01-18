pro ps_difference_plots, folder_names, obs_info, ps_foldernames = ps_foldernames, $
    cube_types, pols, all_type_pol = all_type_pol, $
    uvf_options = uvf_options, ps_options = ps_options, $
    plot_options = plot_options, plot_2d_options = plot_2d_options, $
    refresh_diff = refresh_diff, freq_ch_range = freq_ch_range, $
    plot_slices = plot_slices, slice_type = slice_type, $
    plot_filebase = plot_filebase, save_path = save_path, savefilebase = savefilebase, $
    data_min_abs = data_min_abs, plot_1d = plot_1d, axis_type_1d=axis_type_1d, $
    invert_colorbar = invert_colorbar, $
    quiet = qiet, window_num = window_num

  compare_plot_prep, folder_names, obs_info,  cube_types, pols, 'diff', compare_files, $
    ps_foldernames = ps_foldernames, $
    uvf_options = uvf_options, ps_options = ps_options, $
    plot_options = plot_options, plot_2d_options = plot_2d_options, $
    plot_slices = plot_slices, slice_type = slice_type, $
    freq_ch_range = freq_ch_range, plot_filebase = plot_filebase, $
    save_path = save_path, savefilebase = savefilebase, $
    axis_type_1d = axis_type_1d, full_compare = all_type_pol

  for slice_i=0, compare_files.n_slices-1 do begin
    for cube_i=0, compare_files.n_cubes-1 do begin

      test_save = file_valid(compare_files.mid_savefile_2d[slice_i, cube_i])

      if tag_exist(compare_files, 'savefiles_1d') then begin
        test_save_1d = file_valid(reform(compare_files.savefiles_1d[cube_i,*]))
        test_save = min([test_save, test_save_1d])
      endif

      if test_save eq 0 or keyword_set(refresh_diff) then begin
        if keyword_set(plot_slices) then begin

          ps_slice_differences, compare_files.input_savefile1[slice_i, cube_i], compare_files.input_savefile2[slice_i, cube_i], $
            savefile_diff = compare_files.mid_savefile_2d[slice_i, cube_i]
        endif else begin

          ps_differences, compare_files.input_savefile1[slice_i, cube_i], $
            compare_files.input_savefile2[slice_i, cube_i], refresh = refresh_diff, $
            savefile_3d = compare_files.mid_savefile_3d[slice_i, cube_i], $
            savefile_2d = compare_files.mid_savefile_2d[slice_i, cube_i], $
            savefiles_1d = compare_files.savefiles_1d[cube_i,*], $
            wedge_amp = compare_files.wedge_amp, $
            wt_cutoffs = ps_options.wt_cutoffs, $
            wt_measures = ps_options.wt_measures, plot_options = plot_options

        endelse
      endif

      if plot_options.pub then font = 1 else font = -1

      if compare_files.n_cubes gt 1 then begin
        if cube_i eq 0 then begin
          if compare_files.n_cubes eq 6 then begin
            if plot_2d_options.kperp_linear_axis then begin
              ;; aspect ratio doesn't work out for kperp_linear with multiple rows
              ncol = 6
              nrow = 1
            endif else begin
              ncol = 3
              nrow = 2
            endelse
          endif else begin
            if plot_2d_options.kperp_linear_axis then begin
              ;; aspect ratio doesn't work out for kperp_linear with multiple rows
              ncol = compare_files.n_cubes
              nrow = 1
            endif else begin
              nrow = 2
              ncol = ceil(compare_files.n_cubes/nrow)
            endelse
          endelse
          start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}

          if n_elements(window_num) eq 0 then window_num = 1
        endif else begin
          pos_use = positions[*,cube_i]

        endelse
      endif

      if cube_i eq compare_files.n_cubes-1 and tag_exist(plot_options, 'note') then begin
        note_use = plot_options.note + ', ' + compare_files.kperp_density_names
      endif else begin
        undefine, note_use
      endelse

      if plot_options.pub then plotfile = compare_files.plotfiles_2d[slice_i]

      if keyword_set(plot_slices) then begin
        kpower_slice_plot, compare_files.mid_savefile_2d[slice_i, cube_i], $
          multi_pos = pos_use, start_multi_params = start_multi_params, $
          plot_options = plot_options, plot_2d_options = plot_2d_options, $
          note = note_use, data_min_abs = data_min_abs, $
          plotfile = plotfile, full_title=compare_files.titles[cube_i], $
          window_num = window_num, color_profile = 'sym_log', invert_colorbar = invert_colorbar, $
          wedge_amp = compare_files.wedge_amp

      endif else begin
        kpower_2d_plots, compare_files.mid_savefile_2d[slice_i, cube_i], $
          multi_pos = pos_use, start_multi_params = start_multi_params, $
          plot_options = plot_options, plot_2d_options = plot_2d_options, $
          note = note_use, data_min_abs = data_min_abs, $
          plotfile = plotfile, full_title=compare_files.titles[cube_i], $
          window_num = window_num, color_profile = 'sym_log', invert_colorbar = invert_colorbar, $
          wedge_amp = compare_files.wedge_amp, no_units=no_units
      endelse

      if compare_files.n_cubes gt 1 and cube_i eq 0 then begin
        positions = pos_use
        undefine, start_multi_params
      endif

    endfor
    undefine, positions, pos_use

    if plot_options.pub and compare_files.n_cubes gt 1 then begin
      cgps_close, png = plot_options.png, pdf = plot_options.pdf, $
        delete_ps = plot_options.delete_ps, density=600
    endif
    window_num += 1
  endfor

  if keyword_set(plot_1d) and not keyword_set(plot_slices) then begin
    for i=0, n_elements(wedge_amp) do begin

      if plot_options.pub then plotfiles_use = compare_files.plotfiles_1d[cube_i]

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

            window_num = 2+i
          endif else begin
            pos_use = positions[*,j]

          endelse
        endif

        if i eq compare_files.n_cubes-1 and tag_exist(plot_options, 'note') then begin
          note_use = plot_options.note
        endif else begin
          undefine, note_use
        endelse

        kpower_1d_plots, compare_files.savefiles_1d[j,i], window_num=window_num, $
          start_multi_params = start_multi_params, multi_pos = pos_use, $
          names=compare_files.titles[j], plot_options = plot_options, $
          plotfile = plotfiles_use, title = note_use, yaxis_type = axis_type_1d

        if j eq 0 then begin
          positions = pos_use
          undefine, start_multi_params
        endif

      endfor

    endfor
  endif

end
