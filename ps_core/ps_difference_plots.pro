pro ps_difference_plots, folder_names, obs_info, cube_types, pols, all_type_pol = all_type_pol, $
    refresh_diff = refresh_diff, freq_ch_range = freq_ch_range, $
    plot_slices = plot_slices, slice_type = slice_type, $
    plot_path = plot_path, plot_filebase = plot_filebase, save_path = save_path, savefilebase = savefilebase, $
    note = note, spec_window_types = spec_window_types, ave_removal = ave_removal, $
    data_range = data_range, data_min_abs = data_min_abs, $
    baseline_axis = baseline_axis, delay_axis = delay_axis, $
    kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, plot_1d = plot_1d, axis_type_1d=axis_type_1d, $
    wt_cutoffs = wt_cutoffs, wt_measures = wt_measures, invert_colorbar = invert_colorbar, $
    plot_wedge_line = plot_wedge_line, hinv = hinv, quiet = qiet, png = png, eps = eps, pdf = pdf, window_num = window_num
        
  if keyword_set(png) or keyword_set(eps) or keyword_set(pdf) then pub = 1 else pub = 0
  if pub then begin
    if keyword_set(png) and keyword_set(eps) and keyword_set(pdf) then begin
      print, 'only one of eps, pdf and png can be set, using png'
      eps = 0
    endif
    
    if keyword_set(png) then begin
      plot_exten = '.png'
      delete_ps = 1
    endif else if keyword_set(pdf) then begin
      plot_exten = '.pdf'
      delete_ps = 1
    endif else if keyword_set(eps) then begin
      plot_exten = '.eps'
      delete_ps = 0
    endif
  endif
  
  compare_plot_prep, folder_names, obs_info, cube_types, pols, 'diff', compare_files, $
    plot_slices = plot_slices, slice_type = slice_type, $
    spec_window_types = spec_window_types, freq_ch_range = freq_ch_range, ave_removal = ave_removal, $
    baseline_axis = baseline_axis, delay_axis = delay_axis, plot_wedge_line = plot_wedge_line, hinv = hinv, $
    plot_path = plot_path, plot_filebase = plot_filebase, save_path = save_path, savefilebase = savefilebase, $
    axis_type_1d = axis_type_1d, note = note, wt_cutoffs = wt_cutoffs, wt_measures = wt_measures, $
    pub = pub, plot_exten = plot_exten, full_compare = all_type_pol
    
    
  for slice_i=0, compare_files.n_slices-1 do begin
    for cube_i=0, compare_files.n_cubes-1 do begin
    
      test_save = file_test(compare_files.mid_savefile_2d[slice_i, cube_i]) *  (1 - file_test(compare_files.mid_savefile_2d[slice_i, cube_i], /zero_length))
      
      if tag_exist(compare_files, 'savefiles_1d') then begin
        test_save_1d = file_test(reform(compare_files.savefiles_1d[cube_i,*])) * (1 - file_test(reform(compare_files.savefiles_1d[cube_i,*]) , /zero_length))
        test_save = min([test_save, test_save_1d])
      endif
      
      if test_save eq 0 or keyword_set(refresh_diff) then begin
        if keyword_set(plot_slices) then begin
        
          ps_slice_differences, compare_files.input_savefile1[slice_i, cube_i], compare_files.input_savefile2[slice_i, cube_i], $
            savefile_diff = compare_files.mid_savefile_2d[slice_i, cube_i]
        endif else begin
        
          ps_differences, compare_files.input_savefile1[slice_i, cube_i], compare_files.input_savefile2[slice_i, cube_i], refresh = refresh_diff, $
            savefile_3d = compare_files.mid_savefile_3d[slice_i, cube_i], savefile_2d = compare_files.mid_savefile_3d[slice_i, cube_i], $
            savefiles_1d = compare_files.savefiles_1d[cube_i,*], $
            wedge_amp = compare_files.wedge_amp, wt_cutoffs = wt_cutoffs, wt_measures = wt_measures
            
        endelse
      endif
      
      if keyword_set(pub) then font = 1 else font = -1
      
      if compare_files.n_cubes gt 1 then begin
        if cube_i eq 0 then begin
          if compare_files.n_cubes eq 6 then begin
            if keyword_set(kperp_linear_axis) then begin
              ;; aspect ratio doesn't work out for kperp_linear with multiple rows
              ncol = 6
              nrow = 1
            endif else begin
              ncol = 3
              nrow = 2
            endelse
          endif else begin
            if keyword_set(kperp_linear_axis) then begin
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
      
      if cube_i eq compare_files.n_cubes-1 and n_elements(note) gt 0 then note_use = note + ', ' + compare_files.kperp_density_names else undefine, note_use
      
      if keyword_set(pub) then plotfile = compare_files.plotfiles_2d[slice_i]
      
      if keyword_set(plot_slices) then begin
        case compare_files.slice_tags[slice_i] of
          'xz': begin
            plot_xrange = compare_files.kperp_plot_range
            if n_elements(kpar_plot_range) gt 0 then plot_yrange = kpar_plot_range
          end
          'yz': begin
            plot_xrange = compare_files.kperp_plot_range
            if n_elements(kpar_plot_range) gt 0 then plot_yrange = kpar_plot_range
          end
          'xy': begin
            plot_xrange = compare_files.kperp_plot_range
            plot_yrange = compare_files.kperp_plot_range
          end
        endcase
        if keyword_set(kperp_linear_axis) or keyword_set(kpar_linear_axis) then linear_axes = 1        
        
        kpower_slice_plot, compare_files.mid_savefile_2d[slice_i, cube_i], multi_pos = pos_use, start_multi_params = start_multi_params, $
          plot_xrange = plot_xrange, plot_yrange = plot_yrange, note = note_use, $
          data_range = data_range, data_min_abs = data_min_abs, png = png, eps = eps, pdf = pdf, $
          plotfile = plotfile, full_title=compare_files.titles[cube_i], $
          window_num = window_num, color_profile = 'sym_log', invert_colorbar = invert_colorbar, $
          linear_axes = linear_axes, baseline_axis = baseline_axis, delay_axis = delay_axis, $
          wedge_amp = compare_files.wedge_amp, plot_wedge_line = plot_wedge_line, hinv = hinv

          if compare_files.slice_tags[slice_i] eq 'xz' or compare_files.slice_tags[slice_i] eq 'yz' $
            and n_elements(kpar_plot_range) eq 0 then kpar_plot_range = plot_yrange
      endif else begin
        kpower_2d_plots, compare_files.mid_savefile_2d[slice_i, cube_i], multi_pos = pos_use, start_multi_params = start_multi_params, $
          kperp_plot_range = compare_files.kperp_plot_range, kpar_plot_range = kpar_plot_range, note = note_use, $
          data_range = data_range, data_min_abs = data_min_abs, png = png, eps = eps, pdf = pdf, $
          plotfile = plotfile, full_title=compare_files.titles[cube_i], $
          window_num = window_num, color_profile = 'sym_log', invert_colorbar = invert_colorbar, $
          kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, baseline_axis = baseline_axis, delay_axis = delay_axis, $
          wedge_amp = compare_files.wedge_amp, plot_wedge_line = plot_wedge_line, hinv = hinv, no_units=no_units
      endelse
      
      if compare_files.n_cubes gt 1 and cube_i eq 0 then begin
        positions = pos_use
        undefine, start_multi_params
      endif
      
    endfor
    undefine, positions, pos_use
    
    if keyword_set(pub) and compare_files.n_cubes gt 1 then begin
      cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density=600
    endif
    window_num += 1
  endfor
  
  if keyword_set(plot_1d) and not keyword_set(plot_slices) then begin
    for i=0, n_elements(wedge_amp) do begin
    
      if keyword_set(pub) then plotfiles_use = plotfiles_1d[cube_i]
      
      for j=0, n_cubes-1 do begin
      
        if n_cubes gt 1 then begin
          if j eq 0 then begin
            if n_cubes eq 6 then begin
              ncol = 3
              nrow = 2
            endif else begin
              nrow = 2
              ncol = ceil(n_cubes/nrow)
            endelse
            start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}
            undefine, positions, pos_use
            
            window_num = 2+i
          endif else begin
            pos_use = positions[*,j]
            
          endelse
        endif
        
        if i eq n_cubes-1 and n_elements(note) gt 0 then note_use = note else undefine, note_use
        
        kpower_1d_plots, savefiles_1d[j,i], window_num=window_num, start_multi_params = start_multi_params, multi_pos = pos_use, $
          names=titles[j], hinv = hinv, png = png, eps = eps, pdf = pdf, plotfile = plotfiles_use, note = note_use, yaxis_type = axis_type_1d
          
        if j eq 0 then begin
          positions = pos_use
          undefine, start_multi_params
        endif
        
      endfor
      
    endfor
  endif
  
end