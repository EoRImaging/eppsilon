pro hellebore_diff_wrapper, folder_names, obs_names_in, cube_types = cube_types, pols = pols, all_type_pol = all_type_pol, $
    png = png, eps = eps, pdf = pdf, $
    kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis
    
    
  if n_elements(folder_names) eq 0 then folder_names = base_path('data') + 'fhd_ps_data/128T_cubes/aug23_3hr_first/'
  
  if (n_elements(folder_names) eq 2 or n_elements(obs_names_in) eq 2) $
    and n_elements(cube_types) eq 0 and n_elements(pols) eq 0 and n_elements(all_type_pol) eq 0 then all_type_pol = 1
    
  if keyword_set(all_type_pol) and n_elements(folder_names) lt 2 and n_elements(obs_names_in) lt 2 then $
    message, 'all_type_pol can only be set with 2 folder names and/or 2 obs names'
    
  if not keyword_set(all_type_pol) then begin
    if n_elements(cube_types) eq 0 then if n_elements(folder_names) eq 1 then cube_types = ['dirty', 'res'] else cube_types = 'res'
    if n_elements(pols) eq 0 then pols = 'xx'
  endif else begin
    undefine, cube_types, pols
  endelse
  
  obs_info = hellebore_filenames(folder_names, obs_names_in)
  
  wh_noinfo = where(obs_info.info_files eq '', count_noinfo)
  if count_noinfo gt 0 then message, 'Info files are not all present'
  
  if n_elements(obs_info.folder_names) eq 2 then begin
    save_path = obs_info.diff_save_path
    note = obs_info.diff_note
    plot_path = obs_info.diff_plot_path
  endif else begin
    save_path = obs_info.folder_names[0] + path_sep()
    note = obs_info.fhd_types[0]
    plot_path = obs_info.plot_paths[0]
  endelse
  
  if n_elements(folder_names) eq 1 then begin
    if n_elements(obs_info.obs_names) gt 1 then begin
      plot_filebase = obs_info.fhd_types[0] + '_' + obs_info.obs_names[0] + '_' + cube_types[0] + '_' + pols[0] + $
        '_minus_' + obs_info.obs_names[0] + '_' + cube_types[max_type] + '_' + pols[max_pol]
    endif else begin
      if obs_info.integrated[0] eq 0 then plot_start = obs_info.fhd_types[0] + '_' + obs_info.obs_names[0] else plot_start = obs_info.fhd_types[0]
      
      plot_filebase = plot_start + '_' + cube_types[0] + '_' + pols[0] + $
        '_minus_' + cube_types[max_type] + '_' + pols[max_pol]
    endelse
  endif else begin
    if keyword_set(all_type_pol) then begin
      plot_filebase = obs_info.name_same_parts + '__' + obs_info.name_diff_parts[0]  + $
        '_minus_' + obs_info.name_diff_parts[1]
                
    endif else plot_filebase = obs_info.name_same_parts + '__' + strjoin([obs_info.name_diff_parts[0], cube_types[0], pols[0]], '_')  + $
      '_minus_' + strjoin([obs_info.name_diff_parts[1], cube_types[n_elements(cube_types)-1], pols[n_elements(pols)-1]], '_')
  endelse
  
  if not file_test(plot_path, /directory) then file_mkdir, plot_path
  if not file_test(save_path, /directory) then file_mkdir, save_path
  
  ps_difference_plots, obs_info.info_files, cube_types, pols, all_type_pol = all_type_pol, $
    plot_path = plot_path, plot_filebase = plot_filebase, save_path = save_path, savefilebase = savefilebase, $
    note = note, kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, $
    quiet = quiet, png = png, eps = eps, pdf = pdf
    
end
