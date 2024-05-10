pro ps_cube_images, folder_names_in, obs_names_in, exact_obsnames = exact_obsnames, $
    data_subdirs=data_subdirs, cube_types = cube_types, pols = pols, evenodd = evenodd, $
    rts = rts, casa = casa, png = png, eps = eps, pdf = pdf, $
    slice_range = slice_range, sr2 = sr2, nvis_norm = nvis_norm, ratio = ratio, $
    diff_ratio = diff_ratio, diff_frac = diff_frac, log = log, data_range = data_range, $
    data_min_abs = data_min_abs, color_profile = color_profile, sym_color = sym_color, $
    window_num = window_num, plot_as_map = plot_as_map, plot_path = plot_path, refresh_info=refresh_info

  ;; default to normalizing by nvis
  if n_elements(nvis_norm) eq 0 then nvis_norm = 1

  if n_elements(folder_names_in) gt 2 then message, 'No more than 2 folder_names can be supplied'
  if n_elements(evenodd) eq 0 then evenodd = 'even'
  if n_elements(evenodd) gt 2 then message, 'No more than 2 evenodd values can be supplied'
  if n_elements(obs_names_in) gt 2 then message, 'No more than 2 obs_names can be supplied'

  folder_names = get_folder(folder_names_in, loc_name = loc_name,  rts = rts, $
    dirty_folder = dirty_folder)

  obs_info = ps_filenames(folder_names, obs_names_in, dirty_folder = dirty_folder, $
    exact_obsnames = exact_obsnames, rts = rts, uvf_input = uvf_input, casa = casa, $
    data_subdirs = data_subdirs, save_paths = save_paths, plot_paths = plot_paths, $
    refresh_info = refresh_info, no_wtvar_rts = no_wtvar_rts)

  if n_elements(plot_path) eq 0 then begin
    if tag_exist(obs_info, 'diff_plot_path') then begin
      plot_path = obs_info.diff_plot_path
    endif else begin
      plot_path = obs_info.plot_paths[0]
    endelse
  endif

  if keyword_set(rts) then begin

    if obs_info.info_files[0] ne '' then begin
      datafile = obs_info.info_files[0]
    endif else begin
      if obs_info.cube_files.(0)[0] ne '' then begin
         datafile = obs_info.cube_files.(0)
      endif else begin
         datafile = rts_fits2idlcube(obs_info.datafiles.(0), obs_info.weightfiles.(0), $
         obs_info.variancefiles.(0), pol_inc, save_path = obs_info.folder_names[0]+path_sep(), $
         refresh = refresh_dft, no_wtvar = no_wtvar_rts)
      endelse
    endelse

    if keyword_set(refresh_rtscube) then begin
      datafile = rts_fits2idlcube(obs_info.datafiles.(0), obs_info.weightfiles.(0), $
        obs_info.variancefiles.(0), pol_inc, save_path = obs_info.folder_names[0]+path_sep(), $
        /refresh, no_wtvar = no_wtvar_rts)
    endif

    if keyword_set(no_wtvar_rts) then stop

  endif

  cube_images, folder_names, obs_info, nvis_norm = nvis_norm, pols = pols, $
    cube_types = cube_types, evenodd = evenodd, rts = rts, $
    png = png, eps = eps, pdf = pdf, slice_range = slice_range, sr2 = sr2, $
    ratio = ratio, diff_ratio = diff_ratio, diff_frac = diff_frac, $
    log = log, data_range = data_range, data_min_abs = data_min_abs, $
    color_profile = color_profile, sym_color = sym_color, $
    window_num = window_num, plot_as_map = plot_as_map, plot_path = plot_pat
end
