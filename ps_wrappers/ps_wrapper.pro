pro ps_wrapper, folder_name_in, obs_name, data_subdirs=data_subdirs, $
    exact_obsnames = exact_obsnames, ps_foldername = ps_foldername, $
    version_test = version_test, $
    ps_folder_branch = ps_folder_branch, copy_master_uvf = copy_master_uvf, $
    no_evenodd = no_evenodd, no_wtvar_rts = no_wtvar_rts, $
    set_data_ranges = set_data_ranges, beamfiles = beamfiles, rts = rts, $
    casa = casa, sim = sim, save_slices = save_slices, no_binning = no_binning, $
    refresh_dft = refresh_dft, refresh_ps = refresh_ps, $
    refresh_binning = refresh_binning, refresh_info = refresh_info, $
    refresh_beam = refresh_beam, dft_fchunk = dft_fchunk, require_radec = require_radec, $
    delta_uv_lambda = delta_uv_lambda, max_uv_lambda = max_uv_lambda, $
    full_image = full_image, image_clip = image_clip, $
    pol_inc = pol_inc, type_inc = type_inc, freq_ch_range = freq_ch_range, $
    freq_flags = freq_flags, freq_flag_name = freq_flag_name, $
    allow_beam_approx = allow_beam_approx, uvf_input = uvf_input, uv_avg = uv_avg, $
    uv_img_clip = uv_img_clip, dft_z_use = dft_z_use, std_power = std_power, $
    inverse_covar_weight = inverse_covar_weight, ave_removal = ave_removal, $
    no_wtd_avg = no_wtd_avg, norm_rts_with_fhd = norm_rts_with_fhd, $
    use_weight_cutoff_sim = use_weight_cutoff_sim, $
    wt_cutoffs = wt_cutoffs, wt_measures = wt_measures, fix_sim_input = fix_sim_input, $
    no_spec_window = no_spec_window, spec_window_type = spec_window_type, $
    no_kzero = no_kzero, plot_slices = plot_slices, slice_type = slice_type, $
    uvf_plot_type = uvf_plot_type, plot_stdset = plot_stdset, plot_1to2d = plot_1to2d, $
    plot_2d_masked = plot_2d_masked, plot_kpar_power = plot_kpar_power, $
    plot_kperp_power = plot_kperp_power, plot_k0_power = plot_k0_power, $
    plot_noise_1d = plot_noise_1d, plot_sim_noise = plot_sim_noise, $
    data_range = data_range, sigma_range = sigma_range, nev_range = nev_range, $
    slice_range = slice_range, snr_range = snr_range, noise_range = noise_range, $
    nnr_range = nnr_range, range_1d = range_1d, color_type = color_type, $
    log_kpar = log_kpar, log_kperp = log_kperp, $
    kpar_bin = kpar_bin, kperp_bin = kperp_bin, $
    log_k1d = log_k1d, k1d_bin = k1d_bin, plot_1d_delta = plot_1d_delta, $
    plot_1d_error_bars = plot_1d_error_bars, plot_1d_nsigma = plot_1d_nsigma, $
    set_krange_1dave = set_krange_1dave, kpar_range_1dave = kpar_range_1dave, $
    bin_arr_3d = bin_arr_3d, kperp_range_1dave = kperp_range_1dave, $
    kperp_range_lambda_1dave = kperp_range_lambda_1dave, kx_range_1dave = kx_range_1dave, $
    kx_range_lambda_1dave = kx_range_lambda_1dave, ky_range_1dave = ky_range_1dave, $
    ky_range_lambda_1dave = ky_range_lambda_1dave, $
    kperp_range_lambda_kparpower = kperp_range_lambda_kparpower, $
    kpar_range_kperppower = kpar_range_kperppower, $
    kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, $
    kperp_plot_range = kperp_plot_range, kperp_lambda_plot_range = kperp_lambda_plot_range, $
    kpar_plot_range = kpar_plot_range, baseline_axis = baseline_axis, $
    delay_axis = delay_axis, cable_length_axis = cable_length_axis, hinv = hinv, $
    plot_wedge_line = plot_wedge_line, wedge_angles = wedge_angles, wedge_names = wedge_names, $
    coarse_harm_width = coarse_harm_width, plot_eor_1d = plot_eor_1d, $
    plot_flat_1d = plot_flat_1d, no_text_1d = no_text_1d, save_path = save_path, $
    savefilebase = savefilebase, plot_path = plot_path, plot_filebase = plot_filebase, $
    individual_plots = individual_plots, plot_binning_hist = plot_binning_hist, $
    note = note, png = png, eps = eps, pdf = pdf, cube_power_info = cube_power_info, $
    no_dft_progress = no_dft_progress, loc_name = loc_name

  if n_elements(folder_name_in) ne 1 then message, 'one folder_name must be supplied.'

  folder_name = get_folder(folder_name_in, loc_name = loc_name,  rts = rts, $
    dirty_folder = dirty_folder)

  if keyword_set(version_test) then begin
    if n_elements(ps_folder_branch) eq 0 then ps_folder_branch = 1
    if n_elements(copy_master_uvf) eq 0 then copy_master_uvf = 1
  endif

  if keyword_set(ps_folder_branch) and n_elements(ps_foldername) eq 0 then begin
    git_info = git_info(ps_repository_dir())
    ps_foldername = 'ps_' + git_info.branch

    if git_info.branch eq 'master' and keyword_set(copy_master_uvf) then copy_master_uvf=0
  endif

  if keyword_set(copy_master_uvf) then begin
    ;; copy over initial files from the ps_master folder if they exist and don't exist for this run
    master_folder_test = file_test(folder_name + path_sep() + 'ps_master', /directory)
    if master_folder_test gt 0 then begin
      ;; info files
      master_ps_folder = folder_name + path_sep() + 'ps_master'
      master_info_files = file_search(master_ps_folder + path_sep() + '*_info.idlsave')

      new_ps_folder = folder_name + path_sep() + ps_foldername
      if file_test(new_ps_folder, /directory) eq 0 then file_mkdir, new_ps_folder

      filebases = file_basename(master_info_files)
      new_files_exist = file_test(new_ps_folder + path_sep() + filebases)

      if min(new_files_exist) eq 0 then begin
        wh_copy = where(new_files_exist eq 0)
        file_copy, master_info_files[wh_copy], new_ps_folder, /require_directory
      endif

      ;; uvf files
      master_uvf_folder = master_ps_folder + path_sep() + 'data' + path_sep() + 'uvf_cubes'
      master_uvf_files = [file_search(master_uvf_folder + path_sep() + '*_uvf.idlsave'), $
                          file_search(master_uvf_folder + path_sep() + '*_radec.idlsave')]


      new_uvf_folder = new_ps_folder + path_sep() + 'data' + path_sep() + 'uvf_cubes'
      if file_test(new_uvf_folder, /directory) eq 0 then file_mkdir, new_uvf_folder

      filebases = file_basename(master_uvf_files)
      new_files_exist = file_test(new_uvf_folder + path_sep() + filebases)

      if min(new_files_exist) eq 0 then begin
        wh_copy = where(new_files_exist eq 0)
        file_copy, master_uvf_files[wh_copy], new_uvf_folder, /require_directory
      endif

      ;; beam files
      master_beam_folder = master_ps_folder + path_sep() + 'data' + path_sep() + 'beam_cubes'
      master_beam_files = file_search(master_beam_folder + path_sep() + '*_beam2.idlsave')

      new_beam_folder = new_ps_folder + path_sep() + 'data' + path_sep() + 'beam_cubes'
      if file_test(new_beam_folder, /directory) eq 0 then file_mkdir, new_beam_folder

      filebases = file_basename(master_beam_files)
      new_files_exist = file_test(new_beam_folder + path_sep() + filebases)
      if min(new_files_exist) eq 0 then begin
        wh_copy = where(new_files_exist eq 0)
        file_copy, master_beam_files[wh_copy], new_beam_folder, /require_directory
      endif
    endif
  endif

  obs_info = ps_filenames(folder_name, obs_name, dirty_folder = dirty_folder, $
    exact_obsnames = exact_obsnames, rts = rts, uvf_input = uvf_input, $
    data_subdirs = data_subdirs, ps_foldernames = ps_foldername, $
    save_paths = save_path, plot_paths = plot_path, refresh_info = refresh_info)

  if not file_test(save_path, /directory) then file_mkdir, save_path
  if not file_test(plot_path, /directory) then file_mkdir, plot_path

  if keyword_set(rts) then begin

    if tag_exist(obs_info, 'dirtyfiles') then dirtyfiles = obs_info.dirtyfiles.(0)

    if obs_info.info_files[0] ne '' then datafile = obs_info.info_files[0] else begin

      if obs_info.cube_files.(0)[0] ne '' then datafile = obs_info.cube_files.(0) else begin
        datafile = rts_fits2idlcube(obs_info.datafiles.(0), obs_info.weightfiles.(0), $
          obs_info.variancefiles.(0), dirtyfiles = dirtyfiles, pol_inc = pol_inc, $
          save_path = obs_info.folder_names[0]+path_sep(), $
          refresh = refresh_dft, no_wtvar = no_wtvar_rts)
      endelse
    endelse

    if keyword_set(refresh_rtscube) then $
      datafile = rts_fits2idlcube(obs_info.datafiles.(0), obs_info.weightfiles.(0), $
        obs_info.variancefiles.(0), dirtyfiles = dirtyfiles, pol_inc = pol_inc, $
        save_path = obs_info.folder_names[0]+path_sep(), /refresh, no_wtvar = no_wtvar_rts)

    if keyword_set(no_wtvar_rts) then return

    max_uv_lambda = 300

    note = obs_info.rts_types

  endif else begin

    if obs_info.info_files[0] ne '' then begin
      datafile = obs_info.info_files[0]
    endif else begin
      datafile = obs_info.cube_files.(0)
    endelse

    if keyword_set(uvf_input) then plot_filebase = obs_info.fhd_types[0] + '_uvf'$
    else  plot_filebase = obs_info.fhd_types[0]

    if stregex(obs_info.fhd_types[0], obs_info.obs_names[0], /boolean) eq 0 then $
      plot_filebase = plot_filebase + '_' + obs_info.obs_names[0]

    note = obs_info.fhd_types[0]
    if ps_foldername ne 'ps' then note = note + '_' + ps_foldername
    if keyword_set(uvf_input) then note = note + '_uvf'

    if tag_exist(obs_info, 'beam_files') then beamfiles = obs_info.beam_files

    if keyword_set(sim) then begin
      plot_eor_1d=1
      if n_elements(range_1d) eq 0 then range_1d = [1e3, 1e7]

      if n_elements(use_weight_cutoff_sim) eq 0 then use_weight_cutoff_sim = 1
      if keyword_set(use_weight_cutoff_sim) then begin
        wt_cutoffs = [0,1]
        wt_measures = strarr(2)+'min'
      endif else begin
        wt_cutoffs = 0
      endelse
    endif

    if n_elements(set_data_ranges) eq 0 and not keyword_set(sim) then set_data_ranges = 1
    if n_elements(set_krange_1dave) eq 0 and not keyword_set(sim) then set_krange_1dave = 1
  endelse

  ;; default to doing binning
  if n_elements(no_binning) eq 0 then no_binning = 0
  ;; default to saving slices
  if n_elements(save_slices) eq 0 then save_slices = 1
  if save_slices eq 0 and keyword_set(plot_slices) then begin
    message, 'save_slices cannot be set to 0 if plot_slices is set'
  endif

  print,'datafile = ' + datafile

  if keyword_set(set_data_ranges) then begin
    if n_elements(range_1d) eq 0 then begin
      if keyword_set(plot_1d_delta) then begin
        range_1d = [1e0, 1e10]
      endif else begin
        range_1d = [1e4, 1e15]
      endelse
    endif

    if tag_exist(obs_info, 'integrated') then begin
      integrated = obs_info.integrated[0]
    endif else begin
      integrated = 1
    endelse

    if integrated gt 0 then begin
      if n_elements(sigma_range) eq 0 then sigma_range = [2e5, 2e9]
      if n_elements(nev_range) eq 0 then nev_range = [2e6, 2e10]
    endif else begin
      if n_elements(sigma_range) eq 0 then sigma_range = [1e7, 2e11]
      if n_elements(nev_range) eq 0 then nev_range = [2e8, 2e12]
    endelse

    if n_elements(data_range) eq 0 then data_range = [1e3, 1e15]
    if n_elements(nnr_range) eq 0 then nnr_range = [1e-1, 1e1]
    if n_elements(snr_range) eq 0 then snr_range = [1e-5, 1e7]

    noise_range = nev_range
  endif else begin
    set_data_ranges = 0
  endelse

  if keyword_set(set_krange_1dave) then begin

    ;; set some default kperp & kpar ranges for 1d averaging
    if n_elements(kperp_range_lambda_1dave) eq 0 then kperp_range_lambda_1dave = [10,50]

    ;; This is the kpar range used for Danny's paper -- generally good for EoR limits but not
    ;; for evaluating wedges physics. Commented out for now.
    ;; kpar_range_1dave = [0.2,10]

  endif

  dft_fchunk = 1

  if keyword_set(plot_binning_hist) then refresh_binning = 1

  refresh_options = create_refresh_options(refresh_dft = refresh_dft, $
    refresh_beam = refresh_beam, refresh_ps = refresh_ps, refresh_kcube = refresh_ps, $
    refresh_binning = refresh_binning, refresh_info = refresh_info)

  uvf_options = create_uvf_options(delta_uv_lambda = delta_uv_lambda, $
    max_uv_lambda = max_uv_lambda, full_image = full_image, image_clip = image_clip, $
    uv_avg = uv_avg, uv_img_clip = uv_img_clip, require_radec = require_radec, $
    dft_fchunk = dft_fchunk, no_dft_progress = no_dft_progress)

  ps_options = create_ps_options(ave_removal = ave_removal, wt_cutoffs = wt_cutoffs, $
    wt_measures = wt_measures, spec_window_type = spec_window_type, $
    no_spec_window = no_spec_window, allow_beam_approx = allow_beam_approx, $
    dft_z_use = dft_z_use, std_power = std_power, no_wtd_avg = no_wtd_avg, $
    inverse_covar_weight = inverse_covar_weight)

  binning_2d_options = create_binning_2d_options(no_kzero = no_kzero, $
    log_kpar = log_kpar, log_kperp = log_kperp, kpar_bin = kpar_bin, kperp_bin = kperp_bin)

  binning_1d_options = create_binning_1d_options(wedge_angles = wedge_angles, $
    wedge_names = wedge_names, $
    coarse_harm_width = coarse_harm_width, log_k = log_k1d, k_bin = k1d_bin, $
    kpar_range_1dave = kpar_range_1dave, kperp_range_1dave = kperp_range_1dave, $
    kperp_range_lambda_1dave = kperp_range_lambda_1dave, kx_range_1dave = kx_range_1dave, $
    kx_range_lambda_1dave = kx_range_lambda_1dave, ky_range_1dave = ky_range_1dave, $
    ky_range_lambda_1dave = ky_range_lambda_1dave, $
    kperp_range_lambda_kparpower = kperp_range_lambda_kparpower, $
    kpar_range_kperppower = kpar_range_kperppower)

  plot_options = create_plot_options(hinv = hinv, plot_path = plot_path, $
    plot_filebase = plot_filebase, note = note, individual_plots = individual_plots, $
    png = png, eps = eps, pdf = pdf)

  plot_types = create_plot_types(plot_stdset = plot_stdset, plot_slices = plot_slices, $
    slice_type = slice_type, uvf_plot_type = uvf_plot_type, plot_1to2d = plot_1to2d, $
    plot_2d_masked = plot_2d_masked, plot_kpar_power = plot_kpar_power, $
    plot_kperp_power = plot_kperp_power, plot_k0_power = plot_k0_power, $
    plot_noise_1d = plot_noise_1d, plot_sim_noise = plot_sim_noise, $
    plot_binning_hist = plot_binning_hist)

  plot_2d_options = create_plot_2d_options(plot_wedge_line = plot_wedge_line, $
    kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, $
    kperp_plot_range = kperp_plot_range, kperp_lambda_plot_range = kperp_lambda_plot_range, $
    kpar_plot_range = kpar_plot_range, baseline_axis = baseline_axis, $
    delay_axis = delay_axis, cable_length_axis = cable_length_axis, $
    data_range = data_range, sigma_range = sigma_range, nev_range = nev_range, $
    snr_range = snr_range, noise_range = noise_range, nnr_range = nnr_range, $
    slice_range = slice_range, color_type = color_type)

  plot_1d_options = create_plot_1d_options(range_1d = range_1d, $
    plot_1d_delta = plot_1d_delta, plot_1d_error_bars = plot_1d_error_bars, $
    plot_1d_nsigma = plot_1d_nsigma, plot_eor_1d = plot_eor_1d, $
    plot_flat_1d = plot_flat_1d, no_text_1d = no_text_1d)

  ps_main_plots, datafile, beamfiles = beamfiles, pol_inc = pol_inc, $
    type_inc = type_inc, freq_ch_range = freq_ch_range, freq_flags = freq_flags, $
    freq_flag_name = freq_flag_name, uvf_input = uvf_input, no_evenodd = no_evenodd, $
    rts = rts, norm_rts_with_fhd = norm_rts_with_fhd, casa = casa, sim = sim, $
    fix_sim_input = fix_sim_input, save_path = save_path, $
    savefilebase = savefilebase, cube_power_info = cube_power_info, $
    save_slices = save_slices, no_binning = no_binning, $
    refresh_options = refresh_options, uvf_options = uvf_options, $
    ps_options = ps_options, binning_2d_options = binning_2d_options, $
    binning_1d_options = binning_1d_options, plot_options = plot_options, $
    plot_types = plot_types, plot_2d_options = plot_2d_options, $
    plot_1d_options = plot_1d_options

  if not keyword_set(set_data_ranges) then begin
    if n_elements(data_range) ne 0 then $
      print, 'data_range used: ', number_formatter(data_range, format = '(e7.1)')
    if n_elements(sigma_range) ne 0 then $
      print, 'sigma_range used: ', number_formatter(sigma_range, format = '(e7.1)')
    if n_elements(nev_range) ne 0 then $
      print, 'nev_range used: ', number_formatter(nev_range, format = '(e7.1)')
    if n_elements(nnr_range) ne 0 then $
      print, 'nnr_range used: ', number_formatter(nnr_range, format = '(e7.1)')
    if n_elements(snr_range) ne 0 then $
      print, 'snr_range used: ', number_formatter(snr_range, format = '(e7.1)')
    if n_elements(noise_range) ne 0 then $
      print, 'noise_range used: ', number_formatter(noise_range, format = '(e7.1)')
  endif

end
