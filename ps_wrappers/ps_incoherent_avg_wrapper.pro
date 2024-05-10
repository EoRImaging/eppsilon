pro ps_incoherent_avg_wrapper, folder_names_in, obs_names_in, $
  exact_obsnames = exact_obsnames, ps_foldernames = ps_foldernames, $
  version_test = version_test, loc_name = loc_name, data_subdirs = data_subdirs, $
  new_obsname = new_obsname, save_path = avg_save_path, $
  set_data_ranges = set_data_ranges, sim = sim, save_slices = save_slices, save_sum_cube = save_sum_cube, $
  no_binning = no_binning, refresh_dft = refresh_dft, refresh_ps = refresh_ps, $
  refresh_binning = refresh_binning, refresh_info = refresh_info, $
  refresh_beam = refresh_beam, dft_fchunk = dft_fchunk, require_radec = require_radec, $
  delta_uv_lambda = delta_uv_lambda, max_uv_lambda = max_uv_lambda, $
  full_image = full_image, image_clip = image_clip, $
  pol_inc = pol_inc, type_inc = type_inc, freq_ch_range = freq_ch_range, $
  freq_flags = freq_flags, freq_flag_name = freq_flag_name, $
  freq_flag_repeat = freq_flag_repeat, $
  freq_avg_factor = freq_avg_factor, force_even_freqs = force_even_freqs, $
  freq_avg_bins = freq_avg_bins, freq_bin_name = freq_bin_name, $
  allow_beam_approx = allow_beam_approx, uvf_input = uvf_input, uv_avg = uv_avg, $
  uv_img_clip = uv_img_clip, freq_dft = freq_dft, dft_z_use = dft_z_use, $
  kz_use = kz_use, kzuse_name = kzuse_name, std_power = std_power, $
  inverse_covar_weight = inverse_covar_weight, ave_removal = ave_removal, $
  no_wtd_avg = no_wtd_avg, use_weight_cutoff_sim = use_weight_cutoff_sim, $
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
  plot_flat_1d = plot_flat_1d, no_text_1d = no_text_1d, $
  savefilebase = savefilebase, plot_path = avg_plot_path, plot_filebase = plot_filebase, $
  individual_plots = individual_plots, plot_binning_hist = plot_binning_hist, $
  note = note, png = png, eps = eps, pdf = pdf, cube_power_info = cube_power_info

  compile_opt idl2, logical_predicate

  if n_elements(folder_names_in) eq 0 then message, 'at least 1 folder name must be specified'
  folder_names = get_folder(folder_names_in, loc_name = loc_name)

  obs_info = ps_filenames(folder_names, obs_names_in, $
    exact_obsnames = exact_obsnames, uvf_input = uvf_input, $
    data_subdirs = data_subdirs, ps_foldernames = ps_foldernames)

  if keyword_set(version_test) and obs_info.combined_ps_folder eq "ps" $
      and n_elements(folder_names_in) eq 1 then begin
    git_info = git_info(ps_repository_dir())
    avg_ps_foldername = 'ps_' + git_info.branch
  endif else begin
    avg_ps_foldername = obs_info.combined_ps_folder
  endelse

  _ = where(obs_info.info_files eq '', count_noinfo)
  if count_noinfo gt 0 then message, 'Info files are not all present'

  if n_elements(save_path) eq 0 then begin
    avg_save_path = obs_info.diff_save_path + avg_ps_foldername + path_sep()
  endif

  if n_elements(avg_plot_path) eq 0 then begin
    avg_plot_path = obs_info.diff_plot_path + avg_ps_foldername + path_sep()
  endif

  if n_elements(new_obsname) eq 0 then begin
    new_obsname = obs_info.combined_obsname
  endif

  if not file_test(avg_save_path, /directory) then file_mkdir, avg_save_path
  if not file_test(avg_plot_path, /directory) then file_mkdir, avg_plot_path

  ;; default to doing binning
  if n_elements(no_binning) eq 0 then no_binning = 0
  ;; default to saving slices
  if n_elements(save_slices) eq 0 then save_slices = 1
  if save_slices eq 0 and keyword_set(plot_slices) then begin
    message, 'save_slices cannot be set to 0 if plot_slices is set'
  endif

  if keyword_set(uvf_input) then plot_filebase = obs_info.fhd_types[0] + '_uvf'$
    else plot_filebase = obs_info.fhd_types[0]

  if stregex(obs_info.fhd_types[0], obs_info.obs_names[0], /boolean) eq 0 then $
    plot_filebase = plot_filebase + '_' + obs_info.obs_names[0]

  note = obs_info.fhd_types[0]

  if obs_info.combined_ps_folder ne 'ps' then note = note + '_' + obs_info.combined_ps_folder
  if keyword_set(uvf_input) then note = note + '_uvf'

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

  if keyword_set(plot_binning_hist) then begin
    refresh_binning = 1
  endif

  refresh_options = create_refresh_options(refresh_dft = refresh_dft, $
    refresh_beam = refresh_beam, refresh_ps = refresh_ps, refresh_kcube = refresh_ps, $
    refresh_freq_select_avg = refresh_ps, $
    refresh_binning = refresh_binning, refresh_info = refresh_info)

  uvf_options = create_uvf_options(delta_uv_lambda = delta_uv_lambda, $
    max_uv_lambda = max_uv_lambda, full_image = full_image, image_clip = image_clip, $
    uv_avg = uv_avg, uv_img_clip = uv_img_clip, require_radec = require_radec, $
    dft_fchunk = dft_fchunk, no_dft_progress = no_dft_progress)

  freq_options = create_freq_options( $
    freq_ch_range = freq_ch_range, $
    freq_flags = freq_flags, $
    freq_flag_name = freq_flag_name, $
    freq_flag_repeat = freq_flag_repeat, $
    freq_avg_factor = freq_avg_factor, $
    force_even_freqs = force_even_freqs, $
    freq_avg_bins = freq_avg_bins, $
    freq_bin_name = freq_bin_name)

  ps_options = create_ps_options(ave_removal = ave_removal, wt_cutoffs = wt_cutoffs, $
    wt_measures = wt_measures, spec_window_type = spec_window_type, $
    no_spec_window = no_spec_window, allow_beam_approx = allow_beam_approx, $
    save_sum_cube = save_sum_cube, freq_dft = freq_dft, $
    dft_z_use = dft_z_use, kz_use = kz_use, kzuse_name = kzuse_name, $
    std_power = std_power, no_wtd_avg = no_wtd_avg, $
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

  ; do the integration
  ps_power_incoherent_avg, obs_info.info_files, avg_file_struct=avg_file_struct, $
    pol_inc = pol_inc, freq_options = freq_options, uvf_input = uvf_input, $
    new_obsname=new_obsname, avg_save_path=avg_save_path, $
    savefilebase = savefilebase, save_slices = save_slices, $
    refresh_options = refresh_options, uvf_options = uvf_options, $
    ps_options = ps_options

  ; do the binning and make the plots
  ps_main_plots, avg_file_struct, incoherent_avg=1, pol_inc = pol_inc, $
    type_inc = type_inc, freq_options = freq_options, $
    uvf_input = uvf_input, no_evenodd = no_evenodd, $
    rts = rts, norm_rts_with_fhd = norm_rts_with_fhd, casa = casa, sim = sim, $
    fix_sim_input = fix_sim_input, save_path = avg_save_path, $
    savefilebase = savefilebase, cube_power_info = cube_power_info, $
    save_slices = save_slices, no_binning = no_binning, $
    refresh_options = refresh_options, uvf_options = uvf_options, $
    ps_options = ps_options, binning_2d_options = binning_2d_options, $
    binning_1d_options = binning_1d_options, plot_options = plot_options, $
    plot_types = plot_types, plot_2d_options = plot_2d_options, $
    plot_1d_options = plot_1d_options


end