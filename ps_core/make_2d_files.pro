pro make_2d_files, nfiles, savefile_2d, savefile_k0, power_3D, sim_noise_3D, $
    new_noise_3d, noise_expval_3d, weights_3d, kx_mpc, ky_mpc, kz_mpc, $
    delay_params, hubble_param, hinv, kperp_lambda_conv, freq_mask, vs_name, $
    vs_mean, t_sys_meas, window_int, git_hashes, wt_ave_power, ave_power, $
    ave_weights, wt_ave_power_freq, ave_power_freq, wt_ave_power_uvf, $
    ave_power_uvf, uv_pix_area, uv_area, masked_save_items = masked_save_items, $
    noise_3D = noise_3D, sim_noise_diff_3D = sim_noise_diff_3D, $
    wt_measure = wt_measure, wt_cutoff = wt_cutoff, log_kpar = log_kpar, $
    log_kperp = log_kperp, kperp_bin = kperp_bin, kpar_bin = kpar_bin, $
    fill_holes = fill_holes, freq_flags = freq_flags

  power_rebin = kspace_rebinning_2d(power_3D, kx_mpc, ky_mpc, kz_mpc, $
    kperp_edges_mpc, kpar_edges_mpc, log_kpar = log_kpar, log_kperp = log_kperp, $
    kperp_bin = kperp_bin, kpar_bin = kpar_bin, noise_expval = noise_expval_3d, $
    binned_noise_expval = binned_noise_expval, weights = weights_3d, $
    nbins_2d = nbins_2d, binned_weights = binned_weights, fill_holes = fill_holes, $
    kperp_density_measure = wt_measure, kperp_density_cutoff = wt_cutoff)

  sim_noise_rebin = kspace_rebinning_2d(sim_noise_3D, kx_mpc, ky_mpc, kz_mpc, $
    kperp_edges_mpc, kpar_edges_mpc, log_kpar = log_kpar, log_kperp = log_kperp, $
    kperp_bin = kperp_bin, kpar_bin = kpar_bin, noise_expval = noise_expval_3d, $
    binned_noise_expval = binned_noise_expval, weights = weights_3d, $
    var_power_2d = sim_noise_var_2d, mean_var_2d = mean_var_2d, $
    binned_weights = binned_weights, fill_holes = fill_holes, $
    kperp_density_measure = wt_measure, kperp_density_cutoff = wt_cutoff)

  new_noise_rebin = kspace_rebinning_2d(new_noise_3d, kx_mpc, ky_mpc, kz_mpc, $
    kperp_edges_mpc, kpar_edges_mpc, log_kpar = log_kpar, log_kperp = log_kperp, $
    kperp_bin = kperp_bin, kpar_bin = kpar_bin, noise_expval = noise_expval_3d, $
    binned_noise_expval = binned_noise_expval, weights = weights_3d, $
    var_power_2d = new_noise_var_2d, mean_var_2d = new_mean_var_2d, $
    binned_weights = binned_weights, fill_holes = fill_holes, $
    kperp_density_measure = wt_measure, kperp_density_cutoff = wt_cutoff)

  if nfiles eq 2 then begin
    noise_rebin = kspace_rebinning_2d(noise_3D, kx_mpc, ky_mpc, kz_mpc, $
      kperp_edges_mpc, kpar_edges_mpc, log_kpar = log_kpar, log_kperp = log_kperp, $
      kperp_bin = kperp_bin, kpar_bin = kpar_bin, noise_expval = noise_expval_3d, $
      binned_noise_expval = binned_noise_expval, weights = weights_3d, $
      binned_weights = binned_weights, fill_holes = fill_holes, $
      kperp_density_measure = wt_measure, kperp_density_cutoff = wt_cutoff)

    sim_noise_diff_rebin = kspace_rebinning_2d(sim_noise_diff_3D, kx_mpc, ky_mpc, kz_mpc, $
      kperp_edges_mpc, kpar_edges_mpc, log_kpar = log_kpar, log_kperp = log_kperp, $
      kperp_bin = kperp_bin, kpar_bin = kpar_bin, noise_expval = noise_expval_3d, $
      binned_noise_expval = binned_noise_expval, weights = weights_3d, $
      binned_weights = binned_weights, fill_holes = fill_holes, $
      kperp_density_measure = wt_measure, kperp_density_cutoff = wt_cutoff)
  endif

  power = power_rebin
  sim_noise = sim_noise_rebin
  if nfiles eq 2 then begin
    noise = noise_rebin
    sim_noise_diff = sim_noise_diff_rebin
  endif
  kperp_edges = kperp_edges_mpc
  kpar_edges = kpar_edges_mpc
  weights = binned_weights
  noise_expval = binned_noise_expval

  wh_good_kperp = where(total(weights, 2) gt 0, count)
  if count eq 0 then begin
    print, '2d weights appear to be entirely zero'
    return
  endif

  if n_elements(masked_save_items) then begin

    kperp_range = masked_save_items.kperp_range
    kperp_range_lambda = masked_save_items.kperp_range_lambda
    kpar_range = masked_save_items.kpar_range
    kx_range = masked_save_items.kx_range
    kx_range_lambda = masked_save_items.kx_range_lambda
    ky_range = masked_save_items.ky_range
    ky_range_lambda = masked_save_items.ky_range_lambda
    if tag_exist(masked_save_items, 'coarse_harm0') then begin
      coarse_harm0 = masked_save_items.coarse_harm0
      coarse_width = masked_save_items.coarse_width
    endif

    if n_elements(freq_flags) ne 0 then begin
      save, file = savefile_2d, power, noise, sim_noise, sim_noise_diff, weights, $
        noise_expval, kperp_edges, kpar_edges, kperp_bin, kpar_bin, $
        kperp_lambda_conv, delay_params, hubble_param, freq_mask, vs_name, vs_mean, $
        t_sys_meas, window_int, git_hashes, wt_ave_power, ave_power, ave_weights, $
        wt_ave_power_freq, ave_power_freq, wt_ave_power_uvf, ave_power_uvf, $
        uv_pix_area, uv_area, kperp_range, kperp_range_lambda, kpar_range, $
        kx_range, kx_range_lambda, ky_range, ky_range_lambda, coarse_harm0, coarse_width
    endif else begin
      save, file = savefile_2d, power, noise, sim_noise, sim_noise_diff, weights, $
        noise_expval, kperp_edges, kpar_edges, kperp_bin, kpar_bin, $
        kperp_lambda_conv, delay_params, hubble_param, vs_name, vs_mean, $
        t_sys_meas, window_int, git_hashes, wt_ave_power, ave_power, ave_weights, $
        wt_ave_power_freq, ave_power_freq, wt_ave_power_uvf, ave_power_uvf, $
        uv_pix_area, uv_area, kperp_range, kperp_range_lambda, kpar_range, $
        kx_range, kx_range_lambda, ky_range, ky_range_lambda, coarse_harm0, coarse_width
    endelse

  endif else begin
    if n_elements(freq_flags) ne 0 then begin
      save, file = savefile_2d, power, noise, sim_noise, sim_noise_diff, weights, $
        noise_expval, kperp_edges, kpar_edges, kperp_bin, kpar_bin, $
        kperp_lambda_conv, delay_params, hubble_param, freq_mask, vs_name, vs_mean, $
        t_sys_meas, window_int, git_hashes, wt_ave_power, ave_power, ave_weights, $
        wt_ave_power_freq, ave_power_freq, wt_ave_power_uvf, ave_power_uvf, $
        uv_pix_area, uv_area
    endif else begin
      save, file = savefile_2d, power, noise, sim_noise, sim_noise_diff, weights, $
        noise_expval, kperp_edges, kpar_edges, kperp_bin, kpar_bin, $
        kperp_lambda_conv, delay_params, hubble_param, vs_name, vs_mean, $
        t_sys_meas, window_int, git_hashes, wt_ave_power, ave_power, ave_weights, $
        wt_ave_power_freq, ave_power_freq, wt_ave_power_uvf, ave_power_uvf, $
        uv_pix_area, uv_area
    endelse
  endelse

  ;; save just k0 line for plotting purposes
  if not keyword_set(no_kzero) then begin
    power = power[*,0]
    sim_noise = sim_noise[*,0]
    if n_elements(noise) gt 0 then noise = noise[*,0]
    if nfiles eq 2 and n_elements(sim_noise) gt 0 then sim_noise_diff = sim_noise_diff[*,0]
    weights = weights[*,0]
    noise_expval = noise_expval[*,0]

    k_edges = kperp_edges
    k_bin = kperp_bin

    textfile = strmid(savefile_k0, 0, stregex(savefile_k0, '.idlsave')) + '.txt'
    print, 'saving kpar=0 power to ' + textfile
    save_1D_text, textfile, k_edges, power, weights, noise_expval, hubble_param, $
      noise, sim_noise_power = sim_noise, sim_noise_diff = sim_noise_diff, $
      nfiles = nfiles, hinv = hinv

    if n_elements(freq_flags) ne 0 then begin
      save, file = savefile_k0, power, noise, sim_noise, sim_noise_diff, weights, $
        noise_expval, k_edges, k_bin, kperp_lambda_conv, delay_params, $
        hubble_param, freq_mask, window_int, wt_ave_power, ave_power, ave_weights, $
        uv_pix_area, uv_area, git_hashes
    endif else begin
      save, file = savefile_k0, power, noise, sim_noise, sim_noise_diff, weights, $
        noise_expval, k_edges, k_bin, kperp_lambda_conv, delay_params, hubble_param, $
        window_int, wt_ave_power, ave_power, ave_weights, uv_pix_area, uv_area, git_hashes
    endelse
  endif
end
