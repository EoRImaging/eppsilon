pro fhd_3dps, file_struct, refresh = refresh, kcube_refresh = kcube_refresh, dft_refresh_data = dft_refresh_data, $
    dft_refresh_weight = dft_refresh_weight, dft_ian = dft_ian, cut_image = cut_image, image = image, $
    dft_fchunk = dft_fchunk, freq_ch_range = freq_ch_range, freq_flags = freq_flags, $
    spec_window_type = spec_window_type, delta_uv_lambda = delta_uv_lambda, $
    std_power = std_power, no_kzero = no_kzero, log_kpar = log_kpar, $
    log_kperp = log_kperp, kperp_bin = kperp_bin, kpar_bin = kpar_bin, log_k1d = log_k1d, k1d_bin = k1d_bin, $
    input_units = input_units, fill_holes = fill_holes, quiet = quiet
    
  if tag_exist(file_struct, 'nside') ne 0 then healpix = 1 else healpix = 0
  
  nfiles = n_elements(file_struct.datafile)
  
  if healpix and (keyword_set(dft_refresh_data) or keyword_set(dft_refresh_weight)) then kcube_refresh=1
  if keyword_set(kcube_refresh) then refresh = 1
  
  if n_elements(fill_holes) eq 0 then fill_holes = 0
  
  if n_elements(freq_flags) ne 0 then freq_mask = file_struct.freq_mask
  
  test_powersave = file_test(file_struct.power_savefile) *  (1 - file_test(file_struct.power_savefile, /zero_length))
  
  if test_powersave eq 1 and n_elements(freq_flags) ne 0 then begin
    old_freq_mask = getvar_savefile(file_struct.power_savefile, 'freq_mask')
    if total(abs(old_freq_mask - freq_mask)) ne 0 then test_powersave = 0
  endif
  
  if test_powersave eq 0 or keyword_set(refresh) then begin
  
    test_kcube = file_test(file_struct.kcube_savefile) *  (1 - file_test(file_struct.kcube_savefile, /zero_length))
    
    if test_kcube eq 1 and n_elements(freq_flags) ne 0 then begin
      old_freq_mask = getvar_savefile(file_struct.kcube_savefile, 'freq_mask')
      if total(abs(old_freq_mask - freq_mask)) ne 0 then test_kcube = 0
    endif
    
    if test_kcube eq 0 or keyword_set(kcube_refresh) then $
      fhd_kcube, file_struct, dft_refresh_data = dft_refresh_data, dft_refresh_weight = dft_refresh_weight, dft_ian = dft_ian, $
      dft_fchunk = dft_fchunk, freq_ch_range = freq_ch_range, freq_flags = freq_flags, $
      cut_image = cut_image, image = image, delta_uv_lambda = delta_uv_lambda, $
      spec_window_type = spec_window_type, std_power = std_power, input_units = input_units, /quiet
      
    if nfiles eq 1 then begin
      restore, file_struct.kcube_savefile
      
      n_kx = n_elements(kx_mpc)
      n_ky = n_elements(ky_mpc)
      n_kz = n_elements(kz_mpc)
      
      ;; now construct weights for power (mag. squared) = 1/power variance
      power_weights1 = 1d/(4*(sigma2_1)^2d)
      wh_sig1_0 = where(sigma2_1^2d eq 0, count_sig1_0)
      if count_sig1_0 ne 0 then power_weights1[wh_sig1_0] = 0
      term1 = abs(data_sum_1)^2.*power_weights1
      undefine, data_sum_1
      
      power_weights2 = 1d/(4*(sigma2_2)^2d) ;; inverse variance
      wh_sig2_0 = where(sigma2_2^2d eq 0, count_sig2_0)
      if count_sig2_0 ne 0 then power_weights2[wh_sig2_0] = 0
      term2 = abs(data_sum_2)^2.*power_weights2
      undefine, data_sum_2
      
      noise_expval_3d = sqrt(power_weights1) + sqrt(power_weights2)
      
      weights_3d = (power_weights1 + power_weights2) ;; variance_3d = 1/weights_3d
      undefine, power_weights1, power_weights2
      
      power_3d = (term1 + term2) / weights_3d
      noise_expval_3d = noise_expval_3d / weights_3d
      wh_wt0 = where(weights_3d eq 0, count_wt0)
      if count_wt0 ne 0 then begin
        power_3d[wh_wt0] = 0
        noise_expval_3d[wh_wt0] = 0
      endif
      
      undefine, term1, term2
      
    endif else begin
      ;; nfiles=2
      restore, file_struct.kcube_savefile
      n_kx = n_elements(kx_mpc)
      n_ky = n_elements(ky_mpc)
      n_kz = n_elements(kz_mpc)
      
      ;; now construct weights for power (mag. squared) = 1/power variance
      power_weights1 = 1d/(4*(sigma2_1)^2d)
      wh_sig1_0 = where(sigma2_1^2d eq 0, count_sig1_0)
      if count_sig1_0 ne 0 then power_weights1[wh_sig1_0] = 0
      undefine, sigma2_1
      
      power_weights2 = 1d/(4*(sigma2_2)^2d) ;; inverse variance
      wh_sig2_0 = where(sigma2_2^2d eq 0, count_sig2_0)
      if count_sig2_0 ne 0 then power_weights2[wh_sig2_0] = 0
      undefine, sigma2_1
      
      term1 = (abs(data_sum_1)^2. - abs(data_diff_1)^2.) * power_weights1
      term2 = (abs(data_sum_2)^2. - abs(data_diff_2)^2.) * power_weights2
      undefine, data_sum_1, data_sum_2
      
      noise_t1 = abs(data_diff_1)^2. * power_weights1
      noise_t2 = abs(data_diff_2)^2. * power_weights2
      undefine, data_diff_1, data_diff_2
      
      noise_expval_3d = sqrt(power_weights1) + sqrt(power_weights2)
      
      weights_3d = power_weights1 + power_weights2 ;; variance_3d = 1/weights_3d
      undefine, power_weights1, power_weights2
      
      ;; divide by 4 on power b/c otherwise it would be 4*Re(even-odd crosspower)
      power_3d = (term1 + term2) / (4. * weights_3d)
      noise_3d = (noise_t1 + noise_t2) / weights_3d
      noise_expval_3d= noise_expval_3d / weights_3d
      undefine, term1, term2, noise_t1, noise_t2
      
      wh_wt0 = where(weights_3d eq 0, count_wt0)
      if count_wt0 ne 0 then begin
        power_3d[wh_wt0] = 0
        noise_expval_3d[wh_wt0] = 0
        noise_3d[wh_wt0] = 0
      endif
      
    endelse
    
    save, file = file_struct.power_savefile, power_3d, noise_3d, noise_expval_3d, weights_3d, $
      kx_mpc, ky_mpc, kz_mpc, kperp_lambda_conv, delay_params, hubble_param, n_freq_contrib, freq_mask, vs_name, vs_mean
      
    write_ps_fits, file_struct.fits_power_savefile, power_3d, weights_3d, noise_expval_3d, noise_3d = noise_3d, $
      kx_mpc, ky_mpc, kz_mpc, kperp_lambda_conv, delay_params, hubble_param
      
  endif else restore, file_struct.power_savefile
  
  print, 'power integral:', total(power_3d)
  
  n_kx = n_elements(kx_mpc)
  n_ky = n_elements(ky_mpc)
  n_kz = n_elements(kz_mpc)
  
  if keyword_set (no_kzero) then begin
    ;; leave out kz=0 -- full of foregrounds
    kz_mpc = kz_mpc[1:*]
    power_3d = temporary(power_3d[*, *, 1:*])
    weights_3d = temporary(weights_3d[*,*,1:*])
    noise_expval_3d = temporary(noise_expval_3d[*,*,1:*])
    if nfiles eq 2 then noise_3d = temporary(noise_3d[*,*,1:*])
    n_kz = n_elements(kz_mpc)
  endif
  
  power_tag = file_struct.power_tag
  
  fadd_2dbin = ''
  ;;if keyword_set(fill_holes) then fadd_2dbin = fadd_2dbin + '_nohole'
  if keyword_set(no_kzero) then fadd_2dbin = fadd_2dbin + '_nok0'
  if keyword_set(log_kpar) then fadd_2dbin = fadd_2dbin + '_logkpar'
  if keyword_set(log_kperp) then fadd_2dbin = fadd_2dbin + '_logkperp'
  
  fadd_1dbin = ''
  if keyword_set(log_k) then fadd_1dbin = fadd_1dbin + '_logk'
  
  savefile = file_struct.savefile_froot + file_struct.savefilebase + power_tag + fadd_2dbin + '_2dkpower.idlsave'
  
  print, 'Binning to 2D power spectrum'
  
  power_rebin = kspace_rebinning_2d(power_3D, kx_mpc, ky_mpc, kz_mpc, kperp_edges_mpc, kpar_edges_mpc, log_kpar = log_kpar, $
    log_kperp = log_kperp, kperp_bin = kperp_bin, kpar_bin = kpar_bin, $
    noise_expval = noise_expval_3d, binned_noise_expval = binned_noise_expval, weights = weights_3d, $
    binned_weights = binned_weights, fill_holes = fill_holes)
    
    
  if nfiles eq 2 then $
    noise_rebin = kspace_rebinning_2d(noise_3D, kx_mpc, ky_mpc, kz_mpc, kperp_edges_mpc, kpar_edges_mpc, log_kpar = log_kpar, $
    log_kperp = log_kperp, kperp_bin = kperp_bin, kpar_bin = kpar_bin, $
    noise_expval = noise_expval_3d, binned_noise_expval = binned_noise_expval, $
    weights = weights_3d, binned_weights = binned_weights, fill_holes = fill_holes)
    
    
  power = power_rebin
  if nfiles eq 2 then noise = noise_rebin
  kperp_edges = kperp_edges_mpc
  kpar_edges = kpar_edges_mpc
  weights = binned_weights
  noise_expval = binned_noise_expval
  
  wh_good_kperp = where(total(weights, 2) gt 0, count)
  if count eq 0 then stop
  kperp_plot_range = [min(kperp_edges[wh_good_kperp]), max(kperp_edges[wh_good_kperp+1])]
  
  save, file = savefile, power, noise, weights, noise_expval, kperp_edges, kpar_edges, kperp_bin, kpar_bin, $
    kperp_lambda_conv, delay_params, hubble_param, freq_mask, vs_name, vs_mean
    
  if not keyword_set(quiet) then begin
    kpower_2d_plots, savefile, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
      data_range = data_range
    kpower_2d_plots, savefile, /plot_weights, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
      window_num = 2, title = 'Weights'
  endif
  
  ;; now do slices
  yslice_savefile = file_struct.savefile_froot + file_struct.savefilebase + power_tag + '_xz_plane.idlsave'
  yslice_power = kpower_slice(power_3d, kx_mpc, ky_mpc, kz_mpc, kperp_lambda_conv, delay_params, hubble_param, noise_3d = noise_3d, $
    noise_expval_3d = noise_expval_3d, weights_3d = weights_3d, slice_axis = 1, slice_inds = 0, $
    slice_savefile = yslice_savefile)
    
  xslice_savefile = file_struct.savefile_froot + file_struct.savefilebase + power_tag + '_yz_plane.idlsave'
  xslice_power = kpower_slice(power_3d, kx_mpc, ky_mpc, kz_mpc, kperp_lambda_conv, delay_params, hubble_param, noise_3d = noise_3d, $
    noise_expval_3d = noise_expval_3d, weights_3d = weights_3d, slice_axis = 0, slice_inds = n_kx/2, $
    slice_savefile = xslice_savefile)
  if max(xslice_power) eq 0 then begin
    nloop = 0
    while max(xslice_power) eq 0 do begin
      nloop = nloop+1
      xslice_power = kpower_slice(power_3d, kx_mpc, ky_mpc, kz_mpc, kperp_lambda_conv, delay_params, hubble_param, $
        noise_3d = noise_3d, noise_expval_3d = noise_expval_3d, weights_3d = weights_3d, slice_axis = 0, $
        slice_inds = n_kx/2+nloop, slice_savefile = xslice_savefile)
    endwhile
  endif
  
  zslice_savefile = file_struct.savefile_froot + file_struct.savefilebase + power_tag + '_xy_plane.idlsave'
  zslice_power = kpower_slice(power_3d, kx_mpc, ky_mpc, kz_mpc, kperp_lambda_conv, delay_params, hubble_param, noise_3d = noise_3d, $
    noise_expval_3d = noise_expval_3d, weights_3d = weights_3d, slice_axis = 2, slice_inds = 1, $
    slice_savefile = zslice_savefile)
    
    
  print, 'Binning to 1D power spectrum'
  
  power_1d = kspace_rebinning_1d(power_3d, kx_mpc, ky_mpc, kz_mpc, k_edges_mpc, k_bin = k1d_bin, log_k = log_k1d, $
    noise_expval = noise_expval_3d, binned_noise_expval = noise_expval_1d, weights = weights_3d, $
    binned_weights = weights_1d, mask = mask, pixelwise_mask = pixelwise_mask, k1_mask = k1_mask, $
    k2_mask = k2_mask,  k3_mask = k3_mask)
    
  if nfiles eq 2 then $
    noise_1d = kspace_rebinning_1d(noise_3d, kx_mpc, ky_mpc, kz_mpc, k_edges_mpc, k_bin = k1d_bin, log_k = log_k1d, $
    noise_expval = noise_expval_3d, binned_noise_expval = noise_expval_1d, weights = weights_3d, $
    binned_weights = weights_1d, mask = mask, pixelwise_mask = pixelwise_mask, k1_mask = k1_mask, $
    k2_mask = k2_mask,  k3_mask = k3_mask)
    
  power = power_1d
  if nfiles eq 2 then noise = noise_1d
  weights = weights_1d
  k_edges = k_edges_mpc
  k_bin = k1d_bin
  noise_expval = noise_expval_1d
  
  fadd_1d = ''
  if keyword_set(log_k) then fadd_1d = fadd_1d + '_logk'
  
  savefile = file_struct.savefile_froot + file_struct.savefilebase + power_tag + fadd_1dbin + '_1dkpower.idlsave'
  save, file = savefile, power, noise, weights, noise_expval, k_edges, k_bin, hubble_param, freq_mask
  
  if not keyword_set(quiet) then begin
    kpower_1d_plots, savefile, window_num = 5
  endif
  
  ;; eor_file_1d = base_path() + 'power_spectrum/eor_data/eor_power_1d.idlsave'
  ;; file_arr = [savefile, eor_file_1d]
  ;; if keyword_set(eor_only) then begin
  ;;    if keyword_set(eor_test) then names_arr = 'Input EoR' else names_arr = 'Simulated EoR'
  ;; endif else names_arr = 'Simulation PS'
  ;; names_arr = [names_arr, 'EoR signal']
  ;; colors_arr = [0, 254]
  
  if not keyword_set(quiet) then begin
    kpower_1d_plots, file_arr, window_num = 5, names = names_arr, colors = colors_arr
  endif
end
