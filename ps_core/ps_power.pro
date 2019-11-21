pro ps_power, file_struct, sim = sim, fix_sim_input = fix_sim_input, $
    uvf_input = uvf_input, freq_ch_range = freq_ch_range, freq_flags = freq_flags, $
    input_units = input_units, save_slices = save_slices, $
    refresh_options = refresh_options, uvf_options = uvf_options, $
    ps_options = ps_options

  nfiles = n_elements(file_struct.datafile)

  if n_elements(freq_flags) ne 0 then freq_mask = file_struct.freq_mask

  test_kcube = file_valid(file_struct.kcube_savefile)
  if not test_kcube then test_kcube = check_old_path(file_struct, 'kcube_savefile')

  if test_kcube eq 1 and n_elements(freq_flags) ne 0 then begin
    old_freq_mask = getvar_savefile(file_struct.kcube_savefile, 'freq_mask')
    if total(abs(old_freq_mask - freq_mask)) ne 0 then test_kcube = 0
  endif

  if test_kcube eq 0 or refresh_options.refresh_kcube then begin
    ps_kcube, file_struct, sim = sim, fix_sim_input = fix_sim_input, $
      uvf_input = uvf_input, freq_ch_range = freq_ch_range, freq_flags = freq_flags, $
      input_units = input_units, save_slices = save_slices, $
      refresh_options = refresh_options, uvf_options = uvf_options, $
      ps_options = ps_options
  endif

  if nfiles eq 1 then begin
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
    undefine, sigma2_2

    if keyword_set(no_wtd_avg) then begin
      term1 = abs(data_sum_1)^2.
      term2 = abs(data_sum_2)^2.
      undefine, data_sum_1, data_sum_2

      sim_noise_t1 = abs(sim_noise_sum_1)^2.
      sim_noise_t1 = abs(sim_noise_sum_2)^2.
      undefine, data_sum_1, data_sum_2

      weights_3d = (power_weights1 + power_weights2)
      undefine, power_weights1, power_weights2

      ;; expected noise is just sqrt(variance)
      noise_expval_3d = 1./sqrt(weights_3d)

      ;; Add the 2 terms
      power_3d = (term1 + term2)
      sim_noise_3d = (sim_noise_t1 + sim_noise_t2)
      undefine, term1, term2, sim_noise_t1, sim_noise_t2

      wh_wt0 = where(weights_3d eq 0, count_wt0)
      if count_wt0 ne 0 then begin
        power_3d[wh_wt0] = 0
        noise_expval_3d[wh_wt0] = 0
        sim_noise_3d[wh_wt0] = 0
      endif

    endif else begin
      term1 = abs(data_sum_1)^2.*power_weights1
      term2 = abs(data_sum_2)^2.*power_weights2
      undefine, data_sum_1, data_sum_2

      sim_noise_t1 = abs(sim_noise_sum_1)^2. * power_weights1
      sim_noise_t2 = abs(sim_noise_sum_2)^2. * power_weights2
      undefine, data_sum_1, data_sum_2

      ;; Factor of 2 because we're adding the cosine & sine terms
      noise_expval_3d = (sqrt(power_weights1 + power_weights2))*2.
      ;; except for kparallel=0 b/c there's only one term
      noise_expval_3d[*,*,0] = noise_expval_3d[*,*,0]/2.

      weights_3d = (power_weights1 + power_weights2)
      undefine, power_weights1, power_weights2

      ;; multiply by 2 because power is generally the SUM of the cosine & sine powers
      power_3d = (term1 + term2)*2.
      sim_noise_3d = (sim_noise_t1 + sim_noise_t2)*2.
      ;; except for kparallel=0 b/c there's only one term
      power_3d[*,*,0] = power_3d[*,*,0]/2.
      sim_noise_3d[*,*,0] = sim_noise_3d[*,*,0]/2.

      power_3d = power_3d / weights_3d
      sim_noise_3d = sim_noise_3d / weights_3d
      noise_expval_3d = noise_expval_3d / weights_3d
      undefine, term1, term2

      wh_wt0 = where(weights_3d eq 0, count_wt0)
      if count_wt0 ne 0 then begin
        power_3d[wh_wt0] = 0
        sim_noise_3d[wh_wt0] = 0
        noise_expval_3d[wh_wt0] = 0
      endif

      ;; variance_3d = 4/weights_3d b/c of factors of 2 in power
      ;; in later code variance is taken to be 1/weights so divide by 4 now
      weights_3d = weights_3d/4.
      ;; except for kparallel=0 b/c there's only one term
      weights_3d[*,*,0] = weights_3d[*,*,0]*4.
    endelse


  endif else begin
    ;; nfiles=2
    restore, file_struct.kcube_savefile

    n_kx = n_elements(kx_mpc)
    n_ky = n_elements(ky_mpc)
    n_kz = n_elements(kz_mpc)

    ;; now construct weights for power & diff power cubes
    ;; variance of power is a factor of 2 higher than variance of diff power
    ;; because of sum & difference contributions
    ;; weights = 1/variance
    power_weights1 = 1d/(8*(sigma2_1)^2d)
    wh_sig1_0 = where(sigma2_1^2d eq 0, count_sig1_0)
    if count_sig1_0 ne 0 then power_weights1[wh_sig1_0] = 0

    power_weights2 = 1d/(8*(sigma2_2)^2d)
    wh_sig2_0 = where(sigma2_2^2d eq 0, count_sig2_0)
    if count_sig2_0 ne 0 then power_weights2[wh_sig2_0] = 0

    diff_power_weights1 = 1d/(4*(sigma2_1)^2d)
    if count_sig1_0 ne 0 then diff_power_weights1[wh_sig1_0] = 0
    undefine, sigma2_1

    diff_power_weights2 = 1d/(4*(sigma2_2)^2d)
    if count_sig2_0 ne 0 then diff_power_weights2[wh_sig2_0] = 0
    undefine, sigma2_1

    if keyword_set(no_wtd_avg) then begin
      term1 = (abs(data_sum_1)^2. - abs(data_diff_1)^2.)
      term2 = (abs(data_sum_2)^2. - abs(data_diff_2)^2.)
      undefine, data_sum_1, data_sum_2

      sim_noise_t1 = (abs(sim_noise_sum_1)^2. - abs(sim_noise_diff_1)^2.)
      sim_noise_t2 = (abs(sim_noise_sum_2)^2. - abs(sim_noise_diff_2)^2.)
      undefine, data_sum_1, data_sum_2

      noise_t1 = abs(data_diff_1)^2.
      noise_t2 = abs(data_diff_2)^2.
      undefine, data_diff_1, data_diff_2

      sim_noise_diff_t1 = abs(sim_noise_diff_1)^2.
      sim_noise_diff_t2 = abs(sim_noise_diff_2)^2.
      undefine, sim_noise_diff_1, sim_noise_diff_2

      weights_3d = power_weights1 + power_weights2
      undefine, power_weights1, power_weights2

      ;; expected noise is just sqrt(variance)
      noise_expval_3d = 1./sqrt(diff_power_weights1 + diff_power_weights2)
      diff_weights_3d = diff_power_weights1 + diff_power_weights2
      undefine, diff_power_weights1, diff_power_weights2

      ;; Add the 2 terms
      power_3d = (term1 + term2)
      noise_3d = (noise_t1 + noise_t2)
      sim_noise_3d = (sim_noise_t1 + sim_noise_t2)
      sim_noise_diff_3d = (sim_noise_diff_t1 + sim_noise_diff_t2)
      undefine, term1, term2, noise_t1, noise_t2, sim_noise_t1, sim_noise_t2, $
        sim_noise_diff_t1, sim_noise_diff_t2

      wh_wt0 = where(weights_3d eq 0, count_wt0)
      if count_wt0 ne 0 then begin
        power_3d[wh_wt0] = 0
        noise_expval_3d[wh_wt0] = 0
        noise_3d[wh_wt0] = 0
        sim_noise_3d[wh_wt0] = 0
        sim_noise_diff_3d[wh_wt0] = 0
      endif

    endif else begin
      term1 = (abs(data_sum_1)^2. - abs(data_diff_1)^2.) * power_weights1
      term2 = (abs(data_sum_2)^2. - abs(data_diff_2)^2.) * power_weights2
      undefine, data_sum_1, data_sum_2

      sim_noise_t1 = (abs(sim_noise_sum_1)^2. - abs(sim_noise_diff_1)^2.) * power_weights1
      sim_noise_t2 = (abs(sim_noise_sum_2)^2. - abs(sim_noise_diff_2)^2.) * power_weights2
      undefine, data_sum_1, data_sum_2

      noise_t1 = abs(data_diff_1)^2. * diff_power_weights1
      noise_t2 = abs(data_diff_2)^2. * diff_power_weights2
      undefine, data_diff_1, data_diff_2

      sim_noise_diff_t1 = abs(sim_noise_diff_1)^2. * diff_power_weights1
      sim_noise_diff_t2 = abs(sim_noise_diff_2)^2. * diff_power_weights2
      undefine, sim_noise_diff_1, sim_noise_diff_2

      ;; Factor of 2 because we're adding the cosine & sine terms,
      ;; sqrt(2) because sum of exponentials gives erlang with mean = 2*sigma_exponential
      noise_expval_3d = sqrt(diff_power_weights1 + diff_power_weights2)*2*sqrt(2)
      ;; except for kparallel=0 b/c there's only one term
      noise_expval_3d[*,*,0] = noise_expval_3d[*,*,0]/(2*sqrt(2))


      weights_3d = power_weights1 + power_weights2
      diff_weights_3d = diff_power_weights1 + diff_power_weights2
      undefine, power_weights1, power_weights2, diff_power_weights1, diff_power_weights2

      ;; multiply by 2 because power is generally the SUM of the cosine & sine powers
      power_3d = (term1 + term2)*2.
      noise_3d = (noise_t1 + noise_t2)*2
      sim_noise_3d = (sim_noise_t1 + sim_noise_t2)*2.
      sim_noise_diff_3d = (sim_noise_diff_t1 + sim_noise_diff_t2)*2.
      ;; except for kparallel=0 b/c there's only one term
      power_3d[*,*,0] = power_3d[*,*,0]/2.
      noise_3d[*,*,0] = noise_3d[*,*,0]/2.
      sim_noise_3d[*,*,0] = sim_noise_3d[*,*,0]/2.
      sim_noise_diff_3d[*,*,0] = sim_noise_diff_3d[*,*,0]/2.

      power_3d = power_3d / weights_3d
      noise_3d = noise_3d / diff_weights_3d
      sim_noise_3d = sim_noise_3d / weights_3d
      sim_noise_diff_3d = sim_noise_diff_3d / diff_weights_3d
      noise_expval_3d = noise_expval_3d / diff_weights_3d
      undefine, term1, term2, noise_t1, noise_t2, sim_noise_t1, sim_noise_t2, $
        sim_noise_diff_t1, sim_noise_diff_t2

      wh_wt0 = where(weights_3d eq 0, count_wt0)
      if count_wt0 ne 0 then begin
        power_3d[wh_wt0] = 0
        noise_expval_3d[wh_wt0] = 0
        noise_3d[wh_wt0] = 0
        sim_noise_3d[wh_wt0] = 0
        sim_noise_diff_3d[wh_wt0] = 0
      endif

      ;; variance_3d = 4/weights_3d b/c of factors of 2 in power
      ;; in later code variance is taken to be 1/weights so divide by 4 now
      weights_3d = weights_3d/4.
      diff_weights_3d = diff_weights_3d/4.
      ;; except for kparallel=0 b/c there's only one term
      weights_3d[*,*,0] = weights_3d[*,*,0]*4.
      diff_weights_3d[*,*,0] = diff_weights_3d[*,*,0]*4.
    endelse

  endelse

  git, repo_path = ps_repository_dir(), result=ps_git_hash
  if n_elements(git_hashes) gt 0 then begin
    git_hashes = create_struct(git_hashes, 'ps', ps_git_hash)
  endif else begin
    git_hashes = {uvf:strarr(nfiles), uvf_wt:strarr(nfiles), beam:strarr(nfiles), $
      kcube:'', ps:ps_git_hash}
  endelse

  if n_elements(freq_flags) ne 0 then begin
    save, file = file_struct.power_savefile, power_3d, noise_3d, sim_noise_3d, $
      sim_noise_diff_3d, noise_expval_3d, weights_3d, diff_weights_3d, $
      kx_mpc, ky_mpc, kz_mpc, kperp_lambda_conv, delay_params, hubble_param, $
      n_freq_contrib, freq_mask, vs_name, vs_mean, t_sys_meas, window_int, $
      git_hashes, wt_meas_ave, wt_meas_min, ave_weights, wt_ave_power_freq, $
      ave_power_freq, wt_ave_power_uvf, ave_power_uvf
  endif else begin
    save, file = file_struct.power_savefile, power_3d, noise_3d, sim_noise_3d, $
      sim_noise_diff_3d, noise_expval_3d, weights_3d, diff_weights_3d, $
      kx_mpc, ky_mpc, kz_mpc, kperp_lambda_conv, delay_params, hubble_param, $
      n_freq_contrib, vs_name, vs_mean, t_sys_meas, window_int, git_hashes, $
      wt_meas_ave, wt_meas_min, ave_weights, wt_ave_power_freq, ave_power_freq, $
      wt_ave_power_uvf, ave_power_uvf
  endelse

  write_ps_fits, file_struct.fits_power_savefile, power_3d, weights_3d, $
    noise_expval_3d, noise_3d = noise_3d, kx_mpc, ky_mpc, kz_mpc, $
    kperp_lambda_conv, delay_params, hubble_param

  print, 'power integral:', total(power_3d)

  wt_ave_power = total(weights_3d * power_3d)/total(weights_3d)
  ave_power = mean(power_3d[where(weights_3d ne 0)])

  n_kx = n_elements(kx_mpc)
  n_ky = n_elements(ky_mpc)
  n_kz = n_elements(kz_mpc)

  uv_pix_area = (kx_mpc[1]-kx_mpc[0])*(ky_mpc[1]-ky_mpc[0])*kperp_lambda_conv^2.
  uv_area = uv_pix_area*n_kx*n_ky

  if save_slices then begin
    ;; now do slices
    make_slices, file_struct, type='kspace', data_cube = power_3d, $
                weights_cube = weights_3d, $
                noise_3d = noise_3d, noise_expval_3d = noise_expval_3d, $
                kx_mpc = kx_mpc, ky_mpc = ky_mpc, kz_mpc = kz_mpc, $
                kperp_lambda_conv = kperp_lambda_conv, delay_params = delay_params, $
                hubble_param = hubble_param

  endif

end
