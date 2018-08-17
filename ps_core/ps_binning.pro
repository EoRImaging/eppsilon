pro ps_binning, file_struct, sim = sim, freq_flags = freq_flags, $
    ps_options = ps_options, binning_2d_options = binning_2d_options, $
    binning_1d_options = binning_1d_options, plot_options = plot_options, $
    plot_types = plot_types, $
    savefile_2d = savefile_2d, savefile_1d = savefile_1d, $
    savefile_1to2d_bin = savefile_1to2d_bin, savefile_masked_2d = savefile_masked_2d, $
    savefile_kpar_power = savefile_kpar_power, savefile_kperp_power = savefile_kperp_power, $
    savefile_k0 = savefile_k0, savefile_masked_k0 = savefile_masked_k0, $
    bin_arr_3d = bin_arr_3d

  nfiles = n_elements(file_struct.datafile)

  if n_elements(freq_flags) ne 0 then freq_mask = file_struct.freq_mask

  restore, file_struct.power_savefile

  wt_ave_power = total(weights_3d * power_3d)/total(weights_3d)
  ave_power = mean(power_3d[where(weights_3d ne 0)])

  n_kx = n_elements(kx_mpc)
  n_ky = n_elements(ky_mpc)
  n_kz = n_elements(kz_mpc)

  uv_pix_area = (kx_mpc[1]-kx_mpc[0])*(ky_mpc[1]-ky_mpc[0])*kperp_lambda_conv^2.
  uv_area = uv_pix_area*n_kx*n_ky

  if binning_2d_options.no_kzero then begin
    ;; leave out kz=0 -- full of foregrounds
    kz_mpc = kz_mpc[1:*]
    power_3d = temporary(power_3d[*, *, 1:*])
    weights_3d = temporary(weights_3d[*,*,1:*])
    diff_weights_3d = temporary(diff_weights_3d[*,*,1:*])
    noise_expval_3d = temporary(noise_expval_3d[*,*,1:*])
    if nfiles eq 2 then begin
      noise_3d = temporary(noise_3d[*,*,1:*])
      sim_noise_3d = temporary(sim_noise_3d[*,*,1:*])
      sim_noise_diff_3d = temporary(sim_noise_diff_3d[*,*,1:*])
    endif
    n_kz = n_elements(kz_mpc)
  endif

  power_tag = file_struct.power_tag

  fadd_2dbin = ''
  if binning_2d_options.no_kzero then fadd_2dbin = fadd_2dbin + '_nok0'
  if binning_2d_options.log_kpar then fadd_2dbin = fadd_2dbin + '_logkpar'
  if binning_2d_options.log_kperp then fadd_2dbin = fadd_2dbin + '_logkperp'

  git, repo_path = ps_repository_dir(), result=binning_git_hash
  if n_elements(git_hashes) gt 0 then begin
    git_hashes = create_struct(git_hashes, 'binning', binning_git_hash)
  endif else begin
    git_hashes = {uvf:strarr(nfiles), uvf_wt:strarr(nfiles), beam:strarr(nfiles), $
      kcube:'', ps:'', binning:binning_git_hash}
  endelse

  sigma_3d = sqrt(1./weights_3d)
  sigma_3d[where(weights_3d eq 0)] = 0
  temp_uniform = randomu(seed, n_kx, n_ky, n_kz) - 0.5
  temp_sign = (temp_uniform gt 0) - 1.*(temp_uniform lt 0)
  new_noise_3d = -1 * sigma_3d/sqrt(2.) * temp_sign * alog(1. - 2.*abs(temp_uniform))
  wh_0 = where(temp_uniform eq 0, count_0)
  if count_0 gt 0 then new_noise_3d[wh_0] = 0
  undefine, temp_uniform, temp_sign, sigma_3d


  print, 'Binning to 2D power spectrum'

  for j=0, n_elements(ps_options.wt_cutoffs)-1 do begin

    if ps_options.wt_cutoffs[j] gt 0 then begin
      case ps_options.wt_measures[j] of
        'ave': wt_meas_use = wt_meas_ave
        'min': wt_meas_use = wt_meas_min
      endcase

      wt_cutoff_use = ps_options.wt_cutoffs[j]
    endif else undefine, wt_cutoff_use, wt_meas_use

    make_2d_files, nfiles, savefile_2d[j], savefile_k0[j], power_3D, sim_noise_3D, $
      new_noise_3d, noise_expval_3d, weights_3d, kx_mpc, ky_mpc, kz_mpc, $
      delay_params, hubble_param, plot_options.hinv, kperp_lambda_conv, freq_mask, $
      vs_name, vs_mean, t_sys_meas, window_int, git_hashes, wt_ave_power, ave_power, $
      ave_weights, wt_ave_power_freq, ave_power_freq, wt_ave_power_uvf, ave_power_uvf, $
      uv_pix_area, uv_area, noise_3D = noise_3D, sim_noise_diff_3D = sim_noise_diff_3D, $
      wt_measure = wt_meas_use, wt_cutoff = wt_cutoff_use, freq_flags = freq_flags, $
      binning_2d_options = binning_2d_options

  endfor

  print, 'Binning to 1D power spectrum'


  n_wt_cuts = n_elements(ps_options.wt_cutoffs)

  if n_elements(savefile_1d) ne $
    (n_elements(binning_1d_options.wedge_amps)+1)*(n_wt_cuts) then begin
    message, 'number of elements in savefile_1d is wrong'
  endif

  if tag_exist(binning_1d_options, 'coarse_harm0') then begin
    coarse_harm0 = binning_1d_options.coarse_harm0
    coarse_width = binning_1d_options.coarse_harm_width
  endif

  for i=0, n_elements(binning_1d_options.wedge_amps) do begin
    for j=0, n_wt_cuts-1 do begin
      if i gt 0 then begin
        wedge_amp_use = binning_1d_options.wedge_amps[i-1]
      endif

      if ps_options.wt_cutoffs[j] gt 0 then begin
        case ps_options.wt_measures[j] of
          'ave': wt_meas_use = wt_meas_ave
          'min': wt_meas_use = wt_meas_min
        endcase

        wt_cutoff_use = ps_options.wt_cutoffs[j]

      endif else undefine, wt_cutoff_use, wt_meas_use

      power_1d = kspace_rebinning_1d(power_3d, kx_mpc, ky_mpc, kz_mpc, k_edges_mpc, $
        noise_expval = noise_expval_3d, weights = weights_3d, $
        binned_noise_expval = noise_expval_1d, binned_weights = weights_1d, $
        bin_arr_3d = bin_arr_3d, noise_frac_3d = noise_frac_3d, $
        wedge_amp = wedge_amp_use, kperp_density_measure = wt_meas_use, $
        kperp_density_cutoff = wt_cutoff_use, $
        binning_1d_options = binning_1d_options, plot_options = plot_options, $
        hubble_param = hubble_param, kperp_lambda_conv)

      sim_noise_1d = kspace_rebinning_1d(sim_noise_3d, kx_mpc, ky_mpc, kz_mpc, $
        k_edges_mpc, noise_expval = noise_expval_3d, weights = weights_3d, $
        binned_noise_expval = noise_expval_1d, binned_weights = weights_1d, $
        bin_arr_3d = bin_arr_3d, var_power_1d = var_power_1d, mean_var_1d = mean_var_1d, $
        noise_frac_3d = noise_frac_3d, wedge_amp = wedge_amp_use, $
        kperp_density_measure = wt_meas_use, kperp_density_cutoff = wt_cutoff_use, $
        binning_1d_options = binning_1d_options, plot_options = plot_options, $
        hubble_param = hubble_param, kperp_lambda_conv)

      new_noise_1d = kspace_rebinning_1d(new_noise_3d, kx_mpc, ky_mpc, kz_mpc, $
        k_edges_mpc, noise_expval = noise_expval_3d, weights = weights_3d, $
        binned_noise_expval = noise_expval_1d, binned_weights = weights_1d, $
        bin_arr_3d = bin_arr_3d, var_power_1d = new_var_power_1d, $
        mean_var_1d = new_mean_var_1d, noise_frac_3d = noise_frac_3d, $
        wedge_amp = wedge_amp_use, $
        kperp_density_measure = wt_meas_use, kperp_density_cutoff = wt_cutoff_use, $
        binning_1d_options = binning_1d_options, plot_options = plot_options, $
        hubble_param = hubble_param, kperp_lambda_conv)

      if nfiles eq 2 then begin
        noise_1d = kspace_rebinning_1d(noise_3d, kx_mpc, ky_mpc, kz_mpc, k_edges_mpc, $
          noise_expval = noise_expval_3d, weights = weights_3d, $
          binned_noise_expval = noise_expval_1d, binned_weights = weights_1d, $
          bin_arr_3d = bin_arr_3d, wedge_amp = wedge_amp_use, $
          kperp_density_measure = wt_meas_use, kperp_density_cutoff = wt_cutoff_use, $
          binning_1d_options = binning_1d_options, plot_options = plot_options, $
          hubble_param = hubble_param, kperp_lambda_conv)

        sim_noise_diff_1d = kspace_rebinning_1d(sim_noise_diff_3d, kx_mpc, ky_mpc, $
          kz_mpc, k_edges_mpc, noise_expval = noise_expval_3d, weights = weights_3d, $
          binned_noise_expval = noise_expval_1d, binned_weights = weights_1d, $
          bin_arr_3d = bin_arr_3d, wedge_amp = wedge_amp_use, $
          kperp_density_measure = wt_meas_use, kperp_density_cutoff = wt_cutoff_use, $
          binning_1d_options = binning_1d_options, plot_options = plot_options, $
          hubble_param = hubble_param, kperp_lambda_conv)
      endif
      power = power_1d
      sim_noise = sim_noise_1d
      if nfiles eq 2 then begin
        noise = noise_1d
        sim_noise_diff = sim_noise_diff_1d
      endif
      weights = weights_1d
      k_edges = k_edges_mpc
      k_bin = binning_1d_options.k_bin
      noise_expval = noise_expval_1d
      kperp_range = binning_1d_options.kperp_range_1dave
      kx_range = binning_1d_options.kx_range_1dave
      ky_range = binning_1d_options.ky_range_1dave
      kpar_range = binning_1d_options.kpar_range_1dave
      if plot_options.hinv then begin
        kperp_range_lambda = binning_1d_options.kperp_range_1dave * hubble_param * kperp_lambda_conv
        kx_range_lambda = binning_1d_options.kx_range_1dave * hubble_param * kperp_lambda_conv
        ky_range_lambda = binning_1d_options.ky_range_1dave * hubble_param * kperp_lambda_conv
      endif else begin
        kperp_range_lambda = binning_1d_options.kperp_range_1dave * kperp_lambda_conv
        kx_range_lambda = binning_1d_options.kx_range_1dave * kperp_lambda_conv
        ky_range_lambda = binning_1d_options.ky_range_1dave * kperp_lambda_conv
      endelse
      if i gt 0 then wedge_amp = wedge_amp_use

      ;; This generates warnings about freq_mask, coarse_harm0 and
      ;; coarse_harm_width if they aren't defined. There's not a way to
      ;; eliminate those without repeating this code a bunch of times
      save, file = savefile_1d[j,i], power, noise, sim_noise, sim_noise_diff, $
        weights, noise_expval, k_edges, k_bin, kperp_lambda_conv, delay_params, $
        hubble_param, freq_mask, wedge_amp, kperp_range, kperp_range_lambda, $
        kx_range, kx_range_lambda, ky_range, ky_range_lambda, kpar_range, $
        coarse_harm0, coarse_harm_width, window_int, git_hashes, wt_ave_power, $
        ave_power, ave_weights, wt_ave_power_freq, ave_power_freq, wt_ave_power_uvf, $
        ave_power_uvf, uv_pix_area, uv_area

      textfile = strmid(savefile_1d[j,i], 0, stregex(savefile_1d[j,i], '.idlsave')) + '.txt'
      print, 'saving 1d power to ' + textfile
      save_1D_text, textfile, k_edges, power, weights, noise_expval, hubble_param, $
        noise, sim_noise_power = sim_noise, sim_noise_diff = sim_noise_diff, $
        nfiles = nfiles, hinv = plot_options.hinv

      mask_weights = long(bin_arr_3d gt 0)

      bin_1to2d_ave = kspace_rebinning_2d(bin_arr_3d, kx_mpc, ky_mpc, kz_mpc, $
        kperp_edges_mpc, kpar_edges_mpc, $
        binning_2d_options = binning_2d_options)

      bin_1to2d = kspace_rebinning_2d(bin_arr_3d, kx_mpc, ky_mpc, kz_mpc, kperp_edges_mpc, $
        kpar_edges_mpc, weights = mask_weights, $
        binning_2d_options = binning_2d_options)

      if plot_types.plot_binning_hist then begin
        if n_elements(plotfile_binning_hist) gt 0 then begin
          plotfilebase_use1 = plotfile_binning_hist[j,i] + '_3to1d'
          plotfilebase_use2 = plotfile_binning_hist[j,i] + '_2to1d'
        endif
        binning_hist_plots, power_3d, sim_noise_3d, weights_3d, bin_arr_3d, power_1d, $
          weights_1d, sim_noise_1d, window_start = 1, plotfilebase = plotfilebase_use1, $
          png = plot_options.png, eps = plot_options.eps, pdf = plot_options.pdf
        binning_hist_plots, power_rebin, sim_noise_rebin, binned_weights, bin_1to2d, $
          power_1d, weights_1d, sim_noise_1d, window_start = !d.window+1, $
          plotfilebase = plotfilebase_use2, png = plot_options.png, $
          eps = plot_options.eps, pdf = plot_options.pdf
      endif

      noise_frac_1to2d = kspace_rebinning_2d(noise_frac_3d, kx_mpc, ky_mpc, kz_mpc, $
        kperp_edges_mpc, kpar_edges_mpc, weights = mask_weights, $
        binning_2d_options = binning_2d_options)

      kperp_edges = kperp_edges_mpc
      kpar_edges = kpar_edges_mpc
      kpar_bin = binning_2d_options.kpar_bin
      kperp_bin = binning_2d_options.kperp_bin

      ;; This generates warnings about freq_mask, coarse_harm0 and
      ;; coarse_harm_width if they aren't defined. There's not a way to
      ;; eliminate those without repeating this code a bunch of times
      save, file = savefile_1to2d_bin[j,i], bin_arr_3d, bin_1to2d, bin_1to2d_ave, $
        noise_frac_1to2d, kx_mpc, ky_mpc, kz_mpc, kperp_edges, kpar_edges, $
        kpar_bin, kperp_bin, k_edges, k_bin, kperp_range, kperp_range_lambda, $
        kx_range, kx_range_lambda, ky_range, ky_range_lambda, kpar_range, $
        coarse_harm0, coarse_harm_width, kperp_lambda_conv, delay_params, hubble_param, $
        freq_mask, window_int, git_hashes, wt_ave_power, ave_power, ave_weights, $
        wt_ave_power_freq, ave_power_freq, wt_ave_power_uvf, ave_power_uvf, $
        uv_pix_area, uv_area

      masked_save_items = {kperp_range: kperp_range, kperp_range_lambda: kperp_range_lambda, $
        kpar_range:kpar_range, kx_range: kx_range, kx_range_lambda: kx_range_lambda, $
        ky_range: ky_range, ky_range_lambda: ky_range_lambda}

      if tag_exist(binning_1d_options, 'coarse_harm0') then begin
        masked_save_items = create_struct(masked_save_items, 'coarse_harm0', $
          binning_1d_options.coarse_harm0, 'coarse_width', $
          binning_1d_options.coarse_harm_width)
      endif

      ;; make new 2d files with mask_weights applied so that the 2d ps contain
      ;; the same voxels as the 1d ps
      make_2d_files, nfiles, savefile_masked_2d[j,i], savefile_masked_k0[j,i], $
        power_3D, sim_noise_3D, new_noise_3d, noise_expval_3d, mask_weights*weights_3d, $
        kx_mpc, ky_mpc, kz_mpc, delay_params, hubble_param, plot_options.hinv, kperp_lambda_conv, $
        freq_mask, vs_name, vs_mean, t_sys_meas, window_int, git_hashes, $
        wt_ave_power, ave_power, ave_weights, wt_ave_power_freq, ave_power_freq, $
        wt_ave_power_uvf, ave_power_uvf, uv_pix_area, uv_area, $
        masked_save_items = masked_save_items, noise_3D = noise_3D, $
        sim_noise_diff_3D = sim_noise_diff_3D, wt_measure = wt_meas_use, $
        wt_cutoff = wt_cutoff_use, freq_flags = freq_flags, $
        binning_2d_options = binning_2d_options

      ;; must undefine bin_arr_3d so that a new binning is calculated on next loops.
      undefine, bin_arr_3d

    endfor
  endfor

  ;; bin just in kpar for diagnostic plot
  for j=0, n_wt_cuts-1 do begin
    if ps_options.wt_cutoffs[j] gt 0 then begin
      case ps_options.wt_measures[j] of
        'ave': wt_meas_use = wt_meas_ave
        'min': wt_meas_use = wt_meas_min
      endcase

      wt_cutoff_use = ps_options.wt_cutoffs[j]
    endif else undefine, wt_cutoff_use, wt_meas_use

    kpar_binning_1d_options = create_binning_1d_options(binning_1d_options = binning_1d_options, $
      k_bin = binning_2d_options.kpar_bin, log_k = binning_2d_options.log_kpar, $
      /return_new)

    power_kpar = kspace_rebinning_1d(power_3d, kx_mpc, ky_mpc, kz_mpc, kpar_edges_mpc, $
      noise_expval = noise_expval_3d, weights = weights_3d, $
      binned_noise_expval = noise_expval_kpar, binned_weights = weights_1d, $
      bin_arr_3d = bin_arr_3d, /kpar_power, $
      kperp_density_measure = wt_meas_use, kperp_density_cutoff = wt_cutoff_use, $
      binning_1d_options = kpar_binning_1d_options, plot_options = plot_options, $
      hubble_param = hubble_param, kperp_lambda_conv)

    sim_noise_kpar = kspace_rebinning_1d(sim_noise_3d, kx_mpc, ky_mpc, kz_mpc, $
      kpar_edges_mpc, noise_expval = noise_expval_3d, weights = weights_3d, $
      binned_noise_expval = noise_expval_kpar, binned_weights = weights_1d, $
      var_power_1d = var_power_1d, mean_var_1d = mean_var_1d, $
      bin_arr_3d = bin_arr_3d, /kpar_power, $
      kperp_density_measure = wt_meas_use, kperp_density_cutoff = wt_cutoff_use, $
      binning_1d_options = kpar_binning_1d_options, plot_options = plot_options, $
      hubble_param = hubble_param, kperp_lambda_conv)

    if nfiles eq 2 then begin
      noise_kpar = kspace_rebinning_1d(noise_3d, kx_mpc, ky_mpc, kz_mpc, kpar_edges_mpc, $
        noise_expval = noise_expval_3d, weights = weights_3d, $
        binned_noise_expval = noise_expval_kpar, binned_weights = weights_1d, $
        bin_arr_3d = bin_arr_3d, /kpar_power, $
        kperp_density_measure = wt_meas_use, kperp_density_cutoff = wt_cutoff_use, $
        binning_1d_options = kpar_binning_1d_options, plot_options = plot_options, $
        hubble_param = hubble_param, kperp_lambda_conv)

      sim_noise_diff_kpar = kspace_rebinning_1d(sim_noise_diff_3d, kx_mpc, ky_mpc, $
        kz_mpc, kpar_edges_mpc, noise_expval = noise_expval_3d, weights = weights_3d, $
        binned_noise_expval = noise_expval_kpar, binned_weights = weights_1d, $
        bin_arr_3d = bin_arr_3d, /kpar_power, $
        kperp_density_measure = wt_meas_use, kperp_density_cutoff = wt_cutoff_use, $
        binning_1d_options = kpar_binning_1d_options, plot_options = plot_options, $
        hubble_param = hubble_param, kperp_lambda_conv)
    endif

    power = power_kpar
    sim_noise = sim_noise_kpar
    if nfiles eq 2 then begin
      noise = noise_kpar
      sim_noise_diff = sim_noise_diff_kpar
    endif
    weights = weights_1d
    k_edges = kpar_edges_mpc
    k_bin = kpar_binning_1d_options.k_bin
    noise_expval = noise_expval_kpar
    kperp_range = kpar_binning_1d_options.kperp_range_1dave
    kx_range = kpar_binning_1d_options.kx_range_1dave
    ky_range = kpar_binning_1d_options.ky_range_1dave
    kpar_range = kpar_binning_1d_options.kpar_range_1dave
    if plot_options.hinv then begin
      kperp_range_lambda = kpar_binning_1d_options.kperp_range_1dave * hubble_param * kperp_lambda_conv
      kx_range_lambda = kpar_binning_1d_options.kx_range_1dave * hubble_param * kperp_lambda_conv
      ky_range_lambda = kpar_binning_1d_options.ky_range_1dave * hubble_param * kperp_lambda_conv
    endif else begin
      kperp_range_lambda = kpar_binning_1d_options.kperp_range_1dave * kperp_lambda_conv
      kx_range_lambda = kpar_binning_1d_options.kx_range_1dave * kperp_lambda_conv
      ky_range_lambda = kpar_binning_1d_options.ky_range_1dave * kperp_lambda_conv
    endelse
    if n_elements(freq_flags) ne 0 then begin
      save, file = savefile_kpar_power[j], power, noise, sim_noise, sim_noise_diff, $
        weights, noise_expval, k_edges, k_bin, kperp_lambda_conv, delay_params, $
        hubble_param, freq_mask, kperp_range, kperp_range_lambda, kx_range, $
        kx_range_lambda, ky_range, ky_range_lambda, kpar_range, coarse_harm0, $
        coarse_harm_width, window_int, git_hashes, wt_ave_power, ave_power, ave_weights, $
        wt_ave_power_freq, ave_power_freq, wt_ave_power_uvf, ave_power_uvf, $
        uv_pix_area, uv_area
    endif else begin
      save, file = savefile_kpar_power[j], power, noise, sim_noise, sim_noise_diff, $
        weights, noise_expval, k_edges, k_bin, kperp_lambda_conv, delay_params, $
        hubble_param, kperp_range, kperp_range_lambda, kx_range, kx_range_lambda, $
        ky_range, ky_range_lambda, kpar_range, coarse_harm0, coarse_harm_width, $
        window_int, git_hashes, wt_ave_power, ave_power, ave_weights, $
        wt_ave_power_freq, ave_power_freq, wt_ave_power_uvf, ave_power_uvf, $
        uv_pix_area, uv_area
    endelse

    ;; must undefine bin_arr_3d so that a new binning is calculated on next loops.
    undefine, bin_arr_3d
  endfor

  ;; bin just in kperp for diagnostic plot
  for j=0, n_wt_cuts-1 do begin
    if ps_options.wt_cutoffs[j] gt 0 then begin
      case ps_options.wt_measures[j] of
        'ave': wt_meas_use = wt_meas_ave
        'min': wt_meas_use = wt_meas_min
      endcase

      wt_cutoff_use = ps_options.wt_cutoffs[j]
    endif else undefine, wt_cutoff_use, wt_meas_use

    kperp_binning_1d_options = create_binning_1d_options(binning_1d_options = binning_1d_options, $
      k_bin = binning_2d_options.kperp_bin, log_k = binning_2d_options.log_kperp, $
      /return_new)

    power_kperp = kspace_rebinning_1d(power_3d, kx_mpc, ky_mpc, kz_mpc, $
      kperp_edges_mpc, noise_expval = noise_expval_3d, weights = weights_3d, $
      binned_noise_expval = noise_expval_kperp, binned_weights = weights_1d, $
      bin_arr_3d = bin_arr_3d, /kperp_power, $
      kperp_density_measure = wt_meas_use, kperp_density_cutoff = wt_cutoff_use, $
      binning_1d_options = kperp_binning_1d_options, plot_options = plot_options, $
      hubble_param = hubble_param, kperp_lambda_conv)

    sim_noise_kperp = kspace_rebinning_1d(sim_noise_3d, kx_mpc, ky_mpc, kz_mpc, $
      kperp_edges_mpc, noise_expval = noise_expval_3d, weights = weights_3d, $
      binned_noise_expval = noise_expval_kperp, binned_weights = weights_1d, $
      var_power_1d = var_power_1d, mean_var_1d = mean_var_1d, $
      bin_arr_3d = bin_arr_3d, /kperp_power, $
      kperp_density_measure = wt_meas_use, kperp_density_cutoff = wt_cutoff_use, $
      binning_1d_options = kperp_binning_1d_options, plot_options = plot_options, $
      hubble_param = hubble_param, kperp_lambda_conv)

    if nfiles eq 2 then begin
      noise_kperp = kspace_rebinning_1d(noise_3d, kx_mpc, ky_mpc, kz_mpc, $
        kperp_edges_mpc, noise_expval = noise_expval_3d, weights = weights_3d, $
        binned_noise_expval = noise_expval_kperp, binned_weights = weights_1d, $
        bin_arr_3d = bin_arr_3d, /kperp_power, $
        kperp_density_measure = wt_meas_use, kperp_density_cutoff = wt_cutoff_use, $
        binning_1d_options = kperp_binning_1d_options, plot_options = plot_options, $
        hubble_param = hubble_param, kperp_lambda_conv)

      sim_noise_diff_kperp = kspace_rebinning_1d(sim_noise_diff_3d, kx_mpc, ky_mpc, $
        kz_mpc, kperp_edges_mpc, noise_expval = noise_expval_3d, weights = weights_3d, $
        binned_noise_expval = noise_expval_kperp, binned_weights = weights_1d, $
        bin_arr_3d = bin_arr_3d, /kperp_power, $
        kperp_density_measure = wt_meas_use, kperp_density_cutoff = wt_cutoff_use, $
        binning_1d_options = kperp_binning_1d_options, plot_options = plot_options, $
        hubble_param = hubble_param, kperp_lambda_conv)
    endif

    power = power_kperp
    sim_noise = sim_noise_kperp
    if nfiles eq 2 then begin
      noise = noise_kperp
      sim_noise_diff = sim_noise_diff_kperp
    endif
    weights = weights_1d
    k_edges = kperp_edges_mpc
    k_bin = kperp_binning_1d_options.k_bin
    noise_expval = noise_expval_kperp
    kperp_range = kperp_binning_1d_options.kperp_range_1dave
    kx_range = kperp_binning_1d_options.kx_range_1dave
    ky_range = kperp_binning_1d_options.ky_range_1dave
    kpar_range = kperp_binning_1d_options.kpar_range_1dave
    if plot_options.hinv then begin
      kperp_range_lambda = kperp_binning_1d_options.kperp_range_1dave * hubble_param * kperp_lambda_conv
      kx_range_lambda = kperp_binning_1d_options.kx_range_1dave * hubble_param * kperp_lambda_conv
      ky_range_lambda = kperp_binning_1d_options.ky_range_1dave * hubble_param * kperp_lambda_conv
    endif else begin
      kperp_range_lambda = kperp_binning_1d_options.kperp_range_1dave * kperp_lambda_conv
      kx_range_lambda = kperp_binning_1d_options.kx_range_1dave * kperp_lambda_conv
      ky_range_lambda = kperp_binning_1d_options.ky_range_1dave * kperp_lambda_conv
    endelse

    if n_elements(freq_flags) ne 0 then begin
      save, file = savefile_kperp_power[j], power, noise, sim_noise, sim_noise_diff, $
        weights, noise_expval, k_edges, k_bin, kperp_lambda_conv, delay_params, $
        hubble_param, freq_mask, kperp_range, kperp_range_lambda, kx_range, $
        kx_range_lambda, ky_range, ky_range_lambda, kpar_range, coarse_harm0, $
        coarse_harm_width, window_int, git_hashes, wt_ave_power, ave_power, ave_weights, $
        wt_ave_power_freq, ave_power_freq, wt_ave_power_uvf, ave_power_uvf, $
        uv_pix_area, uv_area
    endif else begin
      save, file = savefile_kperp_power[j], power, noise, sim_noise, sim_noise_diff, $
        weights, noise_expval, k_edges, k_bin, kperp_lambda_conv, delay_params, $
        hubble_param, kperp_range, kperp_range_lambda, kx_range, kx_range_lambda, $
        ky_range, ky_range_lambda, kpar_range, coarse_harm0, coarse_harm_width, window_int, $
        git_hashes, wt_ave_power, ave_power, ave_weights, wt_ave_power_freq, $
        ave_power_freq, wt_ave_power_uvf, ave_power_uvf, uv_pix_area, uv_area
    endelse

    ;; must undefine bin_arr_3d so that a new binning is calculated on next loops.
    undefine, bin_arr_3d
  endfor

end
