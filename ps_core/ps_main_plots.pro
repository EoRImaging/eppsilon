  ;; dft_fchunk applies only to Healpix datasets (it's ignored otherwise) and it specifies how many frequencies to process
  ;;   through the dft at once. This keyword allows for trade-offs between memory use and speed.
  ;; The optimum setting varies by computer and the speed is NOT a linear function of this parameter
  ;;   (it's not even monotonic) so some experimenting is required. The memory required is approximately linear --
  ;;   the higher this value the more memory required.
  ;; The default is 1 (to minimize memory use)
  ;; The maximum value of this parameter is the number of frequency slices in the cube
  ;;   (if its set too large it will be reduced to the maximum)
  
  ;; save_path specifies a location to save the power spectrum files.
  ;; This is also where the code looks for intermediate save files to avoid re-running code.
  ;; If this is parameter is not set, the files will be saved in the same directory as the datafile.
  
  ;; savefilebase specifies a base name to use for the save files
  
  ;; plot_path specifies a location to save plot files.
  ;; If this parameter is not set, the plots will be saved in the same directory as the datafile.
  
  ;; plot_filebase specifies a base name to use for the plot files
  
  
  ;; freq_ch_range specifies which frequency channels to include in the power spectrum.
  ;; Fewer number of channels makes the dfts faster
  
  ;; pol_inc specifies which polarizations to generate the power spectra for.
  ;; The default is ['xx,'yy']
  
  ;; cut_image keyword only applies to Healpix datasets. It allows for limiting the field of view in the
  ;; image plane by only using Healpix pixels inside a 30 degree diameter circle centered in the middle of the field.
  ;; Currently defaults to on. Set equal to 0 to turn it off, 1 to turn it on
  
  
  ;; There are 3 refresh flags to indicate that various stages should be re-calculated
  ;;   (rather than using previous save files if they exist).
  ;; If an early stage is recalculated, all subsequent stages will also be recalculated
  ;; The earliest stage is refresh_dft, which is only used for Healpix datasets (it's ignored otherwise)
  ;; The next stage is refresh_ps and the last stage is refresh_binning.
  ;; To set any of these flags, set them equal to 1 (true)
  
  ;; options for binning:
  ;; log_kperp, log_kpar and log_k1d are flags: set to 1 (true) for logarithmic bins
  ;; kperp_bin, kpar_bin and k1d_bin take scalar values to control bin sizes.
  ;;   (The actual binsize for linear binning and the log binsize for log binning -- bins per decade = 1/binsize)
  
  
  ;; options for plotting:
  ;; kperp_linear_axis is a flag, set to 1 to use a linear kperp axis (default is log axis)
  ;; kpar_linear_axis is a flag, set to 1 to use a linear kpar axis (default is log axis)
  ;; data_range specifies the min & max value of the signal colorbar (values outside that range are clipped to those values)
  ;; sigma_range, nev_range, snr_range, noise_range, nnr_range control the other colorbar ranges
  ;; baseline_axis is a flag (defaulted to true) to mark baseline; length along top axis of 2d plots (set to 0 to turn off)
  ;; delay_axis is a flag (defaulted to true) to mark delay time along right axis of 2d plots (set to 0 to turn off)
  ;; hinv is a flag (defaulted to true) to use h^-1 Mpc rather than physical Mpc in plot units (set to 0 to turn off)
  ;; plot_wedge_line is a flag (defaulted to true) to plot a line marking the wedge (both horizon & FoV) (set to 0 to turn off)
  ;; png & eps are flags to make save plots as png or eps files rather than displaying to the screen
  
  
  
  pro ps_main_plots, datafile, beamfiles = beamfiles, rts = rts, casa = casa, sim = sim, $
      refresh_dft = refresh_dft, refresh_ps = refresh_ps, refresh_binning = refresh_binning, refresh_info = refresh_info, $
      refresh_beam = refresh_beam, dft_fchunk = dft_fchunk, dft_ian = dft_ian, $
      delta_uv_lambda = delta_uv_lambda, max_uv_lambda = max_uv_lambda, cut_image = cut_image, $
      pol_inc = pol_inc, type_inc = type_inc, $
      freq_ch_range = freq_ch_range, freq_flags = freq_flags, freq_flag_name = freq_flag_name, $
      uvf_input = uvf_input, uv_avg = uv_avg, uv_img_clip = uv_img_clip, $
      std_power = std_power, inverse_covar_weight = inverse_covar_weight, no_wtd_avg = no_wtd_avg, norm_rts_with_fhd = norm_rts_with_fhd, $
      wt_cutoffs = wt_cutoffs, wt_measures = wt_measures, fix_sim_input = fix_sim_input, $
      no_spec_window = no_spec_window, spec_window_type = spec_window_type, $
      no_kzero = no_kzero, plot_slices = plot_slices, slice_type = slice_type, $
      uvf_plot_type = uvf_plot_type, plot_stdset = plot_stdset, plot_1to2d = plot_1to2d, $
      plot_kpar_power = plot_kpar_power, plot_kperp_power = plot_kperp_power, plot_k0_power = plot_k0_power, $
      plot_noise_1d = plot_noise_1d, plot_sim_noise = plot_sim_noise, $
      data_range = data_range, sigma_range = sigma_range, nev_range = nev_range, slice_range = slice_range, $
      snr_range = snr_range, noise_range = noise_range, nnr_range = nnr_range, range_1d = range_1d, $
      log_kpar = log_kpar, log_kperp = log_kperp, kpar_bin = kpar_bin, kperp_bin = kperp_bin, $
      log_k1d = log_k1d, k1d_bin = k1d_bin, plot_1d_delta = plot_1d_delta, $
      plot_1d_error_bars = plot_1d_error_bars, plot_1d_nsigma = plot_1d_nsigma, $
      kperp_range_1dave = kperp_range_1dave, kperp_range_lambda_1dave = kperp_range_lambda_1dave, kpar_range_1dave = kpar_range_1dave, $
      kperp_range_lambda_kparpower = kperp_range_lambda_kparpower, kpar_range_kperppower = kpar_range_kperppower, $
      kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, kperp_plot_range = kperp_plot_range, $
      kperp_lambda_plot_range = kperp_lambda_plot_range, kpar_plot_range = kpar_plot_range, $
      baseline_axis = baseline_axis, delay_axis = delay_axis, cable_length_axis = cable_length_axis, hinv = hinv, $
      plot_wedge_line = plot_wedge_line, wedge_angles = wedge_angles, coarse_harm_width = coarse_harm_width, $
      plot_eor_1d = plot_eor_1d, plot_flat_1d = plot_flat_1d, no_text_1d = no_text_1d, $
      save_path = save_path, savefilebase = savefilebase, plot_path = plot_path, plot_filebase = plot_filebase, $
      individual_plots = individual_plots, plot_binning_hist = plot_binning_hist, $
      note = note, png = png, eps = eps, pdf = pdf, cube_power_info = cube_power_info, no_dft_progress = no_dft_progress
      
      
    if keyword_set(refresh_dft) then refresh_beam = 1
    if keyword_set(refresh_beam) then refresh_ps = 1
    if keyword_set(refresh_ps) then refresh_binning = 1
    if keyword_set(plot_binning_hist) then refresh_binning = 1
    
    ;; default to making standard plot set if plot_slices isn't set
    if not keyword_set(plot_slices) then if n_elements(plot_stdset) eq 0 then plot_stdset = 1
    
    ;; default to including baseline axis & delay axis
    if n_elements(baseline_axis) eq 0 then baseline_axis = 1
    if n_elements(delay_axis) eq 0 and n_elements(cable_length_axis) eq 0 then delay_axis = 1
    
    ;; default to hinv
    if n_elements(hinv) eq 0 then hinv = 1
    
    ;; default to plot wedge line
    if n_elements(plot_wedge_line) eq 0 then plot_wedge_line=1
    
    ;; if inverse covariance weighted don't use spectral window
    if keyword_set(inverse_covar_weight) then no_spec_window = 1
    
    ;; default to blackman-harris spectral window
    if not keyword_set(no_spec_window) then begin
      if n_elements(spec_window_type) eq 0 then spec_window_type = 'Blackman-Harris'
    endif else undefine, spec_window_type
    
    ;; default to cutting in image space (down to 30 degree diameter circle)
    if n_elements(cut_image) eq 0 then cut_image = 1
    
    if keyword_set(no_kzero) and keyword_set(plot_k0_power) then message, 'plot_k0_power cannot be set if no_kzero keyword is set.'
    
    if keyword_set(png) or keyword_set(eps) or keyword_set(pdf) then pub = 1 else pub = 0
    if pub eq 1 then begin
    
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
    endif else plot_exten = ''
    
    
    if n_elements(savefilebase) gt 1 then message, 'savefilebase must be a scalar'
    
    time0 = systime(1)
    if keyword_set(rts) then begin
      file_struct_arr = rts_file_setup(datafile, savefilebase = savefilebase, save_path = save_path, $
        spec_window_type = spec_window_type, delta_uv_lambda = delta_uv_lambda, max_uv_lambda = max_uv_lambda, $
        use_fhd_norm = norm_rts_with_fhd, refresh_info = refresh_info)
    endif else if keyword_set(casa) then begin
      file_struct_arr = casa_file_setup(datafile, savefilebase = savefilebase, save_path = save_path, $
        spec_window_type = spec_window_type, delta_uv_lambda = delta_uv_lambda, max_uv_lambda = max_uv_lambda, refresh_info = refresh_info)
    endif else begin
      file_struct_arr = fhd_file_setup(datafile, beamfile = beamfiles, uvf_input = uvf_input, uv_avg = uv_avg, uv_img_clip = uv_img_clip, dft_ian = dft_ian, $
        savefilebase = savefilebase, save_path = save_path, freq_ch_range = freq_ch_range, freq_flags = freq_flags, freq_flag_name = freq_flag_name, $
        spec_window_type = spec_window_type, delta_uv_lambda = delta_uv_lambda, max_uv_lambda = max_uv_lambda, sim = sim, $
        std_power = std_power, inverse_covar_weight = inverse_covar_weight, no_wtd_avg = no_wtd_avg, refresh_info = refresh_info)
    endelse
    time1 = systime(1)
    print, 'file setup time: ' + number_formatter(time1-time0)
    
    if pub then begin
      if n_elements(plot_path) ne 0 then plotfile_path = plot_path $
      else if n_elements(save_path) ne 0 then plotfile_path = save_path else plotfile_path = file_struct_arr.savefile_froot
    endif
    
    if not tag_exist(file_struct_arr, 'beam_int') and keyword_set(refresh_ps) then refresh_beam = 1
    
    nfiles = file_struct_arr[0].nfiles
    
    if nfiles gt 1 then print, 'n_vis % difference between even & odd cubes: ' + $
      number_formatter((file_struct_arr[0].n_vis[1]-file_struct_arr[0].n_vis[0])*100/mean(file_struct_arr[0].n_vis))
      
    if tag_exist(file_struct_arr, 'n_obs') then begin
      print, 'n_obs: ', file_struct_arr[0].n_obs
      if n_elements(note) eq 0 then note = '(' + number_formatter(round(mean(file_struct_arr[0].n_obs))) + ')' $
      else note = note + ' (' + number_formatter(round(mean(file_struct_arr[0].n_obs))) + ')'
    endif
    
    if keyword_set(plot_wedge_line) then begin
      z0_freq = 1420.40 ;; MHz
      redshifts = z0_freq/file_struct_arr[0].frequencies - 1
      mean_redshift = mean(redshifts)
      
      cosmology_measures, mean_redshift, wedge_factor = wedge_factor, hubble_param = hubble_param
      ;; assume 20 degrees from pointing center to first null
      source_dist = 20d * !dpi / 180d
      fov_amp = wedge_factor * source_dist
      
      ;; calculate angular distance to horizon
      horizon_amp = wedge_factor * ((file_struct_arr[0].max_theta+90d) * !dpi / 180d)
      
      wedge_amp = [fov_amp, horizon_amp]
      wedge_names = ['fov', 'horizon']
      if n_elements(wedge_angles) gt 0 then begin
        if min(wedge_angles) le 0 or max(wedge_angles) ge 180 then message, 'wedge_angles must be in degrees and between 0 & 180'
        wedge_amp = [wedge_factor * (wedge_angles*!dpi / 180d)]
        wedge_names = [number_formatter(wedge_angles) + 'deg']
      endif
      
      if keyword_set(coarse_harm_width) then begin
        harm_freq = 1.28
        if n_elements(freq_ch_range) gt 0 then $
          freqs_use = file_struct_arr[0].frequencies[freq_ch_range[0]:freq_ch_range[1]] $
        else freqs_use = file_struct_arr[0].frequencies
        
        bandwidth = max(freqs_use) - min(freqs_use) + freqs_use[1] - freqs_use[0]
        coarse_harm0 = round(bandwidth / harm_freq)
      endif
      
    endif
    
    
    if n_elements(pol_inc) gt 0 or n_elements(type_inc) gt 0 then begin
      if n_elements(pol_inc) gt 0 then begin
        for i=0,n_elements(pol_inc)-1 do begin
          wh_this_pol = where(file_struct_arr.pol eq pol_inc[i], count_this_pol)
          if count_this_pol gt 0 then begin
            if n_elements(wh_pol) eq 0 then begin
              wh_pol = wh_this_pol
              pol_inds = intarr(count_this_pol) + i
            endif else begin
              wh_pol = [wh_pol, wh_this_pol]
              pol_inds = [pol_inds, intarr(count_this_pol) + i]
            endelse
          endif else message, 'desired pol ' + pol_inc[i] + ' is not available'
        endfor
        
        file_struct_arr = file_struct_arr[wh_pol]
        file_struct_arr.pol_index = pol_inds
      endif
      
      if n_elements(type_inc) gt 0 then begin
        for i=0,n_elements(type_inc)-1 do begin
          wh_this_type = where(file_struct_arr.type eq type_inc[i], count_this_type)
          if count_this_type gt 0 then begin
            if n_elements(wh_type) eq 0 then begin
              wh_type = wh_this_type
              type_inds = intarr(count_this_type) + i
            endif else begin
              wh_type = [wh_type, wh_this_type]
              type_inds = [type_inds, intarr(count_this_type) + i]
            endelse
          endif else message, 'desired type ' + type_inc[i] + ' is not available'
        endfor
        
        file_struct_arr = file_struct_arr[wh_type]
        file_struct_arr.type_index = type_inds
      endif
    endif
    
    npol = max(file_struct_arr.pol_index) + 1
    ntype = max(file_struct_arr.type_index) + 1
    n_cubes = n_elements(file_struct_arr)
    if n_cubes ne ntype * npol then message, 'number of cubes does not match expected value'
    
    file_labels = file_struct_arr.file_label
    wt_file_labels = file_struct_arr.wt_file_label
    titles = strarr(n_cubes)
    for i=0, n_cubes-1 do titles[i] = strjoin(strsplit(file_labels[i], '_', /extract), ' ')
    
    weight_ind = file_struct_arr.pol_index
    weight_labels = strupcase(file_struct_arr.pol)
    
    
    power_tag = file_struct_arr[0].power_tag
    if tag_exist(file_struct_arr[0], 'uvf_tag') then uvf_tag = file_struct_arr[0].uvf_tag else uvf_tag = ''
    
    
    fadd_2dbin = ''
    ;;if keyword_set(fill_holes) then fadd_2dbin = fadd_2dbin + '_nohole'
    if keyword_set(no_kzero) then fadd_2dbin = fadd_2dbin + '_nok0'
    if keyword_set(log_kpar) then fadd_2dbin = fadd_2dbin + '_logkpar'
    if keyword_set(log_kperp) then fadd_2dbin = fadd_2dbin + '_logkperp'
    
    fadd_1dbin = ''
    if keyword_set(log_k1d) then fadd_1dbin = fadd_1dbin + '_logk'
    if n_elements(kperp_range_1dave) gt 1 and n_elements(kperp_range_lambda_1dave) gt 1 then $
      message, 'both kperp_range_1dave and kperp_range_lambda_1dave cannot be set'
      
    fadd_kpar_1d = fadd_1dbin
    fadd_kperp_1d = fadd_1dbin
    
    if n_elements(kperp_range_1dave) gt 1 then begin
      fadd_1dbin = fadd_1dbin + '_kperp' + number_formatter(kperp_range_1dave[0]) + '-' + $
        number_formatter(kperp_range_1dave[1])
      note_1d = 'kperp: [' + number_formatter(kperp_range_1dave[0]) + ',' + $
        number_formatter(kperp_range_1dave[1]) + ']'
        
      ;; if we're plotting in [k]=h/Mpc then assume these numbers are in h/Mpc. Convert to 1/Mpc for internal code usage
      if keyword_set(hinv) then kperp_range_1d_use = kperp_range_1dave * hubble_param else kperp_range_1d_use = kperp_range_1dave
      
    endif else begin
      if n_elements(kperp_range_1dave) gt 0 then print, '2 values must be specified for kperp_range_1dave, defaulting to normal range'
      if n_elements(kperp_range_lambda_1dave) gt 1 then begin
        fadd_1dbin = fadd_1dbin + '_kperplambda' + number_formatter(kperp_range_lambda_1dave[0]) + '-' + $
          number_formatter(kperp_range_lambda_1dave[1])
        note_1d = 'kperp: [' + number_formatter(kperp_range_lambda_1dave[0]) + ',' + $
          number_formatter(kperp_range_lambda_1dave[1]) + ']'
      endif else begin
        if n_elements(kperp_range_lambda_1dave) gt 0 then print, '2 values must be specified for kperp_range_lambda_1dave, defaulting to normal range'
        ;; if no range set default to same range as is used in 2D plots
        kperp_range_lambda_1dave = [5., min([file_struct_arr.kspan/2.,file_struct_arr.max_baseline_lambda])]
      endelse
    endelse
    
    if n_elements(kperp_range_lambda_kparpower) gt 0 then begin
      fadd_kpar_1d = fadd_kpar_1d + '_kperplambda' + number_formatter(kperp_range_lambda_kparpower[0]) + '-' + $
        number_formatter(kperp_range_lambda_kparpower[1])
      note_kpar_1d = 'kperp: [' + number_formatter(kperp_range_lambda_kparpower[0]) + ',' + $
        number_formatter(kperp_range_lambda_kparpower[1]) + ']'
    endif else begin
      fadd_kpar_1d = fadd_1dbin
      note_kpar_1d = note_1d
      
      if n_elements(kperp_range_lambda_1dave) gt 1 then kperp_range_lambda_kparpower = kperp_range_lambda_1dave
    endelse
    
    if n_elements(kpar_range_1dave) gt 1 then begin
      fadd_1dbin = fadd_1dbin + '_kpar' + number_formatter(kpar_range_1dave[0]) + '-' + $
        number_formatter(kpar_range_1dave[1])
      if n_elements(note_1d) eq 0 then note_1d='' else note_1d = note_1d + '; '
      note_1d = note_1d + 'kpar: [' + number_formatter(kpar_range_1dave[0]) + ',' + $
        number_formatter(kpar_range_1dave[1]) + ']'
        
      ;; if we're plotting in [k]=h/Mpc then assume these numbers are in h/Mpc. Convert to 1/Mpc for internal code usage
      if keyword_set(hinv) then kpar_range_1d_use = kpar_range_1dave * hubble_param else kpar_range_1d_use = kpar_range_1dave
      
    endif else if n_elements(kpar_range_1dave) gt 0 then print, '2 values must be specified for kpar_range_1dave, defaulting to full range'
    
    if n_elements(kpar_range_kperppower) eq 0 and n_elements(kpar_range_1dave) gt 1 then $
      kpar_range_kperppower = kpar_range_1dave
      
    if n_elements(kpar_range_kperppower) gt 0 then begin
      fadd_kperp_1d = fadd_kperp_1d + '_kpar' + number_formatter(kpar_range_kperppower[0]) + '-' + $
        number_formatter(kpar_range_kperppower[1])
      note_kperp_1d = 'kpar: [' + number_formatter(kpar_range_kperppower[0]) + ',' + $
        number_formatter(kpar_range_kperppower[1]) + ']'
        
      ;; if we're plotting in [k]=h/Mpc then assume these numbers are in h/Mpc. Convert to 1/Mpc for internal code usage
      if keyword_set(hinv) then kpar_range_kperppower_use = kpar_range_kperppower * hubble_param else $
        kpar_range_kperppower_use = kpar_range_kperppower
    endif
    
    
    if keyword_set(plot_wedge_line) then begin
      if keyword_set(coarse_harm_width) then cb_width_name = '_cbw' + number_formatter(coarse_harm_width) else cb_width_name = ''
      wedge_1dbin_names = ['', '_no_' + wedge_names + '_wedge' + cb_width_name]
    endif else wedge_1dbin_names = ''
    
    ;; density correction defaults & file naming for 2D & 1D files
    if n_elements(wt_cutoffs) eq 0 then begin
      ;; default to wt_cutoffs = 1, wt_measures = 'min'
      wt_cutoffs = 1
      wt_measures = 'min'
    endif else if n_elements(wt_measures) eq 0 then begin
      print, 'wt_cutoffs is specified but wt_measures is not. Defaulting wt_measures to "min".'
      wt_measures = strarr(n_elements(wt_cutoffs)) + 'min'
    endif
    
    kperp_density_names = strarr(n_elements(wt_cutoffs))
    wh_cutoff0 = where(wt_cutoffs eq 0, count_cutoff0, complement = wh_cutoff_n0, ncomplement = count_cutoff_n0)
    wh_std = where(wt_cutoffs eq 1 and wt_measures eq 'min', count_std)
    
    if count_cutoff0 gt 0 then kperp_density_names[wh_cutoff0] = '_nodensitycorr'
    if count_cutoff_n0 gt 0 then kperp_density_names[wh_cutoff_n0] = '_kperp_density_' + wt_measures[wh_cutoff_n0] + '_gt' + number_formatter(wt_cutoffs[wh_cutoff_n0])
    
    if count_std gt 0 then kperp_density_names[wh_std] = '_dencorr'
    
    ;; need general_filebase for 1D plotfiles, make sure it doesn't have a full path
    general_filebase = file_struct_arr(0).general_filebase
    for i=0, n_cubes-1 do if file_struct_arr(i).general_filebase ne general_filebase then $
      message, 'general_filebase does not match between 1d savefiles'
      
      
    savefiles_2d = strarr(n_cubes, n_elements(kperp_density_names))
    for j=0, n_elements(kperp_density_names)-1 do begin
      savefiles_2d[*,j] = file_struct_arr.savefile_froot + file_struct_arr.savefilebase + power_tag + $
        fadd_2dbin + kperp_density_names[j] + '_2dkpower.idlsave'
    endfor
    test_save_2d = file_test(savefiles_2d) *  (1 - file_test(savefiles_2d, /zero_length))
    
    savefiles_1d = strarr(n_cubes, n_elements(kperp_density_names), n_elements(wedge_1dbin_names))
    savefiles_1to2d_mask = strarr(n_cubes, n_elements(kperp_density_names), n_elements(wedge_1dbin_names))
    for i=0, n_elements(wedge_1dbin_names)-1 do begin
      for j=0, n_elements(kperp_density_names)-1 do begin
        savefiles_1d[*,j,i] = file_struct_arr.savefile_froot + file_struct_arr.savefilebase + power_tag + $
          kperp_density_names[j] + wedge_1dbin_names[i] + fadd_1dbin + '_1dkpower.idlsave'
        savefiles_1to2d_mask[*,j,i] = file_struct_arr.savefile_froot + file_struct_arr.savefilebase + power_tag + fadd_2dbin + $
          kperp_density_names[j] + wedge_1dbin_names[i] + fadd_1dbin + '_1to2d_mask.idlsave'
      endfor
    endfor
    test_save_1d = file_test(savefiles_1d) *  (1 - file_test(savefiles_1d, /zero_length))
    test_save_mask = file_test(savefiles_1to2d_mask) *  (1 - file_test(savefiles_1to2d_mask, /zero_length))
    
    savefiles_kpar_1d = strarr(n_cubes, n_elements(kperp_density_names))
    for j=0, n_elements(kperp_density_names)-1 do begin
      savefiles_kpar_1d[*,j] = file_struct_arr.savefile_froot + file_struct_arr.savefilebase + power_tag + $
        kperp_density_names[j] + fadd_kpar_1d + '_kpar_power.idlsave'
    endfor
    if keyword_set(plot_kpar_power) then test_save_kpar = file_test(savefiles_kpar_1d) *  (1 - file_test(savefiles_kpar_1d, /zero_length))
    
    savefiles_kperp_1d = strarr(n_cubes, n_elements(kperp_density_names))
    for j=0, n_elements(kperp_density_names)-1 do begin
      savefiles_kperp_1d[*,j] = file_struct_arr.savefile_froot + file_struct_arr.savefilebase + power_tag + $
        kperp_density_names[j] + fadd_kperp_1d + '_kperp_power.idlsave'
    endfor
    if keyword_set(plot_kperp_power) then test_save_kperp = file_test(savefiles_kperp_1d) *  (1 - file_test(savefiles_kperp_1d, /zero_length))
    
    savefiles_k0 = strarr(n_cubes, n_elements(kperp_density_names))
    for j=0, n_elements(kperp_density_names)-1 do begin
      savefiles_k0[*,j] = file_struct_arr.savefile_froot + file_struct_arr.savefilebase + power_tag + fadd_2dbin + $
        kperp_density_names[j]+ '_k0power.idlsave'
    endfor
    if keyword_set(plot_k0_power) then test_save_k0 = file_test(savefiles_k0) *  (1 - file_test(savefiles_k0, /zero_length))
    
    if keyword_set(refresh_binning) then begin
      test_save_2d = test_save_2d*0
      test_save_1d = test_save_1d*0
      if keyword_set(plot_kpar_power) then test_save_kpar = test_save_kpar*0
      if keyword_set(plot_kperp_power) then test_save_kperp = test_save_kperp*0
      if keyword_set(plot_k0_power) then test_save_k0 = test_save_k0*0
    endif
    
    if tag_exist(file_struct_arr[0], 'nside') ne 0 then healpix = 1 else healpix = 0
    
    if tag_exist(file_struct_arr, 'uvf_savefile') then uvf_input = 0 else uvf_input = 1
    
    n_freq = n_elements(file_struct_arr[0].frequencies)
    if n_elements(freq_ch_range) ne 0 then begin
      if max(freq_ch_range) gt n_freq-1 then message, 'invalid freq_ch_range'
      n_freq = freq_ch_range[1]-freq_ch_range[0]+1
    endif
    
    if healpix and n_elements(dft_fchunk) ne 0 then if dft_fchunk gt n_freq then begin
      print, 'dft_fchunk is larger than the number of frequency slices, setting it to the number of slices -- ' + $
        number_formatter(n_freq)
      dft_fchunk = n_freq
    endif
    
    if keyword_set(plot_binning_hist) and pub eq 1 then begin
      if n_elements(plot_filebase) eq 0 then begin
        plotfile_binning_hist = strarr(n_cubes, n_elements(kperp_density_names), n_elements(wedge_1dbin_names))
        
        for i=0, n_elements(wedge_1dbin_names)-1 do begin
          for j=0, n_elements(kperp_density_names)-1 do begin
            plotfile_binning_hist[*,j,i] = plotfile_path + file_struct_arr.savefilebase + power_tag + $
              kperp_density_names[j] + wedge_1dbin_names[i] + fadd_1dbin
          endfor
        endfor
      endif else begin
        plotfile_binning_hist = strarr(n_cubes, n_elements(kperp_density_names), n_elements(wedge_1dbin_names))
        
        for i=0, n_elements(wedge_1dbin_names)-1 do begin
          for j=0, n_elements(kperp_density_names)-1 do begin
            plotfile_binning_hist[*,j,i] = plotfile_path + plot_filebase + uvf_tag + file_struct_arr.file_label + power_tag + $
              kperp_density_names[j] + wedge_1dbin_names[i] + fadd_1dbin
          endfor
        endfor
      endelse
      plotfile_binning_hist = plotfile_binning_hist + '_binning_hist'
    endif
    
    for i=0, n_cubes-1 do begin
    
      savefile_2d_use = reform(savefiles_2d[i,*], n_elements(kperp_density_names))
      test_2d = min(test_save_2d[i,*])
      
      savefile_kpar_use = savefiles_kpar_1d[i,*]
      if keyword_set(plot_kpar_power) then begin
        test_kpar = min(test_save_kpar[i,*])
      endif
      
      savefile_kperp_use = savefiles_kperp_1d[i,*]
      if keyword_set(plot_kperp_power) then begin
        test_kperp = min(test_save_kperp[i,*])
      endif
      
      savefile_k0_use = savefiles_k0[i,*]
      if keyword_set(plot_k0_power) then test_k0 = min(test_save_k0[i,*])
      
      savefile_1d_use = reform(savefiles_1d[i,*,*], n_elements(kperp_density_names), n_elements(wedge_1dbin_names))
      savefiles_mask_use = reform(savefiles_1to2d_mask[i,*,*], n_elements(kperp_density_names), n_elements(wedge_1dbin_names))
      test_1d = min(test_save_1d[i,*,*])
      test_mask = min(test_save_mask[i,*,*])
      
      ;; if binsizes are specified, check that binsize is right
      if (n_elements(kperp_bin) ne 0 or n_elements(kpar_bin) ne 0) and test_2d gt 0 then begin
        if n_elements(kpar_bin) ne 0 then begin
          kpar_bin_file = fltarr(n_elements(savefile_2d_use))
          for j=0, n_elements(savefile_2d_use)-1 do kpar_bin_file[j] = getvar_savefile(savefile_2d_use[j], 'kpar_bin')
          if max(abs(kpar_bin - kpar_bin_file)) gt 0. then test_2d=0
          
          kpar_bin_mask = fltarr(n_elements(savefiles_mask_use))
          for j=0, n_elements(savefiles_mask_use)-1 do kpar_bin_mask[j] = getvar_savefile(savefiles_mask_use[j], 'kpar_bin')
          if max(abs(kpar_bin - kpar_bin_mask)) gt 0. then test_mask=0
        endif
        if n_elements(kperp_bin) ne 0 then begin
          kperp_bin_file = fltarr(n_elements(savefile_2d_use))
          for j=0, n_elements(savefile_2d_use)-1 do kperp_bin_file[j] = getvar_savefile(savefile_2d_use[j], 'kperp_bin')
          if max(abs(kperp_bin - kperp_bin_file)) gt 0. then test_2d=0
          
          kperp_bin_mask = fltarr(n_elements(savefiles_mask_use))
          for j=0, n_elements(savefiles_mask_use)-1 do kperp_bin_mask[j] = getvar_savefile(savefiles_mask_use[j], 'kperp_bin')
          if max(abs(kperp_bin - kperp_bin_mask)) gt 0. then test_mask=0
        endif
      endif
      
      if keyword_set(plot_k0_power) then begin
        if test_k0 gt 0 and n_elements(kperp_bin) ne 0 then begin
          kperp_bin_file = fltarr(n_elements(savefile_k0_use))
          for j=0, n_elements(savefile_k0_use) do kperp_bin_file[j] = getvar_savefile(savefile_k0_use[j], 'k_bin')
          if max(abs(kperp_bin - kperp_bin_file)) gt 0. then test_k0=0
        endif
      endif
      
      if n_elements(k1d_bin) ne 0 and test_1d gt 0 then begin
        k_bin_file = fltarr(n_elements(savefile_1d_use))
        for j=0, n_elements(savefile_1d_use)-1 do k_bin_file[j] = getvar_savefile(savefile_1d_use[j], 'k_bin')
        if max(abs(k1d_bin - k_bin_file)) gt 0. then test_1d=0
        
        k_bin_mask = fltarr(n_elements(savefiles_mask_use))
        for j=0, n_elements(savefiles_mask_use)-1 do k_bin_mask[j] = getvar_savefile(savefiles_mask_use[j], 'k_bin')
        if max(abs(k1d_bin - k_bin_mask)) gt 0. then test_mask=0
      endif
      
      if test_2d gt 0 and n_elements(freq_flags) ne 0 then begin
        for j=0, n_elements(savefile_2d_use)-1 do begin
          old_freq_mask = getvar_savefile(savefile_2d_use[j], 'freq_mask')
          if total(abs(old_freq_mask - file_struct_arr[i].freq_mask)) ne 0 then test_2d = 0
        endfor
        
        for j=0, n_elements(savefiles_mask_use)-1 do begin
          old_freq_mask = getvar_savefile(savefiles_mask_use[j], 'freq_mask')
          if total(abs(old_freq_mask - file_struct_arr[i].freq_mask)) ne 0 then test_mask = 0
        endfor
      endif
      
      if test_1d gt 0 and n_elements(freq_flags) ne 0 then begin
        for j=0, n_elements(savefile_1d_use)-1 do begin
          old_freq_mask = getvar_savefile(savefile_1d_use[j], 'freq_mask')
          if total(abs(old_freq_mask - file_struct_arr[i].freq_mask)) ne 0 then test_1d = 0
        endfor
      endif
      
      if test_1d gt 0 and (n_elements(kperp_range_1d_use) gt 0 or n_elements(kperp_range_lambda_1dave) gt 0 $
        or n_elements(kpar_range_1d_use) gt 0) then begin
        ;; check that 1d binning was over correct ranges
        kperp_range_used = fltarr(2, n_elements(savefile_1d_use))
        kperp_range_lambda_used = fltarr(2, n_elements(savefile_1d_use))
        kpar_range_used = fltarr(2, n_elements(savefile_1d_use))
        for j=0, n_elements(savefile_1d_use)-1 do begin
          kperp_range_used[*,j] = getvar_savefile(savefile_1d_use[j], 'kperp_range')
          kperp_range_lambda_used[*,j] = getvar_savefile(savefile_1d_use[j], 'kperp_range_lambda')
          kpar_range_used[*,j] = getvar_savefile(savefile_1d_use[j], 'kpar_range')
        endfor
        if n_elements(kperp_range_1d_use) gt 0 then if max(abs(1.-kperp_range_used/rebin(kperp_range_1d_use,2,n_elements(savefile_1d_use)))) gt 1e-6 then test_1d = 0
        if n_elements(kperp_range_lambda_1dave) gt 0 then if max(abs(1.-kperp_range_lambda_used/rebin(kperp_range_lambda_1dave,2,n_elements(savefile_1d_use)))) gt 1e-6 then test_1d = 0
        if n_elements(kpar_range_1d_use) gt 0 then if max(abs(1.-kpar_range_used/rebin(kpar_range_1d_use,2,n_elements(savefile_1d_use)))) gt 1e-6 then test_1d = 0
        
        kperp_range_used = fltarr(2, n_elements(savefiles_mask_use))
        kperp_range_lambda_used = fltarr(2, n_elements(savefiles_mask_use))
        kpar_range_used = fltarr(2, n_elements(savefiles_mask_use))
        for j=0, n_elements(savefiles_mask_use)-1 do begin
          kperp_range_used[*,j] = getvar_savefile(savefiles_mask_use[j], 'kperp_range')
          kperp_range_lambda_used[*,j] = getvar_savefile(savefiles_mask_use[j], 'kperp_range_lambda')
          kpar_range_used[*,j] = getvar_savefile(savefiles_mask_use[j], 'kpar_range')
        endfor
        if n_elements(kperp_range_1d_use) gt 0 then if max(abs(1.-kperp_range_used/rebin(kperp_range_1d_use,2,n_elements(savefile_1d_use)))) gt 1e-6 then test_mask = 0
        if n_elements(kperp_range_lambda_1dave) gt 0 then if max(abs(1.-kperp_range_lambda_used/rebin(kperp_range_lambda_1dave,2,n_elements(savefile_1d_use)))) gt 1e-6 then test_mask = 0
        if n_elements(kpar_range_1d_use) gt 0 then if max(abs(1.-kpar_range_used/rebin(kpar_range_1d_use,2,n_elements(savefile_1d_use)))) gt 1e-6 then test_mask = 0
      endif
      
      if keyword_set(plot_kpar_power) then begin
        if test_kpar gt 0 and (n_elements(kperp_range_1d_use) gt 0 or n_elements(kperp_range_lambda_kparpower) gt 0) then begin
          ;; check that 1d binning was over correct ranges
          kperp_range_used = getvar_savefile(savefile_kpar_use, 'kperp_range')
          kperp_range_lambda_used = getvar_savefile(savefile_kpar_use, 'kperp_range_lambda')
          if n_elements(kperp_range_1d_use) gt 0 then if max(abs(1.-kperp_range_used/kperp_range_1d_use)) gt 1e-6 then test_kpar = 0
          if n_elements(kperp_range_lambda_kparpower) gt 0 then if max(abs(1.-kperp_range_lambda_kparpower/rebin(kperp_range_lambda_kparpower,2,n_elements(savefile_1d_use)))) gt 1e-6 then test_kpar = 0
        endif
      endif
      
      if keyword_set(plot_kperp_power) then begin
        if test_kperp gt 0 and n_elements(kpar_range_) gt 0 then begin
          ;; check that 1d binning was over correct ranges
          kperp_range_used = getvar_savefile(savefile_kperp_use, 'kperp_range')
          kperp_range_lambda_used = getvar_savefile(savefile_kperp_use, 'kperp_range_lambda')
          kpar_range_used = getvar_savefile(savefile_kperp_use, 'kpar_range')
          if n_elements(kperp_range_1d_use) gt 0 then if max(abs(1-kperp_range_used/kperp_range_1d_use)) gt 1e-6 then test_kperp = 0
          if n_elements(kperp_range_lambda_1dave) gt 0 then if max(abs(1-kperp_range_lambda_used/rebin(kperp_range_lambda_1dave,2,n_elements(savefile_1d_use)))) gt 1e-6 then test_kperp = 0
          if n_elements(kpar_range_1d_use) gt 0 then if max(abs(1-kpar_range_used/kpar_range_1d_use)) gt 1e-6 then test_kperp = 0
        endif
      endif
      
      test = test_2d * test_1d * test_mask
      if keyword_set(plot_kpar_power) then test = test * test_kpar
      if keyword_set(plot_kperp_power) then test = test * test_kperp
      if keyword_set(plot_k0_power) then test = test * test_k0
      
      if test eq 0 then begin
        if n_elements(plotfile_binning_hist) gt 0 then $
          plotfile_binning_hist_use = reform(plotfile_binning_hist[i,*,*], n_elements(kperp_density_names), n_elements(wedge_1dbin_names))
          
        if healpix or not keyword_set(uvf_input) then begin
          weight_refresh = intarr(n_cubes)
          if keyword_set(refresh_dft) then begin
            temp = weight_ind[uniq(weight_ind, sort(weight_ind))]
            for j=0, n_elements(temp)-1 do weight_refresh[(where(weight_ind eq temp[j]))[0]] = 1
          endif
          
          ps_power, file_struct_arr[i], kcube_refresh = refresh_ps, dft_refresh_data = refresh_dft, $
            dft_refresh_weight = weight_refresh[i], refresh_beam = refresh_beam, $
            dft_ian = dft_ian, cut_image = cut_image, dft_fchunk = dft_fchunk, freq_ch_range = freq_ch_range, freq_flags = freq_flags, $
            spec_window_type = spec_window_type, delta_uv_lambda = delta_uv_lambda, max_uv_lambda = max_uv_lambda, $
            savefile_2d = savefile_2d_use, savefile_1d = savefile_1d_use, savefile_1to2d_mask = savefiles_mask_use, hinv = hinv, $
            savefile_kpar_power = savefile_kpar_use, savefile_kperp_power = savefile_kperp_use, savefile_k0 = savefile_k0_use, $
            std_power = std_power, inverse_covar_weight = inverse_covar_weight, no_wtd_avg = no_wtd_avg, no_kzero = no_kzero, sim=sim, $
            log_k1d = log_k1d, k1d_bin = k1d_bin, log_kpar = log_kpar, log_kperp = log_kperp, kpar_bin = kpar_bin, kperp_bin = kperp_bin, $
            kperp_range_1dave = kperp_range_1d_use, kperp_range_lambda_1dave = kperp_range_lambda_1dave, kpar_range_1dave = kpar_range_1d_use, $
            kperp_range_lambda_kparpower = kperp_range_lambda_kparpower, kpar_range_kperppower = kpar_range_kperppower_use, $
            wt_measures = wt_measures, wt_cutoffs = wt_cutoffs, fix_sim_input = fix_sim_input, $
            wedge_amps = wedge_amp, coarse_harm0 = coarse_harm0, coarse_width = coarse_harm_width, no_dft_progress = no_dft_progress, $
            plot_binning_hist = plot_binning_hist, plotfile_binning_hist = plotfile_binning_hist_use, png = png, eps = eps, pdf = pdf
        endif else $
          ps_power, file_struct_arr[i], kcube_refresh = refresh_ps, refresh_beam = refresh_beam, freq_ch_range = freq_ch_range, $
          freq_flags = freq_flags, spec_window_type = spec_window_type, $
          savefile_2d = savefile_2d_use, savefile_1d = savefile_1d_use, savefile_1to2d_mask = savefiles_mask_use, hinv = hinv, $
          savefile_kpar_power = savefile_kpar_use, savefile_kperp_power = savefile_kperp_use, savefile_k0 = savefile_k0_use, $
          std_power = std_power, inverse_covar_weight = inverse_covar_weight, no_wtd_avg = no_wtd_avg, no_kzero = no_kzero, $
          uvf_input = uvf_input, uv_avg = uv_avg, uv_img_clip = uv_img_clip, sim=sim, $
          log_k1d = log_k1d, k1d_bin = k1d_bin, log_kpar = log_kpar, log_kperp = log_kperp, kpar_bin = kpar_bin, kperp_bin = kperp_bin, $
          kperp_range_1dave = kperp_range_1d_use, kperp_range_lambda_1dave = kperp_range_lambda_1dave, kpar_range_1dave = kpar_range_1d_use, $
          kperp_range_lambda_kparpower = kperp_range_lambda_kparpower, kpar_range_kperppower = kpar_range_kperppower_use, $
          wt_measures = wt_measures, wt_cutoffs = wt_cutoffs, fix_sim_input = fix_sim_input, $
          wedge_amps = wedge_amp, coarse_harm0 = coarse_harm0, coarse_width = coarse_harm_width, no_dft_progress = no_dft_progress, $
          plot_binning_hist = plot_binning_hist, plotfile_binning_hist = plotfile_binning_hist_use, png = png, eps = eps, pdf = pdf
      endif
    endfor
    
    
    restore, savefiles_2d[0]
    if n_elements(window_int) gt 0 then print, 'window integral: ', window_int
    if n_elements(vs_name) ne 0 then begin
      vs_note = vs_name + ': ~' + number_formatter(vs_mean, format = '(f10.2)')
      print, vs_note
    endif
    if n_elements(t_sys_meas) ne 0 then print, 'Tsys range: ', number_formatter(minmax(t_sys_meas))
    if n_elements(t_sys_meas) ne 0 then print, 'Tsys mean: ', number_formatter(mean(t_sys_meas))
    
    if n_elements(git_hashes) ne 0 then print, 'kcube hash: ' + git_hashes.kcube
    
    ;;wh_good_kperp = where(total(power, 2) gt 0, count)
    ;;if count eq 0 then message, '2d power appears to be entirely zero'
    
    ;;kperp_plot_range = [min(kperp_edges[wh_good_kperp]), max(kperp_edges[wh_good_kperp+1])]
    
    ;;kperp_plot_range = [6e-3, min([max(kperp_edges[wh_good_kperp+1]),1.1e-1])]
    ;;kperp_plot_range = [5./kperp_lambda_conv, min([max(kperp_edges[wh_good_kperp+1]),1.1e-1])]
    
    if n_elements(kperp_lambda_plot_range) gt 0 then begin
      kperp_plot_range = kperp_lambda_plot_range / kperp_lambda_conv
      
      ;; if we're plotting in [k]=h/Mpc then need to convert from 1/Mpc
      if keyword_set(hinv) then kperp_plot_range = kperp_plot_range / hubble_param
    endif
    
    if n_elements(kperp_plot_range) eq 0 then begin
      if n_elements(max_uv_lambda) gt 0 then $
        max_kperp_lambda = min([max_uv_lambda, min(file_struct_arr.kspan/2.),min(file_struct_arr.max_baseline_lambda)])$
      else max_kperp_lambda = min([file_struct_arr.kspan/2.,file_struct_arr.max_baseline_lambda])
      kperp_plot_range = [5./kperp_lambda_conv, max_kperp_lambda/kperp_lambda_conv]
      
      ;; if we're plotting in [k]=h/Mpc then need to convert from 1/Mpc
      if keyword_set(hinv) then kperp_plot_range = kperp_plot_range / hubble_param
    endif
    
    if pub then begin
      if n_elements(plot_path) ne 0 then plotfile_path = plot_path $
      else if n_elements(save_path) ne 0 then plotfile_path = save_path else plotfile_path = file_struct_arr.savefile_froot
      
      plot_fadd = ''
      if keyword_set(kperp_linear_axis) and keyword_set(kpar_linear_axis) then plot_fadd = plot_fadd + '_linaxes' else begin
        if keyword_set(kperp_linear_axis) then plot_fadd = plot_fadd + '_kperplinaxis'
        if keyword_set(kpar_linear_axis) then plot_fadd = plot_fadd + '_kparlinaxis'
      endelse
      
      if keyword_set(individual_plots) then begin
        if n_elements(plot_filebase) eq 0 then begin
          plotfile_base = plotfile_path + file_struct_arr.savefilebase + power_tag
          plotfile_base_wt = plotfile_path + general_filebase + '_' + file_struct_arr[uniq(weight_ind, sort(weight_ind))].pol + power_tag
        endif else begin
          plotfile_base = plotfile_path + plot_filebase + uvf_tag + file_struct_arr.file_label + power_tag
          plotfile_base_wt = plotfile_path + plot_filebase + uvf_tag + '_' + file_struct_arr[uniq(weight_ind, sort(weight_ind))].pol + power_tag
        endelse
      endif else begin
        if n_elements(plot_filebase) eq 0 then plotfile_base = plotfile_path + general_filebase + power_tag $
        else plotfile_base = plotfile_path + plot_filebase + uvf_tag + power_tag
        plotfile_base_wt = plotfile_base
      endelse
      
      plotfiles_2d = strarr(n_elements(plotfile_base), n_elements(kperp_density_names))
      for j=0, n_elements(kperp_density_names)-1 do begin
        plotfiles_2d[*,j] = plotfile_base + fadd_2dbin + kperp_density_names[j] + '_2dkpower' + plot_fadd + plot_exten
      endfor
      ;plotfiles_2d = plotfile_base + fadd_2dbin + '_2dkpower' + plot_fadd + plot_exten
      plotfiles_2d_error = plotfile_base_wt + fadd_2dbin + '_2derror' + plot_fadd + plot_exten
      if keyword_set(pub) and keyword_set(individual_plots) then $
        plotfiles_2d_noise_expval = plotfile_base_wt + fadd_2dbin + '_2dnoise_expval' + plot_fadd + plot_exten
      plotfiles_2d_noise = plotfile_base + fadd_2dbin + '_2dnoise' + plot_fadd + plot_exten
      plotfiles_2d_sim_noise = plotfile_base + fadd_2dbin + '_2dsimnoise' + plot_fadd + plot_exten
      plotfiles_2d_sim_noise_diff = plotfile_base + fadd_2dbin + '_2dsimnoisediff' + plot_fadd + plot_exten
      
      plotfiles_2d_snr = strarr(n_elements(plotfile_base), n_elements(kperp_density_names))
      for j=0, n_elements(kperp_density_names)-1 do begin
        plotfiles_2d_snr[*,j] = plotfile_base + fadd_2dbin + kperp_density_names[j] + '_2dsnr' + plot_fadd + plot_exten
      endfor
      
      plotfiles_2d_nnr = plotfile_base + fadd_2dbin + '_2dnnr' + plot_fadd + plot_exten
      plotfiles_2d_sim_snr = plotfile_base + fadd_2dbin + '_2dsimsnr' + plot_fadd + plot_exten
      plotfiles_2d_sim_nnr = plotfile_base + fadd_2dbin + '_2dsimnnr' + plot_fadd + plot_exten
      
      
      if n_elements(plot_filebase) eq 0 then plotfile_1d_base = plotfile_path + general_filebase else $
        plotfile_1d_base = plotfile_path + plot_filebase + uvf_tag
      plotfile_1d = plotfile_1d_base + power_tag + wedge_1dbin_names + fadd_1dbin + '_1dkpower' + plot_exten
      plotfile_1d_noise = plotfile_1d_base + power_tag + wedge_1dbin_names + fadd_1dbin + '_1dnoise' + plot_exten
      plotfile_1d_sim_noise_diff = plotfile_1d_base + power_tag + wedge_1dbin_names + fadd_1dbin + '_1dsimnoisediff' + plot_exten
      plotfile_1d_sim_noise = plotfile_1d_base + power_tag + wedge_1dbin_names + fadd_1dbin + '_1dsimnoise' + plot_exten
      plotfile_kpar_power = plotfile_1d_base + power_tag + fadd_kpar_1d + '_kpar_power' + plot_exten
      plotfile_kperp_power = plotfile_1d_base + power_tag + fadd_kperp_1d + '_kperp_power' + plot_exten
      plotfile_k0_power = plotfile_1d_base + power_tag + fadd_2dbin + '_k0power' + plot_exten
      
    endif
    
    if keyword_set(plot_1to2d) then begin
      if keyword_set(pub) then begin
        plotfile_1to2d_heatmap = plotfile_base + fadd_2dbin + fadd_1dbin + '_1dheatmap' + plot_exten
        plotfile_1to2d_contours = plotfile_base + fadd_2dbin + fadd_1dbin + '_1dcontours' + plot_exten
        plotfile_1to2d_noisefrac = plotfile_base + fadd_2dbin + fadd_1dbin + '_1dnoisefrac' + plot_exten
        plotfile_1to2d_contour_zoom = plotfile_base + fadd_2dbin + fadd_1dbin + '_1dcontour_zoom' + plot_exten
      endif
      type_1to2d_use = 'res'
      pol_1to2d_use = 'xx'
      wh_2d_use = where(file_struct_arr.type eq type_1to2d_use and file_struct_arr.pol eq pol_1to2d_use, count_type_pol)
      if count_type_pol eq 0 then wh_2d_use = 0
    endif
    
    if keyword_set(plot_slices) then begin
      if n_elements(slice_type) eq 0 then slice_type = 'sumdiff'
      slice_type_enum = ['raw', 'divided', 'kspace', 'sumdiff', 'weights']
      
      wh_slice_type = where(slice_type_enum eq slice_type, count_slice_type)
      if count_slice_type eq 0 then begin
        print, 'slice_type not recognized, using default'
        slice_type = 'sumdiff'
      endif
      
      if slice_type ne 'kspace' then begin
        ;if slice_type eq 'weights' then uvf_plot_type='weights'
      
        uvf_type_enum = ['abs', 'phase', 'real', 'imaginary', 'weights']
        if n_elements(uvf_plot_type) eq 0 then uvf_plot_type = 'abs'
        wh = where(uvf_type_enum eq uvf_plot_type, count_uvf_type)
        if count_uvf_type eq 0 then message, 'unknown uvf_plot_type. Use one of: ' + print, strjoin(uvf_type_enum, ', ')
        
        if uvf_plot_type eq 'phase' then uvf_log = 0 else uvf_log=1
      endif
      
      if pub then begin
      
        if keyword_set(individual_plots) then begin
        
          case slice_type of
            'raw': begin
              uf_slice_plotfile = transpose([[plotfile_base +'_even_'],[plotfile_base +'_odd_']]) + '_' + slice_type + '_uf_plane' + plot_fadd + plot_exten
              vf_slice_plotfile = transpose([[plotfile_base +'_even_'],[plotfile_base +'_odd_']]) + '_' + slice_type + '_vf_plane' + plot_fadd + plot_exten
              uv_slice_plotfile = transpose([[plotfile_base +'_even_'],[plotfile_base +'_odd_']]) + '_' + slice_type + '_uv_plane' + plot_fadd + plot_exten
            end
            'divided': begin
              uf_slice_plotfile = transpose([[plotfile_base +'_even_'],[plotfile_base +'_odd_']]) + '_' + slice_type + '_uf_plane' + plot_fadd + plot_exten
              vf_slice_plotfile = transpose([[plotfile_base +'_even_'],[plotfile_base +'_odd_']]) + '_' + slice_type + '_vf_plane' + plot_fadd + plot_exten
              uv_slice_plotfile = transpose([[plotfile_base +'_even_'],[plotfile_base +'_odd_']]) + '_' + slice_type + '_uv_plane' + plot_fadd + plot_exten
            end
            'sumdiff': begin
              uf_slice_plotfile = [plotfile_base +'_sum_',plotfile_base +'_diff_'] + '_' + slice_type + '_uf_plane' + plot_fadd + plot_exten
              vf_slice_plotfile = [plotfile_base +'_sum_',plotfile_base +'_diff_'] + '_' + slice_type + '_vf_plane' + plot_fadd + plot_exten
              uv_slice_plotfile = [plotfile_base +'_sum_',plotfile_base +'_diff_'] + '_' + slice_type + '_uv_plane' + plot_fadd + plot_exten
            end
            'weights': begin
              uf_slice_plotfile = transpose([[plotfile_base +'_even_'],[plotfile_base +'_odd_']]) + '_' + slice_type + '_uf_plane' + plot_fadd + plot_exten
              vf_slice_plotfile = transpose([[plotfile_base +'_even_'],[plotfile_base +'_odd_']]) + '_' + slice_type + '_vf_plane' + plot_fadd + plot_exten
              uv_slice_plotfile = transpose([[plotfile_base +'_even_'],[plotfile_base +'_odd_']]) + '_' + slice_type + '_uv_plane' + plot_fadd + plot_exten
            end
            'kspace': begin
              uf_slice_plotfile = plotfile_base + '_' + slice_type + '_uf_plane' + plot_fadd + plot_exten
              vf_slice_plotfile = plotfile_base + '_' + slice_type + '_vf_plane' + plot_fadd + plot_exten
              uv_slice_plotfile = plotfile_base + '_' + slice_type + '_uv_plane' + plot_fadd + plot_exten
            end
          endcase
          
        endif else begin
          if slice_type ne 'kspace' then begin
            uf_slice_plotfile = plotfile_base + '_' + slice_type + '_uf_plane' + plot_fadd + plot_exten
            vf_slice_plotfile = plotfile_base + '_' + slice_type + '_vf_plane' + plot_fadd + plot_exten
            uv_slice_plotfile = plotfile_base + '_' + slice_type + '_uv_plane' + plot_fadd + plot_exten
          endif else begin
            uf_slice_plotfile = plotfile_base + '_' + slice_type + '_xz_plane' + plot_fadd + plot_exten
            vf_slice_plotfile = plotfile_base + '_' + slice_type + '_yz_plane' + plot_fadd + plot_exten
            uv_slice_plotfile = plotfile_base + '_' + slice_type + '_xy_plane' + plot_fadd + plot_exten
          endelse
        endelse
        
      endif else begin
        case slice_type of
          'raw': slice_titles = file_struct_arr.uvf_label
          'divided': slice_titles = file_struct_arr.uvf_label
          'sumdiff': slice_titles = ['sum ' + file_struct_arr.file_label, 'diff ' + file_struct_arr.file_label]
          'weights': slice_titles = file_struct_arr.uvf_label
          'kspace': slice_titles = file_struct_arr.file_label
        endcase
      endelse
      
      case slice_type of
        'raw': begin
          uf_slice_savefile = file_struct_arr.uf_raw_savefile
          vf_slice_savefile = file_struct_arr.vf_raw_savefile
          uv_slice_savefile = file_struct_arr.uv_raw_savefile
          slice_titles = file_struct_arr.uvf_label
        end
        'divided': begin
          uf_slice_savefile = file_struct_arr.uf_savefile
          vf_slice_savefile = file_struct_arr.vf_savefile
          uv_slice_savefile = file_struct_arr.uv_savefile
          slice_titles = file_struct_arr.uvf_label
        end
        'sumdiff': begin
          uf_slice_savefile = [file_struct_arr.uf_sum_savefile, file_struct_arr.uf_diff_savefile]
          vf_slice_savefile = [file_struct_arr.vf_sum_savefile, file_struct_arr.vf_diff_savefile]
          uv_slice_savefile = [file_struct_arr.uv_sum_savefile, file_struct_arr.uv_diff_savefile]
          slice_titles = ['sum ' + file_struct_arr.file_label, 'diff ' + file_struct_arr.file_label]
        end
        'weights': begin
          uf_slice_savefile = file_struct_arr.uf_weight_savefile
          vf_slice_savefile = file_struct_arr.vf_weight_savefile
          uv_slice_savefile = file_struct_arr.uv_weight_savefile
          slice_titles = file_struct_arr.uvf_label
        end
        'kspace': begin
          uf_slice_savefile = file_struct_arr.savefile_froot + file_struct_arr.savefilebase + power_tag + '_xz_plane.idlsave'
          vf_slice_savefile = file_struct_arr.savefile_froot + file_struct_arr.savefilebase + power_tag + '_yz_plane.idlsave'
          uv_slice_savefile = file_struct_arr.savefile_froot + file_struct_arr.savefilebase + power_tag + '_xy_plane.idlsave'
          slice_titles = file_struct_arr.file_label
        end
      endcase
      
    endif
    
    if keyword_set(pub) then font = 1 else font = -1
    
    window_num=0
    if keyword_set(plot_stdset) then begin
      if keyword_set(pub) and keyword_set(individual_plots) then begin
        for j=0, n_elements(kperp_density_names)-1 do begin
          note_use = note + ', ' + kperp_density_names[j]
          
          for i=0, n_cubes-1 do kpower_2d_plots, savefiles_2d[i,j], png = png, eps = eps, pdf = pdf, plotfile = plotfiles_2d[i,j], $
            kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
            data_range = data_range, title_prefix = titles[i], note = note_use, $
            plot_wedge_line = plot_wedge_line, hinv = hinv, $
            wedge_amp = wedge_amp, baseline_axis = baseline_axis, delay_axis = delay_axis, cable_length_axis = cable_length_axis, $
            kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis
        endfor
        
        for i=0, npol -1 do kpower_2d_plots, savefiles_2d[i,0], /plot_sigma, png = png, eps = eps, pdf = pdf, plotfile = plotfiles_2d_error[i],$
          kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
          data_range = sigma_range, title_prefix = file_struct_arr[i].pol, note = note + ' ' + vs_note, $
          plot_wedge_line = plot_wedge_line, hinv = hinv, $
          wedge_amp = wedge_amp, baseline_axis = baseline_axis, delay_axis = delay_axis, cable_length_axis = cable_length_axis, $
          kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis
          
        for i=0, npol -1 do kpower_2d_plots, savefiles_2d[i,0], /plot_exp_noise, png = png, eps = eps, pdf = pdf, plotfile = plotfiles_2d_noise_expval[i],$
          kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
          data_range = nev_range, title_prefix = file_struct_arr[i].pol, note = note + ' ' + vs_note, $
          plot_wedge_line = plot_wedge_line, hinv = hinv, $
          wedge_amp = wedge_amp, baseline_axis = baseline_axis, delay_axis = delay_axis, cable_length_axis = cable_length_axis, $
          kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis
          
        for j=0, n_elements(kperp_density_names)-1 do begin
          note_use = note + ', ' + kperp_density_names[j]
          
          for i=0, n_cubes-1 do kpower_2d_plots, savefiles_2d[i,j], png = png, eps = eps, pdf = pdf, /snr, plotfile = plotfiles_2d_snr[i,j], $
            kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
            data_range = snr_range, title_prefix = titles[i], note = note, $
            plot_wedge_line = plot_wedge_line, hinv = hinv, $
            wedge_amp = wedge_amp, baseline_axis = baseline_axis, delay_axis = delay_axis, cable_length_axis = cable_length_axis, $
            kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis
        endfor
        
        if nfiles eq 2 then begin
          for i=0, n_cubes-1 do kpower_2d_plots, savefiles_2d[i,0], png = png, eps = eps, pdf = pdf, /plot_noise, $
            plotfile = plotfiles_2d_noise[i], $
            kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
            data_range = noise_range, title_prefix = titles[i], note = note + ' ' + vs_note, $
            plot_wedge_line = plot_wedge_line, hinv = hinv, $
            wedge_amp = wedge_amp, baseline_axis = baseline_axis, delay_axis = delay_axis, cable_length_axis = cable_length_axis, $
            kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis
            
          if keyword_set(plot_sim_noise) then $
            for i=0, n_cubes-1 do kpower_2d_plots, savefiles_2d[i,0], png = png, eps = eps, pdf = pdf, /plot_sim_noise, $
            plotfile = plotfiles_2d_sim_noise[i], $
            kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
            data_range = noise_range, title_prefix = titles[i], note = note + ' ' + vs_note, $
            plot_wedge_line = plot_wedge_line, hinv = hinv, $
            wedge_amp = wedge_amp, baseline_axis = baseline_axis, delay_axis = delay_axis, cable_length_axis = cable_length_axis, $
            kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis
            
          if keyword_set(plot_sim_noise) then $
            for i=0, n_cubes-1 do kpower_2d_plots, savefiles_2d[i,0], png = png, eps = eps, pdf = pdf, /plot_simnoise_diff, $
            plotfile = plotfiles_2d_sim_noise_diff[i], $
            kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
            data_range = noise_range, title_prefix = titles[i], note = note + ' ' + vs_note, $
            plot_wedge_line = plot_wedge_line, hinv = hinv, $
            wedge_amp = wedge_amp, baseline_axis = baseline_axis, delay_axis = delay_axis, cable_length_axis = cable_length_axis, $
            kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis
            
          for i=0, n_cubes-1 do kpower_2d_plots, savefiles_2d[i,0], png = png, eps = eps, pdf = pdf, /nnr, plotfile = plotfiles_2d_nnr[i], $
            kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
            data_range = nnr_range, title_prefix = titles[i], note = note + ' ' + vs_note, $
            plot_wedge_line = plot_wedge_line, hinv = hinv, $
            wedge_amp = wedge_amp, baseline_axis = baseline_axis, delay_axis = delay_axis, cable_length_axis = cable_length_axis, $
            kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis
            
          for i=0, n_cubes-1 do kpower_2d_plots, savefiles_2d[i,0], png = png, eps = eps, pdf = pdf, /sim_nnr, plotfile = plotfiles_2d_sim_nnr[i], $
            kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
            data_range = nnr_range, title_prefix = titles[i], note = note + ' ' + vs_note, $
            plot_wedge_line = plot_wedge_line, hinv = hinv, $
            wedge_amp = wedge_amp, baseline_axis = baseline_axis, delay_axis = delay_axis, cable_length_axis = cable_length_axis, $
            kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis
        endif
      endif else begin
        if ntype gt 1 then begin
          cube_inds = indgen(ntype,npol)
          plot_cube_order = intarr(ntype,npol)
          for pol_i=0, npol-1 do begin
            wh_pol = where(file_struct_arr.pol_index eq pol_i, count_wh_pol)
            if count_wh_pol eq 0 then message, 'no cubes for pol_index = ' + pol_i
            plot_type_order = sort(file_struct_arr[wh_pol].type)
            plot_cube_order[*,pol_i] = cube_inds[plot_type_order,pol_i]
          endfor
        endif else plot_cube_order = indgen(npol)
        
        if keyword_set(kperp_linear_axis) then begin
          ;; aspect ratio doesn't work out for kperp_linear with multiple rows
          ncol = ntype*npol
          nrow = 1
        endif else begin
          ncol = ntype
          nrow = npol
        endelse
        for j=0, n_elements(kperp_density_names)-1 do begin
          start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}
          
          window_num = window_num+1
          for i=0, n_cubes-1 do begin
            cube_i = plot_cube_order[i]
            if i gt 0 then  pos_use = positions[*,i]
            if i eq n_cubes-1 and n_elements(note) gt 0 then note_use = note + ', ' + kperp_density_names[j] else undefine, note_use
            if keyword_set(pub) then plotfile_use = plotfiles_2d[j] else undefine, plotfile_use
            
            kpower_2d_plots, savefiles_2d[cube_i,j], multi_pos = pos_use, start_multi_params = start_multi_params, png = png, eps = eps, pdf = pdf, $
              plotfile = plotfile_use, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
              data_range = data_range, title_prefix = titles[cube_i], note = note_use, $
              plot_wedge_line = plot_wedge_line, hinv = hinv, $
              wedge_amp = wedge_amp, baseline_axis = baseline_axis, delay_axis = delay_axis, cable_length_axis = cable_length_axis, $
              kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, window_num = window_num
              
            if i eq 0 then begin
              positions = pos_use
              undefine, start_multi_params
            endif
          endfor
          undefine, positions, pos_use
          if keyword_set(pub) then begin
            cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density = 600
          endif
        endfor
        if keyword_set(kperp_linear_axis) then begin
          ;; aspect ratio doesn't work out for kperp_linear with multiple rows
          ncol = 2*npol
          nrow = 1
        endif else begin
          ncol = 2
          nrow = npol
        endelse
        start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}
        
        window_num = window_num+1
        for i=0, npol*2-1 do begin
          if i gt 0 then pos_use = positions[*,i]
          if i eq npol*2-1 and n_elements(note) gt 0 then note_use = note + ' ' + vs_note else undefine, note_use
          
          pol_ind = i / 2
          
          if i mod 2 eq 0 then $
            kpower_2d_plots, savefiles_2d[pol_ind,0], multi_pos = pos_use, start_multi_params = start_multi_params, png = png, eps = eps, pdf = pdf, $
            plotfile = plotfiles_2d_error, /plot_sigma, data_range = sigma_range, $
            kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
            title_prefix = file_struct_arr[pol_ind].pol, $
            plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp, hinv = hinv, baseline_axis = baseline_axis, $
            delay_axis = delay_axis, cable_length_axis = cable_length_axis, kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, $
            window_num = window_num else $
            kpower_2d_plots, savefiles_2d[pol_ind,0], multi_pos = pos_use, start_multi_params = start_multi_params, $
            png = png, eps = eps, pdf = pdf, /plot_exp_noise, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
            data_range = nev_range, title_prefix = file_struct_arr[pol_ind].pol, note = note_use, $
            plot_wedge_line = plot_wedge_line, hinv = hinv, $
            wedge_amp = wedge_amp, baseline_axis = baseline_axis, delay_axis = delay_axis, cable_length_axis = cable_length_axis, $
            kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis
          if i eq 0 then begin
            positions = pos_use
            undefine, start_multi_params
          endif
        endfor
        undefine, positions, pos_use
        if keyword_set(pub) then begin
          cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density = 600
        endif
        
        
        ;; now plot SNR -- no separate sigma plots
        if keyword_set(kperp_linear_axis) then begin
          ;; aspect ratio doesn't work out for kperp_linear with multiple rows
          ncol = ntype*npol
          nrow = 1
        endif else begin
          ncol = ntype
          nrow = npol
        endelse
        for j=0, n_elements(kperp_density_names)-1 do begin
          start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}
          
          window_num = window_num+1
          ;;snr_range = [1e0, 1e6]
          for i=0, n_cubes-1 do begin
            cube_i = plot_cube_order[i]
            if i gt 0 then  pos_use = positions[*,i]
            if i eq n_cubes-1 and n_elements(note) gt 0 then note_use = note + ', ' + kperp_density_names[j] else undefine, note_use
            if keyword_set(pub) then plotfile_use = plotfiles_2d_snr[j] else undefine, plotfile_use
            
            kpower_2d_plots, savefiles_2d[cube_i,j], /snr, multi_pos = pos_use, start_multi_params = start_multi_params, png = png, eps = eps, pdf = pdf, $
              plotfile = plotfile_use, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
              data_range = snr_range, title_prefix = titles[cube_i], note = note_use, $
              plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp, $
              hinv = hinv, baseline_axis = baseline_axis, delay_axis = delay_axis, cable_length_axis = cable_length_axis, $
              kperp_linear_axis = kperp_linear_axis, $
              kpar_linear_axis = kpar_linear_axis, window_num = window_num
            if i eq 0 then begin
              positions = pos_use
              undefine, start_multi_params
            endif
          endfor
          undefine, positions, pos_use
          if keyword_set(pub) then begin
            cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density = 600
          endif
        endfor
        
        if keyword_set(plot_sim_noise) then begin
          window_num = window_num+1
          start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}
          undefine, pos_use
          
          for i=0, n_cubes-1 do begin
            cube_i = plot_cube_order[i]
            if i gt 0 then  pos_use = positions[*,i]
            if i eq n_cubes-1 and n_elements(note) gt 0 then note_use = note + ' ' + vs_note else undefine, note_use
            
            kpower_2d_plots, savefiles_2d[cube_i,0], /plot_sim_noise, multi_pos = pos_use, start_multi_params = start_multi_params, png = png, eps = eps, pdf = pdf, $
              plotfile = plotfiles_2d_sim_noise, $
              kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, data_range = sigma_range, $
              title_prefix = titles[cube_i], note = note_use, $
              plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp, hinv = hinv, $
              baseline_axis = baseline_axis, delay_axis = delay_axis, cable_length_axis = cable_length_axis, $
              kperp_linear_axis = kperp_linear_axis, $
              kpar_linear_axis = kpar_linear_axis, window_num = window_num
            if i eq 0 then begin
              positions = pos_use
              undefine, start_multi_params
            endif
          endfor
          if keyword_set(pub) then begin
            cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density = 600
          endif
          
          window_num = window_num+1
          start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}
          undefine, pos_use
          
          for i=0, n_cubes-1 do begin
            cube_i = plot_cube_order[i]
            if i gt 0 then  pos_use = positions[*,i]
            if i eq n_cubes-1 and n_elements(note) gt 0 then note_use = note + ' ' + vs_note else undefine, note_use
            
            kpower_2d_plots, savefiles_2d[cube_i,0], /sim_snr, multi_pos = pos_use, start_multi_params = start_multi_params, $
              png = png, eps = eps, pdf = pdf, plotfile = plotfiles_2d_sim_snr, $
              kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, data_range = nnr_range, $
              title_prefix = titles[cube_i], note = note_use, $
              plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp, hinv = hinv, $
              baseline_axis = baseline_axis, delay_axis = delay_axis, cable_length_axis = cable_length_axis, $
              kperp_linear_axis = kperp_linear_axis, $
              kpar_linear_axis = kpar_linear_axis, window_num = window_num
            if i eq 0 then begin
              positions = pos_use
              undefine, start_multi_params
            endif
          endfor
          if keyword_set(pub) then begin
            cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density = 600
          endif
        endif
        
        if nfiles eq 2 then begin
        
          window_num = window_num+1
          start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}
          undefine, pos_use
          
          for i=0, n_cubes-1 do begin
            cube_i = plot_cube_order[i]
            if i gt 0 then  pos_use = positions[*,i]
            if i eq n_cubes-1 and n_elements(note) gt 0 then note_use = note + ' ' + vs_note else undefine, note_use
            
            kpower_2d_plots, savefiles_2d[cube_i,0], /plot_noise, multi_pos = pos_use, start_multi_params = start_multi_params, png = png, eps = eps, pdf = pdf, $
              plotfile = plotfiles_2d_noise, $
              kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, data_range = noise_range, $
              title_prefix = titles[cube_i], note = note_use, $
              plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp, hinv = hinv, $
              baseline_axis = baseline_axis, delay_axis = delay_axis, cable_length_axis = cable_length_axis, $
              kperp_linear_axis = kperp_linear_axis, $
              kpar_linear_axis = kpar_linear_axis, window_num = window_num
            if i eq 0 then begin
              positions = pos_use
              undefine, start_multi_params
            endif
          endfor
          if keyword_set(pub) then begin
            cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density = 600
          endif
          
          if keyword_set(plot_sim_noise) then begin
            window_num = window_num+1
            start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}
            undefine, pos_use
            
            for i=0, n_cubes-1 do begin
              cube_i = plot_cube_order[i]
              if i gt 0 then  pos_use = positions[*,i]
              if i eq n_cubes-1 and n_elements(note) gt 0 then note_use = note + ' ' + vs_note else undefine, note_use
              
              kpower_2d_plots, savefiles_2d[cube_i,0], /plot_simnoise_diff, multi_pos = pos_use, start_multi_params = start_multi_params, png = png, eps = eps, pdf = pdf, $
                plotfile = plotfiles_2d_sim_noise_diff, $
                kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, data_range = noise_range, $
                title_prefix = titles[cube_i], note = note_use, $
                plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp, hinv = hinv, $
                baseline_axis = baseline_axis, delay_axis = delay_axis, cable_length_axis = cable_length_axis, $
                kperp_linear_axis = kperp_linear_axis, $
                kpar_linear_axis = kpar_linear_axis, window_num = window_num
              if i eq 0 then begin
                positions = pos_use
                undefine, start_multi_params
              endif
            endfor
            if keyword_set(pub) then begin
              cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density = 600
            endif
          endif
          
          window_num = window_num+1
          start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}
          undefine, pos_use
          
          for i=0, n_cubes-1 do begin
            cube_i = plot_cube_order[i]
            if i gt 0 then  pos_use = positions[*,i]
            if i eq n_cubes-1 and n_elements(note) gt 0 then note_use = note + ' ' + vs_note else undefine, note_use
            
            kpower_2d_plots, savefiles_2d[cube_i,0], /nnr, multi_pos = pos_use, start_multi_params = start_multi_params, png = png, eps = eps, pdf = pdf, $
              plotfile = plotfiles_2d_nnr, $
              kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, data_range = nnr_range, $
              title_prefix = titles[cube_i,0], note = note_use, $
              plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp, hinv = hinv, $
              baseline_axis = baseline_axis, delay_axis = delay_axis, cable_length_axis = cable_length_axis, $
              kperp_linear_axis = kperp_linear_axis, $
              kpar_linear_axis = kpar_linear_axis, window_num = window_num
            if i eq 0 then begin
              positions = pos_use
              undefine, start_multi_params
            endif
          endfor
          if keyword_set(pub) then begin
            cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density = 600
          endif
          
          window_num = window_num+1
          start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}
          undefine, pos_use
          
          if keyword_set(plot_sim_noise) then begin
            for i=0, n_cubes-1 do begin
              cube_i = plot_cube_order[i]
              if i gt 0 then  pos_use = positions[*,i]
              if i eq n_cubes-1 and n_elements(note) gt 0 then note_use = note + ' ' + vs_note else undefine, note_use
              
              kpower_2d_plots, savefiles_2d[cube_i,0], /sim_nnr, multi_pos = pos_use, start_multi_params = start_multi_params, png = png, eps = eps, pdf = pdf, $
                plotfile = plotfiles_2d_sim_nnr, $
                kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, data_range = nnr_range, $
                title_prefix = titles[cube_i], note = note_use, $
                plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp, hinv = hinv, $
                baseline_axis = baseline_axis, delay_axis = delay_axis, cable_length_axis = cable_length_axis, $
                kperp_linear_axis = kperp_linear_axis, $
                kpar_linear_axis = kpar_linear_axis, window_num = window_num
              if i eq 0 then begin
                positions = pos_use
                undefine, start_multi_params
              endif
            endfor
            if keyword_set(pub) then begin
              cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density = 600
            endif
          endif
          
        endif
      endelse
    endif
    
    if keyword_set(plot_slices) then begin
    
      if keyword_set(pub) and keyword_set(individual_plots) then begin
        if slice_type eq 'kspace' then nplots = ntype*npol else nplots = ntype*nfiles*npol
        
        for i=0, nplots-1 do begin
        
          if slice_type eq 'kspace' then begin
            kpower_slice_plot, uf_slice_savefile[i], multi_pos = pos_use, start_multi_params = start_multi_params, png = png, eps = eps, pdf = pdf, $
              plotfile = uf_slice_plotfile[i], plot_xrange = kperp_plot_range, plot_yrange = kpar_plot_range, $
              title_prefix = slice_titles[i], note = note_use, data_range = slice_range, $
              plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp, hinv = hinv, $
              baseline_axis = baseline_axis, delay_axis = delay_axis, cable_length_axis = cable_length_axis, linear_axes = kperp_linear_axis, $
              window_num = window_num
              
            kpower_slice_plot, vf_slice_savefile[i], multi_pos = pos_use, start_multi_params = start_multi_params, png = png, eps = eps, pdf = pdf, $
              plotfile = vf_slice_plotfile[i], plot_xrange = kperp_plot_range, plot_yrange = kpar_plot_range, $
              title_prefix = slice_titles[i], note = note_use, data_range = slice_range, $
              plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp, hinv = hinv, $
              baseline_axis = baseline_axis, delay_axis = delay_axis, cable_length_axis = cable_length_axis, linear_axes = kperp_linear_axis, $
              window_num = window_num
              
            kpower_slice_plot, uv_slice_savefile[i], multi_pos = pos_use, start_multi_params = start_multi_params, png = png, eps = eps, pdf = pdf, $
              plotfile = uv_slice_plotfile[i], plot_xrange = kperp_plot_range, plot_yrange = kpar_plot_range, $
              title_prefix = slice_titles[i], note = note_use, data_range = slice_range, $
              plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp, hinv = hinv, $
              baseline_axis = baseline_axis, delay_axis = delay_axis, cable_length_axis = cable_length_axis, linear_axes = kperp_linear_axis, $
              window_num = window_num
              
          endif else begin
          
            uvf_slice_plot, uf_slice_savefile[i], multi_pos = pos_use, start_multi_params = start_multi_params, png = png, eps = eps, pdf = pdf, $
              plotfile = uf_slice_plotfile[i], type = uvf_plot_type, title_prefix = slice_titles[i], note = note_use, $
              hinv = hinv, data_range = slice_range, log = uvf_log, $
              baseline_axis = baseline_axis, window_num = window_num
              
            uvf_slice_plot, vf_slice_savefile[i], multi_pos = pos_use, start_multi_params = start_multi_params, png = png, eps = eps, pdf = pdf, $
              plotfile = vf_slice_plotfile[i], type = uvf_plot_type, title_prefix = slice_titles[i], note = note_use, $
              hinv = hinv, data_range = slice_range, log = uvf_log, $
              baseline_axis = baseline_axis, window_num = window_num
              
            uvf_slice_plot, uv_slice_savefile[i], multi_pos = pos_use, start_multi_params = start_multi_params, png = png, eps = eps, pdf = pdf, $
              plotfile = uv_slice_plotfile[i], type = uvf_plot_type, title_prefix = slice_titles[i], note = note_use, $
              hinv = hinv, data_range = slice_range, log = uvf_log, $
              baseline_axis = baseline_axis, window_num = window_num
              
          endelse
          
        endfor
        
      endif else begin
      
        if keyword_set(kperp_linear_axis) then begin
          ;; aspect ratio doesn't work out for kperp_linear with multiple rows
          nrow = 1
          if slice_type eq 'kspace' then ncol = ntype*npol else ncol = ntype*nfiles*npol
        endif else begin
          if slice_type eq 'kspace' then begin
            ncol=ntype
            nrow = npol
          endif else begin
            if ntype gt 1 then begin
              ncol=ntype
              nrow = npol*nfiles
            endif else begin
              ncol=ntype*nfiles
              nrow = npol
            endelse
          endelse
        endelse
        
        if slice_type eq 'raw' or slice_type eq 'divided' then begin
          slice_titles = transpose(slice_titles)
          uf_slice_savefile = transpose(uf_slice_savefile)
          vf_slice_savefile = transpose(vf_slice_savefile)
          uv_slice_savefile = transpose(uv_slice_savefile)
        endif
        
        window_num = window_num+1
        start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}
        undefine, pos_use
        
        for i=0, (nrow*ncol)-1 do begin
          if i gt 0 then  pos_use = positions[*,i]
          if i eq (nrow*ncol)-1 and n_elements(note) gt 0 then note_use = note else undefine, note_use
          
          if slice_type eq 'kspace' then begin
            kpower_slice_plot, uf_slice_savefile[i], multi_pos = pos_use, start_multi_params = start_multi_params, png = png, eps = eps, pdf = pdf, $
              plotfile = uf_slice_plotfile, plot_xrange = kperp_plot_range, plot_yrange = kpar_plot_range, $
              title_prefix = slice_titles[i], note = note_use, data_range = slice_range, $
              plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp, hinv = hinv, $
              baseline_axis = baseline_axis, delay_axis = delay_axis, cable_length_axis = cable_length_axis, linear_axes = kperp_linear_axis, $
              window_num = window_num
          endif else begin
          
            uvf_slice_plot, uf_slice_savefile[i], multi_pos = pos_use, start_multi_params = start_multi_params, png = png, eps = eps, pdf = pdf, $
              plotfile = uf_slice_plotfile, type = uvf_plot_type, title_prefix = slice_titles[i], note = note_use, $
              hinv = hinv, data_range = slice_range, log = uvf_log, $
              baseline_axis = baseline_axis, window_num = window_num
          endelse
          
          if i eq 0 then begin
            positions = pos_use
            undefine, start_multi_params
          endif
        endfor
        if keyword_set(pub) then begin
          cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density = 600
        endif
        
        window_num = window_num+1
        start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}
        undefine, pos_use
        
        for i=0, (nrow*ncol)-1 do begin
          if i gt 0 then  pos_use = positions[*,i]
          if i eq (nrow*ncol)-1 and n_elements(note) gt 0 then note_use = note else undefine, note_use
          
          if slice_type eq 'kspace' then begin
            kpower_slice_plot, vf_slice_savefile[i], multi_pos = pos_use, start_multi_params = start_multi_params, png = png, eps = eps, pdf = pdf, $
              plotfile = vf_slice_plotfile, $
              plot_xrange = kperp_plot_range, plot_yrange = kpar_plot_range, $
              title_prefix = slice_titles[i], note = note_use, data_range = slice_range, $
              plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp, hinv = hinv, $
              baseline_axis = baseline_axis, delay_axis = delay_axis, cable_length_axis = cable_length_axis, linear_axes = kperp_linear_axis, $
              window_num = window_num
          endif else begin
          
            uvf_slice_plot, vf_slice_savefile[i], multi_pos = pos_use, start_multi_params = start_multi_params, png = png, eps = eps, pdf = pdf, $
              plotfile = vf_slice_plotfile, type = uvf_plot_type, title_prefix = slice_titles[i], note = note_use, $
              hinv = hinv, data_range = slice_range, log = uvf_log, $
              baseline_axis = baseline_axis, window_num = window_num
          endelse
          
          if i eq 0 then begin
            positions = pos_use
            undefine, start_multi_params
          endif
        endfor
        if keyword_set(pub) then begin
          cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density = 600
        endif
        
        window_num = window_num+1
        start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}
        undefine, pos_use
        
        for i=0, (nrow*ncol)-1 do begin
          if i gt 0 then  pos_use = positions[*,i]
          if i eq (nrow*ncol)-1 and n_elements(note) gt 0 then note_use = note else undefine, note_use
          
          if slice_type eq 'kspace' then begin
            kpower_slice_plot, uv_slice_savefile[i], multi_pos = pos_use, start_multi_params = start_multi_params, png = png, eps = eps, pdf = pdf, $
              plotfile = uv_slice_plotfile, $
              plot_xrange = kperp_plot_range, plot_yrange = kperp_plot_range, $
              title_prefix = slice_titles[i], note = note_use, data_range = slice_range, $
              plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp, hinv = hinv, $
              baseline_axis = baseline_axis, delay_axis = delay_axis, cable_length_axis = cable_length_axis, linear_axes = kperp_linear_axis, $
              window_num = window_num
          endif else begin
          
            uvf_slice_plot, uv_slice_savefile[i], multi_pos = pos_use, start_multi_params = start_multi_params, png = png, eps = eps, pdf = pdf, $
              plotfile = uv_slice_plotfile, type = uvf_plot_type, title_prefix = slice_titles[i], note = note_use, $
              hinv = hinv, data_range = slice_range, log = uvf_log, $
              baseline_axis = baseline_axis, window_num = window_num
          endelse
          
          if i eq 0 then begin
            positions = pos_use
            undefine, start_multi_params
          endif
        endfor
        undefine, pos_use
        if keyword_set(pub) then begin
          cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density = 600
        endif
        
      endelse
    endif
    
    ave_power_vals = fltarr(n_cubes)
    wt_ave_power_vals = fltarr(n_cubes)
    ave_weights_vals = fltarr(n_cubes)
    ave_weights_freq_vals = fltarr(n_cubes, n_freq)
    uv_pix_area = fltarr(n_cubes)
    uv_area = fltarr(n_cubes)
    nbsl_lambda2 = fltarr(n_cubes)
    nbsl_lambda2_freq = fltarr(n_cubes, n_freq)
    wt_ave_power_freq_vals = fltarr(n_cubes, n_freq)
    ave_power_freq_vals = fltarr(n_cubes, n_freq)
    ave_power_uvf_vals = fltarr(n_cubes)
    wt_ave_power_uvf_vals = fltarr(n_cubes)
    
    for k=0, n_cubes-1 do begin
    
      void = getvar_savefile(savefiles_1d[k,0,0], names=varnames)
      
      if max(strmatch(varnames, 'ave_power', /fold_case)) gt 0 then $
        ave_power_vals[k] = getvar_savefile(savefiles_1d[k,0,0], 'ave_power') else ave_power_vals[k] = -1
      if max(strmatch(varnames, 'wt_ave_power', /fold_case)) gt 0 then $
        wt_ave_power_vals[k] = getvar_savefile(savefiles_1d[k,0,0], 'wt_ave_power') else wt_ave_power_vals[k] = -1
      if max(strmatch(varnames, 'uv_pix_area', /fold_case)) gt 0 then $
        uv_pix_area[k] = getvar_savefile(savefiles_1d[k,0,0], 'uv_pix_area') else uv_pix_area[k] = -1
      if max(strmatch(varnames, 'uv_area', /fold_case)) gt 0 then begin
        uv_area[k] = getvar_savefile(savefiles_1d[k,0,0], 'uv_area')
        
        nbsl_lambda2[k] = file_struct_arr[0].n_vis[0]/uv_area[k]
        if n_elements(freq_ch_range) gt 0 then $
          nbsl_lambda2_freq[k,*] = file_struct_arr[0].n_vis_freq[0,freq_ch_range[0]:freq_ch_range[1]]/uv_area[k] $
        else nbsl_lambda2_freq[k,*] = file_struct_arr[0].n_vis_freq[0,*]/uv_area[k]
      endif else begin
        uv_area[k] = -1
        nbsl_lambda2[k] = -1
        nbsl_lambda2_freq[k,*] = -1
      endelse
      
      if max(strmatch(varnames, 'ave_weights', /fold_case)) gt 0 and uv_pix_area[k] gt -1 then begin
        ave_weights_vals[k] = mean(getvar_savefile(savefiles_1d[k,0,0], 'ave_weights')/uv_pix_area[k])
        ave_weights_freq_vals[k,*] = (getvar_savefile(savefiles_1d[k,0,0], 'ave_weights')/uv_pix_area[k])[0,*]
      endif else begin
        ave_weights_vals[k] = -1
        ave_weights_freq_vals[k,*] = -1
      endelse
      
      if max(strmatch(varnames, 'wt_ave_power_freq', /fold_case)) gt 0 then $
        wt_ave_power_freq_vals[k,*] = (getvar_savefile(savefiles_1d[k,0,0], 'wt_ave_power_freq'))[0,*] $
      else wt_ave_power_freq_vals[k,*] = -1
      if max(strmatch(varnames, 'ave_power_freq', /fold_case)) gt 0 then $
        ave_power_freq_vals[k,*] = (getvar_savefile(savefiles_1d[k,0,0], 'ave_power_freq'))[0,*] $
      else ave_power_freq_vals[k,*] = -1
      if max(strmatch(varnames, 'ave_power_uvf', /fold_case)) gt 0 then $
        ave_power_uvf_vals[k] = (getvar_savefile(savefiles_1d[k,0,0], 'ave_power_uvf'))[0,*] $
      else ave_power_uvf_vals[k] = -1
      if max(strmatch(varnames, 'wt_ave_power_uvf', /fold_case)) gt 0 then $
        wt_ave_power_uvf_vals[k] = (getvar_savefile(savefiles_1d[k,0,0], 'wt_ave_power_uvf'))[0,*] $
      else wt_ave_power_uvf_vals[k] = -1
    endfor
    cube_power_info = {ave_power:ave_power_vals, wt_ave_power:wt_ave_power_vals, $
      uv_pix_area:uv_pix_area, uv_area:uv_area, $
      ave_weights:ave_weights_vals, ave_weights_freq:ave_weights_freq_vals, wt_ave_power_freq:wt_ave_power_freq_vals, $
      ave_power_freq:ave_power_freq_vals, wt_ave_power_uvf:wt_ave_power_uvf_vals, ave_power_uvf:ave_power_uvf_vals, $
      nbsl_lambda2:nbsl_lambda2, nbsl_lambda2_freq:nbsl_lambda2_freq}
      
    if keyword_set(plot_eor_1d) or keyword_set(plot_flat_1d) then begin
      case strlowcase(!version.os_family) OF
        'windows': split_delim = ';'
        'unix':    split_delim = ':'
      endcase
      path_dirs = strsplit(!path, split_delim, /extract)
      
      fhd_catalog_loc = strpos(path_dirs, 'catalog_data')
      wh_catalog = where(fhd_catalog_loc gt 0, count_catalog)
      if count_catalog gt 0 then begin
        file_path = path_dirs[wh_catalog[0]]
        ;; make sure file_path has a path separator at the end
        pos = strpos(file_path, path_sep(), /reverse_search)
        if pos+1-strlen(file_path) lt 0 then file_path = file_path + path_sep()
        
        eor_file_1d = file_path + 'eor_power_1d.idlsave'
        flat_file_1d = file_path + 'flat_power_1d.idlsave'
        
        flat_power = mean(getvar_savefile(flat_file_1d, 'power'))
        cube_power_info = create_struct(cube_power_info, 'flat_power', flat_power)
        
      endif else print, 'Could not locate catalog_data directory in !path variable'
    endif
    
    if keyword_set(plot_stdset) then begin
      for i=0, n_elements(wedge_1dbin_names)-1 do begin
        file_arr = reform(savefiles_1d[*,*,i], n_cubes*n_elements(kperp_density_names))
        
        if n_elements(plotfile_1d) gt 0 then begin
          plotfile_use = plotfile_1d[i]
          plotfile_noise_use = plotfile_1d_noise[i]
        endif
        
        if i gt 0 then note_use = note + ' ' + strjoin(strsplit(wedge_1dbin_names[i], '_', /extract), ' ') else note_use = note
        
        titles_use = strarr(n_cubes, n_elements(kperp_density_names))
        for j=0, n_elements(kperp_density_names)-1 do titles_use[*,j] = titles + kperp_density_names[j]
        titles_use = reform(titles_use, n_cubes*n_elements(kperp_density_names))
        
        k_range = minmax([kperp_plot_range, kpar_bin, kpar_plot_range[1]])
        
        if keyword_set(plot_noise_1d) then begin
          window_num = window_num+1
          kpower_1d_plots, file_arr, window_num = window_num, names = titles_use, delta = plot_1d_delta, hinv = hinv, $
            png = png, eps = eps, pdf = pdf, plotfile = plotfile_noise_use, k_range = k_range, title = note_use + ' Ob. Noise', $
            note = note_1d, no_text = no_text_1d, data_range = range_1d, $
            plot_error_bars = plot_1d_error_bars, plot_nsigma = plot_1d_nsigma, plot_sim_noise = plot_sim_noise, /plot_noise
        endif
        
        if keyword_set(plot_eor_1d) then begin
          if count_catalog gt 0 then begin
            psyms = [intarr(n_elements(file_arr))+10, -3]
            file_arr = [file_arr, eor_file_1d]
            titles_use = [titles_use, 'EoR signal']
          endif
        endif
        if keyword_set(plot_flat_1d) then begin
          if count_catalog gt 0 then begin
            psyms = [intarr(n_elements(file_arr))+10, -3]
            file_arr = [file_arr, flat_file_1d]
            titles_use = [titles_use, 'input flat power']
          endif
        endif
        
        window_num = window_num+1
        kpower_1d_plots, file_arr, window_num = window_num, colors = colors, names = titles_use, psyms = psyms, delta = plot_1d_delta, $
          hinv = hinv, png = png, eps = eps, pdf = pdf, plotfile = plotfile_use, k_range = k_range, title = note_use, note = note_1d, $
          no_text = no_text_1d, data_range = range_1d, $
          plot_error_bars = plot_1d_error_bars, plot_nsigma = plot_1d_nsigma, plot_sim_noise = plot_sim_noise
          
      endfor
    endif
    
    if keyword_set(plot_kpar_power) then begin
    
      file_arr = reform(savefiles_kpar_1d, n_cubes*n_elements(kperp_density_names))
      
      titles_use = strarr(n_cubes, n_elements(kperp_density_names))
      for j=0, n_elements(kperp_density_names)-1 do titles_use[*,j] = titles + kperp_density_names[j]
      titles_use = reform(titles_use, n_cubes*n_elements(kperp_density_names))
      if keyword_set(plot_eor_1d) then begin
        if count_catalog gt 0 then begin
          psyms = [intarr(n_elements(file_arr))+10, -3]
          file_arr = [file_arr, eor_file_1d]
          titles_use = [titles_use, 'EoR signal']
        endif
      endif
      if keyword_set(plot_flat_1d) then begin
        if count_catalog gt 0 then begin
          psyms = [intarr(n_elements(file_arr))+10, -3]
          file_arr = [file_arr, flat_file_1d]
          titles_use = [titles_use, 'input flat power']
        endif
      endif
      
      k_range = minmax([kperp_plot_range, kpar_bin, kpar_plot_range[1]])
      
      window_num = window_num+1
      kpower_1d_plots, file_arr, window_num = window_num, colors = colors, names = titles_use, psyms = psyms, delta = plot_1d_delta, hinv = hinv, $
        png = png, eps = eps, pdf = pdf, plotfile = plotfile_kpar_power, k_range = k_range, title = note + ' kpar', note = note_kpar_1d, $
        no_text = no_text_1d, plot_error_bars = plot_1d_error_bars, plot_nsigma = plot_1d_nsigma, plot_sim_noise = plot_sim_noise, $
        data_range = range_1d, /kpar_power, delay_axis = delay_axis, cable_length_axis = cable_length_axis
        
    endif
    
    if keyword_set(plot_kperp_power) then begin
    
      file_arr = reform(savefiles_kperp_1d, n_cubes*n_elements(kperp_density_names))
      
      titles_use = strarr(n_cubes, n_elements(kperp_density_names))
      for j=0, n_elements(kperp_density_names)-1 do titles_use[*,j] = titles + kperp_density_names[j]
      titles_use = reform(titles_use, n_cubes*n_elements(kperp_density_names))
      
      if keyword_set(plot_eor_1d) then begin
        if count_catalog gt 0 then begin
          psyms = [intarr(n_elements(file_arr))+10, -3]
          file_arr = [file_arr, eor_file_1d]
          titles_use = [titles_use, 'EoR signal']
        endif else print, 'Could not locate catalog_data directory in !path variable'
      endif
      if keyword_set(plot_flat_1d) then begin
        if count_catalog gt 0 then begin
          psyms = [intarr(n_elements(file_arr))+10, -3]
          file_arr = [file_arr, flat_file_1d]
          titles_use = [titles_use, 'input flat power']
        endif else print, 'Could not locate catalog_data directory in !path variable'
      endif
      
      k_range = minmax([kperp_plot_range, kpar_bin, kpar_plot_range[1]])
      
      window_num = window_num+1
      kpower_1d_plots, file_arr, window_num = window_num, colors = colors, names = titles_use, psyms = psyms, delta = plot_1d_delta, $
        hinv = hinv, png = png, eps = eps, pdf = pdf, plotfile = plotfile_kperp_power, k_range = k_range, title = note + ' kperp', $
        note = note_kperp_1d, no_text = no_text_1d, $
        plot_error_bars = plot_1d_error_bars, plot_nsigma = plot_1d_nsigma, plot_sim_noise = plot_sim_noise, $
        data_range = range_1d, /kperp_power, baseline_axis = baseline_axis
        
    endif
    
    if keyword_set(plot_k0_power) then begin
    
      file_arr = savefiles_k0
      
      titles_use = strarr(n_cubes, n_elements(kperp_density_names))
      for j=0, n_elements(kperp_density_names)-1 do titles_use[*,j] = titles + kperp_density_names[j]
      titles_use = reform(titles_use, n_cubes*n_elements(kperp_density_names))
      
      k_range = minmax([kperp_plot_range, kperp_bin])
      
      window_num = window_num+1
      kpower_1d_plots, file_arr, window_num = window_num, colors = colors, names = titles_use, delta = plot_1d_delta, hinv = hinv, $
        png = png, eps = eps, pdf = pdf, plotfile = plotfile_k0_power, k_range = k_range, title = note + ' kpar=0', note = note_1d, $
        no_text = no_text_1d, plot_error_bars = plot_1d_error_bars, plot_nsigma = plot_1d_nsigma, plot_sim_noise = plot_sim_noise, $
        data_range = range_1d, /kperp_power, baseline_axis = baseline_axis
        
    endif
    
    if keyword_set(plot_1to2d) then begin
    
    
      if keyword_set(kperp_linear_axis) then begin
        ;; aspect ratio doesn't work out for kperp_linear with multiple rows
        ncol = n_elements(wedge_1dbin_names)*n_elements(kperp_density_names)
        nrow = 1
      endif else begin
        ncol = n_elements(wedge_1dbin_names)
        nrow = n_elements(kperp_density_names)
      endelse
      start_multi_params = {ncol:ncol, nrow:nrow, ordering:'col'}
      undefine, pos_use
      
      window_num = window_num+1
      for i=0, nrow*ncol-1 do begin
        if i gt 0 then pos_use = positions[*,i]
        if i eq nrow*ncol-1 and n_elements(note) gt 0 then note_use = note + ',' + type_1to2d_use + '_' + pol_1to2d_use else undefine, note_use
        if keyword_set(pub) then plotfile_use = plotfile_1to2d_heatmap else undefine, plotfile_use
        
        wedge_index = i / n_elements(kperp_density_names)
        density_index = i mod n_elements(kperp_density_names)
        title_use = wedge_1dbin_names[wedge_index] + ' ' + kperp_density_names[density_index]
        
        kpower_2d_plots, savefiles_1to2d_mask[wh_2d_use, density_index, wedge_index], multi_pos = pos_use, start_multi_params = start_multi_params, png = png, eps = eps, pdf = pdf, $
          plotfile = plotfile_use, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
          /plot_mask, /mask_contour, contour_levels = 'all', $
          title_prefix = title_use, note = note_use, $
          plot_wedge_line = plot_wedge_line, hinv = hinv, $
          wedge_amp = wedge_amp, baseline_axis = baseline_axis, delay_axis = delay_axis, cable_length_axis = cable_length_axis, $
          kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, window_num = window_num
          
        if i eq 0 then begin
          positions = pos_use
          undefine, start_multi_params
        endif
      endfor
      undefine, positions, pos_use
      if keyword_set(pub) then begin
        cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density = 600
      endif
      
      
      ;; decide whether to make a zoom plot
      if n_elements(kperp_range_1dave) gt 0 or n_elements(kperp_range_lambda_1dave) gt 0 or n_elements(kpar_range_1dave) gt 0 then begin
        if n_elements(kperp_range_1dave) gt 0 then begin
          if kperp_range_1dave[0] gt kperp_plot_range[0] or kperp_range_1dave[1] lt kperp_plot_range[1] then begin
            zoom = 1
            kperp_range_use = kperp_range_1dave * [.7, 1.1]
            kperp_range_use = [max([kperp_plot_range[0], kperp_range_use[0]]), min([kperp_plot_range[1], kperp_range_use[1]])]
          endif else kperp_range_use = kperp_plot_range
        endif
        if n_elements(kperp_range_lambda_1dave) gt 0 then begin
          if keyword_set(hinv) then kperp_range_use = kperp_range_lambda_1dave / (kperp_lambda_conv * hubble_param) else $
            kperp_range_use = kperp_range_lambda_1dave / kperp_lambda_conv
            
          if kperp_range_use[0] gt kperp_plot_range[0] or kperp_range_use[1] lt kperp_plot_range[1] then begin
            zoom = 1
            kperp_range_use = kperp_range_use * [.7, 1.1]
            kperp_range_use = [max([kperp_plot_range[0], kperp_range_use[0]]), min([kperp_plot_range[1], kperp_range_use[1]])]
          endif else kperp_range_use = kperp_plot_range
        endif
        if n_elements(kpar_range_1dave) gt 0 then begin
          if kpar_range_1dave[0] gt kpar_plot_range[0] or kpar_range_1dave[1] lt kpar_plot_range[1] then begin
            zoom = 1
            kpar_range_use = kpar_range_1dave * [.7, 1.1]
            kpar_range_use = [max([kpar_plot_range[0], kpar_range_use[0]]), min([kpar_plot_range[1], kpar_range_use[1]])]
          endif else kpar_range_use = kpar_plot_range
        endif
      endif
      
      start_multi_params = {ncol:ncol, nrow:nrow, ordering:'col'}
      
      window_num = window_num+1
      for i=0, nrow*ncol-1 do begin
        if i gt 0 then pos_use = positions[*,i]
        if i eq nrow*ncol-1 and n_elements(note) gt 0 then note_use = note + ',' + type_1to2d_use + '_' + pol_1to2d_use else undefine, note_use
        if keyword_set(pub) then plotfile_use = plotfile_1to2d_contours else undefine, plotfile_use
        
        wedge_index = i / n_elements(kperp_density_names)
        density_index = i mod n_elements(kperp_density_names)
        title_use = wedge_1dbin_names[wedge_index] + ' ' + kperp_density_names[density_index]
        
        if keyword_set(zoom) then levels = 1 else levels = 'all'
        ;levels = 'all'
        
        kpower_2d_plots, savefiles_2d[wh_2d_use, density_index], multi_pos = pos_use, start_multi_params = start_multi_params, png = png, eps = eps, pdf = pdf, $
          plotfile = plotfile_use, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
          mask_savefile = savefiles_1to2d_mask[wh_2d_use, density_index, wedge_index], /mask_contour, contour_levels = levels, $
          data_range = data_range, title_prefix = title_use, note = note_use, $
          plot_wedge_line = 0, hinv = hinv, $
          wedge_amp = wedge_amp, baseline_axis = baseline_axis, delay_axis = delay_axis, cable_length_axis = cable_length_axis, $
          kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, window_num = window_num
          
        if i eq 0 then begin
          positions = pos_use
          undefine, start_multi_params
        endif
      endfor
      undefine, positions, pos_use
      if keyword_set(pub) then begin
        cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density = 600
      endif
      
      start_multi_params = {ncol:ncol, nrow:nrow, ordering:'col'}
      
      window_num = window_num+1
      for i=0, nrow*ncol-1 do begin
        if i gt 0 then pos_use = positions[*,i]
        if i eq nrow*ncol-1 and n_elements(note) gt 0 then note_use = note + ',' + type_1to2d_use + '_' + pol_1to2d_use else undefine, note_use
        if keyword_set(pub) then plotfile_use = plotfile_1to2d_noisefrac else undefine, plotfile_use
        
        wedge_index = i / n_elements(kperp_density_names)
        density_index = i mod n_elements(kperp_density_names)
        title_use = wedge_1dbin_names[wedge_index] + ' ' + kperp_density_names[density_index]
        
        levels = 'all'
        
        kpower_2d_plots, savefiles_1to2d_mask[wh_2d_use, density_index, wedge_index], multi_pos = pos_use, start_multi_params = start_multi_params, png = png, eps = eps, pdf = pdf, $
          plotfile = plotfile_use, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
          mask_savefile = savefiles_1to2d_mask[wh_2d_use, density_index, wedge_index], /mask_contour, contour_levels = levels, /plot_1d_noisefrac, $
          data_range = [0,1], title_prefix = title_use, note = note_use, $
          plot_wedge_line = 0, hinv = hinv, $
          wedge_amp = wedge_amp, baseline_axis = baseline_axis, delay_axis = delay_axis, cable_length_axis = cable_length_axis, $
          kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, window_num = window_num
          
        if i eq 0 then begin
          positions = pos_use
          undefine, start_multi_params
        endif
      endfor
      undefine, positions, pos_use
      if keyword_set(pub) then begin
        cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density = 600
      endif
      
      if keyword_set(zoom) then begin
        start_multi_params = {ncol:ncol, nrow:nrow, ordering:'col'}
        
        window_num = window_num+1
        for i=0, nrow*ncol-1 do begin
          if i gt 0 then pos_use = positions[*,i]
          if i eq nrow*ncol-1 and n_elements(note) gt 0 then note_use = note + ',' + type_1to2d_use + '_' + pol_1to2d_use else undefine, note_use
          if keyword_set(pub) then plotfile_use = plotfile_1to2d_contour_zoom else undefine, plotfile_use
          
          wedge_index = i / n_elements(kperp_density_names)
          density_index = i mod n_elements(kperp_density_names)
          title_use = wedge_1dbin_names[wedge_index] + ' ' + kperp_density_names[density_index]
          
          kpower_2d_plots, savefiles_2d[wh_2d_use, density_index], multi_pos = pos_use, start_multi_params = start_multi_params, png = png, eps = eps, pdf = pdf, $
            plotfile = plotfile_use, kperp_plot_range = kperp_range_use, kpar_plot_range = kpar_range_use, $
            mask_savefile = savefiles_1to2d_mask[wh_2d_use, density_index, wedge_index], /mask_contour, contour_levels = 'all', $
            data_range = data_range, title_prefix = title_use, note = note_use, $
            plot_wedge_line = 0, hinv = hinv, $
            wedge_amp = wedge_amp, baseline_axis = baseline_axis, delay_axis = delay_axis, cable_length_axis = cable_length_axis, $
            kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, window_num = window_num
            
          if i eq 0 then begin
            positions = pos_use
            undefine, start_multi_params
          endif
        endfor
        undefine, positions, pos_use
        if keyword_set(pub) then begin
          cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density = 600
        endif
      endif
    endif
  end
