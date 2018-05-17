pro ps_main_plots, datafile, beamfiles = beamfiles, pol_inc = pol_inc, $
    type_inc = type_inc, freq_ch_range = freq_ch_range, freq_flags = freq_flags, $
    freq_flag_name = freq_flag_name, uvf_input = uvf_input, no_evenodd = no_evenodd, $
    rts = rts, norm_rts_with_fhd = norm_rts_with_fhd, casa = casa, sim = sim, $
    fix_sim_input = fix_sim_input, save_path = save_path, $
    savefilebase = savefilebase, cube_power_info = cube_power_info, $
    refresh_options = refresh_options, uvf_options = uvf_options, $
    ps_options = ps_options, binning_2d_options = binning_2d_options, $
    binning_1d_options = binning_1d_options, plot_options = plot_options, $
    plot_types = plot_types, plot_2d_options = plot_2d_options, $
    plot_1d_options = plot_1d_options

  if binning_2d_options.no_kzero and plot_types.plot_k0_power then begin
    message, 'plot_k0_power cannot be set if no_kzero keyword is set.'
  endif

  if n_elements(savefilebase) gt 1 then message, 'savefilebase must be a scalar'

  time0 = systime(1)
  if keyword_set(rts) then begin
    file_struct_arr = rts_file_setup(datafile, savefilebase = savefilebase, $
      save_path = save_path, use_fhd_norm = norm_rts_with_fhd, $
      freq_ch_range = freq_ch_range, freq_flags = freq_flags, freq_flag_name = freq_flag_name, $
      refresh_info = refresh_options.refresh_info, uvf_options = uvf_options, $
      ps_options = ps_options)
  endif else if keyword_set(casa) then begin
    file_struct_arr = casa_file_setup(datafile, savefilebase = savefilebase, $
      save_path = save_path, $
      freq_ch_range = freq_ch_range, freq_flags = freq_flags, freq_flag_name = freq_flag_name, $
      refresh_info = refresh_options.refresh_info, uvf_options = uvf_options, $
      ps_options = ps_options)
  endif else begin
    file_struct_arr = fhd_file_setup(datafile, beamfile = beamfiles, sim = sim, $
      uvf_input = uvf_input, savefilebase = savefilebase, save_path = save_path, $
      freq_ch_range = freq_ch_range, freq_flags = freq_flags, freq_flag_name = freq_flag_name, $
      refresh_info = refresh_options.refresh_info, uvf_options = uvf_options, $
      ps_options = ps_options)
  endelse
  time1 = systime(1)
  print, 'file setup time: ' + number_formatter(time1-time0)

  ;; make sure save paths exist
  save_paths = file_struct_arr[0].savefile_froot + file_struct_arr[0].subfolders.data + $
    ['', file_struct_arr[0].subfolders.uvf + ['', file_struct_arr[0].subfolders.slices], $
     file_struct_arr[0].subfolders.kspace + ['', file_struct_arr[0].subfolders.slices], $
     file_struct_arr[0].subfolders.beams, $
     file_struct_arr[0].subfolders.bin_2d + ['', 'from_1d/'], $
     file_struct_arr[0].subfolders.bin_1d]
  for i = 0, n_elements(save_paths) - 1 do begin
    if not file_test(save_paths[i], /directory) then file_mkdir, save_paths[i]
  endfor

  if plot_options.pub then begin
    if tag_exist(plot_options, 'plot_path') ne 0 then begin
      plotfile_path = plot_options.plot_path
    endif else begin
      plotfile_path = file_struct_arr[0].savefile_froot + file_struct_arr[0].subfolders.plots
    endelse

    ;; make sure plot paths exist
    plot_paths = plotfile_path + ['', file_struct_arr[0].subfolders.slices, $
       file_struct_arr[0].subfolders.bin_2d + ['', 'from_1d/'], $
       file_struct_arr[0].subfolders.bin_1d + ['', 'bin_histograms/']]
    for i = 0, n_elements(plot_paths) - 1 do begin
      if not file_test(plot_paths[i], /directory) then file_mkdir, plot_paths[i]
    endfor
  endif

  if not tag_exist(file_struct_arr, 'beam_int') and refresh_options.refresh_ps then begin
    refresh_options.refresh_beam = 1
  endif

  nfiles = file_struct_arr[0].nfiles
  if nfiles lt 2 and not keyword_set(no_evenodd) then begin
    message, 'Even and odd cubes are not all present. If this is expected, set no_evenodd=1'
  endif

  if nfiles gt 1 then begin
    print, 'n_vis difference between even & odd cubes: ' + $
      number_formatter((file_struct_arr[0].n_vis[1]-file_struct_arr[0].n_vis[0]))
    print, 'n_vis % difference between even & odd cubes: ' + $
      number_formatter((file_struct_arr[0].n_vis[1]-file_struct_arr[0].n_vis[0])*100 / $
        mean(file_struct_arr[0].n_vis))
    print, 'n_vis_freq difference between even & odd cubes: ' + $
      number_formatter(total(abs(file_struct_arr[0].n_vis_freq[0, *] - file_struct_arr[0].n_vis_freq[1,*])))
stop
    if total(abs(total(file_struct_arr[0].n_vis_freq, 2) - file_struct_arr[0].n_vis)) ne 0 then begin
      print, 'number of visibilities in n_vis and nf_vis do not match!'
      print, 'nf_vis difference between even & odd cubes: ' + $
        number_formatter(total(abs(file_struct_arr[0].n_vis_freq[1, *]-file_struct_arr[0].n_vis_freq[0, *])))
      print, 'nf_vis % difference between even & odd cubes: ' + $
        number_formatter(total(abs(file_struct_arr[0].n_vis_freq[1, *]-file_struct_arr[0].n_vis_freq[0, *]))*100 / $
          mean(total(file_struct_arr[0].n_vis_freq, 2)))
    endif
  endif

  if tag_exist(file_struct_arr, 'n_obs') then begin
    print, 'n_obs: ', file_struct_arr[0].n_obs
    if tag_exist(plot_options, 'note') then begin
      plot_options.note += '(' + number_formatter(round(mean(file_struct_arr[0].n_obs))) + ')'
    endif else begin
      plot_options = create_plot_options(plot_options = plot_options, $
        note = ' (' + number_formatter(round(mean(file_struct_arr[0].n_obs))) + ')')
    endelse
  endif

  if plot_2d_options.plot_wedge_line then begin
    z0_freq = 1420.40 ;; MHz
    redshifts = z0_freq/file_struct_arr[0].frequencies - 1
    mean_redshift = mean(redshifts)

    cosmology_measures, mean_redshift, wedge_factor = wedge_factor, hubble_param = hubble_param
    ;; assume 20 degrees from pointing center to first null
    source_dist = 20d * !dpi / 180d
    fov_amp = wedge_factor * source_dist

    ;; calculate angular distance to horizon
    horizon_amp = wedge_factor * ((file_struct_arr[0].max_theta+90d) * !dpi / 180d)

    wedge_amps = [fov_amp, horizon_amp]
    wedge_names = ['fov', 'horizon']
    if tag_exist(binning_1d_options, 'wedge_angles') gt 0 then begin
      if min(binning_1d_options.wedge_angles) le 0 or $
          max(binning_1d_options.wedge_angles) ge 180 then begin
        message, 'wedge_angles must be in degrees and between 0 & 180'
      endif
      wedge_amps = [wedge_factor * (binning_1d_options.wedge_angles*!dpi / 180d)]
      wedge_names = [number_formatter(binning_1d_options.wedge_angles) + 'deg']
    endif
    binning_1d_options = create_binning_1d_options(binning_1d_options = binning_1d_options, $
      wedge_amps = wedge_amps, wedge_names = wedge_names)

    if tag_exist(binning_1d_options, 'coarse_harm_width') then begin
      harm_freq = 1.28
      if n_elements(freq_ch_range) gt 0 then begin
        freqs_use = file_struct_arr[0].frequencies[freq_ch_range[0]:freq_ch_range[1]]
      endif else begin
        freqs_use = file_struct_arr[0].frequencies
      endelse

      bandwidth = max(freqs_use) - min(freqs_use) + freqs_use[1] - freqs_use[0]
      coarse_harm0 = round(bandwidth / harm_freq)
      binning_1d_options = create_binning_1d_options(binning_1d_options = binning_1d_options, $
        coarse_harm0 = coarse_harm0)
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
  if n_cubes ne ntype * npol then begin
    message, 'number of cubes does not match expected value'
  endif

  file_labels = file_struct_arr.file_label
  wt_file_labels = file_struct_arr.wt_file_label
  titles = strarr(n_cubes)
  for i=0, n_cubes-1 do begin
    titles[i] = strjoin(strsplit(file_labels[i], '_', /extract), ' ')
  endfor

  weight_ind = file_struct_arr.pol_index
  weight_labels = strupcase(file_struct_arr.pol)

  power_tag = file_struct_arr[0].power_tag
  if tag_exist(file_struct_arr[0], 'uvf_tag') then begin
    uvf_tag = file_struct_arr[0].uvf_tag
  endif else begin
    uvf_tag = ''
  endelse

  fadd_2dbin = ''
  if binning_2d_options.no_kzero then fadd_2dbin = fadd_2dbin + '_nok0'
  if binning_2d_options.log_kpar then fadd_2dbin = fadd_2dbin + '_logkpar'
  if binning_2d_options.log_kperp then fadd_2dbin = fadd_2dbin + '_logkperp'

  fadd_1dbin = ''
  if binning_1d_options.log_k then fadd_1dbin = fadd_1dbin + '_logk'

  fadd_kpar_1d = fadd_1dbin
  fadd_kperp_1d = fadd_1dbin

  fadd_1dmask = ''
  if tag_exist(binning_1d_options, 'kperp_range_1dave') then begin
    fadd_1dmask = fadd_1dmask + '_kperp' + $
      number_formatter(binning_1d_options.kperp_range_1dave[0]) + '-' + $
      number_formatter(binning_1d_options.kperp_range_1dave[1])
    note_1d = 'kperp: [' + $
      number_formatter(binning_1d_options.kperp_range_1dave[0]) + ',' + $
      number_formatter(binning_1d_options.kperp_range_1dave[1]) + ']'

    ;; if we're plotting in [k]=h/Mpc then assume these numbers are in h/Mpc.
    ;; Convert to 1/Mpc for internal code usage
    if plot_options.hinv then begin
      kperp_range_1d_use = binning_1d_options.kperp_range_1dave * hubble_param
    endif else begin
      kperp_range_1d_use = binning_1d_options.kperp_range_1dave
    endelse
  endif else begin
    if tag_exist(binning_1d_options, 'kperp_range_lambda_1dave') then begin
      fadd_1dmask = fadd_1dmask + '_kperplambda' + $
        number_formatter(binning_1d_options.kperp_range_lambda_1dave[0]) + '-' + $
        number_formatter(binning_1d_options.kperp_range_lambda_1dave[1])
      note_1d = 'kperp: [' + $
        number_formatter(binning_1d_options.kperp_range_lambda_1dave[0]) + ',' + $
        number_formatter(binning_1d_options.kperp_range_lambda_1dave[1]) + ']'
    endif else begin
        ;; if no range set default to same range as is used in 2D plots
        kperp_range_lambda_1dave = [5., min([file_struct_arr.kspan/2., $
          file_struct_arr.max_baseline_lambda])]
        binning_1d_options = create_binning_1d_options(binning_1d_options = binning_1d_options, $
          kperp_range_lambda_1dave = kperp_range_lambda_1dave)
    endelse
  endelse

  if tag_exist(binning_1d_options, 'kx_range_1dave') then begin
    fadd_1dmask = fadd_1dmask + '_kx' + $
      number_formatter(binning_1d_options.kx_range_1dave[0]) + '-' + $
      number_formatter(binning_1d_options.kx_range_1dave[1])
    note_1d = 'kx: [' + $
      number_formatter(binning_1d_options.kx_range_1dave[0]) + ',' + $
      number_formatter(binning_1d_options.kx_range_1dave[1]) + ']'

    ;; if we're plotting in [k]=h/Mpc then assume these numbers are in h/Mpc.
    ;; Convert to 1/Mpc for internal code usage
    if plot_options.hinv then begin
      kx_range_1d_use = binning_1d_options.kx_range_1dave * hubble_param
    endif else begin
      kx_range_1d_use = binning_1d_options.kx_range_1dave
    endelse
  endif else begin
    if tag_exist(binning_1d_options, 'kx_range_lambda_1dave') then begin
      fadd_1dmask = fadd_1dmask + '_kxlambda' + $
        number_formatter(binning_1d_options.kx_range_lambda_1dave[0]) + '-' + $
        number_formatter(binning_1d_options.kx_range_lambda_1dave[1])
      note_1d = 'kx: [' + $
        number_formatter(binning_1d_options.kx_range_lambda_1dave[0]) + ',' + $
        number_formatter(binning_1d_options.kx_range_lambda_1dave[1]) + ']'
    endif
  endelse

  if tag_exist(binning_1d_options, 'ky_range_1dave') then begin
    fadd_1dmask = fadd_1dmask + '_ky' + $
      number_formatter(binning_1d_options.ky_range_1dave[0]) + '-' + $
      number_formatter(binning_1d_options.ky_range_1dave[1])
    note_1d = 'ky: [' + $
      number_formatter(binning_1d_options.ky_range_1dave[0]) + ',' + $
      number_formatter(binning_1d_options.ky_range_1dave[1]) + ']'

    ;; if we're plotting in [k]=h/Mpc then assume these numbers are in h/Mpc.
    ;; Convert to 1/Mpc for internal code usage
    if plot_options.hinv then begin
      ky_range_1d_use = binning_1d_options.ky_range_1dave * hubble_param
    endif else begin
      ky_range_1d_use = binning_1d_options.ky_range_1dave
    endelse

  endif else begin
    if tag_exist(binning_1d_options, 'ky_range_lambda_1dave') then begin
      fadd_1dmask = fadd_1dmask + '_kylambda' + $
        number_formatter(binning_1d_options.ky_range_lambda_1dave[0]) + '-' + $
        number_formatter(binning_1d_options.ky_range_lambda_1dave[1])
      note_1d = 'ky: [' + $
        number_formatter(binning_1d_options.ky_range_lambda_1dave[0]) + ',' + $
        number_formatter(binning_1d_options.ky_range_lambda_1dave[1]) + ']'
    endif
  endelse

  if tag_exist(binning_1d_options, 'kperp_range_lambda_kparpower') then begin
    fadd_kpar_1d = fadd_kpar_1d + '_kperplambda' + $
      number_formatter(binning_1d_options.kperp_range_lambda_kparpower[0]) + '-' + $
      number_formatter(binning_1d_options.kperp_range_lambda_kparpower[1])
    note_kpar_1d = 'kperp: [' + $
      number_formatter(binning_1d_options.kperp_range_lambda_kparpower[0]) + ',' + $
      number_formatter(binning_1d_options.kperp_range_lambda_kparpower[1]) + ']'
  endif else begin
    fadd_kpar_1d = fadd_1dbin + fadd_1dmask
    if n_elements(note_1d) gt 0 then note_kpar_1d = note_1d

    if tag_exist(binning_1d_options, 'kperp_range_lambda_1dave') then begin
      binning_1d_options = create_binning_1d_options(binning_1d_options = binning_1d_options, $
        kperp_range_lambda_kparpower = binning_1d_options.kperp_range_lambda_1dave)
    endif
  endelse

  if tag_exist(binning_1d_options, 'kpar_range_1dave') then begin
    fadd_1dmask = fadd_1dmask + '_kpar' + $
      number_formatter(binning_1d_options.kpar_range_1dave[0]) + '-' + $
      number_formatter(binning_1d_options.kpar_range_1dave[1])
    if n_elements(note_1d) eq 0 then begin
      note_1d=''
    endif else begin
      note_1d = note_1d + '; '
    endelse
    note_1d = note_1d + 'kpar: [' + $
      number_formatter(binning_1d_options.kpar_range_1dave[0]) + ',' + $
      number_formatter(binning_1d_options.kpar_range_1dave[1]) + ']'

    ;; if we're plotting in [k]=h/Mpc then assume these numbers are in h/Mpc.
    ;; Convert to 1/Mpc for internal code usage
    if plot_options.hinv then begin
      kpar_range_1d_use = binning_1d_options.kpar_range_1dave * hubble_param
    endif else begin
      kpar_range_1d_use = binning_1d_options.kpar_range_1dave
    endelse
  endif

  if not tag_exist(binning_1d_options, 'kpar_range_kperppower') and $
      tag_exist(binning_1d_options, 'kpar_range_1dave') then begin
      binning_1d_options = create_binning_1d_options(binning_1d_options = binning_1d_options, $
        kpar_range_kperppower = kpar_range_1dave)
  endif

  if tag_exist(binning_1d_options, 'kpar_range_kperppower') then begin
    fadd_kperp_1d = fadd_kperp_1d + '_kpar' + $
      number_formatter(binning_1d_options.kpar_range_kperppower[0]) + '-' + $
      number_formatter(binning_1d_options.kpar_range_kperppower[1])
    note_kperp_1d = 'kpar: [' + $
      number_formatter(binning_1d_options.kpar_range_kperppower[0]) + ',' + $
      number_formatter(binning_1d_options.kpar_range_kperppower[1]) + ']'

    ;; if we're plotting in [k]=h/Mpc then assume these numbers are in h/Mpc.
    ;; Convert to 1/Mpc for internal code usage
    if plot_options.hinv then begin
      kpar_range_kperppower_use = binning_1d_options.kpar_range_kperppower * hubble_param
    endif else begin
      kpar_range_kperppower_use = binning_1d_options.kpar_range_kperppower
    endelse
  endif
  fadd_1dbin = fadd_1dbin + fadd_1dmask

  if plot_2d_options.plot_wedge_line then begin
    if tag_exist(binning_1d_options, 'coarse_harm_width') then begin
      cb_width_name = '_cbw' + number_formatter(binning_1d_options.coarse_harm_width)
    endif else begin
      cb_width_name = ''
    endelse
    wedge_1dbin_names = ['', '_no_' + binning_1d_options.wedge_names + '_wedge' + cb_width_name]
  endif else wedge_1dbin_names = ''

  ;; density correction file naming for 2D & 1D files
  kperp_density_names = strarr(n_elements(ps_options.wt_cutoffs))
  wh_cutoff0 = where(ps_options.wt_cutoffs eq 0, count_cutoff0, complement = wh_cutoff_n0, $
    ncomplement = count_cutoff_n0)
  wh_std = where(ps_options.wt_cutoffs eq 1 and ps_options.wt_measures eq 'min', count_std)

  if count_cutoff0 gt 0 then kperp_density_names[wh_cutoff0] = '_nodensitycorr'
  if count_cutoff_n0 gt 0 then begin
    kperp_density_names[wh_cutoff_n0] = '_kperp_density_' + ps_options.wt_measures[wh_cutoff_n0] + $
      '_gt' + number_formatter(ps_options.wt_cutoffs[wh_cutoff_n0])
  endif

  if count_std gt 0 then kperp_density_names[wh_std] = '_dencorr'

  ;; need general_filebase for 1D plotfiles, make sure it doesn't have a full path
  general_filebase = file_struct_arr[0].general_filebase
  for i=0, n_cubes-1 do begin
    if file_struct_arr(i).general_filebase ne general_filebase then begin
      message, 'general_filebase does not match between 1d savefiles'
    endif
  endfor

  savefiles_2d = strarr(n_cubes, n_elements(kperp_density_names))
  bin_2d_path = file_struct_arr.savefile_froot + file_struct_arr[0].subfolders.data + $
    file_struct_arr[0].subfolders.bin_2d
  for j=0, n_elements(kperp_density_names)-1 do begin
    savefiles_2d[*,j] = bin_2d_path + file_struct_arr.savefilebase + $
      power_tag + fadd_2dbin + kperp_density_names[j] + '_2dkpower.idlsave'
  endfor
  test_save_2d = file_valid(savefiles_2d)

  savefiles_1d = strarr(n_cubes, n_elements(kperp_density_names), n_elements(wedge_1dbin_names))
  bin_1d_path = file_struct_arr.savefile_froot + file_struct_arr[0].subfolders.data + $
    file_struct_arr[0].subfolders.bin_1d
  savefiles_1to2d_bin = strarr(n_cubes, n_elements(kperp_density_names), $
    n_elements(wedge_1dbin_names))
  savefiles_2d_masked = strarr(n_cubes, n_elements(kperp_density_names), $
    n_elements(wedge_1dbin_names))
  for i=0, n_elements(wedge_1dbin_names)-1 do begin
    for j=0, n_elements(kperp_density_names)-1 do begin
      savefiles_1d[*,j,i] = bin_1d_path + file_struct_arr.savefilebase + $
        power_tag + kperp_density_names[j] + wedge_1dbin_names[i] + fadd_1dbin + $
        '_1dkpower.idlsave'
      savefiles_1to2d_bin[*,j,i] = bin_2d_path + 'from_1d/' + $
        file_struct_arr.savefilebase + power_tag + fadd_2dbin + kperp_density_names[j] + $
        wedge_1dbin_names[i] + fadd_1dbin + '_1to2d_bin.idlsave'
      savefiles_2d_masked[*,j, i] = bin_2d_path + 'from_1d/' + $
        file_struct_arr.savefilebase + power_tag + fadd_2dbin + kperp_density_names[j] + $
        wedge_1dbin_names[i] + fadd_1dmask + '_2dkpower_1dmask.idlsave'
    endfor
  endfor
  test_save_1d = file_valid(savefiles_1d)
  test_save_bin = file_valid(savefiles_1to2d_bin)
  test_save_mask = file_valid(savefiles_2d_masked)

  savefiles_kpar_1d = strarr(n_cubes, n_elements(kperp_density_names))
  for j=0, n_elements(kperp_density_names)-1 do begin
    savefiles_kpar_1d[*,j] = bin_1d_path + file_struct_arr.savefilebase + $
      power_tag + kperp_density_names[j] + fadd_kpar_1d + '_kpar_power.idlsave'
  endfor
  if plot_types.plot_kpar_power then test_save_kpar = file_valid(savefiles_kpar_1d)

  savefiles_kperp_1d = strarr(n_cubes, n_elements(kperp_density_names))
  for j=0, n_elements(kperp_density_names)-1 do begin
    savefiles_kperp_1d[*,j] = bin_1d_path + file_struct_arr.savefilebase + $
      power_tag + kperp_density_names[j] + fadd_kperp_1d + '_kperp_power.idlsave'
  endfor
  if plot_types.plot_kperp_power then test_save_kperp = file_valid(savefiles_kperp_1d)

  savefiles_k0 = strarr(n_cubes, n_elements(kperp_density_names))
  savefiles_k0_masked = strarr(n_cubes, n_elements(kperp_density_names), $
    n_elements(wedge_1dbin_names))
  for j=0, n_elements(kperp_density_names)-1 do begin
    savefiles_k0[*,j] = bin_1d_path + file_struct_arr.savefilebase + $
      power_tag + fadd_2dbin + kperp_density_names[j]+ '_k0power.idlsave'
    for i=0, n_elements(wedge_1dbin_names)-1 do begin
      savefiles_k0_masked[*,j,i] = bin_1d_path + $
        file_struct_arr.savefilebase + power_tag + fadd_2dbin + kperp_density_names[j] + $
        wedge_1dbin_names[i] + fadd_1dmask + '_k0power_1dmask.idlsave'
    endfor
  endfor
  if plot_types.plot_k0_power then test_save_k0 = file_valid(savefiles_k0)
  if plot_types.plot_k0_power then test_save_k0_masked = file_valid(savefiles_k0_masked)

  if refresh_options.refresh_binning then begin
    test_save_2d = test_save_2d*0
    test_save_1d = test_save_1d*0
    if plot_types.plot_kpar_power then test_save_kpar = test_save_kpar*0
    if plot_types.plot_kperp_power then test_save_kperp = test_save_kperp*0
    if plot_types.plot_k0_power then begin
      test_save_k0 = test_save_k0*0
      test_save_k0_masked = test_save_k0_masked*0
    endif
  endif

  if tag_exist(file_struct_arr[0], 'nside') ne 0 then healpix = 1 else healpix = 0

  if tag_exist(file_struct_arr, 'uvf_savefile') then uvf_input = 0 else uvf_input = 1

  n_freq = n_elements(file_struct_arr[0].frequencies)
  if n_elements(freq_ch_range) ne 0 then begin
    if max(freq_ch_range) gt n_freq-1 then message, 'invalid freq_ch_range'
    n_freq = freq_ch_range[1]-freq_ch_range[0]+1
  endif

  if healpix and tag_exist(uvf_options, 'dft_fchunk') then begin
    if tag_exist(uvf_options, 'wt_cutoffs') gt n_freq then begin
      print, 'dft_fchunk is larger than the number of frequency slices, setting ' + $
        'it to the number of slices -- ' + number_formatter(n_freq)
      uvf_options.dft_fchunk = n_freq
    endif
  endif

  if plot_types.plot_binning_hist and plot_options.pub eq 1 then begin
    bin_hist_plot_path = plotfile_path + file_struct_arr[0].subfolders.bin_1d + 'bin_histograms/'
    if not tag_exist(plot_options, 'plot_filebase') then begin
      plotfile_binning_hist = strarr(n_cubes, n_elements(kperp_density_names), $
        n_elements(wedge_1dbin_names))

      for i=0, n_elements(wedge_1dbin_names)-1 do begin
        for j=0, n_elements(kperp_density_names)-1 do begin
          plotfile_binning_hist[*,j,i] = bin_hist_plot_path + file_struct_arr.savefilebase + $
            power_tag + kperp_density_names[j] + wedge_1dbin_names[i] + fadd_1dbin
        endfor
      endfor
    endif else begin
      plotfile_binning_hist = strarr(n_cubes, n_elements(kperp_density_names), $
        n_elements(wedge_1dbin_names))

      for i=0, n_elements(wedge_1dbin_names)-1 do begin
        for j=0, n_elements(kperp_density_names)-1 do begin
          plotfile_binning_hist[*,j,i] = bin_hist_plot_path + plot_options.plot_filebase + $
            uvf_tag + file_struct_arr.file_label + power_tag + kperp_density_names[j] + $
            wedge_1dbin_names[i] + fadd_1dbin
        endfor
      endfor
    endelse
    plotfile_binning_hist = plotfile_binning_hist + '_binning_hist'
  endif

  for i=0, n_cubes-1 do begin

    savefile_2d_use = reform(savefiles_2d[i,*], n_elements(kperp_density_names))
    test_2d = min(test_save_2d[i,*])

    savefile_kpar_use = savefiles_kpar_1d[i,*]
    if plot_types.plot_kpar_power then begin
      test_kpar = min(test_save_kpar[i,*])
    endif

    savefile_kperp_use = savefiles_kperp_1d[i,*]
    if plot_types.plot_kperp_power then begin
      test_kperp = min(test_save_kperp[i,*])
    endif

    savefile_k0_use = savefiles_k0[i,*]
    savefile_k0_masked_use = reform(savefiles_k0_masked[i,*,*], $
      n_elements(kperp_density_names), n_elements(wedge_1dbin_names))
    if plot_types.plot_k0_power then begin
      test_k0 = min(test_save_k0[i,*])
      test_k0_masked = min(test_save_k0_masked[i,*,*])
    endif

    savefile_1d_use = reform(savefiles_1d[i,*,*], n_elements(kperp_density_names), $
      n_elements(wedge_1dbin_names))
    savefiles_bin_use = reform(savefiles_1to2d_bin[i,*,*], n_elements(kperp_density_names), $
      n_elements(wedge_1dbin_names))
    savefiles_mask_use = reform(savefiles_2d_masked[i,*,*], n_elements(kperp_density_names), $
      n_elements(wedge_1dbin_names))
    test_1d = min(test_save_1d[i,*,*])
    test_bin = min(test_save_bin[i,*,*])
    test_mask = min(test_save_mask[i,*,*])

    ;; if binsizes are specified, check that binsize is right
    if (tag_exist(binning_2d_options, 'kperp_bin') or $
        tag_exist(binning_2d_options, 'kpar_bin')) and test_2d gt 0 then begin
      if tag_exist(binning_2d_options, 'kpar_bin') ne 0 then begin
        kpar_bin_file = fltarr(n_elements(savefile_2d_use))
        for j=0, n_elements(savefile_2d_use)-1 do begin
          kpar_bin_file[j] = getvar_savefile(savefile_2d_use[j], 'kpar_bin')
        endfor
        if max(abs(binning_2d_options.kpar_bin - kpar_bin_file)) gt 0. then test_2d=0

        kpar_bin_file = fltarr(n_elements(savefiles_bin_use))
        for j=0, n_elements(savefiles_bin_use)-1 do begin
          kpar_bin_file[j] = getvar_savefile(savefiles_bin_use[j], 'kpar_bin')
        endfor
        if max(abs(binning_2d_options.kpar_bin - kpar_bin_file)) gt 0. then test_bin=0

        kpar_bin_file = fltarr(n_elements(savefiles_mask_use))
        for j=0, n_elements(savefiles_mask_use)-1 do begin
          kpar_bin_file[j] = getvar_savefile(savefiles_mask_use[j], 'kpar_bin')
        endfor
        if max(abs(binning_2d_options.kpar_bin - kpar_bin_file)) gt 0. then test_mask=0
      endif
      if tag_exist(binning_2d_options, 'kperp_bin') ne 0 then begin
        kperp_bin_file = fltarr(n_elements(savefile_2d_use))
        for j=0, n_elements(savefile_2d_use)-1 do begin
          kperp_bin_file[j] = getvar_savefile(savefile_2d_use[j], 'kperp_bin')
        endfor
        if max(abs(binning_2d_options.kperp_bin - kperp_bin_file)) gt 0. then test_2d=0

        kperp_bin_file = fltarr(n_elements(savefiles_bin_use))
        for j=0, n_elements(savefiles_bin_use)-1 do begin
          kperp_bin_file[j] = getvar_savefile(savefiles_bin_use[j], 'kperp_bin')
        endfor
        if max(abs(binning_2d_options.kperp_bin - kperp_bin_file)) gt 0. then test_bin=0

        kperp_bin_file = fltarr(n_elements(savefiles_mask_use))
        for j=0, n_elements(savefiles_mask_use)-1 do begin
          kperp_bin_file[j] = getvar_savefile(savefiles_mask_use[j], 'kperp_bin')
        endfor
        if max(abs(binning_2d_options.kperp_bin - kperp_bin_file)) gt 0. then test_mask=0
      endif
    endif

    if plot_types.plot_k0_power then begin
      if test_k0 gt 0 and tag_exist(binning_2d_options, 'kperp_bin') ne 0 then begin
        kperp_bin_file = fltarr(n_elements(savefile_k0_use))
        for j=0, n_elements(savefile_k0_use) do begin
          kperp_bin_file[j] = getvar_savefile(savefile_k0_use[j], 'k_bin')
        endfor
        if max(abs(binning_2d_options.kperp_bin - kperp_bin_file)) gt 0. then test_k0=0
      endif

      if test_k0_masked gt 0 and tag_exist(binning_2d_options, 'kperp_bin') ne 0 then begin
        kperp_bin_file = fltarr(n_elements(savefile_k0_masked_use))
        for j=0, n_elements(savefile_k0_masked_use) do begin
          kperp_bin_file[j] = getvar_savefile(savefile_k0_masked_use[j], 'k_bin')
        endfor
        if max(abs(binning_2d_options.kperp_bin - kperp_bin_file)) gt 0. then test_k0_masked=0
      endif
    endif

    if tag_exist(binning_1d_options, 'k_bin') ne 0 and test_1d gt 0 then begin
      k_bin_file = fltarr(n_elements(savefile_1d_use))
      for j=0, n_elements(savefile_1d_use)-1 do begin
        k_bin_file[j] = getvar_savefile(savefile_1d_use[j], 'k_bin')
      endfor
      if max(abs(binning_1d_options.k_bin - k_bin_file)) gt 0. then test_1d=0

      k_bin_file = fltarr(n_elements(savefiles_bin_use))
      for j=0, n_elements(savefiles_bin_use)-1 do begin
        k_bin_file[j] = getvar_savefile(savefiles_bin_use[j], 'k_bin')
      endfor
      if max(abs(binning_1d_options.k_bin - k_bin_file)) gt 0. then test_bin=0
    endif

    if test_2d gt 0 and n_elements(freq_flags) ne 0 then begin
      for j=0, n_elements(savefile_2d_use)-1 do begin
        old_freq_mask = getvar_savefile(savefile_2d_use[j], 'freq_mask')
        if total(abs(old_freq_mask - file_struct_arr[i].freq_mask)) ne 0 then test_2d = 0
      endfor

      for j=0, n_elements(savefiles_bin_use)-1 do begin
        old_freq_mask = getvar_savefile(savefiles_bin_use[j], 'freq_mask')
        if total(abs(old_freq_mask - file_struct_arr[i].freq_mask)) ne 0 then test_bin = 0
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

    if test_1d gt 0 and (n_elements(kperp_range_1d_use) gt 0 or $
      tag_exist(binning_1d_options, 'kperp_range_lambda_1dave') gt 0 or $
        n_elements(kpar_range_1d_use) gt 0) then begin

      ;; check that 1d binning was over correct ranges
      kperp_range_used = fltarr(2, n_elements(savefile_1d_use))
      kperp_range_lambda_used = fltarr(2, n_elements(savefile_1d_use))
      kpar_range_used = fltarr(2, n_elements(savefile_1d_use))
      for j=0, n_elements(savefile_1d_use)-1 do begin
        kperp_range_used[*,j] = getvar_savefile(savefile_1d_use[j], 'kperp_range')
        kperp_range_lambda_used[*,j] = getvar_savefile(savefile_1d_use[j], 'kperp_range_lambda')
        kpar_range_used[*,j] = getvar_savefile(savefile_1d_use[j], 'kpar_range')
      endfor
      if n_elements(kperp_range_1d_use) gt 0 then begin
        if max(abs(1.-kperp_range_used/rebin(kperp_range_1d_use, 2, $
          n_elements(savefile_1d_use)))) gt 1e-6 then begin
          test_1d = 0
        endif
      endif
      if tag_exist(binning_1d_options, 'kperp_range_lambda_1dave') gt 0 then begin
        if max(abs(1.-kperp_range_lambda_used / $
               rebin(binning_1d_options.kperp_range_lambda_1dave, 2, $
                     n_elements(savefile_1d_use)))) gt 1e-6 then begin
          test_1d = 0
        endif
      endif
      if n_elements(kpar_range_1d_use) gt 0 then begin
        if max(abs(1.-kpar_range_used/rebin(kpar_range_1d_use, 2, $
                                            n_elements(savefile_1d_use)))) gt 1e-6 then begin
          test_1d = 0
        endif
      endif

      kperp_range_used = fltarr(2, n_elements(savefiles_bin_use))
      kperp_range_lambda_used = fltarr(2, n_elements(savefiles_bin_use))
      kpar_range_used = fltarr(2, n_elements(savefiles_bin_use))
      for j=0, n_elements(savefiles_bin_use)-1 do begin
        kperp_range_used[*,j] = getvar_savefile(savefiles_bin_use[j], 'kperp_range')
        kperp_range_lambda_used[*,j] = getvar_savefile(savefiles_bin_use[j], 'kperp_range_lambda')
        kpar_range_used[*,j] = getvar_savefile(savefiles_bin_use[j], 'kpar_range')
      endfor
      if n_elements(kperp_range_1d_use) gt 0 then begin
        if max(abs(1.-kperp_range_used / $
                   rebin(kperp_range_1d_use, 2, $
                         n_elements(savefiles_bin_use)))) gt 1e-6 then begin
          test_bin = 0
        endif
      endif
      if tag_exist(binning_1d_options, 'kperp_range_lambda_1dave') gt 0 then begin
        if max(abs(1.-kperp_range_lambda_used / $
                   rebin(binning_1d_options.kperp_range_lambda_1dave, 2, $
                         n_elements(savefiles_bin_use)))) gt 1e-6 then begin
          test_bin = 0
        endif
      endif
      if n_elements(kpar_range_1d_use) gt 0 then begin
        if max(abs(1.-kpar_range_used / rebin(kpar_range_1d_use, 2, $
                                              n_elements(savefiles_bin_use)))) gt 1e-6 then begin
          test_bin = 0
        endif
      endif

      kperp_range_used = fltarr(2, n_elements(savefiles_mask_use))
      kperp_range_lambda_used = fltarr(2, n_elements(savefiles_mask_use))
      kpar_range_used = fltarr(2, n_elements(savefiles_mask_use))
      for j=0, n_elements(savefiles_mask_use)-1 do begin
        kperp_range_used[*,j] = getvar_savefile(savefiles_mask_use[j], 'kperp_range')
        kperp_range_lambda_used[*,j] = getvar_savefile(savefiles_mask_use[j], 'kperp_range_lambda')
        kpar_range_used[*,j] = getvar_savefile(savefiles_mask_use[j], 'kpar_range')
      endfor
      if n_elements(kperp_range_1d_use) gt 0 then begin
        if max(abs(1.-kperp_range_used / rebin(kperp_range_1d_use, 2, $
                                               n_elements(savefiles_mask_use)))) gt 1e-6 then begin
          test_mask = 0
        endif
      endif
      if tag_exist(binning_1d_options, 'kperp_range_lambda_1dave') gt 0 then begin
        if max(abs(1.-kperp_range_lambda_used / $
                   rebin(binning_1d_options.kperp_range_lambda_1dave, 2, $
                         n_elements(savefiles_mask_use)))) gt 1e-6 then begin
          test_mask = 0
        endif
      endif
      if n_elements(kpar_range_1d_use) gt 0 then begin
        if max(abs(1.-kpar_range_used / $
                   rebin(kpar_range_1d_use, 2, $
                         n_elements(savefiles_mask_use)))) gt 1e-6 then begin
          test_mask = 0
        endif
      endif
    endif

    if plot_types.plot_kpar_power then begin
      if test_kpar gt 0 and (n_elements(kperp_range_1d_use) gt 0 or $
        tag_exist(binning_1d_options, 'kperp_range_lambda_kparpower')) then begin

        ;; check that 1d binning was over correct ranges
        kperp_range_used = getvar_savefile(savefile_kpar_use, 'kperp_range')
        kperp_range_lambda_used = getvar_savefile(savefile_kpar_use, 'kperp_range_lambda')
        if n_elements(kperp_range_1d_use) gt 0 then begin
          if max(abs(1.-kperp_range_used/kperp_range_1d_use)) gt 1e-6 then begin
            test_kpar = 0
          endif
        endif
        if tag_exist(binning_1d_options, 'kperp_range_lambda_kparpower') then begin
          if max(abs(1.-binning_1d_options.kperp_range_lambda_kparpower / $
                     rebin(binning_1d_options.kperp_range_lambda_kparpower, 2, $
                           n_elements(savefile_1d_use)))) gt 1e-6 then begin
            test_kpar = 0
          endif
        endif
      endif
    endif

    if plot_types.plot_kperp_power then begin
      if test_kperp gt 0 and n_elements(kperp_range_1d_use) gt 0 then begin
        ;; check that 1d binning was over correct ranges
        kperp_range_used = getvar_savefile(savefile_kperp_use, 'kperp_range')
        kperp_range_lambda_used = getvar_savefile(savefile_kperp_use, 'kperp_range_lambda')
        kpar_range_used = getvar_savefile(savefile_kperp_use, 'kpar_range')
        if n_elements(kperp_range_1d_use) gt 0 then begin
          if max(abs(1-kperp_range_used/kperp_range_1d_use)) gt 1e-6 then begin
            test_kperp = 0
          endif
        endif
        if tag_exist(binning_1d_options, 'kperp_range_lambda_1dave') gt 0 then begin
          if max(abs(1-kperp_range_lambda_used / $
                     rebin(binning_1d_options.kperp_range_lambda_1dave, 2, $
                     n_elements(savefile_1d_use)))) gt 1e-6 then begin
            test_kperp = 0
          endif
        endif
        if n_elements(kpar_range_1d_use) gt 0 then begin
          if max(abs(1-kpar_range_used/kpar_range_1d_use)) gt 1e-6 then begin
            test_kperp = 0
          endif
        endif
      endif
    endif

    test = test_2d * test_1d * test_bin
    if plot_types.plot_kpar_power then test = test * test_kpar
    if plot_types.plot_kperp_power then test = test * test_kperp
    if plot_types.plot_k0_power then test = test * test_k0 * test_k0_masked

    if test eq 0 then begin
      if n_elements(plotfile_binning_hist) gt 0 then begin
        plotfile_binning_hist_use = reform(plotfile_binning_hist[i,*,*], $
          n_elements(kperp_density_names), n_elements(wedge_1dbin_names))
      endif

      if healpix or not keyword_set(uvf_input) then begin
        weight_refresh = intarr(n_cubes)
        if refresh_options.refresh_dft then begin
          temp = weight_ind[uniq(weight_ind, sort(weight_ind))]
          for j=0, n_elements(temp)-1 do begin
            weight_refresh[(where(weight_ind eq temp[j]))[0]] = 1
          endfor
        endif
        refresh_options = create_refresh_options(refresh_options = refresh_options, $
          refresh_weight_dft = weight_refresh[i])
      endif

      ps_power, file_struct_arr[i], sim = sim, fix_sim_input = fix_sim_input, $
        uvf_input = uvf_input, freq_ch_range = freq_ch_range, freq_flags = freq_flags, $
        refresh_options = refresh_options, uvf_options = uvf_options, $
        ps_options = ps_options, binning_2d_options = binning_2d_options, $
        binning_1d_options = binning_1d_options, plot_options = plot_options, $
        plot_types = plot_types, $
        savefile_2d = savefile_2d_use, savefile_1d = savefile_1d_use, $
        savefile_1to2d_bin = savefiles_bin_use, savefile_masked_2d = savefiles_mask_use, $
        savefile_masked_k0 = savefile_k0_masked_use, $
        savefile_kpar_power = savefile_kpar_use, savefile_kperp_power = savefile_kperp_use, $
        savefile_k0 = savefile_k0_use, $
        bin_arr_3d = bin_arr_3d

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

  if not tag_exist(binning_2d_options, 'kpar_bin') then begin
    binning_2d_options = create_binning_2d_options(binning_2d_options = binning_2d_options, $
      kpar_bin = kpar_bin)
  endif

  if not tag_exist(binning_2d_options, 'kperp_bin') then begin
    binning_2d_options = create_binning_2d_options(binning_2d_options = binning_2d_options, $
      kperp_bin = kperp_bin)
  endif

  if tag_exist(plot_2d_options, 'kperp_lambda_plot_range') and $
    not tag_exist(plot_2d_options, 'kperp_plot_range') then begin

    kperp_plot_range = plot_2d_options.kperp_lambda_plot_range / kperp_lambda_conv

    ;; if we're plotting in [k]=h/Mpc then need to convert from 1/Mpc
    if plot_options.hinv then kperp_plot_range = kperp_plot_range / hubble_param

    plot_2d_options = create_plot_2d_options(plot_2d_options = plot_2d_options, $
      kperp_plot_range = kperp_plot_range)
  endif

  if not tag_exist(plot_2d_options, 'kperp_plot_range') then begin
    if tag_exist(uvf_options, 'max_uv_lambda') gt 0 then begin
      max_kperp_lambda = min([uvf_options.max_uv_lambda, $
                              min(file_struct_arr.kspan/2.), $
                              min(file_struct_arr.max_baseline_lambda)])
    endif else begin
      max_kperp_lambda = min([file_struct_arr.kspan/2.,file_struct_arr.max_baseline_lambda])
    endelse
    kperp_plot_range = [5./kperp_lambda_conv, max_kperp_lambda/kperp_lambda_conv]

    ;; if we're plotting in [k]=h/Mpc then need to convert from 1/Mpc
    if plot_options.hinv then kperp_plot_range = kperp_plot_range / hubble_param

    plot_2d_options = create_plot_2d_options(plot_2d_options = plot_2d_options, $
      kperp_plot_range = kperp_plot_range)
  endif

  if plot_options.pub then begin
    plot_fadd = ''
    if plot_2d_options.kperp_linear_axis and plot_2d_options.kpar_linear_axis then begin
      plot_fadd = plot_fadd + '_linaxes'
    endif else begin
      if plot_2d_options.kperp_linear_axis then plot_fadd = plot_fadd + '_kperplinaxis'
      if plot_2d_options.kpar_linear_axis then plot_fadd = plot_fadd + '_kparlinaxis'
    endelse

    if plot_options.individual_plots then begin
      if not tag_exist(plot_options, 'plot_filebase') then begin
        plotfile_base = file_struct_arr.savefilebase + power_tag
        plotfile_base_wt = general_filebase + '_' + $
          file_struct_arr[uniq(weight_ind, sort(weight_ind))].pol + power_tag
      endif else begin
        plotfile_base = plot_options.plot_filebase + uvf_tag + $
          file_struct_arr.file_label + power_tag
        plotfile_base_wt = plot_options.plot_filebase + uvf_tag + $
          '_' + file_struct_arr[uniq(weight_ind, sort(weight_ind))].pol + power_tag
      endelse
    endif else begin
      if not tag_exist(plot_options, 'plot_filebase') then begin
        plotfile_base = general_filebase + power_tag
      endif else begin
        plotfile_base = plot_options.plot_filebase + uvf_tag + power_tag
      endelse
      plotfile_base_wt = plotfile_base
    endelse

    plotfile_path_2d = plotfile_path + file_struct_arr[0].subfolders.bin_2d
    plotfiles_shape = [n_elements(plotfile_base), n_elements(kperp_density_names)]
    plotfiles_err_shape = [n_elements(plotfile_base_wt), n_elements(kperp_density_names)]
    plotfiles_2d = strarr(plotfiles_shape)
    plotfiles_2d_error = strarr(plotfiles_err_shape)
    plotfiles_2d_noise_expval = strarr(plotfiles_err_shape)
    plotfiles_2d_noise = strarr(plotfiles_shape)
    plotfiles_2d_sim_noise = strarr(plotfiles_shape)
    plotfiles_2d_sim_noise_diff = strarr(plotfiles_shape)
    plotfiles_2d_snr = strarr(plotfiles_shape)
    plotfiles_2d_nnr = strarr(plotfiles_shape)
    plotfiles_2d_sim_snr = strarr(plotfiles_shape)
    plotfiles_2d_sim_nnr = strarr(plotfiles_shape)
    for j=0, n_elements(kperp_density_names)-1 do begin
      this_fadd = fadd_2dbin + kperp_density_names[j]
      this_begin = plotfile_path_2d + plotfile_base + this_fadd
      this_begin_wt = plotfile_path_2d + plotfile_base_wt + this_fadd
      this_end = plot_fadd + plot_options.plot_exten

      plotfiles_2d[*,j] = this_begin + '_2dkpower' +  this_end
      plotfiles_2d_error[*,j] = this_begin_wt + '_2derror' + this_end
      if plot_options.individual_plots then begin
        plotfiles_2d_noise_expval[*,j] = plotfile_path_2d + plotfile_base_wt + $
          this_fadd + '_2dnoise_expval' + this_end
      endif
      plotfiles_2d_noise[*,j] = this_begin + '_2dnoise' + this_end
      plotfiles_2d_sim_noise[*,j] = this_begin + '_2dsimnoise' + this_end
      plotfiles_2d_sim_noise_diff[*,j] = this_begin + '_2dsimnoisediff' + this_end

      plotfiles_2d_snr[*,j] = this_begin + '_2dsnr' + this_end

      plotfiles_2d_nnr[*,j] = this_begin + '_2dnnr' + this_end
      plotfiles_2d_sim_snr[*,j] = this_begin + '_2dsimsnr' + this_end
      plotfiles_2d_sim_nnr[*,j] = this_begin + '_2dsimnnr' + this_end
    endfor


    if not tag_exist(plot_options, 'plot_filebase') then begin
      plotfile_1d_base = general_filebase
    endif else begin
      plotfile_1d_base = plot_options.plot_filebase + uvf_tag
    endelse
    plotfile_path_1d = plotfile_path + file_struct_arr[0].subfolders.bin_1d
    begin_1d = plotfile_path_1d + plotfile_1d_base + power_tag + wedge_1dbin_names + fadd_1dbin
    plotfile_1d = begin_1d + '_1dkpower' + plot_options.plot_exten
    plotfile_1d_noise = begin_1d + '_1dnoise' + plot_options.plot_exten
    plotfile_1d_sim_noise_diff = begin_1d + '_1dsimnoisediff' + plot_options.plot_exten
    plotfile_1d_sim_noise = begin_1d + '_1dsimnoise' + plot_options.plot_exten

    begin_1d_axis = plotfile_path_1d + plotfile_1d_base + power_tag
    plotfile_kpar_power = begin_1d_axis + fadd_kpar_1d + '_kpar_power' + $
      plot_options.plot_exten
    plotfile_kperp_power = begin_1d_axis + fadd_kperp_1d + '_kperp_power' + $
      plot_options.plot_exten
    plotfile_k0_power = begin_1d_axis + fadd_2dbin + '_k0power' + $
      plot_options.plot_exten
  endif

  if plot_types.plot_1to2d then begin

    if keyword_set(pol_inc) then begin
      if n_elements(pol_inc) eq 1 then pol_1to2d_use = pol_inc else pol_1to2d_use = pol_inc[0]
    endif else pol_1to2d_use='xx'
    if keyword_set(type_inc) then begin
      if n_elements(type_inc) eq 1 then type_1to2d_use = type_inc else type_1to2d_use = type_inc[0]
    endif else type_1to2d_use='res'
    wh_2d_use = where(file_struct_arr.type eq type_1to2d_use and $
      file_struct_arr.pol eq pol_1to2d_use, count_type_pol)
    if count_type_pol eq 0 then wh_2d_use = 0

    if plot_options.pub then begin
      begin_1to2d = plotfile_path_2d + 'from_1d/' + plotfile_base + fadd_2dbin + $
        fadd_1dbin + '_' + type_1to2d_use + '_' + pol_1to2d_use
      plotfile_1to2d_heatmap =  begin_1to2d + '_1dheatmap' + plot_options.plot_exten
      plotfile_1to2d_contours = begin_1to2d + '_1dcontours' + plot_options.plot_exten
      plotfile_1to2d_noisefrac = begin_1to2d + '_1dnoisefrac' + plot_options.plot_exten
      plotfile_1to2d_contour_zoom = begin_1to2d + '_1dcontour_zoom' + plot_options.plot_exten
    endif
  endif

  if plot_options.pub and plot_types.plot_2d_masked then begin

    masked_plotfiles_shape = [n_elements(plotfile_base), $
      n_elements(kperp_density_names), n_elements(wedge_1dbin_names)]

    plotfiles_2d_masked = strarr(masked_plotfiles_shape)
    plotfiles_2d_error_masked = strarr(masked_plotfiles_shape)
    plotfiles_2d_noise_expval_masked = strarr(masked_plotfiles_shape)
    plotfiles_2d_noise_masked = strarr(masked_plotfiles_shape)
    plotfiles_2d_sim_noise_masked = strarr(masked_plotfiles_shape)
    plotfiles_2d_sim_noise_diff_masked = strarr(masked_plotfiles_shape)
    plotfiles_2d_snr_masked = strarr(masked_plotfiles_shape)
    plotfiles_2d_nnr_masked = strarr(masked_plotfiles_shape)
    plotfiles_2d_sim_snr_masked = strarr(masked_plotfiles_shape)
    plotfiles_2d_sim_nnr_masked = strarr(masked_plotfiles_shape)
    for i=0, n_elements(wedge_1dbin_names)-1 do begin
      for j=0, n_elements(kperp_density_names)-1 do begin
        this_fadd = fadd_2dbin + kperp_density_names[j] + wedge_1dbin_names[i] + $
          fadd_1dmask
        this_begin = plotfile_path_2d + 'from_1d/' + plotfile_base + this_fadd
        this_end = plot_fadd + plot_options.plot_exten

        plotfiles_2d_masked[*,j,i] = this_begin + '_2dkpower' + this_end
        plotfiles_2d_error_masked[*,j,i] = this_begin +  '_2derror' + plot_fadd + $
          plot_options.plot_exten
        if plot_options.individual_plots then begin
          plotfiles_2d_noise_expval_masked[*,j,i] = plotfile_path_2d + 'from_1d/' + $
            plotfile_base_wt + this_fadd + '_2dnoise_expval' + this_end
        endif
        plotfiles_2d_noise_masked[*,j,i] = this_begin + '_2dnoise' + this_end
        plotfiles_2d_sim_noise_masked[*,j,i] = this_begin + '_2dsimnoise' + this_end
        plotfiles_2d_sim_noise_diff_masked[*,j,i] = this_begin + '_2dsimnoisediff' + $
          this_end

        plotfiles_2d_snr_masked[*,j,i] = this_begin + '_2dsnr' + this_end

        plotfiles_2d_nnr_masked[*,j,i] = this_begin + '_2dnnr' + this_end
        plotfiles_2d_sim_snr_masked[*,j,i] = this_begin + '_2dsimsnr' + this_end
        plotfiles_2d_sim_nnr_masked[*,j,i] = this_begin + '_2dsimnnr' + this_end
      endfor
    endfor
  endif

  if plot_types.plot_slices then begin
    if not tag_exist(plot_types, 'slice_type') then begin
       slice_type = 'sumdiff'
       plot_types = create_plot_types(plot_types = plot_types, slice_type = slice_type)
    endif
    slice_type_enum = ['raw', 'divided', 'sumdiff', 'weights', 'variance', $
                       'power', 'var_power']

    wh_slice_type = where(slice_type_enum eq plot_types.slice_type, count_slice_type)
    if count_slice_type eq 0 then begin
      message, 'slice_type not recognized. options are: ' + strjoin(slice_type_enum, ', ')
    endif

    if plot_types.slice_type eq 'power' or plot_types.slice_type eq 'var_power' then begin
      slice_space = 'kspace'
    endif else begin
      slice_space = 'uvf'
    endelse

    if slice_space eq 'uvf' then begin
      uvf_type_enum = ['abs', 'phase', 'real', 'imaginary', 'normalized']
      if not tag_exist(plot_types, 'uvf_plot_type') then begin
        uvf_plot_type = 'abs'
        plot_types = create_plot_types(plot_types = plot_types, uvf_plot_type = uvf_plot_type)
      endif
      wh = where(uvf_type_enum eq plot_types.uvf_plot_type, count_uvf_type)
      if count_uvf_type eq 0 then begin
        message, 'unknown uvf_plot_type. Use one of: ' + print, strjoin(uvf_type_enum, ', ')
      endif

      if plot_types.uvf_plot_type eq 'phase' then uvf_log = 0 else uvf_log=1
    endif

    if plot_options.pub then begin

      slice_plotfile_base = plotfile_path + file_struct_arr[0].subfolders.slices + $
        plotfile_base + '_' + plot_types.slice_type
      slice_plotfile_end = plot_fadd + plot_options.plot_exten
      if plot_options.individual_plots then begin

        indv_plotfile_base = plotfile_path + file_struct_arr[0].subfolders.slices + $
          transpose([[plotfile_base +'_even_'], [plotfile_base +'_odd_']]) + $
          '_' + plot_types.slice_type
        if slice_space eq 'uvf' then begin
          if plot_types.slice_type eq 'sumdiff' then begin
              sumdiff_plotfile_base = plotfile_path + file_struct_arr[0].subfolders.slices + $
                [plotfile_base +'_sum_',plotfile_base +'_diff_'] + $
                '_' + plot_types.slice_type
              uf_slice_plotfile = sumdiff_plotfile_base + '_uf_plane' + slice_plotfile_end
              vf_slice_plotfile = sumdiff_plotfile_base + '_vf_plane' + slice_plotfile_end
              uv_slice_plotfile = sumdiff_plotfile_base + '_uv_plane' + slice_plotfile_end
          endif else begin
              uf_slice_plotfile = indv_plotfile_base + '_uf_plane' + slice_plotfile_end
              vf_slice_plotfile = indv_plotfile_base + '_vf_plane' + slice_plotfile_end
              uv_slice_plotfile = indv_plotfile_base + '_uv_plane' + slice_plotfile_end
          endelse
        endif else begin
          uf_slice_plotfile = slice_plotfile_base + '_xz_plane' + slice_plotfile_end
          vf_slice_plotfile = slice_plotfile_base + '_yz_plane' + slice_plotfile_end
          uv_slice_plotfile = slice_plotfile_base + '_xy_plane' + slice_plotfile_end
        endelse
      endif else begin
        if slice_space eq 'uvf' then begin
          uf_slice_plotfile = slice_plotfile_base + '_uf_plane' + slice_plotfile_end
          vf_slice_plotfile = slice_plotfile_base + '_vf_plane' + slice_plotfile_end
          uv_slice_plotfile = slice_plotfile_base + '_uv_plane' + slice_plotfile_end
        endif else begin
          uf_slice_plotfile = slice_plotfile_base + '_xz_plane' + slice_plotfile_end
          vf_slice_plotfile = slice_plotfile_base + '_yz_plane' + slice_plotfile_end
          uv_slice_plotfile = slice_plotfile_base + '_xy_plane' + slice_plotfile_end
        endelse
      endelse

    endif else begin
      case plot_types.slice_type of
        'sumdiff': slice_titles = ['sum ' + file_struct_arr.file_label, 'diff ' + $
          file_struct_arr.file_label]
        'power': slice_titles = file_struct_arr.file_label
        'var_power': slice_titles = file_struct_arr.file_label
        else: slice_titles = file_struct_arr.uvf_label
      endcase
    endelse

    case plot_types.slice_type of
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
      'variance': begin
        uf_slice_savefile = file_struct_arr.uf_var_savefile
        vf_slice_savefile = file_struct_arr.vf_var_savefile
        uv_slice_savefile = file_struct_arr.uv_var_savefile
        slice_titles = file_struct_arr.uvf_label
      end
      'power': begin
        uf_slice_savefile = file_struct_arr.xz_savefile
        vf_slice_savefile = file_struct_arr.yz_savefile
        uv_slice_savefile = file_struct_arr.xy_savefile
        slice_titles = file_struct_arr.file_label
      end
      'var_power': begin
        uf_slice_savefile = file_struct_arr.xz_savefile
        vf_slice_savefile = file_struct_arr.yz_savefile
        uv_slice_savefile = file_struct_arr.xy_savefile
        slice_titles = file_struct_arr.file_label
      end
    endcase

    if plot_types.slice_type eq 'weights' or plot_types.slice_type eq 'variance' $
        or plot_types.slice_type eq 'var_power' then begin
      ;; we don't need separate ones for dirty, model, residual
      wh_use = where(strpos(slice_titles, 'dirty')>0)

      uf_slice_savefile = uf_slice_savefile[wh_use]
      vf_slice_savefile = vf_slice_savefile[wh_use]
      uv_slice_savefile = uv_slice_savefile[wh_use]
      slice_titles = slice_titles[wh_use]
    endif
  endif

  if plot_options.pub then font = 1 else font = -1

  window_num=0
  if plot_types.plot_stdset then begin

    for j=0, n_elements(kperp_density_names)-1 do begin
      note_use = plot_options.note + ', ' + kperp_density_names[j]

      if plot_options.pub then begin
        plotfiles_use = {power: plotfiles_2d[*,j], error: plotfiles_2d_error[*,j], $
          noise_expval: plotfiles_2d_noise_expval[*,j], snr: plotfiles_2d_snr[*,j], $
          noise: plotfiles_2d_noise[*,j], sim_noise: plotfiles_2d_sim_noise[*,j], $
          sim_noise_diff: plotfiles_2d_sim_noise_diff[*,j], nnr: plotfiles_2d_nnr[*,j], $
          sim_nnr: plotfiles_2d_sim_nnr[*,j]}
      endif

      make_std_plots, file_struct_arr, savefiles_2d[*,j], titles, $
        plotfile_struct = plotfiles_use, plot_sim_noise = plot_types.plot_sim_noise, $
        window_num = window_num, vs_note = vs_note, wedge_amp = binning_1d_options.wedge_amps, $
        plot_options = plot_options, plot_2d_options = plot_2d_options
    endfor

  endif

  if plot_types.plot_slices then begin

    if plot_options.pub and plot_options.individual_plots then begin
      case plot_types.slice_type of
        'power': nplots = ntype*npol
        'var_power': nplots = npol
        'weights': nplots = nfiles*npol
        'variance': nplots = nfiles*npol
        else: nplots = ntype*nfiles*npol
      endcase

      for i=0, nplots-1 do begin

        if plot_types.slice_type eq 'power' or plot_types.slice_type eq 'var_power' then begin
          if plot_types.slice_type eq 'power' then begin
            kspace_plot_type = 'power'
          endif else begin
            kspace_plot_type = 'variance'
          endelse

          kpower_slice_plot, uf_slice_savefile[i], multi_pos = pos_use, $
            start_multi_params = start_multi_params, plotfile = uf_slice_plotfile[i], $
            plot_options = plot_options, plot_2d_options = plot_2d_options, $
            plot_type = kspace_plot_type, title_prefix = slice_titles[i], note = note_use, $
            wedge_amp = binning_1d_options.wedge_amps, window_num = window_num

          kpower_slice_plot, vf_slice_savefile[i], multi_pos = pos_use, $
            start_multi_params = start_multi_params, plotfile = vf_slice_plotfile[i], $
            plot_options = plot_options, plot_2d_options = plot_2d_options, $
            plot_type = kspace_plot_type, title_prefix = slice_titles[i], note = note_use, $
            wedge_amp = binning_1d_options.wedge_amps, window_num = window_num

          kpower_slice_plot, uv_slice_savefile[i], multi_pos = pos_use, $
            start_multi_params = start_multi_params, plotfile = uv_slice_plotfile[i], $
            plot_options = plot_options, plot_2d_options = plot_2d_options, $
            plot_type = kspace_plot_type, title_prefix = slice_titles[i], note = note_use,  $
            wedge_amp = binning_1d_options.wedge_amps, window_num = window_num

        endif else begin

          uvf_slice_plot, uf_slice_savefile[i], multi_pos = pos_use, $
            start_multi_params = start_multi_params, $
            plotfile = uf_slice_plotfile[i], type = plot_types.uvf_plot_type, $
            plot_options = plot_options, plot_2d_options = plot_2d_options, $
            title_prefix = slice_titles[i], note = note_use, $
            log = uvf_log, window_num = window_num

          uvf_slice_plot, vf_slice_savefile[i], multi_pos = pos_use, $
            start_multi_params = start_multi_params, $
            plotfile = vf_slice_plotfile[i], type = plot_types.uvf_plot_type, $
            plot_options = plot_options, plot_2d_options = plot_2d_options, $
            title_prefix = slice_titles[i], note = note_use, $
            log = uvf_log, window_num = window_num

          uvf_slice_plot, uv_slice_savefile[i], multi_pos = pos_use, $
            start_multi_params = start_multi_params, $
            plotfile = uv_slice_plotfile[i], type = plot_types.uvf_plot_type, $
            plot_options = plot_options, plot_2d_options = plot_2d_options, $
            title_prefix = slice_titles[i], note = note_use, $
            log = uvf_log, window_num = window_num

        endelse

      endfor

    endif else begin

      if plot_2d_options.kperp_linear_axis then begin
        ;; aspect ratio doesn't work out for kperp_linear with multiple rows
        nrow = 1

        case plot_types.slice_type of
          'power': ncol = ntype*npol
          'var_power': ncol = npol
          'weights': ncol = nfiles*npol
          'variance': ncol = nfiles*npol
          else: ncol = ntype*nfiles*npol
        endcase
      endif else begin
        case plot_types.slice_type of
          'power': begin
            ncol = ntype
            nrow = npol
          end
          'var_power': begin
            ncol = 1
            nrow = npol
          end
          'weights': begin
            ncol = nfiles
            nrow = npol
          end
          'variance': begin
            ncol = nfiles
            nrow = npol
          end
          else: begin
            if ntype gt 1 then begin
              ncol=ntype
              nrow = npol*nfiles
            endif else begin
              ncol=ntype*nfiles
              nrow = npol
            endelse
        end
        endcase
      endelse

      if plot_types.slice_type eq 'raw' or plot_types.slice_type eq 'divided' then begin
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
        if i eq (nrow*ncol)-1 and tag_exist(plot_options, 'note') then begin
          note_use = plot_options.note
        endif else begin
          note_use = ''
        endelse

        if plot_types.slice_type eq 'power' or plot_types.slice_type eq 'var_power' then begin
          if plot_types.slice_type eq 'power' then begin
            kspace_plot_type = 'power'
          endif else begin
            kspace_plot_type = 'variance'
          endelse

          kpower_slice_plot, uf_slice_savefile[i], multi_pos = pos_use, $
            start_multi_params = start_multi_params, plotfile = uf_slice_plotfile, $
            plot_options = plot_options, plot_2d_options = plot_2d_options, $
            plot_type = kspace_plot_type, title_prefix = slice_titles[i], note = note_use,  $
            wedge_amp = binning_1d_options.wedge_amps, window_num = window_num
        endif else begin

          uvf_slice_plot, uf_slice_savefile[i], multi_pos = pos_use, $
            start_multi_params = start_multi_params, $
            plotfile = uf_slice_plotfile, type = plot_types.uvf_plot_type, $
            plot_options = plot_options, plot_2d_options = plot_2d_options, $
            title_prefix = slice_titles[i], note = note_use, $
            log = uvf_log, window_num = window_num
        endelse

        if i eq 0 then begin
          positions = pos_use
          undefine, start_multi_params
        endif
      endfor
      if plot_options.pub then begin
        cgps_close, png = plot_options.png, pdf = plot_options.pdf, $
          delete_ps = plot_options.delete_ps, density = 600
      endif

      window_num = window_num+1
      start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}
      undefine, pos_use

      for i=0, (nrow*ncol)-1 do begin
        if i gt 0 then  pos_use = positions[*,i]
        if i eq (nrow*ncol)-1 and tag_exist(plot_options, 'note') then begin
          note_use = plot_options.note
        endif else begin
          note_use = ''
        endelse

        if plot_types.slice_type eq 'power' or plot_types.slice_type eq 'var_power' then begin
          kpower_slice_plot, vf_slice_savefile[i], multi_pos = pos_use, $
            start_multi_params = start_multi_params, plotfile = vf_slice_plotfile, $
            plot_options = plot_options, plot_2d_options = plot_2d_options, $
            plot_type = kspace_plot_type, title_prefix = slice_titles[i], note = note_use,  $
            wedge_amp = binning_1d_options.wedge_amps, window_num = window_num
        endif else begin

          uvf_slice_plot, vf_slice_savefile[i], multi_pos = pos_use, $
            start_multi_params = start_multi_params, $
            plotfile = vf_slice_plotfile, type = plot_types.uvf_plot_type, $
            plot_options = plot_options, plot_2d_options = plot_2d_options, $
            title_prefix = slice_titles[i], note = note_use, $
            log = uvf_log, window_num = window_num
        endelse

        if i eq 0 then begin
          positions = pos_use
          undefine, start_multi_params
        endif
      endfor
      if plot_options.pub then begin
        cgps_close, png = plot_options.png, pdf = plot_options.pdf, $
          delete_ps = plot_options.delete_ps, density = 600
      endif

      window_num = window_num+1
      start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}
      undefine, pos_use

      for i=0, (nrow*ncol)-1 do begin
        if i gt 0 then  pos_use = positions[*,i]
        if i eq (nrow*ncol)-1 and tag_exist(plot_options, 'note') then begin
          note_use = plot_options.note
        endif else begin
          note_use = ''
        endelse

        if plot_types.slice_type eq 'power' or plot_types.slice_type eq 'var_power' then begin
          kpower_slice_plot, uv_slice_savefile[i], multi_pos = pos_use, $
            start_multi_params = start_multi_params, plotfile = uv_slice_plotfile, $
            plot_options = plot_options, plot_2d_options = plot_2d_options, $
            plot_type = kspace_plot_type, title_prefix = slice_titles[i], note = note_use,  $
            wedge_amp = binning_1d_options.wedge_amps, window_num = window_num
      endif else begin

          uvf_slice_plot, uv_slice_savefile[i], multi_pos = pos_use, $
            start_multi_params = start_multi_params, $
            plotfile = uv_slice_plotfile, type = plot_types.uvf_plot_type, $
            plot_options = plot_options, plot_2d_options = plot_2d_options, $
            title_prefix = slice_titles[i], note = note_use, $
            log = uvf_log, window_num = window_num
        endelse

        if i eq 0 then begin
          positions = pos_use
          undefine, start_multi_params
        endif
      endfor
      undefine, pos_use
      if plot_options.pub then begin
        cgps_close, png = plot_options.png, pdf = plot_options.pdf, $
          delete_ps = plot_options.delete_ps, density = 600
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

    if max(strmatch(varnames, 'ave_power', /fold_case)) gt 0 then begin
      ave_power_vals[k] = getvar_savefile(savefiles_1d[k,0,0], 'ave_power')
    endif else begin
      ave_power_vals[k] = -1
    endelse
    if max(strmatch(varnames, 'wt_ave_power', /fold_case)) gt 0 then begin
      wt_ave_power_vals[k] = getvar_savefile(savefiles_1d[k,0,0], 'wt_ave_power')
    endif else begin
      wt_ave_power_vals[k] = -1
    endelse
    if max(strmatch(varnames, 'uv_pix_area', /fold_case)) gt 0 then begin
      uv_pix_area[k] = getvar_savefile(savefiles_1d[k,0,0], 'uv_pix_area')
    endif else begin
      uv_pix_area[k] = -1
    endelse
    if max(strmatch(varnames, 'uv_area', /fold_case)) gt 0 then begin
      uv_area[k] = getvar_savefile(savefiles_1d[k,0,0], 'uv_area')

      nbsl_lambda2[k] = file_struct_arr[0].n_vis[0]/uv_area[k]
      if n_elements(freq_ch_range) gt 0 then begin
        nbsl_lambda2_freq[k,*] = file_struct_arr[0].n_vis_freq[0, $
          freq_ch_range[0]:freq_ch_range[1]]/uv_area[k]
      endif else begin
        nbsl_lambda2_freq[k,*] = file_struct_arr[0].n_vis_freq[0,*]/uv_area[k]
      endelse
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

    if max(strmatch(varnames, 'wt_ave_power_freq', /fold_case)) gt 0 then begin
      wt_ave_power_freq_vals[k,*] = (getvar_savefile(savefiles_1d[k,0,0], 'wt_ave_power_freq'))[0,*]
    endif else begin
      wt_ave_power_freq_vals[k,*] = -1
    endelse
    if max(strmatch(varnames, 'ave_power_freq', /fold_case)) gt 0 then begin
      ave_power_freq_vals[k,*] = (getvar_savefile(savefiles_1d[k,0,0], 'ave_power_freq'))[0,*]
    endif else begin
      ave_power_freq_vals[k,*] = -1
    endelse
    if max(strmatch(varnames, 'ave_power_uvf', /fold_case)) gt 0 then begin
      ave_power_uvf_vals[k] = (getvar_savefile(savefiles_1d[k,0,0], 'ave_power_uvf'))[0,*]
    endif else begin
      ave_power_uvf_vals[k] = -1
    endelse
    if max(strmatch(varnames, 'wt_ave_power_uvf', /fold_case)) gt 0 then begin
      wt_ave_power_uvf_vals[k] = (getvar_savefile(savefiles_1d[k,0,0], 'wt_ave_power_uvf'))[0,*]
    endif else begin
      wt_ave_power_uvf_vals[k] = -1
    endelse
  endfor
  cube_power_info = {ave_power:ave_power_vals, wt_ave_power:wt_ave_power_vals, $
    uv_pix_area:uv_pix_area, uv_area:uv_area, ave_weights:ave_weights_vals, $
    ave_weights_freq:ave_weights_freq_vals, wt_ave_power_freq:wt_ave_power_freq_vals, $
    ave_power_freq:ave_power_freq_vals, wt_ave_power_uvf:wt_ave_power_uvf_vals, $
    ave_power_uvf:ave_power_uvf_vals, nbsl_lambda2:nbsl_lambda2, $
    nbsl_lambda2_freq:nbsl_lambda2_freq}

  if plot_1d_options.plot_eor_1d or plot_1d_options.plot_flat_1d then begin
    case strlowcase(!version.os_family) OF
      'windows': split_delim = ';'
      'unix':    split_delim = ':'
    endcase
    path_dirs = strsplit(!path, split_delim, /extract)

    fhd_catalog_loc = strpos(path_dirs, 'catalog_data')
    wh_catalog = where(fhd_catalog_loc gt 0, count_catalog)
    if count_catalog gt 1 then begin
      void = min(strlen(path_dirs[wh_catalog]), min_loc)
      wh_catalog = wh_catalog[min_loc]
    endif
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

  if plot_types.plot_stdset then begin
    for i=0, n_elements(wedge_1dbin_names)-1 do begin
      file_arr = reform(savefiles_1d[*,*,i], n_cubes*n_elements(kperp_density_names))

      if n_elements(plotfile_1d) gt 0 then begin
        plotfile_use = plotfile_1d[i]
        plotfile_noise_use = plotfile_1d_noise[i]
      endif

      if i gt 0 then begin
        note_use = plot_options.note + ' ' + strjoin(strsplit(wedge_1dbin_names[i], '_', /extract), ' ')
      endif else begin
        note_use = plot_options.note
      endelse

      titles_use = strarr(n_cubes, n_elements(kperp_density_names))
      for j=0, n_elements(kperp_density_names)-1 do begin
        titles_use[*,j] = titles + kperp_density_names[j]
      endfor
      titles_use = reform(titles_use, n_cubes*n_elements(kperp_density_names))

      k_range = minmax([plot_2d_options.kperp_plot_range, $
        binning_2d_options.kpar_bin, plot_2d_options.kpar_plot_range[1]])

      if plot_types.plot_noise_1d then begin
        window_num = window_num+1
        kpower_1d_plots, file_arr, window_num = window_num, names = titles_use, $
          plot_options = plot_options, plot_1d_options = plot_1d_options, $
          plotfile = plotfile_noise_use, k_range = k_range, title = note_use + ' Ob. Noise', $
          note = note_1d, plot_sim_noise = plot_types.plot_sim_noise, /plot_noise
      endif

      if plot_1d_options.plot_eor_1d then begin
        if count_catalog gt 0 then begin
          psyms = [intarr(n_elements(file_arr))+10, -3]
          file_arr = [file_arr, eor_file_1d]
          titles_use = [titles_use, 'EoR signal']
        endif
      endif
      if plot_1d_options.plot_flat_1d then begin
        if count_catalog gt 0 then begin
          psyms = [intarr(n_elements(file_arr))+10, -3]
          file_arr = [file_arr, flat_file_1d]
          titles_use = [titles_use, 'input flat power']
        endif
      endif

      window_num = window_num+1
      kpower_1d_plots, file_arr, window_num = window_num, colors = colors, $
        plot_options = plot_options, plot_1d_options = plot_1d_options, $
        names = titles_use, psyms = psyms, plotfile = plotfile_use, $
        k_range = k_range, title = note_use, note = note_1d, $
        plot_sim_noise = plot_types.plot_sim_noise

    endfor
  endif

  if plot_types.plot_kpar_power then begin

    file_arr = reform(savefiles_kpar_1d, n_cubes*n_elements(kperp_density_names))

    titles_use = strarr(n_cubes, n_elements(kperp_density_names))
    for j=0, n_elements(kperp_density_names)-1 do begin
      titles_use[*,j] = titles + kperp_density_names[j]
    endfor
    titles_use = reform(titles_use, n_cubes*n_elements(kperp_density_names))
    if plot_1d_options.plot_eor_1d then begin
      if count_catalog gt 0 then begin
        psyms = [intarr(n_elements(file_arr))+10, -3]
        file_arr = [file_arr, eor_file_1d]
        titles_use = [titles_use, 'EoR signal']
      endif
    endif
    if plot_1d_options.plot_flat_1d then begin
      if count_catalog gt 0 then begin
        psyms = [intarr(n_elements(file_arr))+10, -3]
        file_arr = [file_arr, flat_file_1d]
        titles_use = [titles_use, 'input flat power']
      endif
    endif

    if tag_exist(plot_2d_options, 'kpar_plot_range') then begin
      k_range = plot_2d_options.kpar_plot_range
    endif

    window_num = window_num + 1
    kpower_1d_plots, file_arr, window_num = window_num, colors = colors, $
      plot_options = plot_options, plot_1d_options = plot_1d_options, $
      names = titles_use, psyms = psyms, plotfile = plotfile_kpar_power, $
      k_range = k_range, title = plot_options.note + ' kpar', note = note_kpar_1d, $
      plot_sim_noise = plot_types.plot_sim_noise, $
      /kpar_power, delay_axis = plot_2d_options.delay_axis, $
      cable_length_axis = plot_2d_options.cable_length_axis

  endif

  if plot_types.plot_kperp_power then begin

    file_arr = reform(savefiles_kperp_1d, n_cubes*n_elements(kperp_density_names))

    titles_use = strarr(n_cubes, n_elements(kperp_density_names))
    for j=0, n_elements(kperp_density_names)-1 do begin
      titles_use[*,j] = titles + kperp_density_names[j]
    endfor
    titles_use = reform(titles_use, n_cubes*n_elements(kperp_density_names))

    if plot_1d_options.plot_eor_1d then begin
      if count_catalog gt 0 then begin
        psyms = [intarr(n_elements(file_arr))+10, -3]
        file_arr = [file_arr, eor_file_1d]
        titles_use = [titles_use, 'EoR signal']
      endif else print, 'Could not locate catalog_data directory in !path variable'
    endif
    if plot_1d_options.plot_flat_1d then begin
      if count_catalog gt 0 then begin
        psyms = [intarr(n_elements(file_arr))+10, -3]
        file_arr = [file_arr, flat_file_1d]
        titles_use = [titles_use, 'input flat power']
      endif else print, 'Could not locate catalog_data directory in !path variable'
    endif

    if tag_exist(plot_2d_options, 'kpar_plot_range') then begin
      k_range = minmax([plot_2d_options.kperp_plot_range, binning_2d_options.kpar_bin, $
        plot_2d_options.kpar_plot_range[1]])
    endif else begin
      k_range = minmax([plot_2d_options.kperp_plot_range, binning_2d_options.kpar_bin])
    endelse

    window_num = window_num+1
    kpower_1d_plots, file_arr, window_num = window_num, colors = colors, $
      plot_options = plot_options, plot_1d_options = plot_1d_options, $
      names = titles_use, psyms = psyms, plotfile = plotfile_kperp_power, $
      k_range = k_range, title = plot_options.note + ' kperp', note = note_kperp_1d, $
      plot_sim_noise = plot_types.plot_sim_noise, $
      /kperp_power, baseline_axis = plot_2d_options.baseline_axis

  endif

  if plot_types.plot_k0_power then begin

    file_arr = savefiles_k0

    titles_use = strarr(n_cubes, n_elements(kperp_density_names))
    for j=0, n_elements(kperp_density_names)-1 do begin
      titles_use[*,j] = titles + kperp_density_names[j]
    endfor
    titles_use = reform(titles_use, n_cubes*n_elements(kperp_density_names))

    k_range = minmax([plot_2d_options.kperp_plot_range, binning_2d_options.kperp_bin])

    window_num = window_num+1
    kpower_1d_plots, file_arr, window_num = window_num, colors = colors, $
      plot_options = plot_options, plot_1d_options = plot_1d_options, $
      names = titles_use, plotfile = plotfile_k0_power, k_range = k_range, $
      title = plot_options.note + ' kpar=0', note = note_1d, $
      plot_sim_noise = plot_types.plot_sim_noise, $
      /kperp_power, baseline_axis = plot_2d_options.baseline_axis

  endif

  if plot_types.plot_1to2d then begin


    if plot_2d_options.kperp_linear_axis then begin
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
      if i eq nrow*ncol-1 and tag_exist(plot_options, 'note') then begin
        note_use = plot_options.note + ',' + type_1to2d_use + '_' + pol_1to2d_use
      endif else begin
        note_use = ''
      endelse
      if plot_options.pub then plotfile_use = plotfile_1to2d_heatmap else undefine, plotfile_use

      wedge_index = i / n_elements(kperp_density_names)
      density_index = i mod n_elements(kperp_density_names)
      title_use = wedge_1dbin_names[wedge_index] + ' ' + kperp_density_names[density_index]

      kpower_2d_plots, savefiles_1to2d_bin[wh_2d_use, density_index, wedge_index], $
        multi_pos = pos_use, start_multi_params = start_multi_params, $
        plotfile = plotfile_use, $
        plot_options = plot_options, plot_2d_options = plot_2d_options, $
        /plot_bin, /bin_contour, contour_levels = 'all', $
        title_prefix = title_use, note = note_use, $
        wedge_amp = binning_1d_options.wedge_amps, window_num = window_num

      if i eq 0 then begin
        positions = pos_use
        undefine, start_multi_params
      endif
    endfor
    undefine, positions, pos_use
    if plot_options.pub then begin
      cgps_close, png = plot_options.png, pdf = plot_options.pdf, $
        delete_ps = plot_options.delete_ps, density = 600
    endif


    ;; decide whether to make a zoom plot
    if tag_exist(binning_1d_options, 'kperp_range_1dave') gt 0 or $
        tag_exist(binning_1d_options, 'kperp_range_lambda_1dave') gt 0 or $
        tag_exist(binning_1d_options, 'kpar_range_1dave') gt 0 then begin
      if tag_exist(binning_1d_options, 'kperp_range_1dave') gt 0 then begin
        if binning_1d_options.kperp_range_1dave[0] gt plot_2d_options.kperp_plot_range[0] or $
            binning_1d_options.kperp_range_1dave[1] lt plot_2d_options.kperp_plot_range[1] then begin
          zoom = 1
          kperp_range_zoom = binning_1d_options.kperp_range_1dave * [.7, 1.1]
          kperp_range_zoom = [max([plot_2d_options.kperp_plot_range[0], kperp_range_zoom[0]]), min([plot_2d_options.kperp_plot_range[1], kperp_range_zoom[1]])]
        endif else kperp_range_zoom = plot_2d_options.kperp_plot_range
      endif
      if tag_exist(binning_1d_options, 'kperp_range_lambda_1dave') gt 0 then begin
        if plot_options.hinv then kperp_range_zoom = binning_1d_options.kperp_range_lambda_1dave / (kperp_lambda_conv * hubble_param) else $
          kperp_range_zoom = binning_1d_options.kperp_range_lambda_1dave / kperp_lambda_conv

        if kperp_range_zoom[0] gt plot_2d_options.kperp_plot_range[0] or kperp_range_zoom[1] lt plot_2d_options.kperp_plot_range[1] then begin
          zoom = 1
          kperp_range_zoom = kperp_range_zoom * [.7, 1.1]
          kperp_range_zoom = [max([plot_2d_options.kperp_plot_range[0], kperp_range_zoom[0]]), min([plot_2d_options.kperp_plot_range[1], kperp_range_zoom[1]])]
        endif else kperp_range_zoom = plot_2d_options.kperp_plot_range
      endif
      if tag_exist(binning_1d_options, 'kpar_range_1dave') gt 0 then begin
        if binning_1d_options.kpar_range_1dave[0] gt plot_2d_options.kpar_plot_range[0] or $
            binning_1d_options.kpar_range_1dave[1] lt plot_2d_options.kpar_plot_range[1] then begin
          zoom = 1
          kpar_range_zoom = binning_1d_options.kpar_range_1dave * [.7, 1.1]
          kpar_range_zoom = [max([plot_2d_options.kpar_plot_range[0], kpar_range_zoom[0]]), $
            min([plot_2d_options.kpar_plot_range[1], kpar_range_zoom[1]])]
        endif else kpar_range_zoom = plot_2d_options.kpar_plot_range
      endif
    endif

    start_multi_params = {ncol:ncol, nrow:nrow, ordering:'col'}

    window_num = window_num+1
    for i=0, nrow*ncol-1 do begin
      if i gt 0 then pos_use = positions[*,i]
      if i eq nrow*ncol-1 and tag_exist(plot_options, 'note') then begin
        note_use = plot_options.note + ',' + type_1to2d_use + '_' + pol_1to2d_use
      endif else begin
        note_use = ''
     endelse
      if plot_options.pub then begin
        plotfile_use = plotfile_1to2d_contours
      endif else begin
        undefine, plotfile_use
      endelse

      wedge_index = i / n_elements(kperp_density_names)
      density_index = i mod n_elements(kperp_density_names)
      title_use = wedge_1dbin_names[wedge_index] + ' ' + kperp_density_names[density_index]

      if keyword_set(zoom) then levels = 1 else levels = 'all'
      ;levels = 'all'

      kpower_2d_plots, savefiles_2d[wh_2d_use, density_index], $
        multi_pos = pos_use, start_multi_params = start_multi_params, $
        plotfile = plotfile_use, $
        plot_options = plot_options, plot_2d_options = plot_2d_options, $
        bin_savefile = savefiles_1to2d_bin[wh_2d_use, density_index, wedge_index], $
        /bin_contour, contour_levels = levels, title_prefix = title_use, note = note_use, $
        plot_wedge_line = 0, wedge_amp = binning_1d_options.wedge_amps, window_num = window_num

      if i eq 0 then begin
        positions = pos_use
        undefine, start_multi_params
      endif
    endfor
    undefine, positions, pos_use
    if plot_options.pub then begin
      cgps_close, png = plot_options.png, pdf = plot_options.pdf, $
        delete_ps = plot_options.delete_ps, density = 600
    endif

    start_multi_params = {ncol:ncol, nrow:nrow, ordering:'col'}

    window_num = window_num+1
    for i=0, nrow*ncol-1 do begin
      if i gt 0 then pos_use = positions[*,i]
      if i eq nrow*ncol-1 and tag_exist(plot_options, 'note') then begin
        note_use = plot_options.note + ',' + type_1to2d_use + '_' + pol_1to2d_use
      endif else begin
        note_use = ''
      endelse
      if plot_options.pub then begin
        plotfile_use = plotfile_1to2d_noisefrac
      endif else begin
        undefine, plotfile_use
      endelse

      wedge_index = i / n_elements(kperp_density_names)
      density_index = i mod n_elements(kperp_density_names)
      title_use = wedge_1dbin_names[wedge_index] + ' ' + kperp_density_names[density_index]

      levels = 'all'

      kpower_2d_plots, savefiles_1to2d_bin[wh_2d_use, density_index, wedge_index], $
        multi_pos = pos_use, start_multi_params = start_multi_params, $
        plotfile = plotfile_use, $
        plot_options = plot_options, plot_2d_options = plot_2d_options, $
        bin_savefile = savefiles_1to2d_bin[wh_2d_use, density_index, wedge_index], $
        /bin_contour, contour_levels = levels, /plot_1d_noisefrac, $
        data_range = [0,1], title_prefix = title_use, note = note_use, $
        plot_wedge_line = 0, wedge_amp = binning_1d_options.wedge_amps, window_num = window_num

      if i eq 0 then begin
        positions = pos_use
        undefine, start_multi_params
      endif
    endfor
    undefine, positions, pos_use
    if plot_options.pub then begin
      cgps_close, png = plot_options.png, pdf = plot_options.pdf, $
        delete_ps = plot_options.delete_ps, density = 600
    endif

    if keyword_set(zoom) then begin
      start_multi_params = {ncol:ncol, nrow:nrow, ordering:'col'}

      window_num = window_num+1
      for i=0, nrow*ncol-1 do begin
        if i gt 0 then pos_use = positions[*,i]
        if i eq nrow*ncol-1 and tag_exist(plot_options, 'note') then begin
          note_use = plot_options.note + ',' + type_1to2d_use + '_' + pol_1to2d_use
        endif else begin
          note_use = ''
        endelse
        if plot_options.pub then begin
          plotfile_use = plotfile_1to2d_contour_zoom
        endif else begin
          undefine, plotfile_use
        endelse

        wedge_index = i / n_elements(kperp_density_names)
        density_index = i mod n_elements(kperp_density_names)
        title_use = wedge_1dbin_names[wedge_index] + ' ' + kperp_density_names[density_index]

        kpower_2d_plots, savefiles_2d[wh_2d_use, density_index], multi_pos = pos_use, $
          start_multi_params = start_multi_params, $
          plotfile = plotfile_use, $
          plot_options = plot_options, plot_2d_options = plot_2d_options, $
          kperp_plot_range = kperp_range_zoom, kpar_plot_range = kpar_range_zoom, $
          bin_savefile = savefiles_1to2d_bin[wh_2d_use, density_index, wedge_index], $
          /bin_contour, contour_levels = 'all', title_prefix = title_use, note = note_use, $
          plot_wedge_line = 0, wedge_amp = binning_1d_options.wedge_amps, window_num = window_num

        if i eq 0 then begin
          positions = pos_use
          undefine, start_multi_params
        endif
      endfor
      undefine, positions, pos_use
      if plot_options.pub then begin
        cgps_close, png = plot_options.png, pdf = plot_options.pdf, $
          delete_ps = plot_options.delete_ps, density = 600
      endif
    endif
  endif


  if plot_types.plot_2d_masked then begin

    for i=0, n_elements(wedge_1dbin_names)-1 do begin
      for j=0, n_elements(kperp_density_names)-1 do begin
        note_use = plot_options.note + ', ' + kperp_density_names[j]

        if plot_options.pub then begin
          plotfiles_use = {power: plotfiles_2d_masked[*,j,i], $
            error: plotfiles_2d_error_masked[*,j,i], $
            noise_expval: plotfiles_2d_noise_expval_masked[*,j,i], $
            snr: plotfiles_2d_snr_masked[*,j,i], noise: plotfiles_2d_noise_masked[*,j,i], $
            sim_noise: plotfiles_2d_sim_noise_masked[*,j,i], $
            sim_noise_diff: plotfiles_2d_sim_noise_diff_masked[*,j,i], $
            nnr: plotfiles_2d_nnr_masked[*,j,i], sim_nnr: plotfiles_2d_sim_nnr_masked[*,j,i]}
        endif

        make_std_plots, file_struct_arr, savefiles_2d_masked[*,j,i], titles, $
          plotfile_struct = plotfiles_use, plot_sim_noise = plot_types.plot_sim_noise, $
          window_num = window_num, vs_note = vs_note, wedge_amp = binning_1d_options.wedge_amps, $
          plot_options = plot_options, plot_2d_options = plot_2d_options
      endfor
    endfor
  endif

end
