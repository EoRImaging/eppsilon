pro single_cube_dft, folder_name_in, obs_name, data_subdirs=data_subdirs, $
    exact_obsnames = exact_obsnames, ps_foldername = ps_foldername, $
    rts = rts, dirty_folder = dirty_folder, save_path = save_path, $
    savefilebase = savefilebase, refresh_info = refresh_info, refresh_dft = refresh_dft, $
    freq_flags = freq_flags, freq_ch_range = freq_ch_range, cube_type = cube_type, $
    pol = pol, evenodd = evenodd, loc_name=loc_name, $
    delta_uv_lambda = delta_uv_lambda, max_uv_lambda = max_uv_lambda, $
    full_image = full_image

  if n_elements(folder_name_in) ne 1 then begin
    message, 'one folder_name must be supplied.'
  endif
  if n_elements(obs_name) gt 1 then begin
    message, 'no more than one obs_name can be supplied.'
  endif

  if n_elements(cube_type) eq 0 then begin
    message, 'cube_type must be specified'
  endif
  if cube_type eq 'variances' then begin
    message, 'weights and variances must be calculated together. ' + $
      'Use cube_type = weights to calculate both.'
  endif

  if n_elements(pol) eq 0 then begin
    message, 'pol must be specified'
  endif
  if n_elements(evenodd) eq 0 then begin
    message, 'evenodd must be specified'
  endif

  folder_name = get_folder(folder_name_in, loc_name = loc_name,  rts = rts, $
    dirty_folder = dirty_folder)

  obs_info = ps_filenames(folder_name, obs_name, dirty_folder = dirty_folder, $
    exact_obsnames = exact_obsnames, rts = rts, $
    data_subdirs = data_subdirs, ps_foldernames = ps_foldername, save_paths = save_path, $
    refresh_info = refresh_info)

  if not file_test(save_path, /directory) then file_mkdir, save_path

  if keyword_set(rts) then begin

    if tag_exist(obs_info, 'dirtyfiles') then dirtyfiles = obs_info.dirtyfiles.(0)

    if obs_info.info_files[0] ne '' then datafile = obs_info.info_files[0] else begin

      if obs_info.cube_files.(0)[0] ne '' then datafile = obs_info.cube_files.(0) else begin
        datafile = rts_fits2idlcube(obs_info.datafiles.(0), obs_info.weightfiles.(0), $
          obs_info.variancefiles.(0), dirtyfiles = dirtyfiles, pol_inc = pol_inc, $
          save_path = obs_info.folder_names[0]+path_sep(), refresh = refresh_dft, $
          no_wtvar = no_wtvar_rts)
      endelse
    endelse

    if keyword_set(refresh_rtscube) then begin
      datafile = rts_fits2idlcube(obs_info.datafiles.(0), obs_info.weightfiles.(0), $
        obs_info.variancefiles.(0), dirtyfiles = dirtyfiles, pol_inc = pol_inc, $
        save_path = obs_info.folder_names[0]+path_sep(), $
        /refresh, no_wtvar = no_wtvar_rts)
    endif

    if keyword_set(no_wtvar_rts) then return

    max_uv_lambda = 300

  endif else begin

    if obs_info.info_files[0] ne '' then begin
      datafile = obs_info.info_files[0]
    endif else begin
      datafile = obs_info.cube_files.(0)
    endelse

  endelse


  refresh_options = create_refresh_options(refresh_dft = refresh_dft)

  uvf_options = create_uvf_options(delta_uv_lambda = delta_uv_lambda, $
    max_uv_lambda = max_uv_lambda, full_image = full_image, $
    dft_fchunk = dft_fchunk, no_dft_progress = no_dft_progress)

  ;; it doesn't matter what's in ps_options because we don't use any of the filenames
  ;; that depend on this
  ps_options = create_ps_options()

  if keyword_set(rts) then begin
    file_struct_arr = rts_file_setup(datafile, savefilebase = savefilebase, $
      save_path = save_path,use_fhd_norm = norm_rts_with_fhd, refresh_info = refresh_info, $
      uvf_options = uvf_options, ps_options = ps_options)
  endif else if keyword_set(casa) then begin
    file_struct_arr = casa_file_setup(datafile, savefilebase = savefilebase, $
      save_path = save_path, refresh_info = refresh_info, $
      uvf_options = uvf_options, ps_options = ps_options)
  endif else begin
    file_struct_arr = fhd_file_setup(datafile, savefilebase = savefilebase, $
      save_path = save_path, freq_ch_range = freq_ch_range, freq_flags = freq_flags, $
      freq_flag_name = freq_flag_name, sim = sim, refresh_info = refresh_info, $
      uvf_options = uvf_options, ps_options = ps_options)
  endelse

  type_pol_str_list = file_struct_arr.type_pol_str
  this_type_pol_str = cube_type + '_' + pol

  if cube_type eq 'weights' then begin
    cube_ind = where(file_struct_arr.pol eq pol, count_pol)
    if count_pol eq 0 then begin
      message, 'Specified pol is not available'
    endif
    if count_pol gt 0 then cube_ind = cube_ind[0]

    evenodd_mask = stregex(file_struct_arr[cube_ind].weightfile, evenodd, /boolean)
    evenodd_ind = where(evenodd_mask eq 1, count_evenodd)
    if count_evenodd eq 0 then begin
      message, 'specified evenodd value is not available'
    endif
    if count_evenodd gt 1 then begin
      message, 'More than one matching evenodd value'
    endif

    file_struct = file_struct_arr[cube_ind]
    filename = file_struct_arr[cube_ind].weightfile[evenodd_ind]

    weight_varname = file_struct_arr[cube_ind].weightvar
    variance_varname = file_struct_arr[cube_ind].variancevar
    uvf_savefile = file_struct_arr[cube_ind].uvf_weight_savefile[evenodd_ind]
  endif else begin
    cube_ind = where(type_pol_str_list eq this_type_pol_str, count_typepol)
    if count_typepol eq 0 then begin
      message, 'Specified cube type and pol are not available.'
    endif
    if count_typepol gt 1 then begin
      message, 'More than one matching type and pol'
    endif

    evenodd_mask = stregex(file_struct_arr[cube_ind].datafile, evenodd, /boolean)
    evenodd_ind = where(evenodd_mask eq 1, count_evenodd)
    if count_evenodd eq 0 then begin
      message, 'specified evenodd value is not available'
    endif
    if count_evenodd gt 1 then begin
      message, 'More than one matching evenodd value'
    endif

    file_struct = file_struct_arr[cube_ind]
    filename = file_struct_arr[cube_ind].datafile[evenodd_ind]
    uvf_savefile = file_struct_arr[cube_ind].uvf_savefile[evenodd_ind]

    void = getvar_savefile(filename, names=varnames)
    wh_name = where(varnames eq file_struct.datavar, count_name)
    if count_name eq 0 then begin
      message, 'specified cube_type and pol is a derived cube, it does not have a healpix cube'
    endif
  endelse

  test_uvf = file_test(uvf_savefile) *  (1 - file_test(uvf_savefile, /zero_length))
  if test_uvf eq 1 and n_elements(freq_flags) ne 0 then begin
    old_freq_mask = getvar_savefile(file_struct.uvf_savefile[i], 'freq_mask')
    if n_elements(freq_ch_range) ne 0 then begin
      if total(abs(old_freq_mask[min(freq_ch_range):max(freq_ch_range)] - freq_mask)) ne 0 $
        then test_uvf = 0
    endif else begin
      if total(abs(old_freq_mask - freq_mask)) ne 0 then test_uvf = 0
    endelse
  endif

  if test_uvf eq 0 or keyword_set(refresh_dft) then begin

    if tag_exist(file_struct, 'no_var') ne 0 then no_var = 1 else no_var = 0
    if tag_exist(file_struct, 'nside') ne 0 then healpix = 1 else healpix = 0

    frequencies = file_struct.frequencies

    if n_elements(freq_flags) ne 0 then freq_mask = file_struct.freq_mask

    if n_elements(freq_ch_range) ne 0 then begin
      n_freq_orig = n_elements(frequencies)
      frequencies = frequencies[min(freq_ch_range):max(freq_ch_range)]
    endif
    n_freq = n_elements(frequencies)

    pixel_nums = getvar_savefile(file_struct.pixelfile[0], file_struct.pixelvar[0])

    pix_ft_struct = choose_pix_ft(file_struct, pixel_nums = pixel_nums, $
      data_dims = data_dims, uvf_options = uvf_options)

    wh_close = pix_ft_struct.wh_close
    x_use = pix_ft_struct.x_use
    y_use = pix_ft_struct.y_use
    kx_rad_vals = pix_ft_struct.kx_rad_vals
    ky_rad_vals = pix_ft_struct.ky_rad_vals

    pix_window = fltarr(n_elements(wh_close), n_freq) + 1.

    git, repo_path = ps_repository_dir(), result = this_run_git_hash

    uvf_git_hash = this_run_git_hash
    uvf_wt_git_hash = this_run_git_hash

    if cube_type eq 'weights' then begin
      print, 'calculating DFT for ' + file_struct.weightvar + ' in ' + filename

      time0 = systime(1)
      arr = getvar_savefile(filename, file_struct.weightvar)
      time1 = systime(1)

      if time1 - time0 gt 60 then begin
        print, 'Time to restore cube was ' + number_formatter((time1 - time0)/60.) + ' minutes'
      endif

      if size(arr,/type) eq 6 or size(arr,/type) eq 9 then begin
        if max(abs(imaginary(arr))) eq 0 then begin
          arr = real_part(arr)
        endif else begin
          message, 'Weights in Healpix cubes is complex.'
        endelse
      endif
      if not healpix then arr = reform(arr, dims[0]*dims[1], n_freq)

      if n_elements(wh_close) ne n_elements(pixel_nums) then begin
        arr = arr[wh_close, *] * pix_window
      endif
      if n_elements(freq_ch_range) ne 0 then begin
        arr = arr[*, min(freq_ch_range):max(freq_ch_range)]
      endif
      if n_elements(freq_flags) ne 0 then begin
        arr = arr * rebin(reform(freq_mask, 1, n_freq), $
          size(arr, /dimension), /sample)
      endif

      transform = discrete_ft_2D_fast(x_use, y_use, arr, kx_rad_vals, ky_rad_vals, $
        timing = ft_time, fchunk = dft_fchunk, no_progress = no_dft_progress)

      weights_cube = temporary(transform)

      if not no_var then begin
        print, 'calculating DFT for ' + file_struct.variancevar + ' in ' + filename

        time0 = systime(1)
        arr = getvar_savefile(filename, file_struct.variancevar)
        time1 = systime(1)

        if time1 - time0 gt 60 then begin
          print, 'Time to restore cube was ' + number_formatter((time1 - time0)/60.) + ' minutes'
        endif

        if size(arr,/type) eq 6 or size(arr,/type) eq 9 then begin
          if max(abs(imaginary(arr))) eq 0 then begin
            arr = real_part(arr)
          endif else begin
            message, 'Variances in Healpix cubes is complex.'
          endelse
        endif
        if not healpix then arr = reform(arr, dims[0]*dims[1], n_freq)

        if n_elements(wh_close) ne n_elements(pixel_nums) then begin
          arr = arr[wh_close, *] * pix_window^2.
        endif
        if n_elements(freq_ch_range) ne 0 then begin
          arr = arr[*, min(freq_ch_range):max(freq_ch_range)]
        endif
        if n_elements(freq_flags) ne 0 then begin
          arr = arr * rebin(reform(freq_mask, 1, n_freq), $
            size(arr, /dimension), /sample)
        endif

        transform = discrete_ft_2D_fast(x_use, y_use, arr, kx_rad_vals, ky_rad_vals, $
          timing = ft_time, fchunk = dft_fchunk, no_progress = no_dft_progress)

        ;; make variances real, positive definite (amplitude)
        variance_cube = abs(temporary(transform))
        undefine, arr
      endif

      if n_elements(freq_flags) then begin
        save, file = uvf_savefile, kx_rad_vals, ky_rad_vals, weights_cube, $
          variance_cube, freq_mask, uvf_wt_git_hash
      endif else begin
        save, file = uvf_savefile, kx_rad_vals, ky_rad_vals, weights_cube, $
          variance_cube, uvf_wt_git_hash
      endelse
      undefine, new_pix_vec, weights_cube, variance_cube


    endif else begin
      print, 'calculating DFT for ' + file_struct.datavar + ' in ' + filename

      time0 = systime(1)
      arr = getvar_savefile(filename, file_struct.datavar)
      time1 = systime(1)

      if time1 - time0 gt 60 then begin
        print, 'Time to restore cube was ' + number_formatter((time1 - time0)/60.) + ' minutes'
      endif

      if size(arr,/type) eq 6 or size(arr,/type) eq 9 then begin
        if max(abs(imaginary(arr))) eq 0 then begin
          arr = real_part(arr)
        endif else begin
          message, 'Data in Healpix cubes is complex.'
        endelse
      endif
      if not healpix then arr = reform(arr, dims[0]*dims[1], n_freq)

      if n_elements(wh_close) ne n_elements(pixel_nums) then begin
        arr = arr[wh_close, *] * pix_window
      endif
      if n_elements(freq_ch_range) ne 0 then begin
        arr = arr[*, min(freq_ch_range):max(freq_ch_range)]
      endif
      if n_elements(freq_flags) ne 0 then begin
        arr = arr * rebin(reform(freq_mask, 1, n_freq), $
        size(arr, /dimension), /sample)
      endif

      transform = discrete_ft_2D_fast(x_use, y_use, arr, kx_rad_vals, ky_rad_vals, $
        timing = ft_time, fchunk = dft_fchunk, no_progress = no_dft_progress)

      data_cube = temporary(transform)
      undefine, arr

      if n_elements(freq_flags) then begin
        save, file = uvf_savefile, kx_rad_vals, ky_rad_vals, data_cube, freq_mask, uvf_git_hash
      endif else begin
        save, file = uvf_savefile, kx_rad_vals, ky_rad_vals, data_cube, uvf_git_hash
      endelse
      undefine, data_cube
    endelse
  endif
end
