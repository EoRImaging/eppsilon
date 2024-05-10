pro ps_power_incoherent_avg, info_files, avg_file_struct=avg_file_struct, $
  pol_inc = pol_inc, freq_options = freq_options, uvf_input = uvf_input, $
  new_obsname=new_obsname, avg_save_path = avg_save_path, $
  savefilebase = savefilebase,  save_slices = save_slices, $
  refresh_options = refresh_options, uvf_options = uvf_options, $
  ps_options = ps_options

  compile_opt idl2, logical_predicate

  ; just assume FHD.
  n_filesets = n_elements(info_files)
  n_cubes_arr = intarr(n_filesets)

  for fi=0, n_elements(info_files)-1 do begin
    this_file_struct_arr = fhd_file_setup(info_files[fi], sim = sim, $
      uvf_input = uvf_input, savefilebase = savefilebase, save_path = save_path, $
      refresh_info = refresh_options.refresh_info, uvf_options = uvf_options, $
      freq_options = freq_options, ps_options = ps_options, pol_inc=pol_inc)

    n_cubes_arr[fi] = n_elements(this_file_struct_arr)
    power_files = this_file_struct_arr.power_savefile
    test_power_3d = file_valid(power_files)

    if min(test_power_3d) eq 1 and tag_exist(freq_options, 'freq_flags') then begin
      freq_mask = freq_options.freq_mask
      for ci=0, n_cubes_arr[fi]-1 do begin
        old_freq_mask = getvar_savefile(power_files[ci], 'freq_mask')
        if total(abs(old_freq_mask - freq_mask)) ne 0 then test_power_3d[ci] = 0
      endfor
    endif
    if min(test_power_3d) eq 0 then message, "power files do not all exist. " $
      + "Missing files are:" + strjoin(power_files[where(test_power_3d eq 0)], ", ")

    if fi eq 0 then begin
      file_struct_arr = this_file_struct_arr
    endif else begin
      if n_cubes_arr[fi] ne n_cubes_arr[0] then message, "different number of cubes in input ps runs"
      file_struct_arr = [[file_struct_arr], [this_file_struct_arr]]
    endelse
  endfor

  ; would have errored earlier if n_cubes varies
  n_cubes = n_cubes_arr[0]
  struct_arr_dims = size(file_struct_arr,/dim)
  if struct_arr_dims[0] ne n_cubes or struct_arr_dims[1] ne n_filesets then message, "a problem happened with concatenating file structs"

  if n_elements(savefilebase) eq 0 then begin
    savefilebase = new_obsname
  endif

  file_pattern = strarr(n_cubes, n_filesets)
  for ci=0, n_cubes-1 do begin
    for fi=0, n_filesets-1 do begin
      savefilebase_parts = strsplit(file_struct_arr[ci, fi].savefilebase, "__", /extract, /regex, count=parts_count)
      file_pattern[ci, fi] = savefilebase_parts[-1]
    endfor
  endfor

  diff_parts = get_strarray_same_diff(file_pattern, "_", same_str=file_pattern_use)
  diff_parts = reform(diff_parts, n_cubes, n_filesets)

  data_full_path = avg_save_path + file_struct_arr[0].subfolders.data
  data_kspace_path = data_full_path + file_struct_arr[0].subfolders.kspace

  incoh_avg_power_general_filebase = savefilebase + "__" + file_pattern_use
  incoh_avg_power_savefilebase = incoh_avg_power_general_filebase + '_' $
    + diff_parts[*, 0]
  incopherent_avg_power_savefiles = data_kspace_path + incoh_avg_power_savefilebase $
    + file_struct_arr[0].power_tag + '_incoherent_avg_power.idlsave'

  slice_full_path = data_kspace_path + file_struct_arr[0].subfolders.slices
  xy_slice_savefiles = slice_full_path + incoh_avg_power_savefilebase $
    + file_struct_arr[0].power_tag + '_incoherent_avg_xy_plane.idlsave'
  xz_slice_savefiles = slice_full_path + incoh_avg_power_savefilebase $
    + file_struct_arr[0].power_tag + '_incoherent_avg_xz_plane.idlsave'
  yz_slice_savefiles = slice_full_path + incoh_avg_power_savefilebase $
    + file_struct_arr[0].power_tag + '_incoherent_avg_yz_plane.idlsave'

  if not file_test(data_full_path, /directory) then file_mkdir, data_full_path
  if not file_test(data_kspace_path, /directory) then file_mkdir, data_kspace_path
  if not file_test(slice_full_path, /directory) then file_mkdir, slice_full_path
; stop
  ; Now do weighted average
  for ci=0, n_cubes-1 do begin
    n_vis_arr = reform(file_struct_arr[0, *].n_vis)
    n_vis_freq_arr = reform(file_struct_arr[0, *].n_vis_freq)
    this_avg_file_struct = create_struct('savefile_froot', avg_save_path, $
      'n_vis', total(n_vis_arr, 2), $
      'n_vis_freq', total(n_vis_freq_arr, 3), $
      'max_baseline_lambda', max(file_struct_arr[0, *].max_baseline_lambda), $
      'general_filebase', incoh_avg_power_general_filebase, $
      'savefilebase', incoh_avg_power_savefilebase[ci], $
      'power_savefile', incopherent_avg_power_savefiles[ci], $
      'xy_savefile', xy_slice_savefiles, $
      'xz_savefile', xz_slice_savefiles, $
      'yz_savefile', yz_slice_savefiles)

    fields_copy_f0 = ['subfolders', 'power_tag', 'uvf_tag', 'freq_tag', 'nfiles', $
      'pol_index', 'type_index', 'pol', 'type', 'frequencies', 'instrument', $
      'kspan', 'file_label']
    for fdi=0, n_elements(fields_copy_f0)-1 do begin
      fd_name = fields_copy_f0[fdi]
      f0_ind = where(strlowcase(tag_names(file_struct_arr[ci, 0])) eq fd_name, count)
      if count eq 0 then message, "field "+ fd_name + " not found on file_struct_arr"
      this_avg_file_struct = create_struct(this_avg_file_struct, fd_name, file_struct_arr[ci, 0].(f0_ind))
    endfor

    if ci eq 0 then begin
      avg_file_struct = this_avg_file_struct
    endif else begin
      avg_file_struct = [avg_file_struct, this_avg_file_struct]
    endelse
    
    for fi=0, n_filesets-1 do begin

      restore, file_struct_arr[ci, fi].power_savefile

      git, repo_path = ps_repository_dir(), result=ps_git_hash

      if fi eq 0 then begin
        kx_mc_use = kx_mpc
        ky_mpc_use = ky_mpc
        kz_mpc_use = kz_mpc
        kperp_lambda_conv_use = kperp_lambda_conv
        delay_params_use = delay_params
        hubble_param_use = hubble_param
        vs_name_use = vs_name
        
        power_3d_num = power_3d * weights_3d
        noise_3d_num = noise_3d * diff_weights_3d
        sim_noise_3d_num = sim_noise_3d * weights_3d
        sim_noise_diff_3d_num = sim_noise_diff_3d * diff_weights_3d
        ; not sure about this one:
        noise_expval_3d_num = noise_expval_3d * diff_weights_3d
        weights_3d_use = weights_3d
        diff_weights_3d_use = diff_weights_3d

        wt_meas_ave_use = wt_meas_ave
        wt_meas_min_use = wt_meas_min
        ave_weights_use = ave_weights
        wt_ave_power_freq_use = wt_ave_power_freq
        ave_power_freq_use = ave_power_freq

        n_freq_contrib_use = n_freq_contrib
        vs_mean_use = vs_mean
        t_sys_meas_use = t_sys_meas
        window_int_use = window_int
        wt_ave_power_uvf_use = wt_ave_power_uvf
        ave_power_uvf_use = ave_power_uvf

        git_hashes_use = git_hashes
        
      endif else begin
        if max(abs(kx_mpc - kx_mc_use)) ne 0 then message, "kx_mpc do not match"
        if max(abs(ky_mpc - ky_mpc_use)) ne 0 then message, "ky_mpc do not match"
        if max(abs(kz_mpc - kz_mpc_use)) ne 0 then message, "kz_mpc do not match"

        if kperp_lambda_conv ne kperp_lambda_conv_use then message, "kperp_lambda_conv do not match"
        if max(abs(delay_params - delay_params_use)) ne 0then message, "delay_params do not match"
        if hubble_param ne hubble_param_use then message, "hubble_param do not match"
        if vs_name ne vs_name_use then message, "vs_name do not match"

        power_3d_num += power_3d * weights_3d
        noise_3d_num += noise_3d * diff_weights_3d
        sim_noise_3d_num += sim_noise_3d * weights_3d
        sim_noise_diff_3d_num += sim_noise_diff_3d * diff_weights_3d
        ; not sure about this one:
        noise_expval_3d_num += noise_expval_3d * diff_weights_3d
        weights_3d_use += weights_3d
        diff_weights_3d_use += diff_weights_3d

        wt_meas_ave_use = [[[wt_meas_ave_use]], [[wt_meas_ave]]]
        wt_meas_min_use = [[[wt_meas_min_use]], [[wt_meas_min]]]
        ave_weights_use = [[[ave_weights_use]], [[ave_weights]]]
        wt_ave_power_freq_use = [[[wt_ave_power_freq_use]], [[wt_ave_power_freq]]]
        ave_power_freq_use = [[[ave_power_freq_use]], [[ave_power_freq]]]

        n_freq_contrib_use = [[[n_freq_contrib_use]], [[n_freq_contrib]]]
        vs_mean_use = [vs_mean_use, vs_mean]
        t_sys_meas_use = [[t_sys_meas_use], [t_sys_meas]]
        window_int_use = [[window_int_use], [window_int]]
        git_hashes_use = [git_hashes_use, git_hashes]
        wt_ave_power_uvf_use = [[wt_ave_power_uvf_use],[wt_ave_power_uvf]]
        ave_power_uvf_use = [[ave_power_uvf_use], [ave_power_uvf]]

      endelse

    endfor

    ; now divide to get the variance weighted values
    power_3d = power_3d_num / weights_3d_use
    noise_3d = noise_3d_num / diff_weights_3d_use
    sim_noise_3d = sim_noise_3d_num / weights_3d_use
    sim_noise_diff_3d = sim_noise_diff_3d_num / diff_weights_3d_use
    noise_expval_3d = noise_expval_3d_num / diff_weights_3d_use

    wh_wt0 = where(weights_3d_use eq 0, count_wt0)
    if count_wt0 ne 0 then begin
      power_3d[wh_wt0] = 0
      noise_expval_3d[wh_wt0] = 0
      noise_3d[wh_wt0] = 0
      sim_noise_3d[wh_wt0] = 0
      sim_noise_diff_3d[wh_wt0] = 0
    endif

    weights_3d = weights_3d_use
    diff_weights_3d = diff_weights_3d_use

    ; fix variable names for the save files
    wt_meas_ave = mean(wt_meas_ave_use, dimension=3)
    wt_meas_min = min(wt_meas_min_use, dimension=3)
    ave_weights = mean(ave_weights_use, dimension=3)
    wt_ave_power_freq = mean(wt_ave_power_freq_use, dimension=3)
    ave_power_freq = mean(ave_power_freq_use, dimension=3)

    n_freq_contrib = total(n_freq_contrib_use, 3)
    vs_mean = vs_mean_use
    t_sys_meas = t_sys_meas_use
    window_int = window_int_use

    git_hashes = create_struct('incoherent_avg', ps_git_hash, $
      'input_hashes', git_hashes_use)

    wt_ave_power_uvf = mean(wt_ave_power_uvf_use, dimension=2)
    ave_power_uvf = mean(ave_power_uvf_use, dimension=2)

    save, file = incopherent_avg_power_savefiles[ci], power_3d, noise_3d, sim_noise_3d, $
      sim_noise_diff_3d, noise_expval_3d, weights_3d, diff_weights_3d, $
      kx_mpc, ky_mpc, kz_mpc, kperp_lambda_conv, delay_params, hubble_param, $
      n_freq_contrib, freq_mask, vs_name, vs_mean, t_sys_meas, window_int, $
      git_hashes, wt_meas_ave, wt_meas_min, ave_weights, wt_ave_power_freq, $
      ave_power_freq, wt_ave_power_uvf, ave_power_uvf

    ; write_ps_fits, file_struct.fits_power_savefile, power_3d, weights_3d, $
    ;   noise_expval_3d, noise_3d = noise_3d, kx_mpc, ky_mpc, kz_mpc, $
    ;   kperp_lambda_conv, delay_params, hubble_param

    if save_slices then begin
      ;; now do slices
      make_slices, this_avg_file_struct, type='kspace', data_cube = power_3d, $
        weights_cube = weights_3d, noise_3d = noise_3d, noise_expval_3d = noise_expval_3d, $
        kx_mpc = kx_mpc, ky_mpc = ky_mpc, kz_mpc = kz_mpc, kperp_lambda_conv = kperp_lambda_conv, $
        delay_params = delay_params, hubble_param = hubble_param
    endif
  
  endfor

end