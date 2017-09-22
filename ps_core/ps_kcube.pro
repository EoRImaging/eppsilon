pro ps_kcube, file_struct, dft_refresh_data = dft_refresh_data, dft_refresh_weight = dft_refresh_weight, $
    refresh_beam = refresh_beam, sim=sim, fix_sim_input = fix_sim_input, allow_beam_approx = allow_beam_approx, $
    dft_fchunk = dft_fchunk, freq_ch_range = freq_ch_range, freq_flags = freq_flags, $
    spec_window_type = spec_window_type, image_window_name = image_window_name, image_window_frac_size = image_window_frac_size, $
    delta_uv_lambda = delta_uv_lambda, max_uv_lambda = max_uv_lambda, $
    std_power = std_power, inverse_covar_weight = inverse_covar_weight, $
    input_units = input_units, uvf_input = uvf_input, $
    uv_avg = uv_avg, uv_img_clip = uv_img_clip, no_dft_progress = no_dft_progress, $
    ave_removal = ave_removal
    
  if tag_exist(file_struct, 'nside') ne 0 then healpix = 1 else healpix = 0
  if keyword_set(uvf_input) or tag_exist(file_struct, 'uvf_savefile') eq 0 then uvf_input = 1 else uvf_input = 0
  nfiles = n_elements(file_struct.datafile)
  if tag_exist(file_struct, 'no_var') ne 0 then no_var = 1 else no_var = 0
  
  if tag_exist(file_struct, 'beam_savefile') eq 0 then refresh_beam = 0
  
  if n_elements(input_units) eq 0 then input_units = 'jansky'
  units_enum = ['jansky', 'mk']
  wh = where(units_enum eq input_units, count)
  if count eq 0 then message, 'input units not recognized, options are: ' + units_enum
  
  git, repo_path = ps_repository_dir(), result=this_run_git_hash
  
  datavar = strupcase(file_struct.datavar)
  if datavar eq '' then begin
    ;; working with a 'derived' cube that is constructed from uvf_savefiles
    input_uvf_files = reform(file_struct.derived_uvf_inputfiles, nfiles, 2)
    input_uvf_varname = reform(file_struct.derived_uvf_varname, nfiles, 2)
    
    if healpix or not keyword_set(uvf_input) then begin
      input_uvf_wtfiles = file_struct.uvf_weight_savefile
    endif
  endif
  
  frequencies = file_struct.frequencies
  
  if n_elements(freq_flags) ne 0 then freq_mask = file_struct.freq_mask
  
  if n_elements(freq_ch_range) ne 0 then begin
    n_freq_orig = n_elements(frequencies)
    frequencies = frequencies[min(freq_ch_range):max(freq_ch_range)]
  endif
  n_freq = n_elements(frequencies)
  
  z_mpc_mean = z_mpc(frequencies, hubble_param = hubble_param, f_delta = f_delta, even_freq = even_freq, $
    redshifts = redshifts, comov_dist_los = comov_dist_los, z_mpc_delta = z_mpc_delta)
    
  kperp_lambda_conv = z_mpc_mean / (2.*!pi)
  delay_delta = 1e9/(n_freq*f_delta*1e6) ;; equivilent delay bin size for kparallel
  delay_max = delay_delta * n_freq/2.    ;; factor of 2 b/c of neg/positive
  delay_params = [delay_delta, delay_max]
  
  z_mpc_length = max(comov_dist_los) - min(comov_dist_los) + z_mpc_delta
  kz_mpc_range =  (2.*!pi) / (z_mpc_delta)
  kz_mpc_delta = (2.*!pi) / z_mpc_length
  kz_mpc_orig = findgen(round(kz_mpc_range / kz_mpc_delta)) * kz_mpc_delta - kz_mpc_range/2.
  if n_elements(n_kz) ne 0 then begin
    if n_elements(kz_mpc_orig) ne n_kz then message, 'something has gone wrong with kz_mpc calculation.'
  endif else n_kz = n_elements(kz_mpc_orig)
  
  if input_units eq 'jansky' then begin
       ;; converting from Jy (in u,v,f) to mK*str (10^-26 * c^2 * 10^3/ (2*f^2*kb))
    conv_factor = float((3e8)^2 / (2. * (frequencies*1e6)^2. * 1.38065))
    
    ;; convert from mk*str to mK*Mpc^2
    conv_factor = conv_factor * z_mpc_mean^2.
    
    ;; for adrian's weighting we have something in Jy/(wavelength^-1) = Jy * wavelength to start
    ;; converting from wavelength to Mpc uses 1/kperp_lambda_conv, so multiply that by Jy -> mK*Mpc^2 conversion
    conv_factor_adrian = conv_factor / kperp_lambda_conv
    
  endif else conv_factor = 1. + fltarr(n_freq)
  
  ;;t_sys = 440. ; K
  ;;t_sys = 280. * sqrt(2.)* ((1+redshifts)/7.5)^2.3 ;; from skew w/ stu + srt(2) for single pol -- Adam says wrong for Ian's normalization
  t_sys = 280. * ((1+redshifts)/7.5)^2.3 / sqrt(2.) ;; from skew w/ stu + srt(2) for single pol
  ;;eff_area = 16. ; m^2
  eff_area = 21. ; m^2 -- from Aaron's memo
  df = file_struct.freq_resolution ; Hz -- native visibility resolution NOT cube resolution
  tau = file_struct.time_resolution ; seconds
  vis_sigma = (2. * (1.38065e-23) * 1e26) * t_sys / (eff_area * sqrt(df * tau)) ;; in Jy
  vis_sigma = float(vis_sigma)
  
  old_vis_sigma = temporary(vis_sigma)
  
  if tag_exist(file_struct, 'vis_noise') then begin
    vis_sigma_ian = file_struct.vis_noise
    ;; do a straight average over even/odd of sigma because we just want the average noise (should actually be identical)
    if nfiles eq 2 then vis_sigma_ian = total(vis_sigma_ian, 1)/2.
  endif
  
  if n_elements(freq_ch_range) ne 0 then vis_sig_tag = number_formatter(384./n_freq_orig) else vis_sig_tag = number_formatter(384./n_freq)
  vis_sigma_file = file_dirname(file_struct.savefile_froot, /mark_directory) + 'vis_sigma/vis_sigma_measured' + vis_sig_tag + '.sav'
  if file_test(vis_sigma_file) then begin
    vis_sigma_adam = getvar_savefile(vis_sigma_file, 'vis_sigma')
    
    if n_elements(freq_ch_range) ne 0 then begin
      if n_elements(vis_sigma_adam) ne n_freq_orig then message, 'vis_sig file has incorrect number of frequency channels'
      vis_sigma_adam = vis_sigma_adam[min(freq_ch_range):max(freq_ch_range)]
    endif else if n_elements(vis_sigma_adam) ne n_freq then message, 'vis_sig file has incorrect number of frequency channels'
    
    wh_nan = where(finite(vis_sigma_adam) eq 0, count_nan)
    if count_nan gt 0 then vis_sigma_adam[wh_nan] = 0
  endif
  
  if n_elements(vis_sigma_ian) gt 0 then begin
    if max(vis_sigma_ian) gt 5. then begin
      if n_elements(freq_ch_range) ne 0 then vis_sigma_ian = vis_sigma_ian[min(freq_ch_range):max(freq_ch_range)]
      vis_sigma = vis_sigma_ian
      vs_name = 'ian'
    endif
  endif
  
  if n_elements(vis_sigma) eq 0 then begin
    if n_elements(vis_sigma_adam) gt 0 then begin
      ;; nothing in file struct, use file if available
      vis_sigma = vis_sigma_adam
      vs_name = 'adam_high'
    endif else begin
      ;; no vis_sigma information, make a flat vis_sigma
      vis_sigma = old_vis_sigma*0 + old_vis_sigma[0]
      vs_name = 'calc_flat'
    endelse
  endif
  
  ;vis_sigma=fltarr(n_freq)+1.
  
  vs_mean = mean(vis_sigma)
  ;vis_sigma[*] = vs_mean
  
  ;; sqrt(2) is because vis_sigma is for the real or imaginary part separately
  t_sys_meas = (eff_area * sqrt(df * tau) * vis_sigma * sqrt(2)) / ((2. * (1.38065e-23) * 1e26))  ;; in K
  
  n_vis = reform(file_struct.n_vis)
  n_vis_freq = reform(file_struct.n_vis_freq)
  if n_elements(freq_ch_range) ne 0 then begin
    n_vis = total(n_vis_freq[*, min(freq_ch_range):max(freq_ch_range)], 2)
    n_vis_freq = n_vis_freq[*, min(freq_ch_range):max(freq_ch_range)]
  endif
  
  if healpix or not uvf_input then begin
  
    ps_uvf, file_struct, refresh_data = dft_refresh_data, refresh_weight = dft_refresh_weight, $
      refresh_beam = refresh_beam, dft_fchunk = dft_fchunk, freq_ch_range = freq_ch_range, freq_flags = freq_flags, $
      image_window_name = image_window_name, image_window_frac_size = image_window_frac_size, $
      delta_uv_lambda = delta_uv_lambda, max_uv_lambda = max_uv_lambda, no_dft_progress = no_dft_progress, $
      n_vis_freq = n_vis_freq, n_freq = n_freq, input_uvf_varname = input_uvf_varname, $
      input_uvf_files = input_uvf_files, this_run_git_hash = this_run_git_hash
      
    if healpix then begin
      ;; Angular resolution is given in Healpix paper in units of arcminutes, need to convert to radians
      ang_resolution = sqrt(3./!pi) * 3600./file_struct.nside * (1./60.) * (!pi/180.)
      pix_area_rad = ang_resolution^2. ;; by definition of ang. resolution in Healpix paper
    endif else pix_area_rad = (abs(file_struct.degpix) * !pi / 180d)^2.
    
  endif else pix_area_rad = (!dtor*file_struct.degpix)^2.
  
  ;; get weights, data and variance cubes
  if healpix or not keyword_set(uvf_input) then begin
  
    weights_cube1 = prep_uvf_cube(file_struct.uvf_weight_savefile[0], 'weights_cube', z_mpc_mean, healpix, uvf_input, $
      git_hash, kx_mpc, ky_mpc, pix_area_rad = pix_area_rad)
    uvf_wt_git_hashes = git_hash
    data_cube1 = prep_uvf_cube(file_struct.uvf_savefile[0], 'data_cube', z_mpc_mean, healpix, uvf_input, $
      git_hash, kx_mpc, ky_mpc, pix_area_rad = pix_area_rad)
    uvf_git_hashes = git_hash
    if nfiles eq 2 then begin
      weights_cube2 = prep_uvf_cube(file_struct.uvf_weight_savefile[1], 'weights_cube', z_mpc_mean, healpix, uvf_input, $
        git_hash, kx_mpc, ky_mpc, pix_area_rad = pix_area_rad)
      uvf_wt_git_hashes = [uvf_wt_git_hashes, git_hash]
      data_cube2 = prep_uvf_cube(file_struct.uvf_savefile[1], 'data_cube', z_mpc_mean, healpix, uvf_input, $
        git_hash, kx_mpc, ky_mpc, pix_area_rad = pix_area_rad)
      uvf_git_hashes = [uvf_git_hashes, git_hash]
    endif
    
    if not no_var then begin
      variance_cube1 = prep_uvf_cube(file_struct.uvf_weight_savefile[0], 'variance_cube', z_mpc_mean, healpix, uvf_input, $
        git_hash, kx_mpc, ky_mpc, pix_area_rad = pix_area_rad)
      if nfiles eq 2 then variance_cube2 = prep_uvf_cube(file_struct.uvf_weight_savefile[1], 'variance_cube', z_mpc_mean, healpix, uvf_input, $
        git_hash, kx_mpc, ky_mpc, pix_area_rad = pix_area_rad)
    endif
    
  endif else begin
    ;; uvf_input
  
    weights_cube1 = prep_uvf_cube(file_struct.weightfile[0], file_struct.weightvar, z_mpc_mean, healpix, uvf_input, $
      git_hash, kx_mpc, ky_mpc, pol_index=file_struct.pol_index, uv_avg = uv_avg, uv_img_clip = uv_img_clip)
    uvf_wt_git_hashes = git_hash
    if nfiles eq 2 then begin
      weights_cube2 = prep_uvf_cube(file_struct.weightfile[1], file_struct.weightvar, z_mpc_mean, healpix, uvf_input, $
        git_hash, kx_mpc, ky_mpc, pol_index=file_struct.pol_index, uv_avg = uv_avg, uv_img_clip = uv_img_clip)
      uvf_wt_git_hashes = [uvf_wt_git_hashes, git_hash]
    endif
    
    if datavar eq '' then begin
      ;; working with a 'derived' cube (ie residual cube) that is constructed from other cubes
    
      dirty_cube1 = prep_uvf_cube(input_uvf_files[0,0], input_uvf_varname[0,0], z_mpc_mean, healpix, uvf_input, $
        git_hash, kx_mpc, ky_mpc, pol_index=file_struct.pol_index, uv_avg = uv_avg, uv_img_clip = uv_img_clip)
      uvf_git_hashes = git_hash
      model_cube1 = prep_uvf_cube(input_uvf_files[0,1], input_uvf_varname[0,1], z_mpc_mean, healpix, uvf_input, $
        git_hash, kx_mpc, ky_mpc, pol_index=file_struct.pol_index, uv_avg = uv_avg, uv_img_clip = uv_img_clip)
      if nfiles eq 2 then begin
        dirty_cube2 = prep_uvf_cube(input_uvf_files[1,0], input_uvf_varname[1,0], z_mpc_mean, healpix, uvf_input, $
          git_hash, kx_mpc, ky_mpc, pol_index=file_struct.pol_index, uv_avg = uv_avg, uv_img_clip = uv_img_clip)
        uvf_git_hashes = [uvf_git_hashes, git_hash]
        model_cube2 = prep_uvf_cube(input_uvf_files[1,1], input_uvf_varname[1,1], z_mpc_mean, healpix, uvf_input, $
          git_hash, kx_mpc, ky_mpc, pol_index=file_struct.pol_index, uv_avg = uv_avg, uv_img_clip = uv_img_clip)
      endif
      
      data_cube1 = temporary(dirty_cube1) - temporary(model_cube1)
      if nfiles eq 2 then data_cube2 = temporary(dirty_cube2) - temporary(model_cube2)
      
      if max(abs(data_cube1)) eq 0 then message, 'data cube is entirely zero.'
      if nfiles eq 2 then if max(abs(data_cube2)) eq 0 then message, 'data cube is entirely zero.'
      
    endif else begin
    
      data_cube1 = prep_uvf_cube(file_struct.datafile[0], file_struct.datavar, z_mpc_mean, healpix, uvf_input, $
        git_hash, kx_mpc, ky_mpc, pol_index=file_struct.pol_index, uv_avg = uv_avg, uv_img_clip = uv_img_clip)
      uvf_git_hashes = git_hash
      if nfiles eq 2 then begin
        data_cube2 = prep_uvf_cube(file_struct.datafile[1], file_struct.datavar, z_mpc_mean, healpix, uvf_input, $
          git_hash, kx_mpc, ky_mpc, pol_index=file_struct.pol_index, uv_avg = uv_avg, uv_img_clip = uv_img_clip)
        uvf_git_hashes = [uvf_git_hashes, git_hash]
      endif
      
    endelse
    
    if not no_var then begin
      variance_cube1 = prep_uvf_cube(file_struct.variancefile[0], file_struct.variancevar, z_mpc_mean, healpix, uvf_input, $
        git_hash, kx_mpc, ky_mpc, pol_index=file_struct.pol_index, uv_avg = uv_avg, uv_img_clip = uv_img_clip)
      if nfiles eq 2 then variance_cube2 = prep_uvf_cube(file_struct.variancefile[1], file_struct.variancevar, z_mpc_mean, healpix, uvf_input, $
        git_hash, kx_mpc, ky_mpc, pol_index=file_struct.pol_index, uv_avg = uv_avg, uv_img_clip = uv_img_clip)
        
      if max(abs(imaginary(variance_cube1))) gt 0 then begin
        print, 'variance_cube1 is not real, using absolute value'
        variance_cube1 = abs(variance_cube1)
      endif else variance_cube1 = real_part(variance_cube1)
      if nfiles eq 2 then if max(abs(imaginary(variance_cube2))) gt 0 then begin
        print, 'variance_cube2 is not real, using absolute value'
        variance_cube2 = abs(variance_cube2)
      endif else variance_cube2 = real_part(variance_cube2)
      
      if max(abs(variance_cube1)) eq 0 then message, 'variance cube is entirely zero.'
      if nfiles eq 2 then if max(abs(variance_cube2)) eq 0 then message, 'variance cube is entirely zero.'
    endif
    
    ;; get beam sorted out
    if tag_exist(file_struct, 'beam_savefile') then begin
      test_beam = file_test(file_struct.beam_savefile) * ( 1- file_test(file_struct.beam_savefile, /zero_length))
      if min(test_beam) eq 0 or keyword_set(refresh_beam) then begin
      
        for i=0, nfiles-1 do begin
          arr = getvar_savefile(file_struct.beamfile[i], file_struct.beamvar)
          void = getvar_savefile(file_struct.beamfile[i], names = beam_varnames)
          wh_obs = where(stregex(strlowcase(beam_varnames), 'obs', /boolean), count_obs)
          if count_obs gt 0 then obs_struct_name = beam_varnames[wh_obs[0]]
          obs_beam = getvar_savefile(file_struct.beamfile[i], obs_struct_name)
          nfvis_beam = obs_beam.nf_vis
          undefine_fhd, obs_beam
          
          if max(arr) le 1.1 then begin
            ;; beam is peak normalized to 1
            temp = arr * rebin(reform(nfvis_beam, 1, 1, n_elements(nfvis_beam)), n_kx, n_ky, n_elements(nfvis_beam), /sample)
          endif else if max(arr) le file_struct.n_obs[i] then begin
            ;; beam is peak normalized to 1 for each obs, then summed over obs so peak is ~ n_obs
            temp = (arr/file_struct.n_obs[i]) * rebin(reform(nfvis_beam, 1, 1, n_freq), n_kx, n_ky, n_freq, /sample)
          endif else begin
            ;; beam is peak normalized to 1 then multiplied by n_vis_freq for each obs & summed
            temp = arr
          endelse
          
          avg_beam = total(temp, 3) / total(nfvis_beam)
          
          git, repo_path = ps_repository_dir(), result=beam_git_hash
          
          save, file=file_struct.beam_savefile[i], avg_beam, beam_git_hash
        endfor
        
      endif
      
    endif
    
    if keyword_set(uv_avg) then begin
      ;; averging reduces the value of total(weights) ~ n_vis needed for the window int calculation
      n_vis = n_vis/(uv_avg)^2.
    endif
    
  endelse
  
  if max(abs(weights_cube1)) eq 0 then message, 'weights cube is entirely zero.'
  if nfiles eq 2 then if max(abs(weights_cube2)) eq 0 then message, 'weights cube is entirely zero.'
  
  if not no_var then begin
  
    ;; calculate integral of window function
    ;; already cut out negative ky, so multiply by 2
    if nfiles eq 2 then window_int = 2*[total(variance_cube1)/n_vis[0], $
      total(variance_cube2)/n_vis[1]] $
    else window_int = 2*total(variance_cube1)/n_vis[0]
    
    if tag_exist(file_struct, 'beam_savefile') then begin
      beam1 = getvar_savefile(file_struct.beam_savefile[0], 'avg_beam')
      if nfiles eq 2 then beam2 = getvar_savefile(file_struct.beam_savefile[1], 'avg_beam')
      
      void = getvar_savefile(file_struct.beam_savefile[0], names = uvf_varnames)
      wh_hash = where(uvf_varnames eq 'beam_git_hash', count_hash)
      if count_hash gt 0 then begin
        beam_git_hashes = getvar_savefile(file_struct.beam_savefile[0], 'beam_git_hash')
        if nfiles eq 2 then beam_git_hashes = [beam_git_hashes, getvar_savefile(file_struct.beam_savefile[1], 'beam_git_hash')]
      endif else beam_git_hashes = strarr(nfiles)
      
      pix_area_mpc = pix_area_rad * z_mpc_mean^2.
      
      if nfiles eq 2 then window_int_beam = [total(beam1), total(beam2)]*pix_area_mpc*(z_mpc_delta * n_freq) $
      else window_int_beam = total(beam1)*pix_area_mpc*(z_mpc_delta * n_freq)
      
      volume_factor = total(beam1*0+1.)*pix_area_mpc*(z_mpc_delta * n_freq)
      if nfiles eq 2 then volume_factor = fltarr(2) + volume_factor
      
      bandwidth_factor = z_mpc_delta * n_freq
      if nfiles eq 2 then bandwidth_factor = fltarr(2) + bandwidth_factor
      
    endif else beam_git_hashes = ''
    
    if tag_exist(file_struct, 'beam_int') then begin
      if nfiles eq 2 then ave_beam_int = total(file_struct.beam_int * file_struct.n_vis_freq, 2) / total(file_struct.n_vis_freq, 2) $
      else ave_beam_int = total(file_struct.beam_int * file_struct.n_vis_freq) / total(file_struct.n_vis_freq)
      
      ;; fix known units bug in some early runs
      if max(ave_beam_int) lt 0.01 then ave_beam_int = ave_beam_int / (file_struct.kpix)^4. $
      else if max(ave_beam_int) lt 0.03 then ave_beam_int = ave_beam_int / (file_struct.kpix)^2.
      
      ;; convert rad -> Mpc^2, multiply by depth in Mpc
      window_int_beam_obs = ave_beam_int * z_mpc_mean^2. * (z_mpc_delta * n_freq)
    endif
    
  endif
    
  n_kx = n_elements(kx_mpc)
  n_ky = n_elements(ky_mpc)
  kx_mpc_delta = kx_mpc[1] - kx_mpc[0]
  ky_mpc_delta = ky_mpc[1] - ky_mpc[0]
  
  if keyword_set(sim) and keyword_set(fix_sim_input) then begin
    ;; fix for some sims that used saved uvf inputs without correcting for factor of 2 loss
    data_cube1 = data_cube1 * 2.
    if nfiles eq 2 then data_cube2 = data_cube2 * 2.
  endif
  
  ;; save some slices of the raw data cube (before dividing by weights) & weights
  for i=0, nfiles-1 do begin
    if i eq 0 then data_cube = data_cube1 else data_cube = data_cube2
    if i eq 0 then weights_cube = weights_cube1 else weights_cube = weights_cube2
    
    uf_tot = total(total(abs(weights_cube),3),1)
    wh_uf_n0 = where(uf_tot gt 0, count_uf_n0)
    if count_uf_n0 eq 0 then message, 'uvf weights appear to be entirely zero'
    min_dist_uf_n0 = min(wh_uf_n0, min_loc)
    uf_slice_ind = wh_uf_n0[min_loc]
    
    uf_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 1, $
      slice_inds = uf_slice_ind, slice_savefile = file_struct.uf_raw_savefile[i])
      
    uf_weight_slice = uvf_slice(weights_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 1, $
      slice_inds = uf_slice_ind, slice_savefile = file_struct.uf_weight_savefile[i])
    undefine, uf_slice, uf_weight_slice
    
    vf_tot = total(total(abs(weights_cube),3),2)
    wh_vf_n0 = where(vf_tot gt 0, count_vf_n0)
    if count_vf_n0 eq 0 then message, 'uvf weights appear to be entirely zero'
    min_dist_vf_n0 = min(abs(n_kx/2-wh_vf_n0), min_loc)
    vf_slice_ind = wh_vf_n0[min_loc]
    
    vf_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 0, $
      slice_inds = vf_slice_ind, slice_savefile = file_struct.vf_raw_savefile[i])
      
    vf_weight_slice = uvf_slice(weights_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 0, $
      slice_inds = vf_slice_ind, slice_savefile = file_struct.vf_weight_savefile[i])
      
    if max(abs(vf_slice)) eq 0 then message, 'vf data slice is entirely zero'
    undefine, vf_slice, vf_weight_slice
    
    uv_tot = total(total(abs(weights_cube),2),1)
    wh_uv_n0 = where(uv_tot gt 0, count_uv_n0)
    if count_uv_n0 eq 0 then message, 'uvf weights appear to be entirely zero'
    min_dist_uv_n0 = min(wh_uv_n0, min_loc)
    uv_slice_ind = wh_uv_n0[min_loc]
    
    uv_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 2, $
      slice_inds = uv_slice_ind, slice_savefile = file_struct.uv_raw_savefile[i])
      
    uv_weight_slice = uvf_slice(weights_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 2, $
      slice_inds = uv_slice_ind, slice_savefile = file_struct.uv_weight_savefile[i])
      
    if max(abs(uv_slice)) eq 0 then message, 'uv data slice is entirely zero'
    undefine, uv_slice, uv_weight_slice
    
    undefine, data_cube, weights_cube
  endfor
  
  ave_weights = total(total(abs(weights_cube1),2),1)/(n_kx*n_ky)
  if nfiles eq 2 then ave_weights = transpose([[ave_weights], [total(total(abs(weights_cube2),2),1)/(n_kx*n_ky)]])
  
  ;; make sigma2 cubes
  if no_var then begin
    sigma2_cube1 = 1./abs(weights_cube1)
    if nfiles eq 2 then sigma2_cube2 = 1./abs(weights_cube2)
  endif else begin
    sigma2_cube1 = temporary(variance_cube1) / (abs(weights_cube1)^2.)
    if nfiles eq 2 then sigma2_cube2 = temporary(variance_cube2) / (abs(weights_cube2)^2.)
  endelse
  
  wh_wt1_0 = where(abs(weights_cube1) eq 0, count_wt1_0)
  if count_wt1_0 ne 0 then sigma2_cube1[wh_wt1_0] = 0
  if nfiles eq 2 then begin
    wh_wt2_0 = where(abs(weights_cube2) eq 0, count_wt2_0)
    if count_wt2_0 ne 0 then sigma2_cube2[wh_wt2_0] = 0
  endif
  
  if nfiles eq 2 then begin
    if min(sigma2_cube1) lt 0 or min(sigma2_cube2) lt 0 then message, 'sigma2 should be positive definite.'
    if total(abs(sigma2_cube1)) le 0 or total(abs(sigma2_cube2)) le 0 then message, 'one or both sigma2 cubes is all zero'
  endif else begin
    if min(sigma2_cube1) lt 0 then message, 'sigma2 should be positive definite.'
    if total(abs(sigma2_cube1)) le 0 then message, 'sigma2 cube is all zero'
  endelse
  
  
  wt_meas_ave = total(abs(weights_cube1), 3)/n_freq
  wt_meas_min = min(abs(weights_cube1), dimension=3)
  
  ;; divide data by weights
  data_cube1 = data_cube1 / weights_cube1
  if count_wt1_0 ne 0 then data_cube1[wh_wt1_0] = 0
  undefine, weights_cube1, wh_wt1_0, count_wt1_0
  
  if nfiles eq 2 then begin
    data_cube2 = data_cube2 / weights_cube2
    if count_wt2_0 ne 0 then data_cube2[wh_wt2_0] = 0
    undefine, weights_cube2, wh_wt2_0, count_wt2_0
  endif
  
  ;; save some slices of the data cube
  for i=0, nfiles-1 do begin
    if i eq 0 then data_cube = data_cube1 else data_cube = data_cube2
    uf_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 1, $
      slice_inds = uf_slice_ind, slice_savefile = file_struct.uf_savefile[i])
    undefine, uf_slice
    
    vf_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 0, $
      slice_inds = vf_slice_ind, slice_savefile = file_struct.vf_savefile[i])
    undefine, vf_slice
    
    uv_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 2, $
      slice_inds = uv_slice_ind, slice_savefile = file_struct.uv_savefile[i])
    undefine, uv_slice
    
    undefine, data_cube
  endfor
  
  if healpix or not keyword_set(uvf_input) then begin
    ;; fix units on window funtion integral -- now they should be Mpc^3
    ;; checked vs Adam's analytic calculation and it matches to within a factor of 4
    ;; We are using total(weight) = Nvis for Ian's uv plane and multiplying by a factor to convert between Ian's and my uv plane
    ;;   the factor is ((2*!pi)^2*(delta_uv)^2) /  (delta_kperp)^2 * Dm^2 which comes in squared in the denominator.
    ;;   see eq 21e from Adam's memo. Also note that delta D is delta z * n_freq
    ;;   note that we can convert both the weights and variance to uvf from kperp,rz and all jacobians will cancel
    window_int_k = window_int * (z_mpc_delta * n_freq) * (kx_mpc_delta * ky_mpc_delta)*z_mpc_mean^4./((2.*!pi)^2.*file_struct.kpix^4.)
    print, 'window integral from variances: ' + number_formatter(window_int_k[0], format='(e10.4)')
    if tag_exist(file_struct, 'beam_savefile') then print, 'window integral from beam cube: ' + number_formatter(window_int_beam[0], format='(e10.4)')
    if tag_exist(file_struct, 'beam_int') then print, 'window integral from obs.beam_integral: ' + number_formatter(window_int_beam_obs[0], format='(e10.4)')
    
    if tag_exist(file_struct, 'beam_int') then begin
      window_int = window_int_beam_obs
    endif
    
    if (n_elements(window_int) eq 0 or min(window_int) eq 0) then begin
      if keyword_set(allow_beam_approx) then begin
        print, 'WARNING: beam integral in obs structure is zero, using a less good approximation'
      endif else begin
        message, 'Beam integral in obs structure is zero. To use a less good approximation instead, set keyword allow_beam_approx=1'
      endelse
      
      if tag_exist(file_struct, 'beam_savefile') then begin
        window_int = window_int_beam
        if min(window_int) eq 0 then print, 'WARNING: beam cube is zero, using a less good approximation'
      endif
    endif
    if (n_elements(window_int) eq 0 or min(window_int) eq 0) then window_int = window_int_k
  ;if keyword_set(sim) then window_int = 2.39e9 + fltarr(nfiles)
  endif else begin
    window_int_k = window_int * (z_mpc_delta * n_freq) * (2.*!pi)^2. / (kx_mpc_delta * ky_mpc_delta)
    print, 'var_cube multiplier: ', (z_mpc_delta * n_freq) * (2.*!pi)^2. / (kx_mpc_delta * ky_mpc_delta)
    print, 'window integral from variances: ' + number_formatter(window_int_k[0], format='(e10.4)')
    if tag_exist(file_struct, 'beam_savefile') then print, 'window integral from beam cube: ' + number_formatter(window_int_beam[0], format='(e10.4)')
    if tag_exist(file_struct, 'beam_int') then print, 'window integral from obs.beam_integral: ' + number_formatter(window_int_beam_obs[0], format='(e10.4)')
    
    if tag_exist(file_struct, 'beam_int') then begin
      window_int = window_int_beam_obs
      if min(window_int) eq 0 then print, 'WARNING: beam integral in obs structure is zero, using a less good approximation'
    endif
    if (n_elements(window_int) eq 0 or min(window_int) eq 0) and tag_exist(file_struct, 'beam_savefile') then begin
      window_int = window_int_beam
      if min(window_int) eq 0 then print, 'WARNING: beam cube is zero, using a less good approximation'
    endif
    if (n_elements(window_int) eq 0 or min(window_int) eq 0) then window_int = window_int_k
    
  endelse
  
  
  ;; get sigma^2 into Jy^2
  sigma2_cube1 = sigma2_cube1 * rebin(reform(vis_sigma, 1, 1, n_freq), n_kx, n_ky, n_freq)^2.
  if nfiles eq 2 then sigma2_cube2 = sigma2_cube2 * rebin(reform(vis_sigma, 1, 1, n_freq), n_kx, n_ky, n_freq)^2.
  
  ;; get data & sigma into mK Mpc^2 and multiply by 2 (4 for variances) to get to estimate of Stokes I rather than instrumental pol
  for i=0, n_freq-1 do begin
    data_cube1[*,*,i] = data_cube1[*,*,i]*conv_factor[i]*2.
    sigma2_cube1[*,*,i] = sigma2_cube1[*,*,i]*(conv_factor[i])^2.*4.
    if nfiles eq 2 then begin
      data_cube2[*,*,i] = data_cube2[*,*,i]*conv_factor[i]*2.
      sigma2_cube2[*,*,i] = sigma2_cube2[*,*,i]*(conv_factor[i])^2.*4.
    endif
  endfor
  
  
  ;; divide data by sqrt(window_int) and sigma2 by window_int
  data_cube1 = data_cube1 / sqrt(window_int[0])
  sigma2_cube1 = sigma2_cube1 / window_int[0]
  if nfiles eq 2 then begin
    data_cube2 = data_cube2 / sqrt(window_int[1])
    sigma2_cube2 = sigma2_cube2  / window_int[1]
  endif
  
  ;; make simulated noise cubes
  sim_noise1 = randomn(seed, n_kx, n_ky, n_kz) * sqrt(sigma2_cube1) + complex(0,1) * randomn(seed, n_kx, n_ky, n_kz) * sqrt(sigma2_cube1)
  if nfiles eq 2 then sim_noise2 = randomn(seed, n_kx, n_ky, n_kz) * sqrt(sigma2_cube2) + complex(0,1) * randomn(seed, n_kx, n_ky, n_kz) * sqrt(sigma2_cube2)
  
  if nfiles eq 2 then begin
    ;; Now construct added & subtracted cubes (weighted by inverse variance) & new variances
    sum_weights1 = 1./sigma2_cube1
    wh_sig1_0 = where(sigma2_cube1 eq 0, count_sig1_0, complement = wh_sig1_n0)
    if count_sig1_0 ne 0 then sum_weights1[wh_sig1_0] = 0
    undefine, sigma2_cube1, wh_sig1_0, count_sig1_0
    
    sum_weights2 = 1./sigma2_cube2
    wh_sig2_0 = where(sigma2_cube2 eq 0, count_sig2_0, complement = wh_sig2_n0)
    if count_sig2_0 ne 0 then sum_weights2[wh_sig2_0] = 0
    undefine, sigma2_cube2, wh_sig2_0, count_sig2_0
    
    wt_ave_power_freq = fltarr(2, n_freq)
    wt_ave_power_freq[0,*] = total(total(sum_weights1 * abs(data_cube1)^2., 2), 1)/total(total(sum_weights1, 2), 1) * (z_mpc_delta * n_freq)^2.
    wt_ave_power_freq[1,*] = total(total(sum_weights2 * abs(data_cube2)^2., 2), 1)/total(total(sum_weights2, 2), 1) * (z_mpc_delta * n_freq)^2.
    ave_power_freq = fltarr(2, n_freq)
    for i=0, n_freq-1 do ave_power_freq[*, i] = [mean(abs((data_cube1[*,*,i])[where(sum_weights1[*,*,i] ne 0),*])^2.), $
      mean(abs((data_cube2[*,*,i])[where(sum_weights2[*,*,i] ne 0),*])^2.)] * (z_mpc_delta * n_freq)^2.
      
    wt_ave_power_uvf = [total(sum_weights1 * abs(data_cube1)^2.)/total(sum_weights1), $
      total(sum_weights2 * abs(data_cube2)^2.)/total(sum_weights2)] * (z_mpc_delta * n_freq)^2.
    ave_power_uvf = fltarr(2)
    ave_power_uvf[0] = mean(abs(data_cube1[wh_sig1_n0])^2.) * (z_mpc_delta * n_freq)^2.
    ave_power_uvf[1] = mean(abs(data_cube2[wh_sig2_n0])^2.) * (z_mpc_delta * n_freq)^2.
    
    undefine, wh_sig1_n0, wh_sig2_n0
    
    sum_weights_net = sum_weights1 + sum_weights2
    wh_wt0 = where(sum_weights_net eq 0, count_wt0)
    
    data_sum = (sum_weights1 * data_cube1 + sum_weights2 * data_cube2)/sum_weights_net
    data_diff = (sum_weights1 * data_cube1 - sum_weights2 * data_cube2)/sum_weights_net
    sim_noise_sum = (sum_weights1 * sim_noise1 + sum_weights2 * sim_noise2)/sum_weights_net
    sim_noise_diff = (sum_weights1 * sim_noise1 - sum_weights2 * sim_noise2)/sum_weights_net
    undefine, data_cube1, data_cube2, sim_noise1, sim_noise2
    
    if count_wt0 ne 0 then begin
      data_sum[wh_wt0] = 0
      data_diff[wh_wt0] = 0
      sim_noise_sum[wh_wt0] = 0
      sim_noise_diff[wh_wt0] = 0
    endif
    undefine, sum_weights1, sum_weights2
    
    sum_sigma2 = 1./temporary(sum_weights_net)
    if count_wt0 ne 0 then sum_sigma2[wh_wt0] = 0
    
  ;sim_noise = randomn(seed, n_kx, n_ky, n_kz) * sqrt(sum_sigma2) + complex(0,1) * randomn(seed, n_kx, n_ky, n_kz) * sqrt(sum_sigma2)
    
  endif else begin
    sum_weights1 = 1./sigma2_cube1
    wh_sig1_0 = where(sigma2_cube1 eq 0, count_sig1_0, complement = wh_sig1_n0)
    if count_sig1_0 ne 0 then sum_weights1[wh_sig1_0] = 0
    
    wt_ave_power_freq = total(total(sum_weights1 * abs(data_cube1)^2., 2), 1)/total(total(sum_weights1, 2), 1)
    ave_power_freq = fltarr(n_freq)
    for i=0, n_freq-1 do ave_power_freq[i] = mean(abs((data_cube1[*,*,i])[where(sum_weights1[*,*,i] ne 0),*])^2.)
    
    wt_ave_power_uvf = total(sum_weights1 * abs(data_cube1)^2.)/total(sum_weights1)
    ave_power_uvf = mean(abs(data_cube1[wh_sig1_n0])^2.)
    undefine, sum_weights1, wh_sig1_n0
    
    data_sum = temporary(data_cube1)
    sum_sigma2 = temporary(sigma2_cube1)
    sim_noise_sum = temporary(sim_noise1)
  endelse
  
  mask = intarr(n_kx, n_ky, n_kz) + 1
  wh_sig0 = where(sum_sigma2 eq 0, count_sig0)
  if count_sig0 gt 0 then mask[wh_sig0] = 0
  ;; n_pix_contrib = total(total(mask, 2), 1)
  n_freq_contrib = total(mask, 3)
  wh_nofreq = where(n_freq_contrib eq 0, count_nofreq)
  undefine, mask
  
  ;; save some slices of the sum & diff cubes
  uf_slice = uvf_slice(data_sum, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 1, $
    slice_inds = uf_slice_ind, slice_savefile = file_struct.uf_sum_savefile)
  if nfiles eq 2 then $
    uf_slice = uvf_slice(data_diff, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 1, $
    slice_inds = uf_slice_ind, slice_savefile = file_struct.uf_diff_savefile)
  undefine, uf_slice
  
  vf_slice = uvf_slice(data_sum, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 0, $
    slice_inds = vf_slice_ind, slice_savefile = file_struct.vf_sum_savefile)
  if nfiles eq 2 then $
    vf_slice2 = uvf_slice(data_diff, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 0, $
    slice_inds = vf_slice_ind, slice_savefile = file_struct.vf_diff_savefile) $
  else vf_slice2 = vf_slice
  undefine, vf_slice
  
  uv_slice = uvf_slice(data_sum, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 2, $
    slice_inds = uv_slice_ind, slice_savefile = file_struct.uv_sum_savefile)
  if nfiles eq 2 then $
    uv_slice2 = uvf_slice(data_diff, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 2, $
    slice_inds = uv_slice_ind, slice_savefile = file_struct.uv_diff_savefile) $
  else uv_slice2 = uv_slice
  undefine, uv_slice
  
  
  if keyword_set(ave_removal) then begin
  
    data_sum_mean = mean(data_sum, dimension=3)
    data_sum = data_sum - (rebin(real_part(data_sum_mean), n_kx, n_ky, n_freq, /sample) + $
      complex(0,1)*rebin(imaginary(data_sum_mean), n_kx, n_ky, n_freq, /sample))
      
    sim_noise_sum_mean = mean(sim_noise_sum, dimension=3)
    sim_noise_sum = sim_noise_sum - (rebin(real_part(sim_noise_sum_mean), n_kx, n_ky, n_freq, /sample) + $
      complex(0,1)*rebin(imaginary(sim_noise_sum_mean), n_kx, n_ky, n_freq, /sample))
    if nfiles eq 2 then begin
      data_diff_mean = mean(data_diff, dimension=3)
      data_diff = data_diff - (rebin(real_part(data_diff_mean), n_kx, n_ky, n_freq, /sample) + $
        complex(0,1)*rebin(imaginary(data_diff_mean), n_kx, n_ky, n_freq, /sample))
        
      sim_noise_diff_mean = mean(sim_noise_diff, dimension=3)
      sim_noise_diff = sim_noise_diff - (rebin(real_part(sim_noise_diff_mean), n_kx, n_ky, n_freq, /sample) + $
        complex(0,1)*rebin(imaginary(sim_noise_diff_mean), n_kx, n_ky, n_freq, /sample))
    endif
  endif
  
  ;; apply spectral windowing function if desired
  if n_elements(spec_window_type) ne 0 then begin
    window = spectral_window(n_freq, type = spec_window_type, /periodic)
    
    norm_factor = sqrt(n_freq/total(window^2.))
    
    window = window * norm_factor
    
    window_expand = rebin(reform(window, 1, 1, n_freq), n_kx, n_ky, n_freq, /sample)
    
    data_sum = data_sum * window_expand
    sim_noise_sum = sim_noise_sum * window_expand
    if nfiles eq 2 then begin
      data_diff = data_diff * window_expand
      sim_noise_diff = sim_noise_diff * window_expand
    endif
    
    sum_sigma2 = sum_sigma2 * temporary(window_expand^2.)
  endif
  
  if keyword_set(inverse_covar_weight) then begin
  
    max_eta_val = max(sum_sigma2)
    eta_var = shift([max_eta_val, fltarr(n_freq-1)], n_freq/2)
    covar_eta_fg = diag_matrix(eta_var)
    
    identity = diag_matrix([fltarr(n_freq)+1d])
    ft_matrix = fft(identity, dimension=1)*n_freq*z_mpc_delta
    inv_ft_matrix = fft(identity, dimension=1, /inverse)*kz_mpc_delta
    undefine, identity
    
    covar_z_fg = matrix_multiply(inv_ft_matrix, matrix_multiply(shift(covar_eta_fg, n_freq/2, n_freq/2), conj(inv_ft_matrix), /btranspose))
    if max(abs(imaginary(covar_z_fg))) eq 0 then covar_z_fg = real_part(covar_z_fg)
    
    wt_data_sum = fltarr(n_kx, n_ky, n_freq)
    if nfiles eq 2 then begin
      wt_data_diff = fltarr(n_kx, n_ky, n_freq)
    endif
    wt_sum_sigma2 = fltarr(n_kx, n_ky, n_freq)
    wt_power_norm = fltarr(n_kx, n_ky, n_freq, n_freq)
    wt_sum_sigma2_norm = fltarr(n_kx, n_ky, n_freq, n_freq)
    for i=0, n_kx-1 do begin
      for j=0, n_ky-1 do begin
        var_z_inst = reform(sum_sigma2[i,j,*])
        if max(abs(var_z_inst)) eq 0 then continue
        
        ;; need to convert zeros which represet lack of measurements to infinities (or large numbers)
        wh_var0 = where(var_z_inst eq 0, count_var0)
        if count_var0 gt 0 then var_z_inst[wh_var0] = 1e6
        
        covar_z = diag_matrix(var_z_inst) + covar_z_fg
        inv_covar_z = la_invert(covar_z)
        inv_covar_kz = shift(matrix_multiply(ft_matrix, matrix_multiply(inv_covar_z, conj(ft_matrix), /btranspose)), n_freq/2, n_freq/2)
        fisher_kz = abs(inv_covar_kz)^2.
        
        inv_var = 1/var_z_inst
        wh_sigma0 = where(var_z_inst eq 0, count_sigma0)
        if count_sigma0 gt 0 then inv_var[wh_sigma0] = 0
        ;; norm2 = total(inv_var)^2./(n_freq*total(inv_var^2.))
        
        wt_data_sum[i,j,*] = matrix_multiply(inv_covar_z, reform(data_sum[i,j,*])) * z_mpc_delta
        
        if nfiles eq 2 then begin
          wt_data_diff[i,j,*] = matrix_multiply(inv_covar_z, reform(data_diff[i,j,*])) * z_mpc_delta
        endif
        wt_sum_sigma2[i,j,*] = diag_matrix(shift(matrix_multiply(inv_covar_z, matrix_multiply(diag_matrix(var_z_inst), conj(inv_covar_z), /btranspose)), n_freq/2, n_freq/2))
        
        m_matrix = la_invert(matrix_sqrt(fisher_kz))
        wt_power_norm[i,j,*,*] = m_matrix
        wt_sum_sigma2_norm[i,j,*,*] = matrix_multiply(m_matrix, matrix_multiply(fisher_kz, m_matrix, /btranspose))
        
      endfor
    endfor
    undefine, covar_z_fg, ft_matrix, inv_ft_matrix, covar_z, inv_covar_z, inv_covar_kz
    
    data_sum = temporary(wt_data_sum)
    if nfiles eq 2 then begin
      data_diff = temporary(wt_data_diff)
    endif
    sum_sigma2 = temporary(wt_sum_sigma2)
    
  endif
  
  
  ;; drop pixels with less than 1/3 of the frequencies (set weights to 0)
  wh_fewfreq = where(n_freq_contrib lt ceil(n_freq/3d), count_fewfreq)
  if count_fewfreq gt 0 then begin
    mask_fewfreq = n_freq_contrib * 0 + 1
    mask_fewfreq[wh_fewfreq] = 0
    mask_fewfreq = rebin(temporary(mask_fewfreq), n_kx, n_ky, n_kz)
    
    data_sum = temporary(data_sum) * mask_fewfreq
    sim_noise_sum = temporary(sim_noise_sum) * mask_fewfreq
    
    if nfiles gt 1 then begin
      data_diff = temporary(data_diff) * mask_fewfreq
      sim_noise_diff = temporary(sim_noise_diff) * mask_fewfreq
    endif
    
    sum_sigma2 = temporary(sum_sigma2) * mask_fewfreq
  endif
  
  ;; now take frequency FT
  if even_freq then begin
    ;; evenly spaced, just use fft
    ;; old ft convention
    ; data_sum_ft = fft(data_sum, dimension=3) * n_freq * z_mpc_delta / (2.*!pi)
    data_sum_ft = fft(data_sum, dimension=3) * n_freq * z_mpc_delta
    sim_noise_sum_ft = fft(sim_noise_sum, dimension=3) * n_freq * z_mpc_delta
    
    ;; put k0 in middle of cube
    data_sum_ft = shift(data_sum_ft, [0,0,n_kz/2])
    sim_noise_sum_ft = shift(sim_noise_sum_ft, [0,0,n_kz/2])
    
    undefine, data_sum, sim_noise_sum
    if nfiles eq 2 then begin
      data_diff_ft = fft(data_diff, dimension=3) * n_freq * z_mpc_delta
      ;; put k0 in middle of cube
      data_diff_ft = shift(data_diff_ft, [0,0,n_kz/2])
      undefine, data_diff
      
      sim_noise_diff_ft = fft(sim_noise_diff, dimension=3) * n_freq * z_mpc_delta
      ;; put k0 in middle of cube
      sim_noise_diff_ft = shift(sim_noise_diff_ft, [0,0,n_kz/2])
      undefine, sim_noise_diff
    endif
  endif else begin
    ;; Not evenly spaced. Do a dft
    z_exp =  exp(-1.*complex(0,1)*matrix_multiply(comov_dist_los, kz_mpc_orig, /btranspose))
    
    
    data_sum_ft = z_mpc_delta * $
      reform(matrix_multiply(reform(temporary(data_sum), n_kx*n_ky, n_freq), z_exp), n_kx, n_ky, n_kz)
    sim_noise_sum_ft = z_mpc_delta * $
      reform(matrix_multiply(reform(temporary(sim_noise_sum), n_kx*n_ky, n_freq), z_exp), n_kx, n_ky, n_kz)
    if nfiles eq 2 then begin
      data_diff_ft = z_mpc_delta * $
        reform(matrix_multiply(reform(temporary(data_diff), n_kx*n_ky, n_freq), z_exp), n_kx, n_ky, n_kz)
      sim_noise_diff_ft  = z_mpc_delta * $
        reform(matrix_multiply(reform(temporary(sim_noise_diff), n_kx*n_ky, n_freq), z_exp), n_kx, n_ky, n_kz)
    endif
  endelse
  
  n_val = round(kz_mpc_orig / kz_mpc_delta)
  kz_mpc_orig[where(n_val eq 0)] = 0
  
  kz_mpc = kz_mpc_orig[where(n_val ge 0)]
  n_kz = n_elements(kz_mpc)
  
  
  ;; these an and bn calculations don't match the standard
  ;; convention (they are down by a factor of 2) but they make more sense
  ;; and remove factors of 2 we'd otherwise have in the power
  ;; and variance calculations
  ;; note that the 0th mode will have higher noise because there's half as many measurements going into it
  data_sum_cos = complex(fltarr(n_kx, n_ky, n_kz))
  data_sum_sin = complex(fltarr(n_kx, n_ky, n_kz))
  data_sum_cos[*, *, 0] = data_sum_ft[*,*,where(n_val eq 0)]
  data_sum_cos[*, *, 1:n_kz-1] = (data_sum_ft[*,*, where(n_val gt 0)] + data_sum_ft[*,*, reverse(where(n_val lt 0))])/2.
  data_sum_sin[*, *, 1:n_kz-1] = complex(0,1) * (data_sum_ft[*,*, where(n_val gt 0)] - data_sum_ft[*,*, reverse(where(n_val lt 0))])/2.
  undefine, data_sum_ft
  
  sim_noise_sum_cos = complex(fltarr(n_kx, n_ky, n_kz))
  sim_noise_sum_sin = complex(fltarr(n_kx, n_ky, n_kz))
  sim_noise_sum_cos[*, *, 0] = sim_noise_sum_ft[*,*,where(n_val eq 0)]
  sim_noise_sum_cos[*, *, 1:n_kz-1] = (sim_noise_sum_ft[*,*, where(n_val gt 0)] + sim_noise_sum_ft[*,*, reverse(where(n_val lt 0))])/2.
  sim_noise_sum_sin[*, *, 1:n_kz-1] = complex(0,1) * (sim_noise_sum_ft[*,*, where(n_val gt 0)] - sim_noise_sum_ft[*,*, reverse(where(n_val lt 0))])/2.
  undefine, sim_noise_sum_ft
  
  if nfiles gt 1 then begin
    data_diff_cos = complex(fltarr(n_kx, n_ky, n_kz))
    data_diff_sin = complex(fltarr(n_kx, n_ky, n_kz))
    data_diff_cos[*, *, 0] = data_diff_ft[*,*,where(n_val eq 0)]
    data_diff_cos[*, *, 1:n_kz-1] = (data_diff_ft[*,*, where(n_val gt 0)] + data_diff_ft[*,*, reverse(where(n_val lt 0))])/2.
    data_diff_sin[*, *, 1:n_kz-1] = complex(0,1) * (data_diff_ft[*,*, where(n_val gt 0)] - data_diff_ft[*,*, reverse(where(n_val lt 0))])/2.
    undefine, data_diff_ft
    
    sim_noise_diff_cos = complex(fltarr(n_kx, n_ky, n_kz))
    sim_noise_diff_sin = complex(fltarr(n_kx, n_ky, n_kz))
    sim_noise_diff_cos[*, *, 0] = sim_noise_diff_ft[*,*,where(n_val eq 0)]
    sim_noise_diff_cos[*, *, 1:n_kz-1] = (sim_noise_diff_ft[*,*, where(n_val gt 0)] + sim_noise_diff_ft[*,*, reverse(where(n_val lt 0))])/2.
    sim_noise_diff_sin[*, *, 1:n_kz-1] = complex(0,1) * (sim_noise_diff_ft[*,*, where(n_val gt 0)] - sim_noise_diff_ft[*,*, reverse(where(n_val lt 0))])/2.
    undefine, sim_noise_diff_ft
  endif
  
  if keyword_set(ave_removal) then begin
  
    data_sum_cos[*,*,0] = data_sum_cos[*,*,0] + data_sum_mean * n_freq * z_mpc_delta
    sim_noise_sum_cos[*,*,0] = sim_noise_sum_cos[*,*,0] + sim_noise_sum_mean * n_freq * z_mpc_delta
    
    if nfiles eq 2 then begin
      data_diff_cos[*,*,0] = data_diff_cos[*,*,0] + data_diff_mean * n_freq * z_mpc_delta
      sim_noise_diff_cos[*,*,0] = sim_noise_diff_cos[*,*,0] + sim_noise_diff_mean * n_freq * z_mpc_delta
    endif
  endif
  
  if keyword_set(inverse_covar) then begin
  
    wt_power_norm_coscos = complex(fltarr(n_kx, n_ky, n_kz, n_kz))
    wt_power_norm_sinsin = complex(fltarr(n_kx, n_ky, n_kz, n_kz))
    wt_power_norm_cossin = complex(fltarr(n_kx, n_ky, n_kz, n_kz))
    wt_power_norm_sincos = complex(fltarr(n_kx, n_ky, n_kz, n_kz))
    
    wt_power_norm_coscos[*, *, 0, 0] = wt_power_norm[*,*,where(n_val eq 0),where(n_val eq 0)]
    wt_power_norm_coscos[*, *, 1:n_kz-1, 1:n_kz-1] = (wt_power_norm[*,*,where(n_val gt 0),where(n_val gt 0)] $
      + wt_power_norm[*,*,where(n_val gt 0),reverse(where(n_val lt 0))] $
      + wt_power_norm[*,*,reverse(where(n_val lt 0)),where(n_val gt 0)] $
      + wt_power_norm[*,*,reverse(where(n_val lt 0)),reverse(where(n_val lt 0))])/4.
      
    wt_power_norm_sinsin[*, *, 1:n_kz-1, 1:n_kz-1] = (-1*wt_power_norm[*,*,where(n_val gt 0),where(n_val gt 0)] $
      + wt_power_norm[*,*,where(n_val gt 0),reverse(where(n_val lt 0))] $
      + wt_power_norm[*,*,reverse(where(n_val lt 0)),where(n_val gt 0)] $
      - wt_power_norm[*,*,reverse(where(n_val lt 0)),reverse(where(n_val lt 0))])/4.
      
    wt_power_norm_cossin[*, *, 1:n_kz-1, 1:n_kz-1] = complex(0,1) * (wt_power_norm[*,*,where(n_val gt 0),where(n_val gt 0)] $
      - wt_power_norm[*,*,where(n_val gt 0),reverse(where(n_val lt 0))] $
      + wt_power_norm[*,*,reverse(where(n_val lt 0)),where(n_val gt 0)] $
      - wt_power_norm[*,*,reverse(where(n_val lt 0)),reverse(where(n_val lt 0))])/4.
      
    wt_power_norm_sincos[*, *, 1:n_kz-1, 1:n_kz-1] = complex(0,1) * (wt_power_norm[*,*,where(n_val gt 0),where(n_val gt 0)] $
      + wt_power_norm[*,*,where(n_val gt 0),reverse(where(n_val lt 0))] $
      - wt_power_norm[*,*,reverse(where(n_val lt 0)),where(n_val gt 0)] $
      - wt_power_norm[*,*,reverse(where(n_val lt 0)),reverse(where(n_val lt 0))])/4.
      
    undefine, wt_power_norm
    
    
    wt_sum_sigma2_norm_coscos = complex(fltarr(n_kx, n_ky, n_kz, n_kz))
    wt_sum_sigma2_norm_sinsin = complex(fltarr(n_kx, n_ky, n_kz, n_kz))
    wt_sum_sigma2_norm_cossin = complex(fltarr(n_kx, n_ky, n_kz, n_kz))
    wt_sum_sigma2_norm_sincos = complex(fltarr(n_kx, n_ky, n_kz, n_kz))
    
    wt_sum_sigma2_norm_coscos[*, *, 0, 0] = wt_sum_sigma2_norm[*,*,where(n_val eq 0),where(n_val eq 0)]
    wt_sum_sigma2_norm_coscos[*, *, 1:n_kz-1, 1:n_kz-1] = (wt_sum_sigma2_norm[*,*,where(n_val gt 0),where(n_val gt 0)] $
      + wt_sum_sigma2_norm[*,*,where(n_val gt 0),reverse(where(n_val lt 0))] $
      + wt_sum_sigma2_norm[*,*,reverse(where(n_val lt 0)),where(n_val gt 0)] $
      + wt_sum_sigma2_norm[*,*,reverse(where(n_val lt 0)),reverse(where(n_val lt 0))])/4.
      
    wt_sum_sigma2_norm_sinsin[*, *, 1:n_kz-1, 1:n_kz-1] = (-1*wt_sum_sigma2_norm[*,*,where(n_val gt 0),where(n_val gt 0)] $
      + wt_sum_sigma2_norm[*,*,where(n_val gt 0),reverse(where(n_val lt 0))] $
      + wt_sum_sigma2_norm[*,*,reverse(where(n_val lt 0)),where(n_val gt 0)] $
      - wt_sum_sigma2_norm[*,*,reverse(where(n_val lt 0)),reverse(where(n_val lt 0))])/4.
      
    wt_sum_sigma2_norm_cossin[*, *, 1:n_kz-1, 1:n_kz-1] = complex(0,1) * (wt_sum_sigma2_norm[*,*,where(n_val gt 0),where(n_val gt 0)] $
      - wt_sum_sigma2_norm[*,*,where(n_val gt 0),reverse(where(n_val lt 0))] $
      + wt_sum_sigma2_norm[*,*,reverse(where(n_val lt 0)),where(n_val gt 0)] $
      - wt_sum_sigma2_norm[*,*,reverse(where(n_val lt 0)),reverse(where(n_val lt 0))])/4.
      
    wt_sum_sigma2_norm_sincos[*, *, 1:n_kz-1, 1:n_kz-1] = complex(0,1) * (wt_sum_sigma2_norm[*,*,where(n_val gt 0),where(n_val gt 0)] $
      + wt_sum_sigma2_norm[*,*,where(n_val gt 0),reverse(where(n_val lt 0))] $
      - wt_sum_sigma2_norm[*,*,reverse(where(n_val lt 0)),where(n_val gt 0)] $
      - wt_sum_sigma2_norm[*,*,reverse(where(n_val lt 0)),reverse(where(n_val lt 0))])/4.
      
    undefine, wt_sum_sigma2_norm
  endif
  
  
  if keyword_set(std_power) then begin
    ;; comov_dist_los goes from large to small z
    z_relative = dindgen(n_freq)*z_mpc_delta
    freq_kz_arr = rebin(reform(kz_mpc_orig, 1, n_elements(kz_mpc_orig)), n_freq, n_elements(kz_mpc_orig)) * rebin(z_relative, n_freq, n_elements(kz_mpc_orig))
    
    ;; construct FT matrix, hit sigma2 with it from both sides
    identity = diag_matrix([fltarr(n_freq)+1.])
    ft_matrix = fft(identity, dimension=1)*n_freq*z_mpc_delta
    
    sigma2_ft = fltarr(n_kx, n_ky, n_freq)
    for i=0, n_kx-1 do begin
      for j=0, n_ky-1 do begin
        if max(abs(sum_sigma2[i,j,*])) eq 0 then continue
        temp = diag_matrix(reform(sum_sigma2[i,j,*]))
        temp2 = shift(matrix_multiply(ft_matrix, matrix_multiply(temp, conj(ft_matrix), /btranspose)), n_freq/2, n_freq/2)
        
        sigma2_ft[i,j,*] = abs(diag_matrix(temp2))
      endfor
    endfor
    undefine, sum_sigma2
    
    sigma2_1 = fltarr(n_kx, n_ky, n_kz)
    sigma2_2 = fltarr(n_kx, n_ky, n_kz)
    sigma2_1[*, *, 0] = sigma2_ft[*,*,where(n_val eq 0)]
    sigma2_1[*, *, 1:n_kz-1] = (sigma2_ft[*,*, where(n_val gt 0)] + sigma2_ft[*,*, reverse(where(n_val lt 0))])/2.
    sigma2_2[*, *, 1:n_kz-1] = complex(0,1) * (sigma2_ft[*,*, where(n_val gt 0)] - sigma2_ft[*,*, reverse(where(n_val lt 0))])/2.
    undefine, sigma2_ft
    
    data_sum_1 = temporary(data_sum_cos)
    data_sum_2 = temporary(data_sum_sin)
    
    sim_noise_sum_1 = temporary(sim_noise_sum_cos)
    sim_noise_sum_2 = temporary(sim_noise_sum_sin)
    
    if nfiles gt 1 then begin
      data_diff_1 = temporary(data_diff_cos)
      data_diff_2 = temporary(data_diff_sin)
      
      sim_noise_diff_1 = temporary(sim_noise_diff_cos)
      sim_noise_diff_2 = temporary(sim_noise_diff_sin)
    endif
    
    
    if keyword_set(inverse_covar) then begin
    
      wt_power_norm_1 = temporary(wt_power_norm_coscos)
      wt_power_norm_2 = temporary(wt_power_norm_sinsin)
      undefine, wt_power_norm_cossin, wt_power_norm_sincos
      
      wt_sum_sigma2_norm_1 = temporary(wt_sum_sigma2_norm_coscos)
      wt_sum_sigma2_norm_2 = temporary(wt_sum_sigma2_norm_sinsin)
      undefine, wt_sum_sigma2_norm_cossin, wt_sum_sigma2_norm_sincos
      
    endif
    
  endif else begin
    ;; for new power calc, need cos2, sin2, cos*sin transforms
    ;; have to do this in a for loop for memory's sake
    covar_cos = fltarr(n_kx, n_ky, n_kz)
    covar_sin = fltarr(n_kx, n_ky, n_kz)
    covar_cross = fltarr(n_kx, n_ky, n_kz)
    
    ;; comov_dist_los goes from large to small z
    z_relative = dindgen(n_freq)*z_mpc_delta
    freq_kz_arr = rebin(reform(kz_mpc, 1, n_kz), n_freq, n_kz) * rebin(z_relative, n_freq, n_kz)
    
    cos_arr = cos(freq_kz_arr)
    sin_arr = sin(freq_kz_arr)
    
    sum_sigma2 = reform(sum_sigma2, n_kx*n_ky, n_freq)
    ;; doing 2 FTs so need 2 factors of z_mpc_delta/(2*!pi).
    ;; No multiplication by N b/c don't need to fix IDL FFT
    covar_cos = matrix_multiply(sum_sigma2, cos_arr^2d) * (z_mpc_delta)^2.
    covar_sin = matrix_multiply(sum_sigma2, sin_arr^2d) * (z_mpc_delta)^2.
    covar_cross = matrix_multiply(sum_sigma2, cos_arr*sin_arr) * (z_mpc_delta)^2.
    
    wh_0f = where(n_freq_contrib eq 0, count_0f)
    if count_0f gt 0 then begin
      covar_cos[wh_0f, *] = 0
      covar_sin[wh_0f, *] = 0
      covar_cross[wh_0f, *] = 0
    endif
    
    ;; reform to get back to n_kx, n_ky, n_kz dimensions
    covar_cos = reform(covar_cos, n_kx, n_ky, n_kz)
    covar_sin = reform(covar_sin, n_kx, n_ky, n_kz)
    covar_cross = reform(covar_cross, n_kx, n_ky, n_kz)
    
    undefine, sum_sigma2, freq_kz_arr, cos_arr, sin_arr
    
    ;; get rotation angle to diagonalize covariance block
    theta = atan(2.*covar_cross, covar_cos - covar_sin)/2.
    
    cos_theta = cos(theta)
    sin_theta = sin(theta)
    undefine, theta
    
    sigma2_1 = covar_cos*cos_theta^2. + 2.*covar_cross*cos_theta*sin_theta + covar_sin*sin_theta^2.
    sigma2_2 = covar_cos*sin_theta^2. - 2.*covar_cross*cos_theta*sin_theta + covar_sin*cos_theta^2.
    
    undefine, covar_cos, covar_sin, covar_cross
    
    data_sum_1 = data_sum_cos*cos_theta + data_sum_sin*sin_theta
    data_sum_2 = (-1d)*data_sum_cos*sin_theta + data_sum_sin*cos_theta
    undefine, data_sum_cos, data_sum_sin
    
    sim_noise_sum_1 = sim_noise_sum_cos*cos_theta + sim_noise_sum_sin*sin_theta
    sim_noise_sum_2 = (-1d)*sim_noise_sum_cos*sin_theta + sim_noise_sum_sin*cos_theta
    undefine, sim_noise_sum_cos, sim_noise_sum_sin
    
    if nfiles eq 2 then begin
      data_diff_1 = data_diff_cos*cos_theta + data_diff_sin*sin_theta
      data_diff_2 = (-1d)*data_diff_cos*sin_theta + data_diff_sin*cos_theta
      undefine, data_diff_cos, data_diff_sin
      
      sim_noise_diff_1 = sim_noise_diff_cos*cos_theta + sim_noise_diff_sin*sin_theta
      sim_noise_diff_2 = (-1d)*sim_noise_diff_cos*sin_theta + sim_noise_diff_sin*cos_theta
      undefine, sim_noise_diff_cos, sim_noise_diff_sin
    endif
    
    if keyword_set(inverse_covar) then begin
    
      wt_power_norm_1 = wt_power_norm_coscos*cos_theta^2. + wt_power_norm_cossin*cos_theta*sin_theta $
        + wt_power_norm_sincos*cos_theta*sin_theta + wt_power_norm_sinsin*sin_theta^2.
      wt_power_norm_2 = wt_power_norm_coscos*sin_theta^2. - wt_power_norm_cossin*cos_theta*sin_theta $
        - wt_power_norm_sincos*cos_theta*sin_theta + wt_power_norm_sinsin*cos_theta^2.
        
      undefine, wt_power_norm_coscos, wt_power_norm_cossin, wt_power_norm_sincos, wt_power_norm_sinsin
      
      wt_sum_sigma2_norm_1 = wt_sum_sigma2_norm_coscos*cos_theta^2. + wt_sum_sigma2_norm_cossin*cos_theta*sin_theta $
        + wt_sum_sigma2_norm_sincos*cos_theta*sin_theta + wt_sum_sigma2_norm_sinsin*sin_theta^2.
      wt_sum_sigma2_norm_2 = wt_sum_sigma2_norm_coscos*sin_theta^2. - wt_sum_sigma2_norm_cossin*cos_theta*sin_theta $
        - wt_sum_sigma2_norm_sincos*cos_theta*sin_theta + wt_sum_sigma2_norm_sinsin*cos_theta^2.
        
      undefine, wt_sum_sigma2_norm_coscos, wt_sum_sigma2_norm_cossin, wt_sum_sigma2_norm_sincos, wt_sum_sigma2_norm_sinsin
      
    endif
    
  endelse
  
  git, repo_path = ps_repository_dir(), result=kcube_git_hash
  git_hashes = {uvf:uvf_git_hashes, uvf_wt:uvf_wt_git_hashes, beam:beam_git_hashes, kcube:kcube_git_hash}
  
  ;; This is very repetative, but avoids getting warnings about not saving non-exisitent variables when they're not present
  if keyword_set(inverse_covar) then begin
    if n_elements(freq_flags) gt 0 then begin
      save, file = file_struct.kcube_savefile, data_sum_1, data_sum_2, data_diff_1, data_diff_2, $
        sim_noise_sum_1, sim_noise_sum_2, sim_noise_diff_1, sim_noise_diff_2, sigma2_1, sigma2_2, $
        wt_power_norm_1, wt_power_norm_2, wt_sum_sigma2_norm_1, wt_sum_sigma2_norm_2, n_val, $
        kx_mpc, ky_mpc, kz_mpc, kperp_lambda_conv, delay_params, hubble_param, n_freq_contrib, freq_mask, $
        vs_name, vs_mean, t_sys_meas, window_int, git_hashes, $
        wt_meas_ave, wt_meas_min, ave_weights, wt_ave_power_freq, ave_power_freq, wt_ave_power_uvf, ave_power_uvf
    endif else begin
      save, file = file_struct.kcube_savefile, data_sum_1, data_sum_2, data_diff_1, data_diff_2, $
        sim_noise_sum_1, sim_noise_sum_2, sim_noise_diff_1, sim_noise_diff_2, sigma2_1, sigma2_2, $
        wt_power_norm_1, wt_power_norm_2, wt_sum_sigma2_norm_1, wt_sum_sigma2_norm_2, n_val, $
        kx_mpc, ky_mpc, kz_mpc, kperp_lambda_conv, delay_params, hubble_param, n_freq_contrib, $
        vs_name, vs_mean, t_sys_meas, window_int, git_hashes, $
        wt_meas_ave, wt_meas_min, ave_weights, wt_ave_power_freq, ave_power_freq, wt_ave_power_uvf, ave_power_uvf
    endelse
  endif else begin
    if n_elements(freq_flags) gt 0 then begin
      save, file = file_struct.kcube_savefile, data_sum_1, data_sum_2, data_diff_1, data_diff_2, $
        sim_noise_sum_1, sim_noise_sum_2, sim_noise_diff_1, sim_noise_diff_2, sigma2_1, sigma2_2, n_val, $
        kx_mpc, ky_mpc, kz_mpc, kperp_lambda_conv, delay_params, hubble_param, n_freq_contrib, freq_mask, $
        vs_name, vs_mean, t_sys_meas, window_int, git_hashes, $
        wt_meas_ave, wt_meas_min, ave_weights, wt_ave_power_freq, ave_power_freq, wt_ave_power_uvf, ave_power_uvf
    endif else begin
      save, file = file_struct.kcube_savefile, data_sum_1, data_sum_2, data_diff_1, data_diff_2, $
        sim_noise_sum_1, sim_noise_sum_2, sim_noise_diff_1, sim_noise_diff_2, sigma2_1, sigma2_2, n_val, $
        kx_mpc, ky_mpc, kz_mpc, kperp_lambda_conv, delay_params, hubble_param, n_freq_contrib, $
        vs_name, vs_mean, t_sys_meas, window_int, git_hashes, $
        wt_meas_ave, wt_meas_min, ave_weights, wt_ave_power_freq, ave_power_freq, wt_ave_power_uvf, ave_power_uvf
    endelse
  endelse
  
end
