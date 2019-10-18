pro ps_kcube, file_struct, sim = sim, fix_sim_input = fix_sim_input, $
    uvf_input = uvf_input, freq_ch_range = freq_ch_range, freq_flags = freq_flags, $
    input_units = input_units, save_slices = save_slices, $
    refresh_options = refresh_options, uvf_options = uvf_options, $
    ps_options = ps_options

  if tag_exist(file_struct, 'nside') ne 0 then healpix = 1 else healpix = 0

  if keyword_set(uvf_input) or tag_exist(file_struct, 'uvf_savefile') eq 0 then begin
    uvf_input = 1
  endif else begin
    uvf_input = 0
  endelse

  nfiles = n_elements(file_struct.datafile)
  if tag_exist(file_struct, 'no_var') ne 0 then no_var = 1 else no_var = 0

  if n_elements(input_units) eq 0 then input_units = 'jansky'
  units_enum = ['jansky', 'mk']
  wh = where(units_enum eq input_units, count)
  if count eq 0 then message, 'input units not recognized, options are: ' + units_enum

  git, repo_path = ps_repository_dir(), result=this_run_git_hash

  frequencies = file_struct.frequencies

  if n_elements(freq_flags) ne 0 then freq_mask = file_struct.freq_mask

  if n_elements(freq_ch_range) ne 0 then begin
    n_freq_orig = n_elements(frequencies)
    frequencies = frequencies[min(freq_ch_range):max(freq_ch_range)]
  endif
  n_freq = n_elements(frequencies)

  z_mpc_mean = z_mpc(frequencies, hubble_param = hubble_param, f_delta = f_delta, $
    redshifts = redshifts, comov_dist_los = comov_dist_los, $
    z_mpc_delta = z_mpc_delta)

  kperp_lambda_conv = z_mpc_mean / (2.*!pi)
  delay_delta = 1e9/(n_freq*f_delta*1e6) ;; equivilent delay bin size for kparallel
  delay_max = delay_delta * n_freq/2.    ;; factor of 2 b/c of neg/positive
  delay_params = [delay_delta, delay_max]

  z_mpc_length = max(comov_dist_los) - min(comov_dist_los) + z_mpc_delta
  kz_mpc_range =  (2.*!pi) / (z_mpc_delta)
  kz_mpc_delta = (2.*!pi) / z_mpc_length

  ;; need to figure out where the zero bin is (for later computations)
  ;; The location of the zero bin also dictates what the most negative kz bin is (min_kz)
  ;; which we use to construct the kz_mpc_orig array
  ;; The calculation depends a bit on whether we have an odd or even number of frequencies.
  ;; The following calculation is taken directly from the IDL FFT documentation
  int_arr = findgen((n_freq - 1)/2) + 1
  is_n_freq_even = (n_freq mod 2) eq 0
  fft_shift_val = -(n_freq/2 + 1)
  if (is_n_freq_even) then begin
    kz_integers = [0.0, int_arr, n_freq/2, -n_freq/2 + int_arr]
  endif else begin
    kz_integers = [0.0, int_arr, fft_shift_val + int_arr]
  endelse
  kz_integers = shift(kz_integers, fft_shift_val)
  if where(kz_integers eq min(kz_integers)) ne 0 then begin
    message, 'something went very wrong with shifting!'
  endif

  min_kz = min(kz_integers)*kz_mpc_delta
  kz_mpc_orig = findgen(round(kz_mpc_range / kz_mpc_delta)) * kz_mpc_delta + min_kz

  n_kz = n_elements(kz_mpc_orig)
  if n_kz ne n_freq then message, 'something has gone wrong with kz_mpc calculation.'

  if input_units eq 'jansky' then begin
    ;; converting from Jy (in u,v,f) to mK*str (10^-26 * c^2 * 10^3/ (2*f^2*kb))
    conv_factor = double((299792458.)^2 / (2. * (frequencies*1e6)^2. * 1.38065))

    ;; convert from mk*str to mK*Mpc^2
    conv_factor = conv_factor * z_mpc_mean^2.

    ;; for adrian's weighting we have something in Jy/(wavelength^-1) = Jy * wavelength to start
    ;; converting from wavelength to Mpc uses 1/kperp_lambda_conv,
    ;; so multiply that by Jy -> mK*Mpc^2 conversion
    conv_factor_adrian = conv_factor / kperp_lambda_conv

  endif else if input_units eq 'mKstr' then begin
    conv_factor = 1. + fltarr(n_freq)
  endif else begin
    message, 'input_units not recognized'
  endelse

  t_sys = 280. * ((1+redshifts)/7.5)^2.3 / sqrt(2.) ;; from skew w/ stu + srt(2) for single pol
  eff_area = 21. ; m^2 -- from Aaron's memo
  df = file_struct.freq_resolution ; Hz -- native visibility resolution NOT cube resolution
  tau = file_struct.time_resolution ; seconds
  vis_sigma = (2. * (1.38065e-23) * 1e26) * t_sys / (eff_area * sqrt(df * tau)) ;; in Jy
  vis_sigma = float(vis_sigma)

  old_vis_sigma = temporary(vis_sigma)

  if tag_exist(file_struct, 'vis_noise') then begin
    vis_sigma_ian = file_struct.vis_noise
    ;; do a straight average over even/odd of sigma because we just want the
    ;; average noise (should actually be identical)
    if nfiles eq 2 then vis_sigma_ian = total(vis_sigma_ian, 1)/2.
  endif

  if n_elements(freq_ch_range) ne 0 then begin
    vis_sig_tag = number_formatter(384./n_freq_orig)
  endif else begin
    vis_sig_tag = number_formatter(384./n_freq)
  endelse
  vis_sigma_file = file_dirname(file_struct.savefile_froot, /mark_directory) + $
    'vis_sigma/vis_sigma_measured' + vis_sig_tag + '.sav'
  if file_valid(vis_sigma_file) then begin
    vis_sigma_adam = getvar_savefile(vis_sigma_file, 'vis_sigma')

    if n_elements(freq_ch_range) ne 0 then begin
      if n_elements(vis_sigma_adam) ne n_freq_orig then message, $
        'vis_sig file has incorrect number of frequency channels'
      vis_sigma_adam = vis_sigma_adam[min(freq_ch_range):max(freq_ch_range)]
    endif else begin
      if n_elements(vis_sigma_adam) ne n_freq then message, $
        'vis_sig file has incorrect number of frequency channels'
    endelse

    wh_nan = where(finite(vis_sigma_adam) eq 0, count_nan)
    if count_nan gt 0 then vis_sigma_adam[wh_nan] = 0
  endif

  if n_elements(vis_sigma_ian) gt 0 then begin
    if max(vis_sigma_ian) gt 5. then begin
      if n_elements(freq_ch_range) ne 0 then begin
        vis_sigma_ian = vis_sigma_ian[min(freq_ch_range):max(freq_ch_range)]
      endif
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

  vs_mean = mean(vis_sigma)

  ;; in K. sqrt(2) is because vis_sigma is for the real or imaginary part separately
  t_sys_meas = (eff_area * sqrt(df * tau) * vis_sigma * sqrt(2)) / $
    ((2. * (1.38065e-23) * 1e26))

  n_vis = reform(file_struct.n_vis)
  n_vis_freq = reform(file_struct.n_vis_freq)
  if n_elements(freq_ch_range) ne 0 then begin
    n_vis = total(n_vis_freq[*, min(freq_ch_range):max(freq_ch_range)], 2)
    n_vis_freq = n_vis_freq[*, min(freq_ch_range):max(freq_ch_range)]
  endif

  if healpix or not uvf_input then begin
    ps_image_to_uvf, file_struct, n_vis_freq, kx_rad_vals, ky_rad_vals, $
      freq_ch_range = freq_ch_range, freq_flags = freq_flags, no_var = no_var, $
      refresh_options = refresh_options, uvf_options = uvf_options

    n_kx = n_elements(kx_rad_vals)
    kx_rad_delta = kx_rad_vals[1] - kx_rad_vals[0]
    kx_mpc = temporary(kx_rad_vals) / z_mpc_mean
    kx_mpc_delta = kx_mpc[1] - kx_mpc[0]

    n_ky = n_elements(ky_rad_vals)
    ky_rad_delta = ky_rad_vals[1] - ky_rad_vals[0]
    ky_mpc = temporary(ky_rad_vals) / z_mpc_mean
    ky_mpc_delta = ky_mpc[1] - ky_mpc[0]

    if healpix then begin
      ;; Angular resolution is given in Healpix paper in units of arcminutes,
      ;; need to convert to radians
      ang_resolution = sqrt(3./!pi) * 3600./file_struct.nside * (1./60.) * (!pi/180.)
      pix_area_rad = ang_resolution^2. ;; by definition of ang. resolution in Healpix paper
    endif else pix_area_rad = (abs(file_struct.degpix) * !pi / 180d)^2.

    pix_area_mpc = pix_area_rad * z_mpc_mean^2.

    weights_cube1 = getvar_savefile(file_struct.uvf_weight_savefile[0], 'weights_cube')
    if nfiles eq 2 then begin
      weights_cube2 = getvar_savefile(file_struct.uvf_weight_savefile[1], 'weights_cube')
    endif

    if max(abs(weights_cube1)) eq 0 then begin
      message, 'weights cube is entirely zero.'
    endif
    if nfiles eq 2 then if max(abs(weights_cube2)) eq 0 then begin
      message, 'weights cube is entirely zero.'
    endif

    void = getvar_savefile(file_struct.uvf_weight_savefile[0], names = uvf_varnames)
    wh_hash = where(uvf_varnames eq 'uvf_wt_git_hash', count_hash)
    if count_hash gt 0 then begin
      uvf_wt_git_hashes = getvar_savefile(file_struct.uvf_weight_savefile[0], $
        'uvf_wt_git_hash')
      if nfiles eq 2 then begin
        uvf_wt_git_hashes = [uvf_wt_git_hashes, getvar_savefile(file_struct.uvf_weight_savefile[1], $
          'uvf_wt_git_hash')]
      endif
    endif else begin
      uvf_wt_git_hashes = strarr(nfiles)
    endelse

    if min(ky_mpc) lt 0 then begin
      ;; negative ky values haven't been cut yet
      ;; need to cut uvf cubes in half because image is real -- we'll cut negative ky
      weights_cube1 = weights_cube1[*, n_ky/2:n_ky-1,*]
      if nfiles eq 2 then weights_cube2 = weights_cube2[*, n_ky/2:n_ky-1,*]
    endif

  endif else begin
    ;; uvf_input

    weights_cube1 = getvar_savefile(file_struct.weightfile[0], file_struct.weightvar)
    if nfiles eq 2 then begin
      weights_cube2 = getvar_savefile(file_struct.weightfile[1], file_struct.weightvar)
    endif

    uvf_wt_git_hashes = strarr(nfiles)

    wt_size = size(weights_cube1)
    if wt_size[n_elements(wt_size)-2] eq 10 then begin
      ;; weights cube is a pointer
      dims2 = size(*weights_cube1[0], /dimension)
      iter=1
      while min(dims2) eq 0 and iter lt wt_size[2] do begin
        print, 'warning: some frequency slices have null pointers'
        dims2 = size(*weights_cube1[iter], /dimension)
        iter = iter+1
      endwhile

      temp = complex(fltarr([dims2, n_freq]))
      if nfiles eq 2 then temp2 = complex(fltarr([dims2, n_freq]))
      for i = 0, n_freq-1 do begin
        foo = *weights_cube1[file_struct.pol_index, i]
        if foo ne !null then begin
          temp[*,*,i] = temporary(foo)
          if nfiles eq 2 then temp2[*,*,i] = *weights_cube2[file_struct.pol_index, i]
        endif
      endfor
      undefine_fhd, weights_cube1, weights_cube2

      weights_cube1 = temporary(temp)
      if nfiles eq 2 then weights_cube2 = temporary(temp2)

    endif else dims2 = size(weights_cube1, /dimension)

    if max(abs(weights_cube1)) eq 0 then begin
      message, 'weights cube is entirely zero.'
    endif
    if nfiles eq 2 then if max(abs(weights_cube2)) eq 0 then begin
      message, 'weights cube is entirely zero.'
    endif

    n_kx = dims2[0]
    if abs(file_struct.kpix-1/(n_kx[0] * (abs(file_struct.degpix) * !pi / 180d)))/file_struct.kpix gt 1e-4 then begin
      message, 'Something has gone wrong with calculating uv pixel size'
    endif
    kx_mpc_delta = (2.*!pi)*file_struct.kpix / z_mpc_mean
    kx_mpc = (dindgen(n_kx)-n_kx/2) * kx_mpc_delta

    n_ky = dims2[1]
    ky_mpc_delta = (2.*!pi)*file_struct.kpix / z_mpc_mean
    ky_mpc = (dindgen(n_ky)-n_ky/2) * kx_mpc_delta

    ave_weights = total(total(abs(weights_cube1),2),1)/(n_kx*n_ky)
    if nfiles eq 2 then begin
      ave_weights = transpose([[ave_weights], [total(total(abs(weights_cube2),2),1)/(n_kx*n_ky)]])
    endif

    pix_area_rad = (!dtor*file_struct.degpix)^2.
    pix_area_mpc = pix_area_rad * z_mpc_mean^2.

    ;; get beam sorted out
    if tag_exist(file_struct, 'beam_savefile') then begin
      test_beam = file_valid(file_struct.beam_savefile)
      if min(test_beam) eq 0 then test_beam = check_old_path(file_struct, 'beam_savefile')

      if min(test_beam) eq 0 or refresh_options.refresh_beam then begin

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
            temp = arr * rebin(reform(nfvis_beam, 1, 1, n_elements(nfvis_beam)), $
              n_kx, n_ky, n_elements(nfvis_beam), /sample)
          endif else if max(arr) le file_struct.n_obs[i] then begin
            ;; beam is peak normalized to 1 for each obs, then summed over obs so peak is ~ n_obs
            temp = (arr/file_struct.n_obs[i]) * rebin(reform(nfvis_beam, 1, 1, n_freq), $
              n_kx, n_ky, n_freq, /sample)
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

    if tag_exist(uvf_options, 'uv_avg') then begin
      nkx_new = floor(n_kx / uvf_options.uv_avg)
      temp = complex(fltarr(nkx_new, n_ky, n_freq))
      if nfiles eq 2 then temp2 = complex(fltarr(nkx_new, n_ky, n_freq))
      temp_kx = fltarr(nkx_new)
      for i=0, nkx_new-1 do begin
        temp[i,*,*] = total(weights_cube1[i*uvf_options.uv_avg:(i+1)*uvf_options.uv_avg-1,*,*], 1) / uvf_options.uv_avg
        if nfiles eq 2 then begin
          temp2[i,*,*] = total(weights_cube2[i*uvf_options.uv_avg:(i+1)*uvf_options.uv_avg-1,*,*], 1) / uvf_options.uv_avg
        endif
        temp_kx[i] = total(kx_mpc[i*uvf_options.uv_avg:(i+1)*uvf_options.uv_avg-1]) / uvf_options.uv_avg
      endfor

      nky_new = floor(n_ky / uvf_options.uv_avg)
      temp3 = complex(fltarr(nkx_new, nky_new, n_freq))
      if nfiles eq 2 then temp4 = complex(fltarr(nkx_new, nky_new, n_freq))
      temp_ky = fltarr(nky_new)
      for i=0, nky_new-1 do begin
        temp3[*,i,*] = total(temp[*,i*uvf_options.uv_avg:(i+1)*uvf_options.uv_avg-1,*], 2) / uvf_options.uv_avg
        if nfiles eq 2 then temp4[*,i,*] = total(temp2[*,i*uvf_options.uv_avg:(i+1)*uvf_options.uv_avg-1,*], 2) / uvf_options.uv_avg
        temp_ky[i] = total(ky_mpc[i*uvf_options.uv_avg:(i+1)*uvf_options.uv_avg-1]) / uvf_options.uv_avg
      endfor
      undefine, temp, temp2

      ;; averging reduces the value of total(weights) ~ n_vis needed for the
      ;; window int calculation
      n_vis = n_vis/(uvf_options.uv_avg)^2.

      weights_cube1 = temporary(temp3)
      if nfiles eq 2 then weights_cube2 = temporary(temp4)

      if max(abs(weights_cube1)) eq 0 then message, 'weights cube is entirely zero.'
      if nfiles eq 2 then if max(abs(weights_cube2)) eq 0 then begin
        message, 'weights cube is entirely zero.'
      endif

      kx_mpc = temporary(temp_kx)
      ky_mpc = temporary(temp_ky)
      n_kx = nkx_new
      n_ky = nky_new
    endif

    if tag_exist(uvf_options, 'uv_img_clip') then begin
      kx_mpc_delta_old = kx_mpc_delta
      ky_mpc_delta_old = ky_mpc_delta
      temp = shift(fft(fft(shift(weights_cube1,dims2[0]/2,dims2[1]/2,0), $
        dimension=1),dimension=2),dims2[0]/2,dims2[1]/2,0)
      temp = temp[(dims2[0]/2)-(dims2[0]/uvf_options.uv_img_clip)/2:(dims2[0]/2)+(dims2[0]/uvf_options.uv_img_clip)/2-1, *, *]
      temp = temp[*, (dims2[1]/2)-(dims2[1]/uvf_options.uv_img_clip)/2:(dims2[1]/2)+(dims2[1]/uvf_options.uv_img_clip)/2-1, *]
      if nfiles eq 2 then begin
        temp2 = shift(fft(fft(shift(weights_cube2,dims2[0]/2,dims2[1]/2,0), dimension=1), $
          dimension=2),dims2[0]/2,dims2[1]/2,0)
        temp2 = temp2[(dims2[0]/2)-(dims2[0]/uvf_options.uv_img_clip)/2:(dims2[0]/2)+(dims2[0]/uvf_options.uv_img_clip)/2-1, *, *]
        temp2 = temp2[*, (dims2[1]/2)-(dims2[1]/uvf_options.uv_img_clip)/2:(dims2[1]/2)+(dims2[1]/uvf_options.uv_img_clip)/2-1, *]
      endif
      temp_dims = size(temp, /dimension)
      n_kx = temp_dims[0]
      n_ky = temp_dims[1]
      kx_mpc_delta = kx_mpc_delta * uvf_options.uv_img_clip
      kx_mpc = (dindgen(n_kx)-n_kx/2) * kx_mpc_delta
      ky_mpc_delta = ky_mpc_delta * uvf_options.uv_img_clip
      ky_mpc = (dindgen(n_kx)-n_kx/2) * ky_mpc_delta

      temp = shift(fft(fft(shift(temp,n_kx/2,n_ky/2,0), dimension=1, /inverse), $
        dimension=2, /inverse),n_kx/2,n_ky/2,0)
      if nfiles eq 2 then begin
        temp2 = shift(fft(fft(shift(temp2,n_kx/2,n_ky/2,0), dimension=1, /inverse), $
          dimension=2, /inverse), n_kx/2, n_ky/2,0)
      endif

      weights_cube1 = temp
      if nfiles eq 2 then weights_cube2 = temp2

      if max(abs(weights_cube1)) eq 0 then message, 'weights cube is entirely zero.'
      if nfiles eq 2 then if max(abs(weights_cube2)) eq 0 then begin
        message, 'weights cube is entirely zero.'
      endif
    endif

    if n_elements(freq_ch_range) ne 0 then begin
      weights_cube1 = weights_cube1[*, *, min(freq_ch_range):max(freq_ch_range)]
      if nfiles eq 2 then begin
        weights_cube2 = weights_cube2[*, *, min(freq_ch_range):max(freq_ch_range)]
      endif
    endif
    if n_elements(freq_flags) ne 0 then begin
      flag_arr = rebin(reform(freq_mask, 1, 1, n_freq), $
        size(weights_cube1,/dimension), /sample)
      weights_cube1 = weights_cube1 * flag_arr
      if nfiles eq 2 then weights_cube2 = weights_cube2 * flag_arr
    endif

    ;; need to cut uvf cubes in half because image is real -- we'll cut negative ky
    weights_cube1 = weights_cube1[*, n_ky/2:n_ky-1,*]
    if nfiles eq 2 then weights_cube2 = weights_cube2[*, n_ky/2:n_ky-1,*]
  endelse

  ;; Also need to drop 1/2 of ky=0 line -- drop ky=0, kx<0
  weights_cube1[0:n_kx/2-1, 0, *] = 0
  if nfiles eq 2 then weights_cube2[0:n_kx/2-1, 0, *] = 0

  if max(abs(weights_cube1)) eq 0 then message, 'weights cube is entirely zero.'
  if nfiles eq 2 then if max(abs(weights_cube2)) eq 0 then begin
    message, 'weights cube is entirely zero.'
  endif

  if not no_var then begin
    if healpix or not keyword_set(uvf_input) then begin
      variance_cube1 = getvar_savefile(file_struct.uvf_weight_savefile[0], 'variance_cube')
      if nfiles eq 2 then begin
        variance_cube2 = getvar_savefile(file_struct.uvf_weight_savefile[1], 'variance_cube')
      endif

      if max(abs(variance_cube1)) eq 0 then message, 'variance cube is entirely zero.'
      if nfiles eq 2 then if max(abs(variance_cube2)) eq 0 then begin
        message, 'variance cube is entirely zero.'
      endif

      if min(ky_mpc) lt 0 then begin
        ;; calculate integral of window function before cut for comparison
        if nfiles eq 2 then window_int_orig = [total(variance_cube1)*pix_area_rad/n_vis[0], $
          total(variance_cube2)*pix_area_rad/n_vis[1]] $
        else window_int_orig = total(variance_cube1)*pix_area_rad/n_vis[0]

        ;; negative ky values haven't been cut yet
        ;; need to cut uvf cubes in half because image is real -- we'll cut negative ky
        variance_cube1 = variance_cube1[*, n_ky/2:n_ky-1,*]
        if nfiles eq 2 then variance_cube2 = variance_cube2[*, n_ky/2:n_ky-1,*]

        if max(abs(variance_cube1)) eq 0 then message, 'variance cube is entirely zero.'
        if nfiles eq 2 then if max(abs(variance_cube2)) eq 0 then begin
          message, 'variance cube is entirely zero.'
        endif
      endif

      ;; Also need to drop 1/2 of ky=0 line -- drop ky=0, kx<0
      variance_cube1[0:n_kx/2-1, 0, *] = 0
      if nfiles eq 2 then variance_cube2[0:n_kx/2-1, 0, *] = 0

      ;; calculate integral of window function (use pix_area_rad for FT normalization)
      ;; already cut out negative ky, so multiply by 2
      if nfiles eq 2 then window_int = 2*[total(variance_cube1)*pix_area_rad/n_vis[0], $
        total(variance_cube2)*pix_area_rad/n_vis[1]] $
      else window_int = 2*total(variance_cube1)*pix_area_rad/n_vis[0]

      if tag_exist(file_struct, 'beam_savefile') then begin
        beam1 = getvar_savefile(file_struct.beam_savefile[0], 'avg_beam')
        if nfiles eq 2 then beam2 = getvar_savefile(file_struct.beam_savefile[1], 'avg_beam')

        void = getvar_savefile(file_struct.beam_savefile[0], names = uvf_varnames)
        wh_hash = where(uvf_varnames eq 'beam_git_hash', count_hash)
        if count_hash gt 0 then begin
          beam_git_hashes = getvar_savefile(file_struct.beam_savefile[0], 'beam_git_hash')
          if nfiles eq 2 then begin
            beam_git_hashes = [beam_git_hashes, $
              getvar_savefile(file_struct.beam_savefile[1], 'beam_git_hash')]
          endif
        endif else beam_git_hashes = strarr(nfiles)

        if nfiles eq 2 then begin
          window_int_beam = [total(beam1), total(beam2)]*pix_area_mpc*(z_mpc_delta * n_freq)
        endif else begin
          window_int_beam = total(beam1)*pix_area_mpc*(z_mpc_delta * n_freq)
        endelse

        volume_factor = total(beam1*0+1.)*pix_area_mpc*(z_mpc_delta * n_freq)
        if nfiles eq 2 then volume_factor = fltarr(2) + volume_factor

        bandwidth_factor = z_mpc_delta * n_freq
        if nfiles eq 2 then bandwidth_factor = fltarr(2) + bandwidth_factor

      endif else beam_git_hashes = ''

      if tag_exist(file_struct, 'beam_int') then begin
        if nfiles eq 2 then begin
          if n_elements(freq_ch_range) ne 0 then begin
            ave_beam_int = total(file_struct.beam_int[*,freq_ch_range] * file_struct.n_vis_freq[*,freq_ch_range], 2) $
              / total(file_struct.n_vis_freq[*,freq_ch_range], 2)
          endif else begin
            ave_beam_int = total(file_struct.beam_int * file_struct.n_vis_freq, 2) $
              / total(file_struct.n_vis_freq, 2)
          endelse
        endif else begin
          if n_elements(freq_ch_range) ne 0 then begin
            ave_beam_int = total(file_struct.beam_int[freq_ch_range] * file_struct.n_vis_freq[freq_ch_range]) $
              / total(file_struct.n_vis_freq[freq_ch_range])
          endif else begin
            ave_beam_int = total(file_struct.beam_int * file_struct.n_vis_freq) $
              / total(file_struct.n_vis_freq)
          endelse
        endelse

        ;; convert rad -> Mpc^2, multiply by depth in Mpc
        window_int_beam_obs = ave_beam_int * z_mpc_mean^2. * (z_mpc_delta * n_freq)
      endif


    endif else begin
      ;; uvf_input
      variance_cube1 = getvar_savefile(file_struct.variancefile[0], file_struct.variancevar)
      if nfiles eq 2 then begin
        variance_cube2 = getvar_savefile(file_struct.variancefile[1], file_struct.variancevar)
      endif

      var_size = size(variance_cube1)
      if var_size[n_elements(var_size)-2] eq 10 then begin
        ;; variance cube is a pointer
        dims2 = size(*variance_cube1[0], /dimension)
        iter=1
        while min(dims2) eq 0 and iter lt var_size[2] do begin
          print, 'warning: some frequency slices have null pointers'
          dims2 = size(*variance_cube1[iter], /dimension)
          iter = iter+1
        endwhile
        temp = complex(fltarr([dims2, n_freq]))
        if nfiles eq 2 then temp2 = complex(fltarr([dims2, n_freq]))
        for i = 0, n_freq-1 do begin
          foo = *variance_cube1[file_struct.pol_index, i]
          if foo ne !null then begin
            temp[*,*,i] = foo
            if nfiles eq 2 then temp2[*,*,i] = *variance_cube2[file_struct.pol_index, i]
          endif
        endfor
        undefine_fhd, variance_cube1, variance_cube2

        variance_cube1 = temporary(temp)
        if nfiles eq 2 then variance_cube2 = temporary(temp2)
      endif

      if max(abs(variance_cube1)) eq 0 then message, 'variance cube is entirely zero.'
      if nfiles eq 2 then if max(abs(variance_cube2)) eq 0 then begin
        message, 'variance cube is entirely zero.'
      endif

      if tag_exist(uvf_options, 'uv_avg') then begin
        temp = complex(fltarr(nkx_new, dims2[1], n_freq))
        if nfiles eq 2 then temp2 = complex(fltarr(nkx_new, dims2[1], n_freq))
        for i=0, nkx_new-1 do begin
          temp[i,*,*] = total(variance_cube1[i*uvf_options.uv_avg:(i+1)*uvf_options.uv_avg-1,*,*], 1) / uvf_options.uv_avg
          if nfiles eq 2 then temp2[i,*,*] = total(variance_cube2[i*uvf_options.uv_avg:(i+1)*uvf_options.uv_avg-1,*,*], 1) / uvf_options.uv_avg
        endfor

        temp3 = complex(fltarr(nkx_new, nky_new, n_freq))
        if nfiles eq 2 then temp4 = complex(fltarr(nkx_new, nky_new, n_freq))
        for i=0, nky_new-1 do begin
          temp3[*,i,*] = total(temp[*,i*uvf_options.uv_avg:(i+1)*uvf_options.uv_avg-1,*], 2) / uvf_options.uv_avg
          if nfiles eq 2 then temp4[*,i,*] = total(temp2[*,i*uvf_options.uv_avg:(i+1)*uvf_options.uv_avg-1,*], 2) / uvf_options.uv_avg
        endfor
        undefine, temp, temp2

        variance_cube1 = temporary(temp3)
        if nfiles eq 2 then variance_cube2 = temporary(temp4)

        if max(abs(variance_cube1)) eq 0 then message, 'variance cube is entirely zero.'
        if nfiles eq 2 then if max(abs(variance_cube2)) eq 0 then begin
          message, 'variance cube is entirely zero.'
        endif
      endif

      if tag_exist(uvf_options, 'uv_img_clip') then begin
        temp = shift(fft(fft(shift(variance_cube1,dims2[0]/2,dims2[1]/2,0), dimension=1), dimension=2),dims2[0]/2,dims2[1]/2,0)
        temp = temp[(dims2[0]/2)-(dims2[0]/uvf_options.uv_img_clip)/2:(dims2[0]/2)+(dims2[0]/uvf_options.uv_img_clip)/2-1, *, *]
        temp = temp[*, (dims2[1]/2)-(dims2[1]/uvf_options.uv_img_clip)/2:(dims2[1]/2)+(dims2[1]/uvf_options.uv_img_clip)/2-1, *]
        temp = shift(fft(fft(shift(temp,n_kx/2,n_ky/2,0), dimension=1, /inverse), dimension=2, /inverse),n_kx/2,n_ky/2,0)
        if nfiles eq 2 then begin
          temp2 = shift(fft(fft(shift(variance_cube2,dims2[0]/2,dims2[1]/2,0), dimension=1), dimension=2),dims2[0]/2,dims2[1]/2,0)
          temp2 = temp2[(dims2[0]/2)-(dims2[0]/uvf_options.uv_img_clip)/2:(dims2[0]/2)+(dims2[0]/uvf_options.uv_img_clip)/2-1, *, *]
          temp2 = temp2[*, (dims2[1]/2)-(dims2[1]/uvf_options.uv_img_clip)/2:(dims2[1]/2)+(dims2[1]/uvf_options.uv_img_clip)/2-1, *]
          temp2 = shift(fft(fft(shift(temp2,n_kx/2,n_ky/2,0), dimension=1, /inverse), dimension=2, /inverse),n_kx/2,n_ky/2,0)
        endif

        variance_cube1 = temp
        if nfiles eq 2 then variance_cube2 = temp2

        if max(abs(variance_cube1)) eq 0 then begin
          message, 'variance cube is entirely zero.'
        endif
        if nfiles eq 2 then if max(abs(variance_cube2)) eq 0 then begin
          message, 'variance cube is entirely zero.'
        endif
      endif

      if max(abs(imaginary(variance_cube1))) gt 0 then begin
        print, 'variance_cube1 is not real, using absolute value'
        variance_cube1 = abs(variance_cube1)
      endif else variance_cube1 = real_part(variance_cube1)
      if nfiles eq 2 then if max(abs(imaginary(variance_cube2))) gt 0 then begin
        print, 'variance_cube2 is not real, using absolute value'
        variance_cube2 = abs(variance_cube2)
      endif else variance_cube2 = real_part(variance_cube2)

      ;; need to cut uvf cubes in half because image is real -- we'll cut negative ky
      variance_cube1 = variance_cube1[*, n_ky/2:n_ky-1,*]
      if nfiles eq 2 then variance_cube2 = variance_cube2[*, n_ky/2:n_ky-1,*]

      if n_elements(freq_ch_range) ne 0 then begin
        variance_cube1 = variance_cube1[*, *, min(freq_ch_range):max(freq_ch_range)]
        if nfiles eq 2 then variance_cube2 = variance_cube2[*, *, min(freq_ch_range):max(freq_ch_range)]
      endif
      if n_elements(freq_flags) ne 0 then begin
        flag_arr = rebin(reform(freq_mask, 1, 1, n_freq), $
          size(variance_cube1,/dimension), /sample)
        variance_cube1 = variance_cube1 * flag_arr
        if nfiles eq 2 then variance_cube2 = variance_cube2 * flag_arr
      endif

      ;; Also need to drop 1/2 of ky=0 line -- drop ky=0, kx<0
      variance_cube1[0:n_kx/2-1, 0, *] = 0
      if nfiles eq 2 then variance_cube2[0:n_kx/2-1, 0, *] = 0

      if max(abs(variance_cube1)) eq 0 then message, 'variance cube is entirely zero.'
      if nfiles eq 2 then if max(abs(variance_cube2)) eq 0 then begin
        message, 'variance cube is entirely zero.'
      endif

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
          if nfiles eq 2 then begin
            beam_git_hashes = [beam_git_hashes, getvar_savefile(file_struct.beam_savefile[1], $
              'beam_git_hash')]
          endif
        endif else beam_git_hashes = strarr(nfiles)

        if nfiles eq 2 then begin
          window_int_beam = [total(beam1), total(beam2)]*pix_area_mpc*(z_mpc_delta * n_freq)
        endif else begin
          window_int_beam = total(beam1)*pix_area_mpc*(z_mpc_delta * n_freq)
        endelse

        volume_factor = total(beam1*0+1.)*pix_area_mpc*(z_mpc_delta * n_freq)
        if nfiles eq 2 then volume_factor = fltarr(2) + volume_factor

        bandwidth_factor = z_mpc_delta * n_freq
        if nfiles eq 2 then bandwidth_factor = fltarr(2) + bandwidth_factor

        undefine, beam1, beam2
      endif else beam_git_hashes = ''

      if tag_exist(file_struct, 'beam_int') then begin
        if nfiles eq 2 then begin
          if N_elements(freq_ch_range) ne 0 then begin
            ave_beam_int = total(file_struct.beam_int[*, freq_ch_range] * file_struct.n_vis_freq[*, freq_ch_range], 2) / $
              total(file_struct.n_vis_freq[*, freq_ch_range], 2)
          endif else begin
            ave_beam_int = total(file_struct.beam_int * file_struct.n_vis_freq, 2) / $
              total(file_struct.n_vis_freq, 2)
          endelse
        endif else begin
          if N_elements(freq_ch_range) ne 0 then begin
            ave_beam_int = total(file_struct.beam_int[freq_ch_range] * file_struct.n_vis_freq[freq_ch_range]) / $
              total(file_struct.n_vis_freq[freq_ch_range])
          endif else begin
            ave_beam_int = total(file_struct.beam_int * file_struct.n_vis_freq) / $
              total(file_struct.n_vis_freq)
          endelse
        endelse

        ;; convert rad -> Mpc^2, multiply by depth in Mpc
        window_int_beam_obs = ave_beam_int * z_mpc_mean^2. * (z_mpc_delta * n_freq)
      endif

    endelse

  endif

  ;; now get data cubes
  if healpix or not keyword_set(uvf_input) then begin
    data_cube1 = getvar_savefile(file_struct.uvf_savefile[0], 'data_cube')
    if nfiles eq 2 then data_cube2 = getvar_savefile(file_struct.uvf_savefile[1], 'data_cube')

    void = getvar_savefile(file_struct.uvf_savefile[0], names = uvf_varnames)
    wh_hash = where(uvf_varnames eq 'uvf_git_hash', count_hash)
    if count_hash gt 0 then begin
      uvf_git_hashes = getvar_savefile(file_struct.uvf_savefile[0], 'uvf_git_hash')
      if nfiles eq 2 then begin
        uvf_git_hashes = [uvf_git_hashes, getvar_savefile(file_struct.uvf_savefile[1], $
          'uvf_git_hash')]
      endif
    endif else uvf_git_hashes = strarr(nfiles)

    if min(ky_mpc) lt 0 then begin
      ;; negative ky values haven't been cut yet
      ;; need to cut uvf cubes in half because image is real -- we'll cut negative ky
      data_cube1 = data_cube1[*, n_ky/2:n_ky-1,*]
      if nfiles eq 2 then data_cube2 = data_cube2[*, n_ky/2:n_ky-1,*]

      ky_mpc = ky_mpc[n_ky/2:n_ky-1]
      n_ky = n_elements(ky_mpc)
    endif

    if max(abs(data_cube1)) eq 0 then message, 'data cube is entirely zero.'
    if nfiles eq 2 then if max(abs(data_cube2)) eq 0 then message, 'data cube is entirely zero.'

  endif else begin
    ;; uvf_input
    datavar = strupcase(file_struct.datavar)
    if datavar eq '' then begin
      ;; working with a 'derived' cube that is constructed from uvf_savefiles
      input_uvf_files = reform(file_struct.derived_uvf_inputfiles, nfiles, 2)
      input_uvf_varname = reform(file_struct.derived_uvf_varname, nfiles, 2)

      if healpix or not keyword_set(uvf_input) then begin
        input_uvf_wtfiles = file_struct.uvf_weight_savefile
      endif
    endif

    if datavar eq '' then begin
      ;; working with a 'derived' cube (ie residual cube) that is constructed from other cubes
      ;; need to cut uvf cubes in half because image is real -- we'll cut negative ky
      dirty_cube1 = getvar_savefile(input_uvf_files[0,0], input_uvf_varname[0,0])
      model_cube1 = getvar_savefile(input_uvf_files[0,1], input_uvf_varname[0,1])
      if nfiles eq 2 then begin
        dirty_cube2 = getvar_savefile(input_uvf_files[1,0], input_uvf_varname[1,0])
        model_cube2 = getvar_savefile(input_uvf_files[1,1], input_uvf_varname[1,1])
      endif

      dirty_size = size(dirty_cube1)
      if dirty_size[n_elements(dirty_size)-2] eq 10 then begin
        ;; dirty cube is a pointer
        ;; handle null frequency pointers carefully
        for ind_i=0, product(dirty_size[1:dirty_size[0]]) - 1 do begin
          dims2 = size(*dirty_cube1[ind_i], /dimension)
          if total(dims2) NE 0 then break
        endfor
        temp = complex(fltarr([dims2, n_freq]))
        temp_m = complex(fltarr([dims2, n_freq]))
        if nfiles eq 2 then begin
          temp2 = complex(fltarr([dims2, n_freq]))
          temp_m2 = complex(fltarr([dims2, n_freq]))
        endif
        for i = 0, n_freq-1 do begin
          if *dirty_cube1[file_struct.pol_index, i] NE !NULL then $
            temp[*,*,i] = *dirty_cube1[file_struct.pol_index, i]
          if *model_cube1[file_struct.pol_index, i] NE !NULL then $
            temp_m[*,*,i] = *model_cube1[file_struct.pol_index, i]
          if nfiles eq 2 then begin
            if *dirty_cube2[file_struct.pol_index, i] NE !NULL then $
              temp2[*,*,i] = *dirty_cube2[file_struct.pol_index, i]
            if *model_cube2[file_struct.pol_index, i] NE !NULL then $
              temp_m2[*,*,i] = *model_cube2[file_struct.pol_index, i]
          endif
        endfor
        undefine_fhd, dirty_cube1, model_cube1, dirty_cube2, model_cube2

        dirty_cube1 = temporary(temp)
        model_cube1 = temporary(temp_m)
        if nfiles eq 2 then begin
          dirty_cube2 = temporary(temp2)
          model_cube2 = temporary(temp_m2)
        endif
      endif

      data_cube1 = temporary(dirty_cube1) - temporary(model_cube1)

      if nfiles eq 2 then begin
        data_cube2 = temporary(dirty_cube2) - temporary(model_cube2)
      endif

      if max(abs(data_cube1)) eq 0 then message, 'data cube is entirely zero.'
      if nfiles eq 2 then if max(abs(data_cube2)) eq 0 then message, 'data cube is entirely zero.'

      void = getvar_savefile(input_uvf_files[0,0], names = uvf_varnames)
      wh_hash = where(uvf_varnames eq 'uvf_git_hash', count_hash)
      if count_hash gt 0 then begin
        uvf_git_hashes = getvar_savefile(input_uvf_files[0,0], 'uvf_git_hash')
        if nfiles eq 2 then begin
          uvf_git_hashes = [uvf_git_hashes, getvar_savefile(input_uvf_files[0,1], $
            'uvf_git_hash')]
        endif
      endif else uvf_git_hashes = strarr(nfiles)

    endif else begin
      data_cube1 = getvar_savefile(file_struct.datafile[0], file_struct.datavar)
      if nfiles eq 2 then begin
        data_cube2 = getvar_savefile(file_struct.datafile[1], file_struct.datavar)
      endif

      uvf_git_hashes = strarr(nfiles)

      data_size = size(data_cube1)
      if data_size[n_elements(data_size)-2] eq 10 then begin
        ;; data cube is a pointer
        dims2 = size(*data_cube1[0], /dimension)
        iter=1
        while min(dims2) eq 0 and iter lt data_size[2] do begin
          print, 'warning: some frequency slices have null pointers'
          dims2 = size(*data_cube1[iter], /dimension)
          iter = iter+1
        endwhile

        temp = complex(fltarr([dims2, n_freq]))
        if nfiles eq 2 then temp2 = complex(fltarr([dims2, n_freq]))
        for i = 0, n_freq-1 do begin
          foo = *data_cube1[file_struct.pol_index, i]
          if foo ne !null then begin
            temp[*,*,i] = foo
            if nfiles eq 2 then temp2[*,*,i] = *data_cube2[file_struct.pol_index, i]
          endif
        endfor
        undefine_fhd, data_cube1, data_cube2

        data_cube1 = temporary(temp)
        if nfiles eq 2 then data_cube2 = temporary(temp2)
      endif

      if max(abs(data_cube1)) eq 0 then message, 'data cube is entirely zero.'
      if nfiles eq 2 then if max(abs(data_cube2)) eq 0 then message, 'data cube is entirely zero.'

    endelse

    if tag_exist(uvf_options, 'uv_avg') then begin
      temp = complex(fltarr(nkx_new, dims2[1], n_freq))
      if nfiles eq 2 then temp2 = complex(fltarr(nkx_new, dims2[1], n_freq))
      for i=0, nkx_new-1 do begin
        temp[i, *, *] = total(data_cube1[i * uvf_options.uv_avg:(i+1) * uvf_options.uv_avg-1, *, *], 1) / $
          uvf_options.uv_avg
        if nfiles eq 2 then begin
          temp2[i, *, *] = total(data_cube2[i * uvf_options.uv_avg:(i+1) * uvf_options.uv_avg-1, *, *], 1) / $
            uvf_options.uv_avg
        endif
      endfor

      temp3 = complex(fltarr(nkx_new, nky_new, n_freq))
      if nfiles eq 2 then temp4 = complex(fltarr(nkx_new, nky_new, n_freq))
      for i=0, nky_new-1 do begin
        temp3[*,i,*] = total(temp[*,i*uvf_options.uv_avg:(i+1)*uvf_options.uv_avg-1,*], 2) / $
          uvf_options.uv_avg
        if nfiles eq 2 then begin
          temp4[*,i,*] = total(temp2[*,i*uvf_options.uv_avg:(i+1)*uvf_options.uv_avg-1,*], 2) / $
            uvf_options.uv_avg
        endif
      endfor
      undefine, temp, temp2

      data_cube1 = temporary(temp3)
      if nfiles eq 2 then data_cube2 = temporary(temp4)

      if max(abs(data_cube1)) eq 0 then message, 'data cube is entirely zero.'
      if nfiles eq 2 then if max(abs(data_cube2)) eq 0 then message, 'data cube is entirely zero.'
    endif

    if tag_exist(uvf_options, 'uv_img_clip') then begin
      temp = shift(fft(fft(shift(data_cube1,dims2[0]/2,dims2[1]/2,0), dimension=1), $
        dimension=2), dims2[0]/2, dims2[1]/2,0)
      temp = temp[(dims2[0]/2)-(dims2[0]/uvf_options.uv_img_clip)/2:(dims2[0]/2) + $
        (dims2[0]/uvf_options.uv_img_clip)/2-1, *, *]
      temp = temp[*, (dims2[1]/2)-(dims2[1]/uvf_options.uv_img_clip)/2:(dims2[1]/2) + $
        (dims2[1]/uvf_options.uv_img_clip)/2-1, *]
      temp = shift(fft(fft(shift(temp,n_kx/2,n_ky/2,0), dimension=1, /inverse), $
        dimension=2, /inverse),n_kx/2,n_ky/2,0)
      if nfiles eq 2 then begin
        temp2 = shift(fft(fft(shift(data_cube2,dims2[0]/2,dims2[1]/2,0), dimension=1), $
          dimension=2),dims2[0]/2,dims2[1]/2,0)
        temp2 = temp2[(dims2[0]/2)-(dims2[0]/uvf_options.uv_img_clip)/2:(dims2[0]/2) + $
          (dims2[0]/uvf_options.uv_img_clip)/2-1, *, *]
        temp2 = temp2[*, (dims2[1]/2)-(dims2[1]/uvf_options.uv_img_clip)/2:(dims2[1]/2) + $
          (dims2[1]/uvf_options.uv_img_clip)/2-1, *]
        temp2 = shift(fft(fft(shift(temp2,n_kx/2,n_ky/2,0), dimension=1, /inverse), $
          dimension=2, /inverse),n_kx/2,n_ky/2,0)
      endif

      data_cube1 = temp
      if nfiles eq 2 then data_cube2 = temp2

      if max(abs(data_cube1)) eq 0 then message, 'data cube is entirely zero.'
      if nfiles eq 2 then if max(abs(data_cube2)) eq 0 then message, 'data cube is entirely zero.'
    endif

    ;; need to cut uvf cubes in half because image is real -- we'll cut negative ky
    data_cube1 = data_cube1[*, n_ky/2:n_ky-1,*]
    if nfiles eq 2 then data_cube2 = data_cube2[*, n_ky/2:n_ky-1,*]

    if n_elements(freq_ch_range) ne 0 then begin
      data_cube1 = data_cube1[*, *, min(freq_ch_range):max(freq_ch_range)]
      if nfiles eq 2 then data_cube2 = data_cube2[*, *, min(freq_ch_range):max(freq_ch_range)]
    endif
    if n_elements(freq_flags) ne 0 then begin
      flag_arr = rebin(reform(freq_mask, 1, 1, n_freq), size(data_cube1,/dimension), /sample)
      data_cube1 = data_cube1 * flag_arr
      if nfiles eq 2 then data_cube2 = data_cube2 * flag_arr
    endif

    ky_mpc = ky_mpc[n_ky/2:n_ky-1]
    n_ky = n_elements(ky_mpc)
  endelse

  ;; Also need to drop 1/2 of ky=0 line -- drop ky=0, kx<0
  data_cube1[0:n_kx/2-1, 0, *] = 0
  if nfiles eq 2 then data_cube2[0:n_kx/2-1, 0, *] = 0

  if keyword_set(sim) and keyword_set(fix_sim_input) then begin
    ;; fix for some sims that used saved uvf inputs without correcting for factor of 2 loss
    data_cube1 = data_cube1 * 2.
    if nfiles eq 2 then data_cube2 = data_cube2 * 2.
  endif

  if max(abs(data_cube1)) eq 0 then message, 'data cube is entirely zero.'
  if nfiles eq 2 then if max(abs(data_cube2)) eq 0 then message, 'data cube is entirely zero.'

  if save_slices then begin
    ;; save some slices of the raw data cube (before dividing by weights) & weights
    for i=0, nfiles-1 do begin
      if i eq 0 then data_cube = data_cube1 else data_cube = data_cube2
      if i eq 0 then weights_cube = weights_cube1 else weights_cube = weights_cube2
      if i eq 0 then var_cube = variance_cube1 else var_cube = variance_cube2

      make_slices, file_struct, type='raw', file_ind = i, data_cube = data_cube, $
                  weights_cube = weights_cube, var_cube = var_cube, $
                  kx_mpc = kx_mpc, ky_mpc = ky_mpc, frequencies = frequencies, $
                  kperp_lambda_conv = kperp_lambda_conv, $
                  delay_params = delay_params, hubble_param = hubble_param

      undefine, data_cube, weights_cube, var_cube
    endfor
  endif

  if healpix or not keyword_set(uvf_input) then begin
    ;; multiply data, weights & variance cubes by pixel_area_rad to get proper units from DFT
    ;; (not squared for variance because they weren't treated as units squared in FHD code)
    data_cube1 = data_cube1 * pix_area_rad
    weights_cube1 = weights_cube1 * pix_area_rad
    variance_cube1 = variance_cube1 * pix_area_rad
    if nfiles eq 2 then begin
      data_cube2 = data_cube2 * pix_area_rad
      weights_cube2 = weights_cube2 * pix_area_rad
      variance_cube2 = variance_cube2 * pix_area_rad
    endif

    ave_weights = total(total(abs(weights_cube1),2),1)/(n_kx*n_ky)
    if nfiles eq 2 then begin
      ave_weights = transpose([[ave_weights], $
        [total(total(abs(weights_cube2),2),1)/(n_kx*n_ky)]])
    endif
  endif

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
    if min(sigma2_cube1) lt 0 or min(sigma2_cube2) lt 0 then begin
      message, 'sigma2 should be positive definite.'
    endif
    if total(abs(sigma2_cube1)) le 0 or total(abs(sigma2_cube2)) le 0 then begin
      message, 'one or both sigma2 cubes is all zero'
    endif
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

  if save_slices then begin
    ;; save some slices of the data cube
    for i=0, nfiles-1 do begin
      if i eq 0 then data_cube = data_cube1 else data_cube = data_cube2

      make_slices, file_struct, type='divided', file_ind = i, data_cube = data_cube, $
                  kx_mpc = kx_mpc, ky_mpc = ky_mpc, frequencies = frequencies, $
                  kperp_lambda_conv = kperp_lambda_conv, $
                  delay_params = delay_params, hubble_param = hubble_param

      undefine, data_cube
    endfor
  endif

  if healpix or not keyword_set(uvf_input) then begin
    ;; fix units on window funtion integral -- now they should be Mpc^3
    ;; checked vs Adam's analytic calculation and it matches to within a factor of 4
    ;; We are using total(weight) = Nvis for Ian's uv plane and multiplying by a
    ;;   factor to convert between Ian's and my uv plane
    ;;   the factor is ((2*!pi)^2*(delta_uv)^2) /  (delta_kperp)^2 * Dm^2 which
    ;;   comes in squared in the denominator.
    ;;   see eq 21e from Adam's memo. Also note that delta D is delta z * n_freq
    ;;   note that we can convert both the weights and variance to uvf from kperp,
    ;;   rz and all jacobians will cancel
    window_int_k = window_int * (z_mpc_delta * n_freq) * (kx_mpc_delta * ky_mpc_delta) * $
      z_mpc_mean^4./((2.*!pi)^2.*file_struct.kpix^4.)
    print, 'window integral from variances: ' + number_formatter(window_int_k[0], $
      format='(e10.4)')
    if tag_exist(file_struct, 'beam_savefile') then begin
      print, 'window integral from beam cube: ' + number_formatter(window_int_beam[0], $
        format='(e10.4)')
    endif
    if tag_exist(file_struct, 'beam_int') then begin
      print, 'window integral from obs.primary_beam_sq_area: ' + $
        number_formatter(window_int_beam_obs[0], format='(e10.4)')
    endif
    if tag_exist(file_struct, 'beam_int') then begin
      window_int = window_int_beam_obs
    endif

    if (n_elements(window_int) eq 0 or min(window_int) eq 0) then begin
      if ps_options.allow_beam_approx then begin
        print, 'WARNING: beam integral in obs structure is zero, using a less good approximation'
      endif else begin
        message, 'Beam integral in obs structure is zero. To use a less good ' + $
          'approximation instead, set keyword allow_beam_approx=1'
      endelse

      if tag_exist(file_struct, 'beam_savefile') then begin
        window_int = window_int_beam
        if min(window_int) eq 0 then print, 'WARNING: beam cube is zero, using a ' + $
          'less good approximation'
      endif
    endif
    if (n_elements(window_int) eq 0 or min(window_int) eq 0) then window_int = window_int_k
  endif else begin

    window_int_k = window_int * (z_mpc_delta * n_freq) * (2.*!pi)^2. / $
      (kx_mpc_delta * ky_mpc_delta)
    print, 'var_cube multiplier: ', (z_mpc_delta * n_freq) * (2.*!pi)^2. / $
      (kx_mpc_delta * ky_mpc_delta)
    print, 'window integral from variances: ' + number_formatter(window_int_k[0], $
      format='(e10.4)')
    if tag_exist(file_struct, 'beam_savefile') then begin
      print, 'window integral from beam cube: ' + number_formatter(window_int_beam[0], $
        format='(e10.4)')
    endif
    if tag_exist(file_struct, 'beam_int') then begin
      print, 'window integral from obs.primary_beam_sq_area: ' + $
        number_formatter(window_int_beam_obs[0], format='(e10.4)')
    endif

    if tag_exist(file_struct, 'beam_int') then begin
      window_int = window_int_beam_obs
      if min(window_int) eq 0 then begin
        print, 'WARNING: beam integral in obs structure is zero, using a less good approximation'
      endif
    endif
    if (n_elements(window_int) eq 0 or min(window_int) eq 0) and $
        tag_exist(file_struct, 'beam_savefile') then begin
      window_int = window_int_beam
      if min(window_int) eq 0 then begin
        print, 'WARNING: beam cube is zero, using a less good approximation'
      endif
    endif
    if (n_elements(window_int) eq 0 or min(window_int) eq 0) then window_int = window_int_k

  endelse


  ;; get sigma^2 into Jy^2
  sigma2_cube1 = sigma2_cube1 * rebin(reform(vis_sigma, 1, 1, n_freq), n_kx, n_ky, n_freq)^2.
  if nfiles eq 2 then begin
    sigma2_cube2 = sigma2_cube2 * rebin(reform(vis_sigma, 1, 1, n_freq), n_kx, n_ky, n_freq)^2.
  endif

  ;; get data & sigma into mK Mpc^2 and multiply by 2 (4 for variances) to get
  ;; to estimate of Stokes I rather than instrumental pol
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
  sim_noise1 = randomn(seed, n_kx, n_ky, n_kz) * sqrt(sigma2_cube1) + $
    complex(0,1) * randomn(seed, n_kx, n_ky, n_kz) * sqrt(sigma2_cube1)
  if nfiles eq 2 then begin
    sim_noise2 = randomn(seed, n_kx, n_ky, n_kz) * sqrt(sigma2_cube2) + $
      complex(0,1) * randomn(seed, n_kx, n_ky, n_kz) * sqrt(sigma2_cube2)
  endif

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
    wt_ave_power_freq[0,*] = total(total(sum_weights1 * abs(data_cube1)^2., 2), 1) / $
      total(total(sum_weights1, 2), 1) * (z_mpc_delta * n_freq)^2.
    wt_ave_power_freq[1,*] = total(total(sum_weights2 * abs(data_cube2)^2., 2), 1) / $
      total(total(sum_weights2, 2), 1) * (z_mpc_delta * n_freq)^2.
    ave_power_freq = fltarr(2, n_freq)
    for i=0, n_freq-1 do begin
      ave_power_freq[*, i] = [mean(abs((data_cube1[*,*,i])[where(sum_weights1[*,*,i] ne 0),*])^2.), $
        mean(abs((data_cube2[*,*,i])[where(sum_weights2[*,*,i] ne 0),*])^2.)] * (z_mpc_delta * n_freq)^2.
    endfor

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

  endif else begin
    sum_weights1 = 1./sigma2_cube1
    wh_sig1_0 = where(sigma2_cube1 eq 0, count_sig1_0, complement = wh_sig1_n0)
    if count_sig1_0 ne 0 then sum_weights1[wh_sig1_0] = 0

    wt_ave_power_freq = total(total(sum_weights1 * abs(data_cube1)^2., 2), 1)/$
      total(total(sum_weights1, 2), 1)
    ave_power_freq = fltarr(n_freq)
    for i=0, n_freq-1 do begin
      ave_power_freq[i] = mean(abs((data_cube1[*,*,i])[where(sum_weights1[*,*,i] ne 0),*])^2.)
    endfor

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
  n_freq_contrib = total(mask, 3)
  wh_nofreq = where(n_freq_contrib eq 0, count_nofreq)
  undefine, mask

  if save_slices then begin
    ;; save some slices of the sum & diff cubes

    make_slices, file_struct, type='sum', data_cube = data_sum, $
                kx_mpc = kx_mpc, ky_mpc = ky_mpc, frequencies = frequencies, $
                kperp_lambda_conv = kperp_lambda_conv, $
                delay_params = delay_params, hubble_param = hubble_param

    if nfiles eq 2 then begin
      make_slices, file_struct, type='diff', data_cube = data_diff, $
                  kx_mpc = kx_mpc, ky_mpc = ky_mpc, frequencies = frequencies, $
                  kperp_lambda_conv = kperp_lambda_conv, $
                  delay_params = delay_params, hubble_param = hubble_param
    endif
  endif

  if ps_options.ave_removal then begin

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
  if tag_exist(ps_options, 'spec_window_type') then begin
    window = spectral_window(n_freq, type = ps_options.spec_window_type, /periodic)

    ;; this is the equivilent bandwidth reduction (approx. 2 for blackman-harris)
    norm_factor = sqrt(n_freq/total(window^2.))

    window = window * norm_factor

    if ps_options.ave_removal then begin
      ;; need to divide by bandwidth factor to account for the
      ;; reduction in the mean that should happen when the window is applied
      data_sum_mean = data_sum_mean / norm_factor
      sim_noise_sum_mean = sim_noise_sum_mean / norm_factor
      data_diff_mean = data_diff_mean / norm_factor
      sim_noise_diff_mean = sim_noise_diff_mean / norm_factor
    endif

    window_expand = rebin(reform(window, 1, 1, n_freq), n_kx, n_ky, n_freq, /sample)

    data_sum = data_sum * window_expand
    sim_noise_sum = sim_noise_sum * window_expand
    if nfiles eq 2 then begin
      data_diff = data_diff * window_expand
      sim_noise_diff = sim_noise_diff * window_expand
    endif

    sum_sigma2 = sum_sigma2 * temporary(window_expand^2.)
  endif

  if ps_options.inverse_covar_weight then begin

    max_eta_val = max(sum_sigma2)
    eta_var = shift([max_eta_val, fltarr(n_freq-1)], fft_shift_val)
    covar_eta_fg = diag_matrix(eta_var)

    identity = diag_matrix([fltarr(n_freq)+1d])
    ft_matrix = fft(identity, dimension=1)*n_freq*z_mpc_delta
    inv_ft_matrix = fft(identity, dimension=1, /inverse)*kz_mpc_delta
    undefine, identity

    covar_z_fg = matrix_multiply(inv_ft_matrix, $
      matrix_multiply(shift(covar_eta_fg, fft_shift_val, fft_shift_val), conj(inv_ft_matrix), /btranspose))

    wt_data_sum = fltarr(n_kx, n_ky, n_freq)
    wt_data_diff = fltarr(n_kx, n_ky, n_freq)
    wt_sum_sigma2 = fltarr(n_kx, n_ky, n_freq)
    wt_power_norm = fltarr(n_kx, n_ky, n_freq)
    for i=0, n_kx-1 do begin
      for j=0, n_ky-1 do begin
        var_z_inst = reform(sum_sigma2[i,j,*])
        if max(abs(var_z_inst)) eq 0 then continue

        ;; need to convert zeros which represet lack of measurements to infinities
        ;; (or large numbers)
        wh_var0 = where(var_z_inst eq 0, count_var0)
        if count_var0 gt 0 then var_z_inst[wh_var0] = 1e6

        covar_z = diag_matrix(var_z_inst) + covar_z_fg
        inv_covar_z = la_invert(covar_z)
        inv_covar_kz = shift(matrix_multiply(ft_matrix, $
          matrix_multiply(inv_covar_z, conj(ft_matrix), /btranspose)), fft_shift_val, fft_shift_val)

        inv_var = 1/var_z_inst
        wh_sigma0 = where(var_z_inst eq 0, count_sigma0)
        if count_sigma0 gt 0 then inv_var[wh_sigma0] = 0
        norm2 = total(inv_var)^2./(n_freq*total(inv_var^2.))

        wt_data_sum[i,j,*] = matrix_multiply(inv_covar_z, reform(data_sum[i,j,*])) * z_mpc_delta
        wt_data_diff[i,j,*] = matrix_multiply(inv_covar_z, reform(data_diff[i,j,*])) * z_mpc_delta
        wt_power_norm[i,j,*] = norm2 * total(abs(inv_covar_kz)^2.,2) * kz_mpc_delta

        wt_sum_sigma2[i,j,*] = diag_matrix(shift(matrix_multiply(inv_covar_z, $
          matrix_multiply(diag_matrix(var_z_inst), conj(inv_covar_z), /btranspose)), $
          fft_shift_val, fft_shift_val))

      endfor
    endfor
    undefine, covar_z_fg, ft_matrix, inv_ft_matrix, covar_z, inv_covar_z, inv_covar_kz

    data_sum = temporary(wt_data_sum)
    data_diff = temporary(wt_data_diff)
    sum_sigma2 = temporary(wt_sum_sigma2)

  endif

  n_val = kz_integers
  wh_pos = where(n_val gt 0, n_pos)
  wh_neg = where(n_val lt 0, n_neg)

  ;; this should be zero by construction, but add a test just in case
  if kz_mpc_orig[where(n_val eq 0)] ne 0 then begin
    message, 'something has gone terribly wrong with calculating the kz_mpc_orig values.'
  endif

  kz_mpc = kz_mpc_orig[where(n_val ge 0)]

  if is_n_freq_even then begin
    ;; if n_freq is even, there's an extra positive mode (because there's always a zero mode.)
    if max(abs(n_val(wh_pos))) le max(abs(n_val(wh_neg))) then begin
      message, 'something has gone terribly wrong with calculating the kz_integers.'
    endif
    ;; remove unmatched max k positive mode:
    trim_max_pos = 1
    n_val = n_val[0:-2]
    kz_mpc = kz_mpc[0:-2]
    kz_mpc_orig_trim = kz_mpc_orig[0:-2]
    wh_pos = wh_pos[0:-2]
    n_pos = n_elements(wh_pos)
  endif else begin
    kz_mpc_orig_trim = kz_mpc_orig
    trim_max_pos = 0
  endelse

  ;; now take frequency FT
  if not ps_options.freq_dft then begin
    ;; evenly spaced and dft not requested, just use fft
    print, "Using FFT for evenly spaced frequencies"
    data_sum_ft = fft(data_sum, dimension=3) * n_freq * z_mpc_delta
    sim_noise_sum_ft = fft(sim_noise_sum, dimension=3) * n_freq * z_mpc_delta

    ;; put k0 in middle of cube
    data_sum_ft = shift(data_sum_ft, [0,0,fft_shift_val])
    sim_noise_sum_ft = shift(sim_noise_sum_ft, [0,0,fft_shift_val])

    if trim_max_pos then begin
      ;; remove unmatched max k positive mode:
      data_sum_ft = data_sum_ft[*, *, 0:-2]
      sim_noise_sum_ft = sim_noise_sum_ft[*, *, 0:-2]
    endif

    undefine, data_sum, sim_noise_sum
    if nfiles eq 2 then begin
      data_diff_ft = fft(data_diff, dimension=3) * n_freq * z_mpc_delta
      ;; put k0 in middle of cube
      data_diff_ft = shift(data_diff_ft, [0,0,fft_shift_val])
      if trim_max_pos then begin
        ;; remove unmatched max k positive mode:
        data_diff_ft = data_diff_ft[*, *, 0:-2]
      endif
      undefine, data_diff

      sim_noise_diff_ft = fft(sim_noise_diff, dimension=3) * n_freq * z_mpc_delta
      ;; put k0 in middle of cube
      sim_noise_diff_ft = shift(sim_noise_diff_ft, [0,0,fft_shift_val])
      if trim_max_pos then begin
        ;; remove unmatched max k positive mode:
        sim_noise_diff_ft = sim_noise_diff_ft[*, *, 0:-2]
      endif
      undefine, sim_noise_diff
    endif
  endif else begin
    ;; Not evenly spaced or DFT requested.
    ;; use the kz_mpc_orig_trim to avoid the possible extra positive mode.
    print, "Uneven frequency structure found"
    print, "Performing Discrete Fourier Transform"

    ;; There are 2 options for the DFT. One will get the (very nearly) the same answer
    ;; as the FFT for regular gridding because it uses fixed conversion between delta_z and delta_f
    ;; in some simple testing, this had a difference ratio vs the FFT (abs(diff)/abs(fft)) of ~10^-3
    z_reg = (frequencies-frequencies[0])*z_mpc_delta/f_delta
    z_exp_zreg = exp(-1.*complex(0,1)*matrix_multiply(z_reg, kz_mpc_orig_trim, /btranspose))

    ;; This one uses the true comov_dist_los values, so it is a little more different compared to the fft
    ;; (which run large to small)
    ;; in some simple testing, this had a difference ratio vs the FFT (abs(diff)/abs(fft)) of ~10^1
    z_exp_ztrue =  exp(-1.*complex(0,1)*matrix_multiply(reverse(comov_dist_los-min(comov_dist_los)), kz_mpc_orig_trim, /btranspose))

    case ps_options.dft_z_use of
      'true': z_exp = z_exp_ztrue
      'regular': z_exp = z_exp_zreg
    endcase

    n_kz_trim = n_elements(kz_mpc_orig_trim)

    data_sum_ft = z_mpc_delta * $
      reform(matrix_multiply(reform(temporary(data_sum), n_kx*n_ky, n_freq), z_exp), $
      n_kx, n_ky, n_kz_trim)
    sim_noise_sum_ft = z_mpc_delta * $
      reform(matrix_multiply(reform(temporary(sim_noise_sum), n_kx*n_ky, n_freq), z_exp), $
      n_kx, n_ky, n_kz_trim)
    if nfiles eq 2 then begin
      data_diff_ft = z_mpc_delta * $
        reform(matrix_multiply(reform(temporary(data_diff), n_kx*n_ky, n_freq), z_exp), $
        n_kx, n_ky, n_kz_trim)
      sim_noise_diff_ft  = z_mpc_delta * $
        reform(matrix_multiply(reform(temporary(sim_noise_diff), n_kx*n_ky, n_freq), z_exp), $
        n_kx, n_ky, n_kz_trim)
    endif
  endelse

  n_kz = n_elements(kz_mpc)

  ;; these an and bn calculations don't match the standard
  ;; convention (they ares down by a factor of 2) but they make more sense
  ;; and remove factors of 2 we'd otherwise have in the power
  ;; and variance calculations
  ;; note that the 0th mode will have higher noise because there's half as many
  ;; measurements going into it
  a1_0 = data_sum_ft[*,*,where(n_val eq 0)]
  a1_n = (data_sum_ft[*, *, wh_pos] + data_sum_ft[*, *, reverse(wh_neg)])/2.
  b1_n = complex(0,1) * (data_sum_ft[*, *, wh_pos] - $
        data_sum_ft[*, *, reverse(wh_neg)])/2.
  undefine, data_sum_ft

  a3_0 = sim_noise_sum_ft[*,*,where(n_val eq 0)]
  a3_n = (sim_noise_sum_ft[*, *, wh_pos] + sim_noise_sum_ft[*, *, reverse(wh_neg)])/2.
  b3_n = complex(0,1) * (sim_noise_sum_ft[*, *, wh_pos] - $
        sim_noise_sum_ft[*, *, reverse(wh_neg)])/2.
  undefine, sim_noise_sum_ft

  if nfiles gt 1 then begin
    a2_0 = data_diff_ft[*,*,where(n_val eq 0)]
    a2_n = (data_diff_ft[*, *, wh_pos] + data_diff_ft[*, *, reverse(wh_neg)])/2.
    b2_n = complex(0,1) * (data_diff_ft[*, *, wh_pos] - $
          data_diff_ft[*, *, reverse(wh_neg)])/2.
    undefine, data_diff_ft

    a4_0 = sim_noise_diff_ft[*,*,where(n_val eq 0)]
    a4_n = (sim_noise_diff_ft[*, *, wh_pos] + sim_noise_diff_ft[*, *, reverse(wh_neg)])/2.
    b4_n = complex(0,1) * (sim_noise_diff_ft[*, *, wh_pos] - $
          sim_noise_diff_ft[*, *, reverse(wh_neg)])/2.
    undefine, sim_noise_diff_ft
  endif

  ;; drop pixels with less than 1/3 of the frequencies (set weights to 0)
  wh_fewfreq = where(n_freq_contrib lt ceil(n_freq/3d), count_fewfreq)
  if count_fewfreq gt 0 then begin
    mask_fewfreq = n_freq_contrib * 0 + 1
    mask_fewfreq[wh_fewfreq] = 0
    mask_fewfreq = rebin(temporary(mask_fewfreq), n_kx, n_ky, n_kz)

    a1_0 = temporary(a1_0) * mask_fewfreq[*,*,0]
    a1_n = temporary(a1_n) * mask_fewfreq[*,*,1:*]
    b1_n = temporary(b1_n) * mask_fewfreq[*,*,1:*]

    a3_0 = temporary(a3_0) * mask_fewfreq[*,*,0]
    a3_n = temporary(a3_n) * mask_fewfreq[*,*,1:*]
    b3_n = temporary(b3_n) * mask_fewfreq[*,*,1:*]
    if nfiles gt 1 then begin
      a2_0 = temporary(a2_0) * mask_fewfreq[*,*,0]
      a2_n = temporary(a2_n) * mask_fewfreq[*,*,1:*]
      b2_n = temporary(b2_n) * mask_fewfreq[*,*,1:*]

      a4_0 = temporary(a4_0) * mask_fewfreq[*,*,0]
      a4_n = temporary(a4_n) * mask_fewfreq[*,*,1:*]
      b4_n = temporary(b4_n) * mask_fewfreq[*,*,1:*]
    endif
  endif

  data_sum_cos = dcomplex(fltarr(n_kx, n_ky, n_kz))
  data_sum_sin = dcomplex(fltarr(n_kx, n_ky, n_kz))
  data_sum_cos[*, *, 0] = a1_0
  data_sum_cos[*, *, 1:n_kz-1] = a1_n
  data_sum_sin[*, *, 1:n_kz-1] = b1_n

  sim_noise_sum_cos = dcomplex(fltarr(n_kx, n_ky, n_kz))
  sim_noise_sum_sin = dcomplex(fltarr(n_kx, n_ky, n_kz))
  sim_noise_sum_cos[*, *, 0] = a3_0
  sim_noise_sum_cos[*, *, 1:n_kz-1] = a3_n
  sim_noise_sum_sin[*, *, 1:n_kz-1] = b3_n

  if nfiles gt 1 then begin
    data_diff_cos = dcomplex(fltarr(n_kx, n_ky, n_kz))
    data_diff_sin = dcomplex(fltarr(n_kx, n_ky, n_kz))
    data_diff_cos[*, *, 0] = a2_0
    data_diff_cos[*, *, 1:n_kz-1] = a2_n
    data_diff_sin[*, *, 1:n_kz-1] = b2_n

    sim_noise_diff_cos = dcomplex(fltarr(n_kx, n_ky, n_kz))
    sim_noise_diff_sin = dcomplex(fltarr(n_kx, n_ky, n_kz))
    sim_noise_diff_cos[*, *, 0] = a4_0
    sim_noise_diff_cos[*, *, 1:n_kz-1] = a4_n
    sim_noise_diff_sin[*, *, 1:n_kz-1] = b4_n
  endif

  if ps_options.ave_removal then begin

    data_sum_cos[*,*,0] = data_sum_cos[*,*,0] + data_sum_mean * n_freq * z_mpc_delta
    sim_noise_sum_cos[*,*,0] = sim_noise_sum_cos[*,*,0] + sim_noise_sum_mean * $
      n_freq * z_mpc_delta

    if nfiles eq 2 then begin
      data_diff_cos[*,*,0] = data_diff_cos[*,*,0] + data_diff_mean * n_freq * z_mpc_delta
      sim_noise_diff_cos[*,*,0] = sim_noise_diff_cos[*,*,0] + sim_noise_diff_mean * $
        n_freq * z_mpc_delta
    endif
  endif

  ;; for new power calc, need cos2, sin2, cos*sin transforms
  covar_cos = fltarr(n_kx, n_ky, n_kz)
  covar_sin = fltarr(n_kx, n_ky, n_kz)
  covar_cross = fltarr(n_kx, n_ky, n_kz)

  if not ps_options.freq_dft then begin
    ;; comov_dist_los goes from large to small z
    z_relative = dindgen(n_freq)*z_mpc_delta
    freq_kz_arr = rebin(reform(kz_mpc, 1, n_kz), n_freq, n_kz) * $
      rebin(z_relative, n_freq, n_kz)
  endif else begin
    ;; unevenly spaced frequencies. Need to match what is done above
    freq_kz_arr_reg = rebin(reform(kz_mpc, 1, n_kz), n_freq, n_kz) * $
      rebin(z_reg, n_freq, n_kz)

    freq_kz_arr_true = rebin(reform(kz_mpc, 1, n_kz), n_freq, n_kz) * $
      rebin(reverse(comov_dist_los-min(comov_dist_los)), n_freq, n_kz)

    case ps_options.dft_z_use of
      'true': freq_kz_arr = freq_kz_arr_true
      'regular': freq_kz_arr = freq_kz_arr_reg
    endcase
  endelse

  cos_arr = cos(freq_kz_arr)
  sin_arr = sin(freq_kz_arr)

  sum_sigma2 = reform(sum_sigma2, n_kx*n_ky, n_freq)
  ;; doing 2 FTs so need 2 factors of z_mpc_delta.
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

  ;; drop pixels with less than 1/3 of the frequencies
  if count_fewfreq gt 0 then begin
    covar_cos = temporary(covar_cos) * mask_fewfreq
    covar_sin = temporary(covar_sin) * mask_fewfreq
    covar_cross = temporary(covar_cross) * mask_fewfreq
    undefine, mask_fewfreq
  endif

  undefine, sum_sigma2, freq_kz_arr, cos_arr, sin_arr

  if ps_options.std_power then begin
    ;; standard power means no rotation between sine & cosine modes
    print, 'WARNING: using the standard power turns off the Lomb-Scargle ' $
      + 'rotation, which can result in power gain or loss.'
    theta = fltarr(n_elements(covar_cos))
  endif else begin
    ;; get rotation angle to diagonalize covariance block
    ;; There is a numerical issue where the denominator is near zero.
    ;; If the numerator is near zero, theta should be zero (no matter what the denom is)
    ;; If the numerator is not near zero and the denominator is near zero, theta should be pi/4
    ;; We test these cases explicitly and set theta to the correct values
    zeroish_thresh = 1e-8
    wh_num_zeroish = where(abs(covar_cross) lt zeroish_thresh, count_num_zeroish)
    denominator = covar_cos - covar_sin
    wh_denom_zeroish = where(abs(denominator) lt zeroish_thresh, count_denom_zeroish)

    ;; use atan with a single argument because it should be between -pi/2 and pi/2 not -pi and pi
    theta = atan(2.*covar_cross/denominator)/2.

    if count_denom_zeroish gt 0 then theta[wh_denom_zeroish] = !dpi / 2.
    if count_num_zeroish gt 0 then theta[wh_num_zeroish] = 0.
  endelse

  cos_theta = cos(theta)
  sin_theta = sin(theta)
  undefine, theta

  sigma2_1 = covar_cos*cos_theta^2. + 2.*covar_cross*cos_theta*sin_theta + $
    covar_sin*sin_theta^2.
  sigma2_2 = covar_cos*sin_theta^2. - 2.*covar_cross*cos_theta*sin_theta + $
    covar_sin*cos_theta^2.

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

  git, repo_path = ps_repository_dir(), result=kcube_git_hash
  git_hashes = {uvf:uvf_git_hashes, uvf_wt:uvf_wt_git_hashes, beam:beam_git_hashes, $
    kcube:kcube_git_hash}

  if n_elements(freq_flags) gt 0 then begin
    save, file = file_struct.kcube_savefile, data_sum_1, data_sum_2, data_diff_1, $
      data_diff_2, sigma2_1, sigma2_2, sim_noise_sum_1, sim_noise_sum_2, $
      sim_noise_diff_1, sim_noise_diff_2, kx_mpc, ky_mpc, kz_mpc, kperp_lambda_conv, $
      delay_params, hubble_param, n_freq_contrib, freq_mask, vs_name, vs_mean, $
      t_sys_meas, window_int, git_hashes, wt_meas_ave, wt_meas_min, ave_weights, $
      wt_ave_power_freq, ave_power_freq, wt_ave_power_uvf, ave_power_uvf
  endif else begin
    save, file = file_struct.kcube_savefile, data_sum_1, data_sum_2, data_diff_1, $
      data_diff_2, sigma2_1, sigma2_2, sim_noise_sum_1, sim_noise_sum_2, $
      sim_noise_diff_1, sim_noise_diff_2, kx_mpc, ky_mpc, kz_mpc, kperp_lambda_conv, $
      delay_params, hubble_param, n_freq_contrib, vs_name, vs_mean, t_sys_meas, $
      window_int, git_hashes, wt_meas_ave, wt_meas_min, ave_weights, $
      wt_ave_power_freq, ave_power_freq, wt_ave_power_uvf, ave_power_uvf
  endelse


end
