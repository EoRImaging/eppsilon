pro fhd_kcube, file_struct, dft_refresh_data = dft_refresh_data, dft_refresh_weight = dft_refresh_weight, dft_ian = dft_ian, $
    dft_fchunk = dft_fchunk, freq_ch_range = freq_ch_range, spec_window_type = spec_window_type, cut_image = cut_image, $
    noise_sim = noise_sim, std_power = std_power, input_units = input_units, image = image, quiet = quiet
    
  if tag_exist(file_struct, 'nside') ne 0 then healpix = 1 else healpix = 0
  if keyword_set(image) or tag_exist(file_struct, 'uvf_savefile') ne 0 then image = 1 else image = 0
  nfiles = n_elements(file_struct.datafile)
  if tag_exist(file_struct, 'no_var') ne 0 then no_var = 1 else no_var = 0

  
  if n_elements(input_units) eq 0 then input_units = 'jansky'
  units_enum = ['jansky', 'mk']
  wh = where(units_enum eq input_units, count)
  if count eq 0 then message, 'input units not recognized, options are: ' + units_enum
  
  
  datavar = strupcase(file_struct.datavar)
  if datavar eq '' then begin
    ;; working with a 'derived' cube (ie residual cube or noise simulation cube) that is constructed from uvf_savefiles
    if keyword_set(noise_sim) then begin
      input_uvf_files = reform(file_struct.res_uvf_inputfiles)
      input_uvf_varname = reform(file_struct.res_uvf_varname)
    endif else begin
      input_uvf_files = reform(file_struct.res_uvf_inputfiles, nfiles, 2)
      input_uvf_varname = reform(file_struct.res_uvf_varname, nfiles, 2)
    endelse
    
    if healpix or keyword_set(image) then begin
      input_uvf_wtfiles = file_struct.uvf_weight_savefile
    endif
  endif
  
  frequencies = file_struct.frequencies
  
  if n_elements(freq_ch_range) ne 0 then begin
    n_freq_orig = n_elements(frequencies)
    frequencies = frequencies[min(freq_ch_range):max(freq_ch_range)]
  endif
  n_freq = n_elements(frequencies)
  
  ;; check whether or not the frequencies are evenly spaced.
  freq_diff = frequencies - shift(frequencies, 1)
  freq_diff = freq_diff[1:*]
  
  z0_freq = 1420.40d ;; MHz
  redshifts = z0_freq/frequencies - 1d
  cosmology_measures, redshifts, comoving_dist_los = comov_dist_los, hubble_param = hubble_param
  
  comov_los_diff = comov_dist_los - shift(comov_dist_los, -1)
  comov_los_diff = comov_los_diff[0:n_elements(comov_dist_los)-2]
  
  if max(freq_diff-freq_diff[0]) gt 1e-12 then begin
    ;; frequencies are not evenly spaced, need to be careful about z_mpc_delta/mean
    even_freq = 0
    
    freq_diff_hist = histogram(freq_diff, binsize = min(freq_diff)*.1, locations=locs, reverse_indices = ri)
    if max(freq_diff_hist)/float(n_freq) lt .9 then stop else begin
      peak_bin = (where(freq_diff_hist eq max(freq_diff_hist), count_peak))[0]
      if count_peak eq 1 then peak_diffs = freq_diff[ri[ri[peak_bin] : ri[peak_bin+1]-1]]
      
      f_delta = mean(peak_diffs)
    endelse
    
    nominal_freqs = findgen(floor(((max(frequencies)-min(frequencies))/f_delta))+1)*f_delta + min(frequencies)
    nominal_z = z0_freq/nominal_freqs - 1
    cosmology_measures, nominal_z, comoving_dist_los = nominal_comov_dist_los
    nominal_comov_diffs = nominal_comov_dist_los - shift(nominal_comov_dist_los, -1)
    nominal_comov_diffs = nominal_comov_diffs[0:n_elements(nominal_comov_diffs)-2]
    
    z_mpc_delta = mean(nominal_comov_diffs)
    z_mpc_mean = mean(nominal_comov_dist_los)
    
  endif else begin
    even_freq = 1
    
    f_delta = freq_diff[0]
    z_mpc_delta = float(mean(comov_los_diff))
    z_mpc_mean = float(mean(comov_dist_los))
    n_kz = n_freq
    
  endelse
  kperp_lambda_conv = z_mpc_mean / (2.*!pi)
  delay_delta = 1e9/(n_freq*f_delta*1e6) ;; equivilent delay bin size for kparallel
  delay_max = delay_delta * n_freq/2.    ;; factor of 2 b/c of neg/positive
  delay_params = [delay_delta, delay_max]
  
  z_mpc_length = max(comov_dist_los) - min(comov_dist_los) + z_mpc_delta
  kz_mpc_range =  (2.*!pi) / (z_mpc_delta)
  kz_mpc_delta = (2.*!pi) / z_mpc_length
  kz_mpc_orig = findgen(round(kz_mpc_range / kz_mpc_delta)) * kz_mpc_delta - kz_mpc_range/2.
  if n_elements(n_kz) ne 0 then begin
    if n_elements(kz_mpc_orig) ne n_kz then stop
  endif else n_kz = n_elements(kz_mpc_orig)
  
  if input_units eq 'jansky' then begin
    ;; beam_diameter_rad = (3d * 10^8d) / (frequencies * 10^6d * max_baseline)
    ;; beam_area_str = !pi * beam_diameter_rad^2d /4d
  
    ;; conv_factor = (10^(double(-26+16+3-12+23)) * 9d) / (beam_area_str * 2d * frequencies^2d * 1.38)
    ;; if max(conv_factor-conv_factor[0]) gt 1e-8 then stop else conv_factor = conv_factor[0]
    ;; conv_factor = float(2. * max_baseline^2. / (!pi * 1.38065))
  
    ;; converting from Jy (in u,v,f) to mK*str (10^-26 * c^2 * 10^-3/ (2*f^2*kb))
    conv_factor = float((3e8)^2 / (2. * (frequencies*1e6)^2. * 1.38065))
    
    ;; convert from mk*str to mK*Mpc^2
    conv_factor = conv_factor * z_mpc_mean^2.
    
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
  
  old_vis_sigma = vis_sigma
  
  if n_elements(freq_ch_range) ne 0 then vis_sig_tag = number_formatter(384./n_freq_orig) else vis_sig_tag = number_formatter(384./n_freq)
  vis_sigma_file = file_dirname(file_struct.savefile_froot, /mark_directory) + 'vis_sigma/vis_sigma_measured' + vis_sig_tag + '.sav'
  
  if file_test(vis_sigma_file) then begin
    restore, vis_sigma_file
    
    if n_elements(freq_ch_range) ne 0 then begin
      if n_elements(vis_sigma) ne n_freq_orig then stop
      vis_sigma = vis_sigma[min(freq_ch_range):max(freq_ch_range)]
    endif else if n_elements(vis_sigma) ne n_freq then stop
    
    wh_nan = where(finite(vis_sigma) eq 0, count_nan)
    if count_nan gt 0 then vis_sigma[wh_nan] = 0
  endif else begin
    ;; make a flat vis_sigma
    vis_sigma = old_vis_sigma*0 + old_vis_sigma[0]
  endelse
  
  n_vis = file_struct.n_vis
  if n_elements(freq_ch_range) ne 0 then n_vis = n_vis * (float(n_freq)/float(n_freq_orig))
  
  if healpix or image then begin
  
    for i=0, nfiles-1 do begin
      test_uvf = file_test(file_struct.uvf_savefile[i]) *  (1 - file_test(file_struct.uvf_savefile[i], /zero_length))
      
      test_wt_uvf = file_test(file_struct.uvf_weight_savefile[i]) * (1 - file_test(file_struct.uvf_weight_savefile[i], /zero_length))
      
      if test_uvf eq 0 and not keyword_set(dft_refresh_data) and n_elements(freq_ch_range) ne 0 then begin
        ;; if this is a limited freq. range cube, check for the full cube to avoid redoing the DFT
        full_uvf_file = strjoin(strsplit(file_struct.uvf_savefile[i], '_ch[0-9]+-[0-9]+', /regex, /extract))
        test_full_uvf = file_test(full_uvf_file) *  (1 - file_test(full_uvf_file, /zero_length))
        if test_full_uvf eq 1 then begin
          restore, full_uvf_file
          
          data_cube = data_cube[*, *, min(freq_ch_range):max(freq_ch_range)]
          
          if keyword_set(dft_ian) then save, file = file_struct.uvf_savefile[i], u_lambda_vals, v_lambda_vals, data_cube $
          else save, file = file_struct.uvf_savefile[i], kx_rad_vals, ky_rad_vals, data_cube
          undefine, data_cube
          
          test_uvf = 1
        endif
      endif
      
      if test_wt_uvf eq 0 and not keyword_set(dft_refresh_weight) and n_elements(freq_ch_range) ne 0 then begin
        ;; if this is a limited freq. range cube, check for the full cube to avoid redoing the DFT
        full_uvf_wt_file = strjoin(strsplit(file_struct.uvf_weight_savefile[i], '_ch[0-9]+-[0-9]+', /regex, /extract))
        test_full_wt_uvf = file_test(full_uvf_wt_file) *  (1 - file_test(full_uvf_wt_file, /zero_length))
        if test_full_wt_uvf eq 1 then begin
          restore, full_uvf_wt_file
          
          weights_cube = weights_cube[*, *, min(freq_ch_range):max(freq_ch_range)]
          variance_cube = variance_cube[*, *, min(freq_ch_range):max(freq_ch_range)]
          
          if keyword_set(dft_ian) then save, file = file_struct.uvf_weight_savefile[i], u_lambda_vals, v_lambda_vals, weights_cube, variance_cube else $
            save, file = file_struct.uvf_weight_savefile[i], kx_rad_vals, ky_rad_vals, weights_cube, variance_cube
          undefine, weights_cube, variance_cube
          
          test_wt_uvf = 1
        endif
      endif
      
      if test_uvf eq 0 or test_wt_uvf eq 0 or keyword_set(dft_refresh_data) or keyword_set(dft_refresh_weight) then begin
        if datavar eq '' then begin
          ;; working with a 'derived' cube (ie residual cube or noise simulation cube) that is constructed from uvf_savefiles
        
          if keyword_set(noise_sim) then begin
            variance_cube = getvar_savefile(file_struct.uvf_weight_savefile[i], 'variance_cube')
            kx_rad_vals = getvar_savefile(file_struct.uvf_weight_savefile[i], 'kx_rad_vals')
            ky_rad_vals = getvar_savefile(file_struct.uvf_weight_savefile[i], 'ky_rad_vals')
            
            seed = systime(1)
            noise = randomn(seed, dims) * sqrt(variance_cube) + $
              complex(0,1)*randomn(seed, dims) * sqrt(variance_cube)
            data_cube = temporary(noise)
            if keyword_set(dft_ian) then save, file = file_struct.uvf_savefile[i], u_lambda_vals, v_lambda_vals, data_cube $
            else save, file = file_struct.uvf_savefile[i], kx_rad_vals, ky_rad_vals, data_cube
            undefine, data_cube, seed, variance_cube
          endif else begin
            dirty_cube = getvar_savefile(input_uvf_files[i,0], input_uvf_varname[i,0])
            kx_dirty = getvar_savefile(input_uvf_files[i,0], 'kx_rad_vals')
            ky_dirty = getvar_savefile(input_uvf_files[i,0], 'ky_rad_vals')
            
            model_cube = getvar_savefile(input_uvf_files[i,1], input_uvf_varname[i,1])
            kx_rad_vals = getvar_savefile(input_uvf_files[i,1], 'kx_rad_vals')
            ky_rad_vals = getvar_savefile(input_uvf_files[i,1], 'ky_rad_vals')
            
            if total(abs(kx_rad_vals - kx_dirty)) ne 0 then message, 'kx_rad_vals for dirty and model cubes must match'
            if total(abs(ky_rad_vals - ky_dirty)) ne 0 then message, 'kx_rad_vals for dirty and model cubes must match'
            undefine, kx_dirty, ky_dirty
          endelse
          
          data_cube = temporary(dirty_cube) - temporary(model_cube)
          if keyword_set(dft_ian) then save, file = file_struct.uvf_savefile[i], u_lambda_vals, v_lambda_vals, data_cube $
          else save, file = file_struct.uvf_savefile[i], kx_rad_vals, ky_rad_vals, data_cube
          undefine, data_cube
        endif else begin
        
          if healpix then begin
            pixel_nums1 = getvar_savefile(file_struct.pixelfile[0], file_struct.pixelvar[0])
            
            pixel_dims = size(pixel_nums1, /dimension)
            if datavar ne '' and total(abs(dims - pixel_dims)) ne 0 then message, 'pixel and data dimensions do not match'
            
            if nfiles eq 2 then begin
              ;; check that they have the same set of healpix pixels
              pixel_nums2 = getvar_savefile(file_struct.pixelfile[1], file_struct.pixelvar[1])
              if n_elements(pixel_nums1) ne n_elements(pixel_nums2) then message, 'Different number of Healpix pixels in cubes'
              
              if total(abs(pixel_nums1-pixel_nums2)) ne 0 then message, 'Pixel numbers are not consistent between cubes'
            endif
            
            ;; get pixel vectors
            pix2vec_ring, file_struct.nside, pixel_nums1, pix_center_vec
            ;; find mid point (work in x/y because of possible jumps in phi)
            vec_mid = [mean(pix_center_vec[*,0]), mean(pix_center_vec[*,1]), mean(pix_center_vec[*,2])]
            theta0 = acos(vec_mid[2])
            phi0 = atan(vec_mid[1], vec_mid[0])
            
            ;; To go to flat sky, rotate patch to zenith and flatten.
            ;; To get to current location, need to first rotate around z by
            ;; phi, then around y by -theta, then around z by -phi
            ;; use inverse to rotate back to zenith
            rot_matrix = get_rot_matrix(theta0, phi0, /inverse)
            new_pix_vec = rot_matrix ## pix_center_vec
            
          endif else begin
            ;; gridded image to dft to parallel Healpix computation
            pix_size_rad = abs(file_struct.degpix) * !pi / 180d
            x_vec = (findgen(dims[0]) - dims[0]/2.) * pix_size_rad
            y_vec = (findgen(dims[1]) - dims[1]/2.) * pix_size_rad
            
            new_pix_vec = fltarr(dims[0]*dims[1], 3)
            new_pix_vec[*,0] = reform(rebin(x_vec, dims[0], dims[1], /sample), dims[0]*dims[1])
            new_pix_vec[*,1] = reform(rebin(reform(y_vec, 1, dims[1]), dims[0], dims[1], /sample), dims[0]*dims[1])
            new_pix_vec[*,2] = 1.
          endelse
          
          ;; figure out k values to calculate dft
          uv_cellsize_m = 5 ;; based on calculations of beam FWHM by Aaron
          if keyword_set(dft_ian) then begin
            ;;delta_u_lambda = file_struct.kpix
            delta_u_lambda = uv_cellsize_m * mean(frequencies*1e6) / (3e8)
            
            ;; go a little beyond max_baseline to account for expansion due to w projection
            ;;max_u_lambda = (file_struct.max_baseline_lambda) * 1.1
            ;; use kspan of Ian's cubes
            max_kperp_rad = min([file_struct.kspan/2.,file_struct.max_baseline_lambda])
            
            if keyword_set(cut_image) then begin
              ;; limit field of view to match calculated k-modes
              xy_len = 1/delta_u_lambda
              wh_close = where(new_pix_vec[*,0] le xy_len/2 and new_pix_vec[*,0] ge -1*xy_len/2 and $
                new_pix_vec[*,1] le xy_len/2 and new_pix_vec[*,1] ge -1*xy_len/2, count_close, $
                ncomplement = count_far)
              if count_far ne 0 then new_pix_vec = new_pix_vec[wh_close, *] else begin
                ;; image may be smaller than expected, may need to adjust delta_kperp_rad
                image_len = max(sqrt(new_pix_vec[*,0]^2. + new_pix_vec[*,1]^2.))*2.
                if image_len/xy_len lt .9 then begin
                  print, 'Image FoV is smaller than expected, increasing delta kperp to match image FoV'
                  delta_kperp_rad = 1/image_len
                endif
              endelse
            endif else count_far = 0
            
            n_u = round(max_u_lambda / delta_u_lambda) * 2 + 1
            u_lambda_vals = (findgen(n_u) - (n_u-1)/2) * delta_u_lambda
            
            ;; need to cut uvf cubes in half because image is real -- we'll cut in v
            ;; drop the unused half before the DFT to save time
            v_lambda_vals = u_lambda_vals[n_u/2:n_u-1]
            
          endif else begin
            ;;delta_kperp_rad = file_struct.kpix * z_mpc_mean / kperp_lambda_conv
          
            delta_kperp_rad = uv_cellsize_m * mean(frequencies*1e6) * z_mpc_mean / (3e8 * kperp_lambda_conv)
            ;; go a little beyond max_baseline to account for expansion due to w projection
            ;; max_kperp_rad = (file_struct.max_baseline_lambda/kperp_lambda_conv) * z_mpc_mean * 1.1
            ;; use kspan of Ian's cubes
            if tag_exist(file_struct, 'kspan') then begin
              max_kperp_rad = (min([file_struct.kspan/2.,file_struct.max_baseline_lambda])/kperp_lambda_conv) * z_mpc_mean
            endif else max_kperp_rad = (min([file_struct.max_baseline_lambda])/kperp_lambda_conv) * z_mpc_mean
            
            if keyword_set(cut_image) then begin
              ;; limit field of view to match calculated k-modes
              xy_len = 2*!pi/delta_kperp_rad
              wh_close = where(new_pix_vec[*,0] le xy_len/2 and new_pix_vec[*,0] ge -1*xy_len/2 and $
                new_pix_vec[*,1] le xy_len/2 and new_pix_vec[*,1] ge -1*xy_len/2, count_close, $
                ncomplement = count_far)
              if count_far ne 0 then new_pix_vec = new_pix_vec[wh_close, *] else begin
                ;; image may be smaller than expected, may need to adjust delta_kperp_rad
                image_len = max(sqrt(new_pix_vec[*,0]^2. + new_pix_vec[*,1]^2.))*2.
                if image_len/xy_len lt .9 then begin
                  print, 'Image FoV is smaller than expected, increasing delta kperp to match image FoV'
                  delta_kperp_rad = 2*!pi/image_len
                endif
              endelse
            endif else count_far = 0
            
            n_kperp = round(max_kperp_rad / delta_kperp_rad) * 2 + 1
            kx_rad_vals = (findgen(n_kperp) - (n_kperp-1)/2) * delta_kperp_rad
            
            ;; need to cut uvf cubes in half because image is real -- we'll cut in v
            ;; drop the unused half before the DFT to save time
            ky_rad_vals = kx_rad_vals[n_kperp/2:n_kperp-1]
            
          endelse
          
          ;; do DFT.
          if test_uvf eq 0 or keyword_set(dft_refresh_data) then begin
            print, 'calculating DFT for ' + file_struct.datavar + ' in ' + file_struct.datafile[i]
            
            arr = getvar_savefile(file_struct.datafile[i], file_struct.datavar)
            
            if size(arr,/type) eq 6 or size(arr,/type) eq 9 then $
              if max(abs(imaginary(arr))) eq 0 then arr = real_part(arr) else stop
            if not healpix then arr = reform(arr, dims[0]*dims[1], n_freq)
            
            if count_far ne 0 then arr = arr[wh_close, *]
            if n_elements(freq_ch_range) ne 0 then arr = arr[*, min(freq_ch_range):max(freq_ch_range)]
            
            if keyword_set(dft_ian) then begin
              transform = discrete_ft_2D_fast(new_pix_vec[*,0], new_pix_vec[*,1], arr, u_lambda_vals, v_lambda_vals, /exp2pi, $
                timing = ft_time, fchunk = dft_fchunk)
            endif else begin
              transform = discrete_ft_2D_fast(new_pix_vec[*,0], new_pix_vec[*,1], arr, kx_rad_vals, ky_rad_vals, $
                timing = ft_time, fchunk = dft_fchunk)
            endelse
            data_cube = temporary(transform)
            undefine, arr
            
            if keyword_set(dft_ian) then save, file = file_struct.uvf_savefile[i], u_lambda_vals, v_lambda_vals, data_cube $
            else save, file = file_struct.uvf_savefile[i], kx_rad_vals, ky_rad_vals, data_cube
            undefine, data_cube
          endif
          
          if test_wt_uvf eq 0 or keyword_set(dft_refresh_weight) then begin
            print, 'calculating DFT for ' + file_struct.weightvar + ' in ' + file_struct.weightfile[i]
            arr = getvar_savefile(file_struct.weightfile[i], file_struct.weightvar)
            
            if size(arr,/type) eq 6 or size(arr,/type) eq 9 then $
              if max(abs(imaginary(arr))) eq 0 then arr = real_part(arr) else stop
            if not healpix then arr = reform(arr, dims[0]*dims[1], n_freq)
            
            if count_far ne 0 then arr = arr[wh_close, *]
            if n_elements(freq_ch_range) ne 0 then arr = arr[*, min(freq_ch_range):max(freq_ch_range)]
            
            if keyword_set(dft_ian) then begin
              transform = discrete_ft_2D_fast(new_pix_vec[*,0], new_pix_vec[*,1], arr, u_lambda_vals, v_lambda_vals, /exp2pi, $
                timing = ft_time, fchunk = dft_fchunk)
            endif else begin
              transform = discrete_ft_2D_fast(new_pix_vec[*,0], new_pix_vec[*,1], arr, kx_rad_vals, ky_rad_vals, $
                timing = ft_time, fchunk = dft_fchunk)
            endelse
            
            weights_cube = temporary(transform)
            
            if not no_var then begin
              print, 'calculating DFT for ' + file_struct.variancevar + ' in ' + file_struct.variancefile[i]
              arr = getvar_savefile(file_struct.variancefile[i], file_struct.variancevar)
              
              if size(arr,/type) eq 6 or size(arr,/type) eq 9 then $
                if max(abs(imaginary(arr))) eq 0 then arr = real_part(arr) else stop
              if not healpix then arr = reform(arr, dims[0]*dims[1], n_freq)
              
              if count_far ne 0 then arr = arr[wh_close, *]
              if n_elements(freq_ch_range) ne 0 then arr = arr[*, min(freq_ch_range):max(freq_ch_range)]
              
              if keyword_set(dft_ian) then begin
                transform = discrete_ft_2D_fast(new_pix_vec[*,0], new_pix_vec[*,1], arr, u_lambda_vals, v_lambda_vals, /exp2pi, $
                  timing = ft_time, fchunk = dft_fchunk)
              endif else begin
                transform = discrete_ft_2D_fast(new_pix_vec[*,0], new_pix_vec[*,1], arr, kx_rad_vals, ky_rad_vals, $
                  timing = ft_time, fchunk = dft_fchunk)
              endelse
              variance_cube = abs(temporary(transform)) ;; make variances real, positive definite (amplitude)
              undefine, arr
            endif
            
            if keyword_set(dft_ian) then $
              save, file = file_struct.uvf_weight_savefile[i], u_lambda_vals, v_lambda_vals, weights_cube, variance_cube else $
              save, file = file_struct.uvf_weight_savefile[i], kx_rad_vals, ky_rad_vals, weights_cube, variance_cube
            undefine, new_pix_vec, weights_cube, variance_cube
          endif
        endelse
      endif else begin
        if keyword_set(dft_ian) then begin
          u_lambda_vals = getvar_savefile(file_struct.uvf_savefile[0], 'u_lambda_vals')
          v_lambda_vals = getvar_savefile(file_struct.uvf_savefile[0], 'v_lambda_vals')
        endif else begin
          kx_rad_vals = getvar_savefile(file_struct.uvf_savefile[0], 'kx_rad_vals')
          ky_rad_vals = getvar_savefile(file_struct.uvf_savefile[0], 'ky_rad_vals')
        endelse
      endelse
    endfor
  endif
  
  if healpix or keyword_set(image) then begin
    if keyword_set(dft_ian) then begin
      n_kx = n_elements(u_lambda_vals)
      kx_mpc = 2.*!pi*temporary(u_lambda_vals) / z_mpc_mean
      kx_mpc_delta = kx_mpc[1] - kx_mpc[0]
      
      n_ky = n_elements(v_lambda_vals)
      ky_mpc = 2.*!pi*temporary(v_lambda_vals) / z_mpc_mean
      ky_mpc_delta = ky_mpc[1] - ky_mpc[0]
    endif else begin
      n_kx = n_elements(kx_rad_vals)
      kx_rad_delta = kx_rad_vals[1] - kx_rad_vals[0]
      kx_mpc = temporary(kx_rad_vals) / z_mpc_mean
      kx_mpc_delta = kx_mpc[1] - kx_mpc[0]
      
      n_ky = n_elements(ky_rad_vals)
      ky_rad_delta = ky_rad_vals[1] - ky_rad_vals[0]
      ky_mpc = temporary(ky_rad_vals) / z_mpc_mean
      ky_mpc_delta = ky_mpc[1] - ky_mpc[0]
    endelse
    
    if healpix then begin
      ;; Angular resolution is given in Healpix paper in units of arcminutes, need to convert to radians
      ang_resolution = sqrt(3./!pi) * 3600./file_struct.nside * (1./60.) * (!pi/180.)
      pix_area_rad = ang_resolution^2. ;; by definition of ang. resolution in Healpix paper
    endif else pix_area_rad = (abs(file_struct.degpix) * !pi / 180d)^2.
    
    pix_area_mpc = pix_area_rad * z_mpc_mean^2.
    
  endif else begin
  
    x_rad_delta = abs(file_struct.degpix) * !pi / 180d
    n_kx = dims[0]
    x_rad_length = dims[0] * x_rad_delta
    if abs(file_struct.kpix-1/x_rad_length)/file_struct.kpix gt 1e-4 then stop
    
    kx_mpc_delta = (2.*!dpi)*file_struct.kpix / z_mpc_mean
    kx_mpc = (dindgen(n_kx)-n_kx/2) * kx_mpc_delta
    
    y_rad_delta = abs(file_struct.degpix) * !pi / 180.
    n_ky = dims[1]
    y_rad_length = dims[1] * y_rad_delta
    
    ky_mpc_delta = (2.*!dpi)*file_struct.kpix / z_mpc_mean
    ky_mpc = (dindgen(n_ky)-n_ky/2) * kx_mpc_delta
    
    pix_area_rad = x_rad_delta * y_rad_delta
    
  endelse
  
  if nfiles eq 2 then begin
    if no_var then begin
      ;; use 1/abs(weights) instead
      if healpix or keyword_set(image) then begin
        weights_cube1 = getvar_savefile(file_struct.uvf_weight_savefile[0], 'weights_cube')
        weights_cube2 = getvar_savefile(file_struct.uvf_weight_savefile[1], 'weights_cube')
        
        if min(ky_mpc) lt 0 then begin
          ;; negative ky values haven't been cut yet
          ;; need to cut uvf cubes in half because image is real -- we'll cut negative ky
          weights_cube1 = weights_cube1[*, n_ky/2:n_ky-1,*]
          weights_cube2 = weights_cube2[*, n_ky/2:n_ky-1,*]
        endif
        
      endif else begin
        weights_cube1 = getvar_savefile(file_struct.weightfile[0], file_struct.weightvar)
        weights_cube2 = getvar_savefile(file_struct.weightfile[1], file_struct.weightvar)
        
        ;; need to cut uvf cubes in half because image is real -- we'll cut negative ky
        weights_cube1 = weights_cube1[*, n_ky/2:n_ky-1,*]
        weights_cube2 = weights_cube2[*, n_ky/2:n_ky-1,*]
      endelse
      
      ;; Also need to drop 1/2 of ky=0 line -- drop ky=0, kx<0
      weights_cube1[0:n_kx/2-1, 0, *] = 0
      weights_cube2[0:n_kx/2-1, 0, *] = 0
      
    endif else begin
      if healpix or keyword_set(image) then begin
        variance_cube1 = getvar_savefile(file_struct.uvf_weight_savefile[0], 'variance_cube')
        variance_cube2 = getvar_savefile(file_struct.uvf_weight_savefile[1], 'variance_cube')
        weights_cube1 = getvar_savefile(file_struct.uvf_weight_savefile[0], 'weights_cube')
        weights_cube2 = getvar_savefile(file_struct.uvf_weight_savefile[1], 'weights_cube')
        
        if min(ky_mpc) lt 0 then begin
          ;; calculate integral of window function before cut for comparison
          window_int_orig = [total(variance_cube1)*pix_area_rad/n_vis[0], $
            total(variance_cube2)*pix_area_rad/n_vis[1]]
            
          ;; negative ky values haven't been cut yet
          ;; need to cut uvf cubes in half because image is real -- we'll cut negative ky
          variance_cube1 = variance_cube1[*, n_ky/2:n_ky-1,*]
          variance_cube2 = variance_cube2[*, n_ky/2:n_ky-1,*]
          weights_cube1 = weights_cube1[*, n_ky/2:n_ky-1,*]
          weights_cube2 = weights_cube2[*, n_ky/2:n_ky-1,*]
        endif
        
      endif else begin
        ;; need to cut uvf cubes in half because image is real -- we'll cut negative ky
        variance_cube1 = getvar_savefile(file_struct.variancefile[0], file_struct.variancevar)
        if max(abs(imaginary(variance_cube1))) gt 0 then stop else variance_cube1 = real_part(variance_cube1)
        variance_cube1 = variance_cube1[*, n_ky/2:n_ky-1,*]
        variance_cube2 = getvar_savefile(file_struct.variancefile[1], file_struct.variancevar)
        if max(abs(imaginary(variance_cube2))) gt 0 then stop else variance_cube2 = real_part(variance_cube2)
        variance_cube2 = variance_cube2[*, n_ky/2:n_ky-1,*]
        
        weights_cube1 = getvar_savefile(file_struct.weightfile[0], file_struct.weightvar)
        weights_cube1 = weights_cube1[*, n_ky/2:n_ky-1,*]
        weights_cube2 = getvar_savefile(file_struct.weightfile[1], file_struct.weightvar)
        weights_cube2 = weights_cube2[*, n_ky/2:n_ky-1,*]
      endelse
      
      ;; Also need to drop 1/2 of ky=0 line -- drop ky=0, kx<0
      weights_cube1[0:n_kx/2-1, 0, *] = 0
      weights_cube2[0:n_kx/2-1, 0, *] = 0
      variance_cube1[0:n_kx/2-1, 0, *] = 0
      variance_cube2[0:n_kx/2-1, 0, *] = 0
      
      ;; calculate integral of window function (use pix_area_rad for FT normalization)
      ;; already cut out negative ky, so multiply by 2
      window_int = 2*[total(variance_cube1)*pix_area_rad/n_vis[0], $
        total(variance_cube2)*pix_area_rad/n_vis[1]]
        
    endelse
    
    ;; now get data cubes
    if healpix or keyword_set(image) then begin
      data_cube1 = getvar_savefile(file_struct.uvf_savefile[0], 'data_cube')
      data_cube2 = getvar_savefile(file_struct.uvf_savefile[1], 'data_cube')
      
      if min(ky_mpc) lt 0 then begin
        ;; negative ky values haven't been cut yet
        ;; need to cut uvf cubes in half because image is real -- we'll cut negative ky
        data_cube1 = data_cube1[*, n_ky/2:n_ky-1,*]
        data_cube2 = data_cube2[*, n_ky/2:n_ky-1,*]
        
        ky_mpc = ky_mpc[n_ky/2:n_ky-1]
        n_ky = n_elements(ky_mpc)
      endif
      
    endif else begin
      if datavar eq '' then begin
        ;; working with a 'derived' cube (ie residual cube or noise simulation cube) that is constructed from other cubes
        if keyword_set(noise_sim) then begin
          ;; need to cut uvf cubes in half because image is real -- we'll cut negative ky
          variance_cube1 = getvar_savefile(file_struct.variancefile[0], file_struct.variancevar)
          variance_cube1 = variance_cube1[*, n_ky/2:n_ky-1,*]
          variance_cube2 = getvar_savefile(file_struct.variancefile[1], file_struct.variancevar)
          variance_cube2 = variance_cube2[*, n_ky/2:n_ky-1,*]
          
          seed = systime(1)
          noise = randomn(seed, dims) * sqrt(variance_cube) + $
            complex(0,1)*randomn(seed, dims) * sqrt(variance_cube)
          data_cube = temporary(noise)
          save, file = file_struct.uvf_savefile[i], kx_rad_vals, ky_rad_vals, data_cube
          undefine, data_cube, seed, variance_cube
        endif else begin
          ;; need to cut uvf cubes in half because image is real -- we'll cut negative ky
          dirty_cube1 = getvar_savefile(input_uvf_files[0,0], input_uvf_varname[0,0])
          dirty_cube1 = dirty_cube1[*, n_ky/2:n_ky-1,*]
          model_cube1 = getvar_savefile(input_uvf_files[0,1], input_uvf_varname[0,1])
          model_cube1 = model_cube1[*, n_ky/2:n_ky-1,*]
          data_cube1 = temporary(dirty_cube1) - temporary(model_cube1)
          
          dirty_cube2 = getvar_savefile(input_uvf_files[1,0], input_uvf_varname[1,0])
          dirty_cube2 = dirty_cube2[*, n_ky/2:n_ky-1,*]
          model_cube2 = getvar_savefile(input_uvf_files[1,1], input_uvf_varname[1,1])
          model_cube2 = model_cube2[*, n_ky/2:n_ky-1,*]
          data_cube2 = temporary(dirty_cube2) - temporary(model_cube2)
        endelse
      endif else begin
        ;; need to cut uvf cubes in half because image is real -- we'll cut negative ky
        data_cube1 = getvar_savefile(file_struct.datafile[0], file_struct.datavar)
        data_cube1 = data_cube1[*, n_ky/2:n_ky-1,*]
        data_cube2 = getvar_savefile(file_struct.datafile[1], file_struct.datavar)
        data_cube2 = data_cube2[*, n_ky/2:n_ky-1,*]
      endelse
      
      ky_mpc = ky_mpc[n_ky/2:n_ky-1]
      n_ky = n_elements(ky_mpc)
    endelse
    
    ;; Also need to drop 1/2 of ky=0 line -- drop ky=0, kx<0
    data_cube1[0:n_kx/2-1, 0, *] = 0
    data_cube2[0:n_kx/2-1, 0, *] = 0
    
    
  endif else begin
    ;; single file mode
    if healpix or keyword_set(image) then begin
      data_cube1 = getvar_savefile(file_struct.uvf_savefile, 'data_cube')
      weights_cube1 = getvar_savefile(file_struct.uvf_weight_savefile, 'weights_cube')
      
      if min(ky_mpc) lt 0 then begin
        ;; negative ky values haven't been cut yet
        ;; need to cut uvf cubes in half because image is real -- we'll cut negative ky
        weights_cube1 = weights_cube1[*, n_ky/2:n_ky-1,*]
        data_cube1 = data_cube1[*, n_ky/2:n_ky-1,*]
      endif
      
      ;; Also need to drop 1/2 of ky=0 line -- drop ky=0, kx<0
      weights_cube1[0:n_kx/2-1, 0, *] = 0
      data_cube1[0:n_kx/2-1, 0, *] = 0
      
      if not no_var then begin
        variance_cube1 = getvar_savefile(file_struct.uvf_weight_savefile, 'variance_cube')
        
        if min(ky_mpc) lt 0 then begin
          ;; negative ky values haven't been cut yet
          ;; need to cut uvf cubes in half because image is real -- we'll cut negative ky
          variance_cube1 = variance_cube1[*, n_ky/2:n_ky-1,*]
        endif
        
        ;; Also need to drop 1/2 of ky=0 line -- drop ky=0, kx<0
        variance_cube1[0:n_kx/2-1, 0, *] = 0
        
        ;; calculate integral of window function
        ;; already cut out negative ky, so multiply by 2
        window_int = 2*total(variance_cube1)*pix_area_rad/n_vis
        
      endif
      
      if min(ky_mpc) lt 0 then begin
        ;; negative ky values haven't been cut yet
        ;; need to cut uvf cubes in half because image is real -- we'll cut negative ky
        ky_mpc = ky_mpc[n_ky/2:n_ky-1]
        n_ky = n_elements(ky_mpc)
      endif
      
    endif else begin
      ;; need to cut uvf cubes in half because image is real -- we'll cut negative ky
      data_cube1 = getvar_savefile(file_struct.datafile, file_struct.datavar)
      data_cube1 = data_cube1[*, n_ky/2:n_ky-1,*]
      weights_cube1 = getvar_savefile(file_struct.weightfile, file_struct.weightvar)
      weights_cube1 = weights_cube1[*, n_ky/2:n_ky-1,*]
      
      ;; Also need to drop 1/2 of ky=0 line -- drop ky=0, kx<0
      weights_cube1[0:n_kx/2-1, 0, *] = 0
      data_cube1[0:n_kx/2-1, 0, *] = 0
      
      if not no_var then begin
        variance_cube1 = getvar_savefile(file_struct.variancefile, file_struct.variance_var)
        ;; need to cut uvf cubes in half because image is real -- we'll cut negative ky
        variance_cube1 = variance_cube1[*, n_ky/2:n_ky-1,*]
        
        ;; Also need to drop 1/2 of ky=0 line -- drop ky=0, kx<0
        variance_cube1[0:n_kx/2-1, 0, *] = 0
        
        ;; calculate integral of window function
        ;; already cut out negative ky, so multiply by 2
        window_int = 2*total(variance_cube1)*pix_area_rad/n_vis
        
        
      endif
      
      ;; need to cut uvf cubes in half because image is real -- we'll cut negative ky
      ky_mpc = ky_mpc[n_ky/2:n_ky-1]
      n_ky = n_elements(ky_mpc)
      
    endelse
  endelse
  
  ;; save some slices of the raw data cube (before dividing by weights)
  for i=0, nfiles-1 do begin
    if i eq 0 then data_cube = data_cube1 else data_cube = data_cube2
    uf_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 1, $
      slice_inds = 0, slice_savefile = file_struct.uf_raw_savefile[i])
      
    vf_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 0, $
      slice_inds = n_kx/2, slice_savefile = file_struct.vf_raw_savefile[i])
      
    if max(abs(vf_slice)) eq 0 then begin
      nloop = 0
      while max(abs(vf_slice)) eq 0 do begin
        nloop = nloop+1
        vf_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, $
          slice_axis = 0, slice_inds = n_kx/2+nloop, slice_savefile = file_struct.vf_raw_savefile[i])
      endwhile
    endif
    
    uv_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 2, $
      slice_inds = 0, slice_savefile = file_struct.uv_raw_savefile[i])
  endfor
  
  if healpix or keyword_set(image) then begin
    if keyword_set (dft_ian) then begin
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
      
    endif else begin
      ;; multiply data, weights & variance cubes by pixel_area_mpc to get proper units from DFT
      ;; divide by (2*!pi)^2d to get right FT convention
      ;; (not squared for variance because they weren't treated as units squared in FHD code)
      data_cube1 = data_cube1 * pix_area_mpc
      weights_cube1 = weights_cube1 * pix_area_mpc
      variance_cube1 = variance_cube1 * pix_area_mpc
      if nfiles eq 2 then begin
        data_cube2 = data_cube2 * pix_area_mpc
        weights_cube2 = weights_cube2 * pix_area_mpc
        variance_cube2 = variance_cube2 * pix_area_mpc
      endif
      
      ;; apply FT jacobian because fft was in Jy/str and dft was in Jy/Mpc^2
      ;; same for variance as data
      data_cube1 = data_cube1 / (z_mpc_mean^2.)
      weights_cube1 = weights_cube1 / (z_mpc_mean^2.)
      variance_cube1 = variance_cube1 / (z_mpc_mean^2.)
      if nfiles eq 2 then begin
        data_cube2 = data_cube2 / (z_mpc_mean^2.)
        weights_cube2 = weights_cube2 / (z_mpc_mean^2.)
        variance_cube2 = variance_cube2 / (z_mpc_mean^2.)
      endif
    endelse
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
    if min(sigma2_cube1) lt 0 or min(sigma2_cube2) lt 0 then message, 'sigma2 should be positive definite.'
    if total(abs(sigma2_cube1)) le 0 or total(abs(sigma2_cube2)) le 0 then message, 'one or both sigma2 cubes is all zero'
  endif else begin
    if min(sigma2_cube1) lt 0 then message, 'sigma2 should be positive definite.'
    if total(abs(sigma2_cube1)) le 0 then message, 'sigma2 cube is all zero'
  endelse
  
  ;; divide data by weights
  data_cube1 = data_cube1 / weights_cube1
  if count_wt1_0 ne 0 then data_cube1[wh_wt1_0] = 0
  undefine, weights_cube1, wh_wt1_0, count_wt1_0
  
  if nfiles eq 2 then begin
    data_cube2 = data_cube2 / weights_cube2
    if count_wt2_0 ne 0 then data_cube2[wh_wt2_0] = 0
    undefine, weights_cube2, wh_wt2_0, count_wt2_0
  endif
  
  if keyword_set(noise_sim) then begin
    noise = randomn(seed, dims) * sqrt(sigma2_cube1) + $
      complex(0,1)*randomn(seed, dims) * sqrt(sigma2_cube1)
    data_cube1 = temporary(noise)
  endif
  
  ;; take care of FT convention for EoR
  data_cube1 = data_cube1 / (2.*!pi)^2.
  sigma2_cube1 = sigma2_cube1 / (2.*!pi)^4.
  if nfiles eq 2 then begin
    data_cube2 = data_cube2 / (2.*!pi)^2.
    sigma2_cube2 = sigma2_cube2 / (2.*!pi)^4.
  endif
  
  ;; get sigma^2 into Jy^2
  sigma2_cube1 = sigma2_cube1 * rebin(reform(vis_sigma, 1, 1, n_freq), n_kx, n_ky, n_freq)^2.
  if nfiles eq 2 then sigma2_cube2 = sigma2_cube2 * rebin(reform(vis_sigma, 1, 1, n_freq), n_kx, n_ky, n_freq)^2.
  
  ;; get data & sigma into mK Mpc^2
  for i=0, n_freq-1 do begin
    data_cube1[*,*,i] = data_cube1[*,*,i]*conv_factor[i]
    sigma2_cube1[*,*,i] = sigma2_cube1[*,*,i]*(conv_factor[i])^2.
    if nfiles eq 2 then begin
      data_cube2[*,*,i] = data_cube2[*,*,i]*conv_factor[i]
      sigma2_cube2[*,*,i] = sigma2_cube2[*,*,i]*(conv_factor[i])^2.
    endif
  endfor
  
  ;; fix units on window funtion integral -- now they should be Mpc^3
  ;; use delta f for native visibility frequency resolution
  ;; calculate integral of window in r^3 to compare w/ Adam
  window_int_r = window_int * (z_mpc_delta * n_freq) * (kx_mpc_delta * ky_mpc_delta)*z_mpc_mean^4./((2.*!dpi)^2.*file_struct.kpix^4.)
  print, 'window integral in r^3: ' + number_formatter(window_int_r[0], format='(e10.4)')
  
  ;; need window integral in k^3
  window_int = window_int_r * (2*!dpi)^3.
  
  ;; divide data by sqrt(window_int) and sigma2 by window_int
  data_cube1 = data_cube1 / sqrt(window_int[0])
  sigma2_cube1 = sigma2_cube1 / window_int[0]
  if nfiles eq 2 then begin
    data_cube2 = data_cube2 / sqrt(window_int[1])
    sigma2_cube2 = sigma2_cube2  / window_int[1]
  endif
  
  ;; temp = data_cube1[*,*,10]
  ;; temp_exp = complex(dblarr(165, 83+82))
  ;; temp_exp[*, 82:*] = temp
  ;; temp_exp[*, 0:81] = conj(reverse(reverse(temp[*,1:*]),2))
  ;; temp_exp[0:81, 82] =  reverse(conj(temp[83:*,0]))
  ;; temp_exp = temp_exp[0:163, 0:163]
  ;; temp_img = fft_shift(fft(fft_shift(temp_exp)))
  ;; temp_2 = fft_shift(fft(fft_shift(temp_img), /inverse))
  ;; quick_image, abs(temp_exp), title='107',/log
  
  ;; save some slices of the data cube
  for i=0, nfiles-1 do begin
    if i eq 0 then data_cube = data_cube1 else data_cube = data_cube2
    uf_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 1, $
      slice_inds = 0, slice_savefile = file_struct.uf_savefile[i])
      
    vf_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 0, $
      slice_inds = n_kx/2, slice_savefile = file_struct.vf_savefile[i])
      
    if max(abs(vf_slice)) eq 0 then begin
      nloop = 0
      while max(abs(vf_slice)) eq 0 and (n_kx/2+nloop) lt (n_kx-1) do begin
        nloop = nloop+1
        vf_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, $
          slice_axis = 0, slice_inds = n_kx/2+nloop, slice_savefile = file_struct.vf_savefile[i])
      endwhile
    endif
    
    uv_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 2, $
      slice_inds = 0, slice_savefile = file_struct.uv_savefile[i])
  endfor
  
  if nfiles eq 2 then begin
    ;; Now construct added & subtracted cubes (weighted by inverse variance) & new variances
    sum_weights1 = 1./sigma2_cube1
    wh_sig1_0 = where(sigma2_cube1 eq 0, count_sig1_0)
    if count_sig1_0 ne 0 then sum_weights1[wh_sig1_0] = 0
    undefine, sigma2_cube1, wh_sig1_0, count_sig1_0
    
    sum_weights2 = 1./sigma2_cube2
    wh_sig2_0 = where(sigma2_cube2 eq 0, count_sig2_0)
    if count_sig2_0 ne 0 then sum_weights2[wh_sig2_0] = 0
    undefine, sigma2_cube2, wh_sig2_0, count_sig2_0
    
    sum_weights_net = sum_weights1 + sum_weights2
    wh_wt0 = where(sum_weights_net eq 0, count_wt0)
    
    data_sum = (sum_weights1 * data_cube1 + sum_weights2 * data_cube2)/sum_weights_net
    data_diff = (sum_weights1 * data_cube1 - sum_weights2 * data_cube2)/sum_weights_net
    undefine, data_cube1, data_cube2
    
    if count_wt0 ne 0 then data_sum[wh_wt0] = 0
    if count_wt0 ne 0 then data_diff[wh_wt0] = 0
    undefine, sum_weights1, sum_weights2
    
    sum_sigma2 = 1./temporary(sum_weights_net)
    if count_wt0 ne 0 then sum_sigma2[wh_wt0] = 0
    
  endif else begin
    data_sum = temporary(data_cube1)
    sum_sigma2 = temporary(sigma2_cube1)
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
    slice_inds = 0, slice_savefile = file_struct.uf_sum_savefile)
  if nfiles eq 2 then $
    uf_slice = uvf_slice(data_diff, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 1, $
    slice_inds = 0, slice_savefile = file_struct.uf_diff_savefile)
    
  vf_slice = uvf_slice(data_sum, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 0, $
    slice_inds = n_kx/2, slice_savefile = file_struct.vf_sum_savefile)
  if nfiles eq 2 then $
    vf_slice2 = uvf_slice(data_diff, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 0, $
    slice_inds = n_kx/2, slice_savefile = file_struct.vf_diff_savefile) $
  else vf_slice2 = vf_slice
  
  test_vf = 1
  if max(abs(vf_slice)) eq 0 then test_vf=0
  if nfiles eq 2 then if max(vf_slice2) eq 0 and max(abs(data_diff)) gt 0 then test_vf=0
  
  if test_vf eq 0 then begin
    nloop = 0
    slice_inds = shift(indgen(n_kx), n_kx/2)
    while test_vf eq 0 do begin
      nloop = nloop+1
      vf_slice = uvf_slice(data_sum, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, $
        slice_axis = 0, slice_inds = slice_inds[nloop], slice_savefile = file_struct.vf_sum_savefile)
      if nfiles eq 2 then $
        vf_slice2 = uvf_slice(data_diff, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, $
        slice_axis = 0, slice_inds = slice_inds[nloop], slice_savefile = file_struct.vf_diff_savefile)
        
      test_vf = 1
      if max(abs(vf_slice)) eq 0 then test_vf=0
      if nfiles eq 2 then if max(vf_slice2) eq 0 and max(abs(data_diff)) gt 0 then test_vf=0
    endwhile
  endif
  
  uv_slice = uvf_slice(data_sum, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 2, $
    slice_inds = 0, slice_savefile = file_struct.uv_sum_savefile)
  if nfiles eq 2 then $
    uv_slice = uvf_slice(data_diff, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 2, $
    slice_inds = 0, slice_savefile = file_struct.uv_diff_savefile)
    
  ;; apply spectral windowing function if desired
  if n_elements(spec_window_type) ne 0 then begin
    window = spectral_window(n_freq, type = spec_window_type, /periodic)
    
    norm_factor = n_freq / total(window)
    window = window * norm_factor
    
    window_expand = rebin(reform(window, 1, 1, n_freq), n_kx, n_ky, n_freq, /sample)
    
    data_sum = data_sum * window_expand
    if nfiles eq 2 then data_diff = data_diff * window_expand
    
    sum_sigma2 = sum_sigma2 * temporary(window_expand^2.)
  endif
  
  ;; now take frequency FT
  if even_freq then begin
    ;; evenly spaced, just use fft
    data_sum_ft = fft(data_sum, dimension=3) * n_freq * z_mpc_delta / (2.*!dpi)
    ;; put k0 in middle of cube
    data_sum_ft = shift(data_sum_ft, [0,0,n_kz/2])
    undefine, data_sum
    if nfiles eq 2 then begin
      data_diff_ft = fft(data_diff, dimension=3) * n_freq * z_mpc_delta / (2.*!dpi)
      ;; put k0 in middle of cube
      data_diff_ft = shift(data_diff_ft, [0,0,n_kz/2])
      undefine, data_diff
    endif
  endif else begin
    ;; Not evenly spaced. Do a dft
    z_exp =  exp(-1.*complex(0,1)*matrix_multiply(comov_dist_los, kz_mpc_orig, /btranspose))
    
    data_sum_ft = z_mpc_delta/(2.*!dpi) * $
      reform(matrix_multiply(reform(temporary(data_sum), n_kx*n_ky, n_freq), z_exp), n_kx, n_ky, n_kz)
    if nfiles eq 2 then $
      data_diff_ft = z_mpc_delta/(2.*!dpi) * $
      reform(matrix_multiply(reform(temporary(data_diff), n_kx*n_ky, n_freq), z_exp), n_kx, n_ky, n_kz)
  endelse
  
  n_val = round(kz_mpc_orig / kz_mpc_delta)
  kz_mpc_orig[where(n_val eq 0)] = 0
  ;; these an and bn calculations don't match the standard
  ;; convention (they down by a factor of 2) but they make more sense
  ;; and remove factors of 2 we'd otherwise have in the power
  ;; and variance calculations
  a1_0 = data_sum_ft[*,*,where(n_val eq 0)]
  a1_n = (data_sum_ft[*,*, where(n_val gt 0)] + data_sum_ft[*,*, reverse(where(n_val lt 0))])/2.
  b1_n = complex(0,1) * (data_sum_ft[*,*, where(n_val gt 0)] - data_sum_ft[*,*, reverse(where(n_val lt 0))])/2.
  undefine, data_sum_ft
  
  if nfiles gt 1 then begin
    a2_0 = data_diff_ft[*,*,where(n_val eq 0)]
    a2_n = (data_diff_ft[*,*, where(n_val gt 0)] + data_diff_ft[*,*, reverse(where(n_val lt 0))])/2.
    b2_n = complex(0,1) * (data_diff_ft[*,*, where(n_val gt 0)] - data_diff_ft[*,*, reverse(where(n_val lt 0))])/2.
    undefine, data_diff_ft
  endif
  
  kz_mpc = kz_mpc_orig[where(n_val ge 0)]
  n_kz = n_elements(kz_mpc)
  
  if keyword_set(std_power) then begin
    ;; for standard power calc. just need ft of sigma2 (sigma has squared units relative to data, so use z_mpc_delta^2d)
    sigma2_ft = fft(sum_sigma2, dimension=3) * n_freq * z_mpc_delta^2. / (2.*!pi)
    sigma2_ft = shift(sigma2_ft, [0,0,n_kz/2])
    undefine, sum_sigma2
    
    sigma_a0 = abs(sigma2_ft[*,*,where(n_val eq 0)])
    sigma_an_bn = sqrt(abs(sigma2_ft[*,*, where(n_val gt 0)])^2. + abs(sigma2_ft[*,*, reverse(where(n_val lt 0))])^2.)/2.
    undefine, sigma2_ft
    
    save, file = file_struct.kcube_savefile, a1_0, a1_n, b1_n, a2_0, a2_n, b2_n, sigma_a0, sigma_an_bn, $
      kx_mpc, ky_mpc, kz_mpc, kperp_lambda_conv, delay_params, hubble_param, n_freq_contrib
      
  endif else begin
  
    ;; drop pixels with less than 1/3 of the frequencies (set weights to 0)
    wh_fewfreq = where(n_freq_contrib lt ceil(n_freq/3d), count_fewfreq)
    if count_fewfreq gt 0 then begin
      mask_fewfreq = n_freq_contrib * 0 + 1
      mask_fewfreq[wh_fewfreq] = 0
      mask_fewfreq = rebin(temporary(mask_fewfreq), n_kx, n_ky, n_kz)
      
      a1_0 = temporary(a1_0) * mask_fewfreq[*,*,0]
      a1_n = temporary(a1_n) * mask_fewfreq[*,*,1:*]
      b1_n = temporary(b1_n) * mask_fewfreq[*,*,1:*]
      if nfiles gt 1 then begin
        a2_0 = temporary(a2_0) * mask_fewfreq[*,*,0]
        a2_n = temporary(a2_n) * mask_fewfreq[*,*,1:*]
        b2_n = temporary(b2_n) * mask_fewfreq[*,*,1:*]
      endif
    endif
    
    data_sum_cos = complex(fltarr(n_kx, n_ky, n_kz))
    data_sum_sin = complex(fltarr(n_kx, n_ky, n_kz))
    data_sum_cos[*, *, 0] = a1_0 /2.
    data_sum_cos[*, *, 1:n_kz-1] = a1_n
    data_sum_sin[*, *, 1:n_kz-1] = b1_n
    if nfiles gt 1 then begin
      data_diff_cos = complex(fltarr(n_kx, n_ky, n_kz))
      data_diff_sin = complex(fltarr(n_kx, n_ky, n_kz))
      data_diff_cos[*, *, 0] = a2_0 /2.
      data_diff_cos[*, *, 1:n_kz-1] = a2_n
      data_diff_sin[*, *, 1:n_kz-1] = b2_n
    endif
    
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
    covar_cos = matrix_multiply(sum_sigma2, cos_arr^2d) * (z_mpc_delta / (2.*!pi))^2.
    covar_sin = matrix_multiply(sum_sigma2, sin_arr^2d) * (z_mpc_delta / (2.*!pi))^2.
    covar_cross = matrix_multiply(sum_sigma2, cos_arr*sin_arr) * (z_mpc_delta / (2.*!pi))^2.
    
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
    
    ;; cos 0 term has different normalization
    covar_cos[*,*,0] = covar_cos[*,*,0]/4d
    undefine, sum_sigma2, freq_kz_arr, cos_arr, sin_arr
    
    ;; factor to go to eor theory FT convention
    ;; I don't think I need these factors in the covariance
    ;; matrix because I've use the FT & inv FT -- should cancel
    ;; covar_cos = factor * temporary(covar_cos2)
    ;; covar_sin = factor * temporary(covar_sin2)
    ;; covar_cross = factor * temporary(covar_cross)
    
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
    if nfiles eq 2 then begin
      data_diff_1 = data_diff_cos*cos_theta + data_diff_sin*sin_theta
      data_diff_2 = (-1d)*data_diff_cos*sin_theta + data_diff_sin*cos_theta
      undefine, data_diff_cos, data_diff_sin
    endif
    
    save, file = file_struct.kcube_savefile, data_sum_1, data_sum_2, data_diff_1, data_diff_2, sigma2_1, sigma2_2, $
      kx_mpc, ky_mpc, kz_mpc, kperp_lambda_conv, delay_params, hubble_param, n_freq_contrib
  endelse
  
end
