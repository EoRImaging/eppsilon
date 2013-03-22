pro fhd_kcube, file_struct, dft_refresh_data = dft_refresh_data, dft_refresh_weight = dft_refresh_weight, $
              dft_fchunk = dft_fchunk, std_power = std_power, input_units = input_units

  if n_elements(file_struct.nside) ne 0 then healpix = 1 else healpix = 0
  nfiles = n_elements(file_struct.datafile)

  if n_elements(input_units) eq 0 then input_units = 'jansky'
  units_enum = ['jansky', 'mk']
  wh = where(units_enum eq input_units, count)
  if count eq 0 then message, 'input units not recognized, options are: ' + units_enum
  
  for i=0, nfiles-1 do begin
     datafile_obj = obj_new('IDL_Savefile', file_struct.datafile[i])
     datafile_names = datafile_obj->names()
     datavar = strupcase(file_struct.datavar)
     wh = where(datafile_names eq datavar[i], count)
     if count eq 0 then message, 'specified datavar is not present in datafile (datafile=' + file_struct.datafile[i] + $
                                 ', datavar=' + file_struct.datavar[i] + ')'
     data_dims = datafile_obj->size(file_struct.datavar[i], /dimensions)
     obj_destroy, datafile_obj
 
     if i gt 0 then if total(abs(data_dims - dims)) ne 0 then message, 'data dimensions in files do not match'
 
     if not keyword_set(no_weighting) then begin
        weightfile_obj = obj_new('IDL_Savefile', file_struct.weightfile[i])
        weightfile_names = weightfile_obj->names()
        weightvar = strupcase(file_struct.weightvar)
        wh = where(weightfile_names eq weightvar[i], count)
        if count eq 0 then message, 'specified weightvar is not present in weightfile (weightfile=' + file_struct.weightfile[i] + $
                                    ', weightvar=' + file_struct.weightvar[i] + ')'
        weight_dims = weightfile_obj->size(weightvar[i], /dimensions)
        obj_destroy, weightfile_obj
        
        if total(abs(data_dims - weight_dims)) ne 0 then message, 'data and weight dimensions do not match'
        undefine, weight_dims
     endif
     dims = data_dims
     undefine, data_dims
  endfor
     
  if healpix then n_freq = dims[1] else n_freq = dims[2]
  frequencies = file_struct.frequencies
  if n_elements(frequencies) ne n_freq then message, 'number of frequencies does not match frequency dimension of data'
     
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
     
     f_delta = freq_resolution
     nominal_freqs = findgen(floor(((max(frequencies)-min(frequencies))/freq_resolution))+1)*freq_resolution + min(frequencies)
     nominal_z = z0_freq/nominal_freqs - 1
     comoving_distance_los, nominal_z, nominal_comov_dist_los
     nominal_comov_diffs = nominal_comov_dist_los - shift(nominal_comov_dist_los, -1)
     nominal_comov_diffs = nominal_comov_diffs[0:n_elements(nominal_comov_diffs)-2]
     
     z_mpc_delta = mean(nominal_comov_diffs)
     z_mpc_mean = mean(nominal_comov_dist_los)
     
  endif else begin
     
     f_delta = freq_diff[0]
     z_mpc_delta = mean(comov_los_diff)
     z_mpc_mean = mean(comov_dist_los)
     n_kz = n_freq
     
  endelse
  kperp_lambda_conv = z_mpc_mean / (2d*!dpi)
  delay_delta = 1e9/(n_freq*f_delta*1e6) ;; equivilent delay bin size for kparallel
  delay_max = delay_delta * n_freq/2d    ;; factor of 2 b/c of neg/positive
  delay_params = [delay_delta, delay_max]     
  
  z_mpc_length = max(comov_dist_los) - min(comov_dist_los) + z_mpc_delta
  kz_mpc_range =  (2.*!pi) / (z_mpc_delta)
  kz_mpc_delta = (2.*!pi) / z_mpc_length
  kz_mpc_orig = findgen(round(kz_mpc_range / kz_mpc_delta)) * kz_mpc_delta - kz_mpc_range/2.
  if n_elements(n_kz) ne 0 then begin
     if n_elements(kz_mpc_orig) ne n_kz then stop
  endif else n_kz = n_elements(kz_mpc_orig)
  
  print, 'z delta: ', z_mpc_delta
  print,  'kz delta: ', kz_mpc_delta
  
  max_baseline = file_struct.max_baseline
  if input_units eq 'jansky' then begin
     ;; beam_diameter_rad = (3d * 10^8d) / (frequencies * 10^6d * max_baseline)
     ;; beam_area_str = !pi * beam_diameter_rad^2d /4d
     
     ;; conv_factor = (10^(double(-26+16+3-12+23)) * 9d) / (beam_area_str * 2d * frequencies^2d * 1.38)
     ;; if max(conv_factor-conv_factor[0]) gt 1e-8 then stop else conv_factor = conv_factor[0]
     conv_factor = float(2. * max_baseline^2d / (!dpi * 1.38))
  endif else conv_factor = 1.
  
  
  if healpix then begin
     for i=0, nfiles-1 do begin     
        test_uvf = file_test(file_struct.uvf_savefile[i]) *  (1 - file_test(file_struct.uvf_savefile[i], /zero_length))
     
        if not keyword_set(no_weighting) then $
           test_wt_uvf = file_test(file_struct.uvf_weight_savefile[i]) * $
                         (1 - file_test(file_struct.uvf_weight_savefile[i], /zero_length))
        
        if test_uvf eq 0 or test_wt_uvf eq 0 or keyword_set(dft_refresh_data) or keyword_set(dft_refresh_weight) then begin          
           pixel_nums = getvar_savefile(file_struct.pixelfile[i], file_struct.pixelvar[i])
           pixel_dims = size(pixel_nums, /dimension)
           if total(abs(dims - pixel_dims)) ne 0 then message, 'pixel and data dimensions do not match'
         
           test_setup = file_test(file_struct.hpx_dftsetup_savefile[i]) * $
                        (1 - file_test(file_struct.hpx_dftsetup_savefile[i], /zero_length))
           if test_setup eq 0 then begin
              ;; figure out k values to calculate dft
              healpix_setup_ft, pixel_nums, file_struct.nside, new_pix_vec, limits, kx_rad_vals, ky_rad_vals, /quiet
              save, file = file_struct.hpx_dftsetup_savefile[i], new_pix_vec, limits, kx_rad_vals, ky_rad_vals
           endif else restore, file_struct.hpx_dftsetup_savefile[i]
           
           ;; do DFT.
           if test_uvf eq 0 or keyword_set(dft_refresh_data) then begin
              arr = getvar_savefile(file_struct.datafile[i], file_struct.datavar[i])
              transform = discrete_ft_2D_fast(new_pix_vec[*,0], new_pix_vec[*,1], arr, kx_rad_vals, ky_rad_vals, timing = ft_time, $
                                              fchunk = dft_fchunk)
              data_cube = temporary(transform)
              undefine, arr

              save, file = file_struct.uvf_savefile[i], kx_rad_vals, ky_rad_vals, data_cube
              undefine, data_cube
           endif
           
           if test_wt_uvf eq 0 or keyword_set(dft_refresh_weight) then begin
              arr = getvar_savefile(file_struct.weightfile[i], file_struct.weightvar[i])
              transform = discrete_ft_2D_fast(new_pix_vec[*,0], new_pix_vec[*,1], arr, kx_rad_vals, ky_rad_vals, timing = ft_time, $
                                              fchunk = dft_fchunk)
              
              weights_cube = abs(transform) ;; make weights real, positive definite (amplitude)
              undefine, arr, transform

              save, file = file_struct.uvf_weight_savefile[i], kx_rad_vals, ky_rad_vals, weights_cube
              undefine, new_pix_vec, weights_cube
           endif
           
        endif else begin
           kx_rad_vals = getvar_savefile(file_struct.uvf_savefile[0], 'kx_rad_vals')
           ky_rad_vals = getvar_savefile(file_struct.uvf_savefile[0], 'ky_rad_vals')
        endelse
     endfor
  endif


  ;; if 2 file mode and we're using weights then check that uvf weights are close enough to proceed.
  if nfiles eq 2 then begin
     if not keyword_set(no_weighting) then begin
        if healpix then begin
           weights_cube1 = getvar_savefile(file_struct.uvf_weight_savefile[0], 'weights_cube')
           weights_cube2 = getvar_savefile(file_struct.uvf_weight_savefile[1], 'weights_cube')
        endif else begin
           weights_cube1 = getvar_savefile(file_struct.weightfile[0], file_struct.weight_var[i])
           weights_cube2 = getvar_savefile(file_struct.weightfile[1], file_struct.weight_var[i])
        endelse

        if min(weights_cube1) lt 0 or min(weights_cube2) lt 0 then message, 'Weights should be positive definite.'
        if total(abs(weights_cube1)) le 0 or total(abs(weights_cube2)) le 0 then message, 'one or both weights cubes is all zero'
         
        ;; freq_channel = 0
        
        ;; if windowavailable(1) then wset, 1 else window, 1
        ;; quick_image, weights_cube1[*,*,freq_channel], kx_rad_vals, ky_rad_vals, /log, title = 'Even cube', $
        ;;              xtitle = 'kx (Mpc!U-1!N)', ytitle = 'ky (Mpc!U-1!N)', data_range = data_range
        ;; if windowavailable(2) then wset, 2 else window, 2
        ;; quick_image, weights_cube2[*,*,freq_channel], kx_rad_vals, ky_rad_vals, /log, title = 'Odd cube', $
        ;;              xtitle = 'kx (Mpc!U-1!N)', ytitle = 'ky (Mpc!U-1!N)', data_range = data_range


        diff = abs(weights_cube1 - weights_cube2)
        ;; quick_image, diff[*,*,freq_channel], kx_rad_vals, ky_rad_vals, /log, title = 'Difference', $
        ;;              xtitle = 'kx (Mpc!U-1!N)', ytitle = 'ky (Mpc!U-1!N)'

        weights_ave = (weights_cube1 + weights_cube2) / 2.
        diff_frac = temporary(diff) / weights_ave

        ;; if windowavailable(3) then wset, 3 else window, 3
        ;; quick_image, diff_frac[*,*,freq_channel], kx_rad_vals, ky_rad_vals, /log, title = 'Difference Fraction', $
        ;;              xtitle = 'kx (Mpc!U-1!N)', ytitle = 'ky (Mpc!U-1!N)'

        ;; diff_frac_mask = fix(diff_frac*0)
        ;; wh_large = where(diff_frac ge 1, count_large)
        ;; if count_large gt 0 then diff_frac_mask[wh_large] = 1
        ;; ;;quick_image, diff_frac_mask[*,*,freq_channel], kx_rad_vals, ky_rad_vals, title = 'Difference Fraction above 1', $
        ;; ;;             xtitle = 'kx (Mpc!U-1!N)', ytitle = 'ky (Mpc!U-1!N)', data_range = [0,1], /grey_scale

        cube_max = max([weights_cube1, weights_cube2])
        cube_min_n0 = min([weights_cube1[where(weights_cube1 ne 0)], weights_cube2[where(weights_cube2 ne 0)]])
        ;; cube_binsize = 10.^(floor(alog10(cube_max))-2)
        ;;cube1_hist = histogram(weights_cube1, binsize = cube_binsize, min = 0, max = cube_max, locations = cube_locs)
        ;;cube2_hist = histogram(weights_cube2, binsize = cube_binsize, min = 0, max = cube_max, locations = cube_locs)
        ;;diff_hist = histogram(diff, binsize = cube_binsize*0.01, min = 0, omax = diff_max, locations = diff_locs)
        diff_frac_hist = histogram(diff_frac, binsize = 0.01, min = 0, omax = df_max, locations = diff_frac_locs, $
                                   reverse_indices = df_ri)

        min_log_wt = alog10(cube_min_n0)
        min_log_df = alog10(min(diff_frac[where(diff_frac gt 0)]))
        ;;wt_log_hist = histogram(alog10(weights_ave), binsize = 0.1, min = min_log_wt, omax = max_log_wt, locations = wt_log_locs)
        ;;df_log_hist = histogram(alog10(diff_frac), binsize = 0.1, min = min_log_df, omax = max_log_df, locations = df_log_locs)

        df_wt_hist = hist_2d(alog10(weights_ave), alog10(diff_frac), bin1=0.1, bin2=0.1, $
                             min1 = min_log_wt, min2=min_log_df, max1=max_log_wt, max2=max_log_df)

        ;; if windowavailable(4) then wset, 4 else window, 4
        ;; quick_image, df_wt_hist, 10^wt_log_locs, 10^df_log_locs,/log, /xlog, /ylog, xtitle = 'average weight', $
        ;;              ytitle = 'difference fraction', title = 'Weight vs difference fraction histogram'


        ;; ;;binarea = matrix_multiply(10^(wt_log_locs + 0.1) - 10^(wt_log_locs), 10^(df_log_locs + 0.1) - 10^(df_log_locs))
        
        ;; ;;quick_image, df_wt_hist/binarea, 10^wt_log_locs, 10^df_log_locs,/log, /xlog, /ylog, xtitle = 'average weight', $
        ;; ;;             ytitle = 'difference fraction', title = 'Weight vs difference fraction density'

        df_cut_level = 0.1
        wh_df_cut = where(diff_frac_locs gt df_cut_level, count_df_cut)
        if count_df_cut gt 0 then begin
           print, 'cutting out pixels with a weight difference greater than ' + number_formatter(df_cut_level*100.) + '%'

           cut_inds = df_ri[df_ri[min(wh_df_cut)]:df_ri[max(wh_df_cut)+1]-1]
           if (n_elements(cut_inds) / total(diff_frac_hist)) gt .1 then stop

           print, 'removed ' + number_formatter(n_elements(cut_inds) / total(diff_frac_hist) * 100, format = '(f7.1)') + '% of pixels'
           weights_ave[cut_inds] = 0.
        endif
        ;;undefine, diff, diff_frac_mask, cube1_hist, cube2_hist, diff_hist, wt_log_hist, df_log_hist
        undefine, weights_cube1, weights_cube2, diff_frac, diff_frac_hist, df_wt_hist, df_ri
     endif
        
     ;; now get 2 data cubes
     if healpix then begin
        data_cube1 = getvar_savefile(file_struct.uvf_savefile[0], 'data_cube') * float(conv_factor)
        data_cube2 = getvar_savefile(file_struct.uvf_savefile[1], 'data_cube') * float(conv_factor)
     endif else begin
        data_cube1 = getvar_savefile(file_struct.datafile[0], file_struct.datavar[i]) * float(conv_factor)
        data_cube2 = getvar_savefile(file_struct.datafile[1], file_struct.datavar[i]) * float(conv_factor)
     endelse

     if keyword_set(no_weighting) then weights_ave = real_part(data_cube1)*0. + 1.
  endif else begin
     ;; single file mode
     if healpix then begin
        data_cube1 = getvar_savefile(file_struct.uvf_savefile, 'data_cube') * float(conv_factor)
        weights_ave = getvar_savefile(file_struct.uvf_weight_savefile, 'weights_cube')
     endif else begin
        data_cube1 = getvar_savefile(file_struct.datafile, file_struct.datavar) * float(conv_factor)
        weights_ave = getvar_savefile(file_struct.weightfile, file_struct.weight_var)
     endelse

     if keyword_set(no_weighting) then weights_ave = real_part(data_cube1)*0. + 1.
  endelse


  if healpix then begin
     dims = size(data_cube1, /dimensions)
     n_kx = dims[0]
     kx_mpc = temporary(kx_rad_vals) / z_mpc_mean
     
     n_ky = dims[1]
     ky_mpc = temporary(ky_rad_vals) / z_mpc_mean
     
  endif else begin
     
     x_rad_delta = abs(file_struct.degpix) * !pi / 180d
     n_kx = dims[0]
     x_rad_length = dims[0] * x_rad_delta
     x_mpc_delta = x_rad_delta * z_mpc_mean
     x_mpc_length = x_rad_length * z_mpc_mean
     kx_mpc_range = (2.*!pi) / x_mpc_delta
     kx_mpc_delta = (2.*!pi) / x_mpc_length
     kx_mpc = (dindgen(n_kx)-n_kx/2) * kx_mpc_delta
     
     y_rad_delta = abs(file_struct.degpix) * !pi / 180.
     n_ky = dims[1]
     y_rad_length = dims[1] * y_rad_delta
     y_mpc_delta = y_rad_delta * z_mpc_mean
     y_mpc_length = y_rad_length * z_mpc_mean
     ky_mpc_range = (2.*!pi) / y_mpc_delta
     ky_mpc_delta = (2.*!pi) / y_mpc_length
     ky_mpc = (dindgen(n_ky)-n_ky/2) * ky_mpc_delta
     
  endelse
  
  sigma2_cube = 1d/(weights_ave)
  wh_wt0 = where(weights_ave eq 0, count_wt0)
  ;; wh = where(weights_cube le 1e-10, count)
  if count_wt0 ne 0 then sigma2_cube[wh_wt0] = 0
  undefine, weights_ave
  
  mask = intarr(n_kx, n_ky, n_kz) + 1
  if count_wt0 gt 0 then mask[wh_wt0] = 0
  ;; n_pix_contrib = total(total(mask, 2), 1)
  n_freq_contrib = total(mask, 3)
  wh_nofreq = where(n_freq_contrib eq 0, count_nofreq)
  undefine, mask
  
  print, 'pre-weighting sum(data_cube^2)*z_delta:', total(abs(data_cube1)^2d)*z_mpc_delta
  if nfiles eq 2 then print, 'pre-weighting sum(data_cube2^2)*z_delta:', total(abs(data_cube2)^2d)*z_mpc_delta

  ;; divide by weights (~array beam) to estimate true sky
  data_cube1 = data_cube1 * sigma2_cube
  if count_wt0 ne 0 then data_cube1[wh_wt0] = 0
    if nfiles eq 2 then begin
     data_cube2 = data_cube2 * sigma2_cube
     if count_wt0 ne 0 then data_cube2[wh_wt0] = 0
  endif

  print, 'sum(data_cube^2)*z_delta (after weighting):', total(abs(data_cube1)^2d)*z_mpc_delta
  if nfiles eq 2 then print, 'sum(data_cube2^2)*z_delta (after weighting):', total(abs(data_cube2)^2d)*z_mpc_delta
  
  ;; save some slices of the data cube
  for i=0, nfiles-1 do begin
     if i eq 0 then data_cube = data_cube1 else data_cube = data_cube2
     uf_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 1, $
                          slice_inds = n_ky/2, slice_savefile = file_struct.uf_savefile[i])
     
     vf_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 0, $
                          slice_inds = n_kx/2, slice_savefile = file_struct.vf_savefile[i])
     
     if max(abs(vf_slice)) eq 0 then begin
        nloop = 0
        while max(abs(vf_slice)) eq 0 do begin
           nloop = nloop+1
           vf_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, $
                                slice_axis = 0, slice_inds = n_kx/2+nloop, slice_savefile = file_struct.vf_savefile[i])
        endwhile
     endif
     
     uv_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 2, $
                          slice_inds = 0, slice_savefile = file_struct.uv_savefile[i])
  endfor     
  undefine, data_cube

  ;; need to cut uvf cubes in half because image is real -- we'll cut in v
  data_cube1 = data_cube1[*, n_ky/2:n_ky-1,*]
  if nfiles eq 2 then data_cube2 = data_cube2[*, n_ky/2:n_ky-1,*]
  sigma2_cube = sigma2_cube[*, n_ky/2:n_ky-1,*]
  n_freq_contrib = n_freq_contrib[*, n_ky/2:n_ky-1]
  
  ky_mpc = ky_mpc[n_ky/2:n_ky-1]
  n_ky = n_elements(ky_mpc)
     
  print, 'sum(data_cube^2)*z_delta (after cut):', total(abs(data_cube1)^2d)*z_mpc_delta
  if nfiles eq 2 then  print, 'sum(data_cube2^2)*z_delta (after cut):', total(abs(data_cube2)^2d)*z_mpc_delta

  ;; now take FFT
  data1_ft = fft(data_cube1, dimension=3) * n_freq * z_mpc_delta / (2.*!pi)
  ;; put k0 in middle of cube
  data1_ft = shift(data1_ft, [0,0,n_kz/2])
  undefine, data_cube1
  if nfiles eq 2 then begin
     data2_ft = fft(data_cube2, dimension=3) * n_freq * z_mpc_delta / (2.*!pi)
     ;; put k0 in middle of cube
     data2_ft = shift(data2_ft, [0,0,n_kz/2])
     undefine, data_cube2
  endif

  print, 'full_ft^2 integral:', total(abs(data1_ft)^2d)
  print, 'full_ft^2 integral * 2pi*delta_k^2:', total(abs(data1_ft)^2d) * kz_mpc_delta * 2. * !pi
  if nfiles eq 2 then begin
     print, 'full2_ft^2 integral :', total(abs(data2_ft)^2d)
     print, 'full2_ft^2 integral * 2pi*delta_k^2:', total(abs(data2_ft)^2d) * kz_mpc_delta * 2. * !pi
  endif

  ;; factor to go to eor theory FT convention
  ;; Only 1 factor of 2pi b/c only transforming along z
  ;; factor = (2d*!pi)
  ;; if not keyword_set(eor_test) then data_ft = factor * temporary(data_ft)
  
  ;; print, 'full_ft^2d integral (after theory factor):', total(abs(data_ft)^2d)
  
  n_val = round(kz_mpc_orig / kz_mpc_delta)
  kz_mpc_orig[where(n_val eq 0)] = 0
  a1_0 = 2. * data1_ft[*,*,where(n_val eq 0)]
  a1_n = data1_ft[*,*, where(n_val gt 0)] + data1_ft[*,*, reverse(where(n_val lt 0))]
  b1_n = complex(0,1) * (data1_ft[*,*, where(n_val gt 0)] - data1_ft[*,*, reverse(where(n_val lt 0))])
  undefine, data1_ft

  if nfiles gt 1 then begin
     a2_0 = 2. * data2_ft[*,*,where(n_val eq 0)]
     a2_n = data2_ft[*,*, where(n_val gt 0)] + data2_ft[*,*, reverse(where(n_val lt 0))]
     b2_n = complex(0,1) * (data2_ft[*,*, where(n_val gt 0)] - data2_ft[*,*, reverse(where(n_val lt 0))])
     undefine, data2_ft
  endif

  kz_mpc = kz_mpc_orig[where(n_val ge 0)]
  n_kz = n_elements(kz_mpc)
  
  if keyword_set(std_power) then begin
     ;; for standard power calc. just need ft of sigma2
     sigma2_ft = fft(sigma2_cube, dimension=3) * n_freq * z_mpc_delta / (2.*!pi)
     sigma2_ft = shift(sigma2_ft, [0,0,n_kz/2])
     undefine, sigma2_cube
     
     sigma_a0 = 2. * abs(sigma2_ft[*,*,where(n_val eq 0)])
     sigma_an_bn = sqrt(abs(sigma2_ft[*,*, where(n_val gt 0)])^2. + abs(sigma2_ft[*,*, reverse(where(n_val lt 0))])^2.)
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
        
        a1_0 = temporary(a_0) * mask_fewfreq[*,*,0]
        a1_n = temporary(a_n) * mask_fewfreq[*,*,1:*]
        b1_n = temporary(b_n) * mask_fewfreq[*,*,1:*]
        if nfiles gt 1 then begin
           a2_0 = temporary(a_0) * mask_fewfreq[*,*,0]
           a2_n = temporary(a_n) * mask_fewfreq[*,*,1:*]
           b2_n = temporary(b_n) * mask_fewfreq[*,*,1:*]
        endif
     endif
     
     data1_cos = complex(fltarr(n_kx, n_ky, n_kz))
     data1_sin = complex(fltarr(n_kx, n_ky, n_kz))
     data1_cos[*, *, 0] = a1_0 /2.
     data1_cos[*, *, 1:n_kz-1] = a1_n
     data1_sin[*, *, 1:n_kz-1] = b1_n
     if nfiles gt 1 then begin
        data2_cos = complex(fltarr(n_kx, n_ky, n_kz))
        data2_sin = complex(fltarr(n_kx, n_ky, n_kz))
        data2_cos[*, *, 0] = a2_0 /2.
        data2_cos[*, *, 1:n_kz-1] = a2_n
        data2_sin[*, *, 1:n_kz-1] = b2_n
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

     sigma2_cube = reform(sigma2_cube, n_kx*n_ky, n_freq)
     covar_cos = matrix_multiply(sigma2_cube, cos_arr^2d)
     covar_sin = matrix_multiply(sigma2_cube, sin_arr^2d)
     covar_cross = matrix_multiply(sigma2_cube, cos_arr*sin_arr)

     wh_divide = where(n_freq_contrib gt 0, count_divide, complement = wh_0f, ncomplement = count_0f)
     n_divide = double(rebin(reform(n_freq_contrib[wh_divide]),count_divide, n_kz))

     if count_divide gt 0 then begin
        covar_cos[wh_divide, *] = covar_cos[wh_divide, *] / n_divide
        covar_sin[wh_divide, *] = covar_sin[wh_divide, *] / n_divide
        covar_cross[wh_divide, *] = covar_cross[wh_divide, *] / n_divide
     endif

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
     undefine, sigma2_cube, freq_kz_arr, cos_arr, sin_arr, sigma2_arr
     
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
     
     sigma1_2 = covar_cos*cos_theta^2. + 2.*covar_cross*cos_theta*sin_theta + covar_sin*sin_theta^2.
     sigma2_2 = covar_cos*sin_theta^2. - 2.*covar_cross*cos_theta*sin_theta + covar_sin*cos_theta^2.
     
     undefine, covar_cos, covar_sin, covar_cross
     
     data1_1 = data1_cos*cos_theta + data1_sin*sin_theta
     data1_2 = (-1d)*data1_cos*sin_theta + data1_sin*cos_theta
     undefine, data1_cos, data1_sin
     if nfiles eq 2 then begin
        data2_1 = data2_cos*cos_theta + data2_sin*sin_theta
        data2_2 = (-1d)*data2_cos*sin_theta + data2_sin*cos_theta
        undefine, data2_cos, data2_sin
     endif

     save, file = file_struct.kcube_savefile, data1_1, data1_2, data2_1, data2_2, sigma1_2, sigma2_2, $
           kx_mpc, ky_mpc, kz_mpc, kperp_lambda_conv, delay_params, hubble_param, n_freq_contrib
  endelse
  
end
