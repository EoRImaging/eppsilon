pro fhd_kcube, file_struct, dft_refresh_data = dft_refresh_data, dft_refresh_weight = dft_refresh_weight, $
              dft_fchunk = dft_fchunk, std_power = std_power, input_units = input_units, quiet = quiet

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

     
     variancefile_obj = obj_new('IDL_Savefile', file_struct.variancefile[i])
     variancefile_names = variancefile_obj->names()
     variancevar = strupcase(file_struct.variancevar)
     wh = where(variancefile_names eq variancevar[i], count)
     if count eq 0 then begin
        print, 'specified variancevar is not present in variancefile (variancefile=' + file_struct.variancefile[i] $
               +  ', variancevar=' + file_struct.variancevar[i] + '). Weights will be used instead' 
        no_var = 1
     endif else begin
        if n_elements(no_var) eq 0 then no_var = 0

        variance_dims = variancefile_obj->size(variancevar[i], /dimensions)
        if total(abs(data_dims - variance_dims)) ne 0 then message, 'data and variance dimensions do not match'
        undefine, variance_dims
     endelse
     obj_destroy, variancefile_obj
     
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
     conv_factor = float(2. * max_baseline^2d / (!dpi * 1.38065))
  endif else conv_factor = 1.
  
  t_sys = 440. ; K
  eff_area = 16. ; m^2
  df = 40.e3 ; Hz
  tau = 1. ; seconds
  vis_sigma = (2. * (1.38065e-23) * 1e26) * t_sys / (eff_area * sqrt(df * tau)) ;; in Jy
  vis_sigma = vis_sigma * conv_factor ;; convert to mK

  if healpix then begin
     if nfiles eq 2 then begin
        ;; check that they have the same set of healpix pixels
        pixel_nums1 = getvar_savefile(file_struct.pixelfile[0], file_struct.pixelvar[0])
        pixel_nums2 = getvar_savefile(file_struct.pixelfile[1], file_struct.pixelvar[1])
        if n_elements(pixel_nums1) ne n_elements(pixel_nums2) then message, 'Different number of Healpix pixels in cubes'
        
        if total(abs(pixel_nums1-pixel_nums2)) ne 0 then message, 'Pixel numbers are not consistent between cubes'
  
        pixel_dims = size(pixel_nums1, /dimension)
        if total(abs(dims - pixel_dims)) ne 0 then message, 'pixel and data dimensions do not match'

     endif

     for i=0, nfiles-1 do begin     
        test_uvf = file_test(file_struct.uvf_savefile[i]) *  (1 - file_test(file_struct.uvf_savefile[i], /zero_length))
     
        test_wt_uvf = file_test(file_struct.uvf_weight_savefile[i]) * (1 - file_test(file_struct.uvf_weight_savefile[i], /zero_length))
        
        if test_uvf eq 0 or test_wt_uvf eq 0 or keyword_set(dft_refresh_data) or keyword_set(dft_refresh_weight) then begin          
           test_setup = file_test(file_struct.hpx_dftsetup_savefile) * $
                        (1 - file_test(file_struct.hpx_dftsetup_savefile, /zero_length))
           if test_setup eq 0 then begin
              ;; figure out k values to calculate dft
              healpix_setup_ft, pixel_nums1, file_struct.nside, new_pix_vec, limits, kx_rad_vals, ky_rad_vals, /quiet
              save, file = file_struct.hpx_dftsetup_savefile, new_pix_vec, limits, kx_rad_vals, ky_rad_vals
           endif else restore, file_struct.hpx_dftsetup_savefile
           
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
              weights_cube = temporary(transform)
              
              if not no_var then begin
                 arr = getvar_savefile(file_struct.variancefile[i], file_struct.variancevar[i])
                 transform = discrete_ft_2D_fast(new_pix_vec[*,0], new_pix_vec[*,1], arr, kx_rad_vals, ky_rad_vals, timing = ft_time, $
                                                 fchunk = dft_fchunk)            
                 variance_cube = abs(temporary(transform)) ;; make variances real, positive definite (amplitude)
                 undefine, arr
              endif
 
              save, file = file_struct.uvf_weight_savefile[i], kx_rad_vals, ky_rad_vals, weights_cube, variance_cube
              undefine, new_pix_vec, weights_cube, variance_cube
           endif
           
        endif else begin
           kx_rad_vals = getvar_savefile(file_struct.uvf_savefile[0], 'kx_rad_vals')
           ky_rad_vals = getvar_savefile(file_struct.uvf_savefile[0], 'ky_rad_vals')
        endelse
     endfor
  endif


  ;; if 2 file mode check that uvf variances are close enough to proceed.
  if nfiles eq 2 then begin
     if no_var then begin
        ;; use 1/abs(weights) instead
        if healpix then begin
           weights_cube1 = getvar_savefile(file_struct.uvf_weight_savefile[0], 'weights_cube')
           weights_cube2 = getvar_savefile(file_struct.uvf_weight_savefile[1], 'weights_cube')
           sigma2_cube1 = 1./abs(weights_cube1)
           sigma2_cube2 = 1./abs(weights_cube2)
        endif else begin
           weights_cube1 = getvar_savefile(file_struct.weightfile[0], file_struct.weightvar[0])
           weights_cube2 = getvar_savefile(file_struct.weightfile[1], file_struct.weightvar[1])
           sigma2_cube1 = 1./abs(weights_cube1)
           sigma2_cube2 = 1./abs(weights_cube2)
        endelse
     endif else begin
        if healpix then begin
           variance_cube1 = getvar_savefile(file_struct.uvf_weight_savefile[0], 'variance_cube')
           variance_cube2 = getvar_savefile(file_struct.uvf_weight_savefile[1], 'variance_cube')
           weights_cube1 = getvar_savefile(file_struct.uvf_weight_savefile[0], 'weights_cube')
           weights_cube2 = getvar_savefile(file_struct.uvf_weight_savefile[1], 'weights_cube')
        endif else begin
           variance_cube1 = getvar_savefile(file_struct.variancefile[0], file_struct.variancevar[i])
           variance_cube2 = getvar_savefile(file_struct.variancefile[1], file_struct.variancevar[i])
           weights_cube1 = getvar_savefile(file_struct.weightfile[0], file_struct.weightvar[0])
           weights_cube2 = getvar_savefile(file_struct.weightfile[1], file_struct.weightvar[1])
       endelse

        sigma2_cube1 = temporary(variance_cube1) / abs(weights_cube1)^2.
        sigma2_cube2 = temporary(variance_cube2) / abs(weights_cube2)^2.
     endelse

     wh_wt1_0 = where(abs(weights_cube1) eq 0, count_wt1_0)
     if count_wt1_0 ne 0 then sigma2_cube1[wh_wt1_0] = 0
     wh_wt2_0 = where(abs(weights_cube2) eq 0, count_wt2_0)
     if count_wt2_0 ne 0 then sigma2_cube2[wh_wt2_0] = 0
 

    if min(sigma2_cube1) lt 0 or min(sigma2_cube2) lt 0 then message, 'sigma2 should be positive definite.'
     if total(abs(sigma2_cube1)) le 0 or total(abs(sigma2_cube2)) le 0 then message, 'one or both sigma2 cubes is all zero'
     
     if not keyword_set(quiet) then begin
        freq_channel = 0
     
        win_num=1
        if windowavailable(win_num) then wset, win_num else window, win_num
        quick_image, sigma2_cube1[*,*,freq_channel], kx_rad_vals, ky_rad_vals, /log, title = 'Even cube variance', $
                     xtitle = 'kx (Mpc!U-1!N)', ytitle = 'ky (Mpc!U-1!N)', data_range = data_range
        win_num=2
        if windowavailable(win_num) then wset, win_num else window, win_num
        quick_image, sigma2_cube2[*,*,freq_channel], kx_rad_vals, ky_rad_vals, /log, title = 'Odd cube variance', $
                     xtitle = 'kx (Mpc!U-1!N)', ytitle = 'ky (Mpc!U-1!N)', data_range = data_range
     end
     
     sigma2_ave = (sigma2_cube1 + sigma2_cube2) / 2.
     diff = abs(1./sigma2_cube1 - 1./sigma2_cube2)
     wh_invsig0 = where(sigma2_cube1 eq 0 or sigma2_cube2 eq 0, count_invsig0)
     if count_invsig0 gt 0 then diff[wh_invsig0] = 0

     if not keyword_set(quiet) then begin
        win_num=3
        if windowavailable(win_num) then wset, win_num else window, win_num
        quick_image, diff[*,*,freq_channel], kx_rad_vals, ky_rad_vals, /log, title = 'Inverse Variance Difference', $
                     xtitle = 'kx (Mpc!U-1!N)', ytitle = 'ky (Mpc!U-1!N)'
     endif

     ave_inv_var = (1./sigma2_cube1 + 1./sigma2_cube2) /2.
     if count_invsig0 gt 0 then ave_inv_var[wh_invsig0] = 0
     diff_frac = diff / ave_inv_var
     if count_invsig0 gt 0 then diff_frac[wh_invsig0] = 0
    

     if not keyword_set(quiet) then begin
        win_num=4 
       if windowavailable(win_num) then wset, win_num else window, win_num
        quick_image, diff_frac[*,*,freq_channel], kx_rad_vals, ky_rad_vals, /log, title = 'Inverse Variance Difference Fraction', $
                     xtitle = 'kx (Mpc!U-1!N)', ytitle = 'ky (Mpc!U-1!N)'
        
        ;;diff_frac_mask = fix(diff_frac*0)
        ;;wh_large = where(diff_frac ge 1, count_large)
        ;;if count_large gt 0 then diff_frac_mask[wh_large] = 1
        ;;quick_image, diff_frac_mask[*,*,freq_channel], kx_rad_vals, ky_rad_vals, title = 'Difference Fraction above 1', $
        ;;             xtitle = 'kx (Mpc!U-1!N)', ytitle = 'ky (Mpc!U-1!N)', data_range = [0,1], /grey_scale
     endif

     cube_max = max([sigma2_cube1, sigma2_cube2])
     cube_min_n0 = min([sigma2_cube1[where(sigma2_cube1 ne 0)], sigma2_cube2[where(sigma2_cube2 ne 0)]])
     diff_frac_hist = histogram(diff_frac, binsize = 0.01, min = 0, omax = df_max, locations = diff_frac_locs, $
                                reverse_indices = df_ri)
     
     min_log_sigma = alog10(cube_min_n0)
     min_log_df = alog10(min(diff_frac[where(diff_frac gt 0)]))
     invvar_log_hist = histogram(alog10(ave_inv_var), binsize = 0.1, min = min_log_invvar, omax = max_log_invvar, $
                                locations = invvar_log_locs)
     df_log_hist = histogram(alog10(diff_frac), binsize = 0.1, min = min_log_df, omax = max_log_df, locations = df_log_locs)
     
     df_invvar_hist = hist_2d(alog10(ave_inv_var), alog10(diff_frac), bin1=0.1, bin2=0.1, $
                          min1 = min_log_invvar, min2=min_log_df, max1=max_log_invvar, max2=max_log_df)
     
     if not keyword_set(quiet) then begin
        win_num=5
        if windowavailable(win_num) then wset, win_num else window, win_num
        quick_image, df_invvar_hist, 10^invvar_log_locs, 10^df_log_locs,/log, /xlog, /ylog, xtitle = 'average inverse variance', $
                     ytitle = 'inverse variance difference fraction', title = 'inverse variance vs difference fraction histogram'
        
        ;;binarea = matrix_multiply(10^(sigma_log_locs + 0.1) - 10^(sigma_log_locs), 10^(df_log_locs + 0.1) - 10^(df_log_locs))
        ;;quick_image, df_wt_hist/binarea, 10^wt_log_locs, 10^df_log_locs,/log, /xlog, /ylog, xtitle = 'average sigma2', $
        ;;             ytitle = 'difference fraction', title = 'sigma2 vs difference fraction density'
        
     endif

     df_cut_level = 0.1
     if not keyword_set(quiet) then cgplot, /overplot, 10^invvar_log_locs, invvar_log_locs*0+df_cut_level

     wh_df_cut = where(diff_frac_locs gt df_cut_level, count_df_cut)
     if count_df_cut gt 0 then begin
        print, 'cutting out pixels with an inverse variance difference greater than ' + number_formatter(df_cut_level*100.) + '%'
        
        cut_inds = df_ri[df_ri[min(wh_df_cut)]:df_ri[max(wh_df_cut)+1]-1]
        if (n_elements(cut_inds) / total(diff_frac_hist)) gt .1 then stop
        
        print, 'removed ' + number_formatter(n_elements(cut_inds) / total(diff_frac_hist) * 100, format = '(f7.1)') + '% of pixels'
        sigma2_ave[cut_inds] = 0.
     endif
     ;;undefine, cube1_hist, cube2_hist, diff_frac_mask
     undefine, sigma2_cube1, sigma2_cube2
     undefine, diff, ave_inv_var, diff_frac, diff_frac_hist, invvar_log_hist, df_log_hist, df_invvar_hist, df_ri

     if not keyword_set(quiet) then stop
     
     ;; now get data cubes
     if healpix then begin
        data_cube1 = getvar_savefile(file_struct.uvf_savefile[0], 'data_cube') * float(conv_factor)
        data_cube2 = getvar_savefile(file_struct.uvf_savefile[1], 'data_cube') * float(conv_factor)
     endif else begin
        data_cube1 = getvar_savefile(file_struct.datafile[0], file_struct.datavar[0]) * float(conv_factor)
        data_cube2 = getvar_savefile(file_struct.datafile[1], file_struct.datavar[1]) * float(conv_factor)
     endelse
     
  endif else begin
     ;; single file mode
     if healpix then begin
        data_cube1 = getvar_savefile(file_struct.uvf_savefile, 'data_cube') * float(conv_factor)
        weights_cube1 = getvar_savefile(file_struct.uvf_weight_savefile, 'weights_cube')
        if not no_var then begin
           variance_cube1 = getvar_savefile(file_struct.uvf_weight_savefile, 'variance_cube') 
           sigma2_ave = temporary(variance_cube1) / abs(weights_cube1)^2.
        endif else sigma2_ave = 1./abs(weights_cube1)

        wh_wt1_0 = where(abs(weights_cube1) eq 0, count_wt1_0)
        if count_wt1_0 ne 0 then sigma2_ave[wh_wt1_0] = 0
     endif else begin
        data_cube1 = getvar_savefile(file_struct.datafile, file_struct.datavar) * float(conv_factor)
        weights_cube1 = getvar_savefile(file_struct.weightfile, file_struct.weightvar)
        if not no_var then begin
           variance_ave = getvar_savefile(file_struct.variancefile, file_struct.variance_var) 
           sigma2_ave = temporary(variance_cube1) / abs(weights_cube1)^2.
        endif else sigma2_ave = 1./abs(weights_cube1)

        wh_wt1_0 = where(abs(weights_cube1) eq 0, count_wt1_0)
        if count_wt1_0 ne 0 then sigma2_ave[wh_wt1_0] = 0
     endelse
  endelse
  
  ;; get sigma into mK
  sigma2_ave = sigma2_ave * vis_sigma^2.

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

  mask = intarr(n_kx, n_ky, n_kz) + 1
  wh_sig0 = where(sigma2_ave eq 0, count_sig0)
  if count_sig0 gt 0 then mask[wh_sig0] = 0
  ;; n_pix_contrib = total(total(mask, 2), 1)
  n_freq_contrib = total(mask, 3)
  wh_nofreq = where(n_freq_contrib eq 0, count_nofreq)
  undefine, mask
  
  print, 'pre-weighting sum(data_cube^2)*z_delta:', total(abs(data_cube1)^2d)*z_mpc_delta
  if nfiles eq 2 then print, 'pre-weighting sum(data_cube2^2)*z_delta:', total(abs(data_cube2)^2d)*z_mpc_delta

  ;; divide by weights (~array beam) to estimate true sky
  data_cube1 = data_cube1 / weights_cube1
  wh_wt1_0 = where(weights_cube1 eq 0, count_wt1_0)
  if count_wt1_0 ne 0 then data_cube1[wh_wt1_0] = 0
  undefine, weights_cube1

  if nfiles eq 2 then begin
     data_cube2 = data_cube2 / weights_cube2
     wh_wt2_0 = where(weights_cube2 eq 0, count_wt2_0)
     if count_wt2_0 ne 0 then data_cube2[wh_wt2_0] = 0
     undefine, weights_cube2
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
  sigma2_ave = sigma2_ave[*, n_ky/2:n_ky-1,*]
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
     sigma2_ft = fft(sigma2_ave, dimension=3) * n_freq * z_mpc_delta / (2.*!pi)
     sigma2_ft = shift(sigma2_ft, [0,0,n_kz/2])
     undefine, sigma2_ave
     
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

     sigma2_ave = reform(sigma2_ave, n_kx*n_ky, n_freq)
     covar_cos = matrix_multiply(sigma2_ave, cos_arr^2d)
     covar_sin = matrix_multiply(sigma2_ave, sin_arr^2d)
     covar_cross = matrix_multiply(sigma2_ave, cos_arr*sin_arr)

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
     undefine, sigma2_ave, freq_kz_arr, cos_arr, sin_arr
     
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
