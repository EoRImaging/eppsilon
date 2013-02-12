pro fhd_3dps, datafile, datavar, weightfile, weightvar, frequencies, max_baseline, degpix = degpix, healpix=healpix, $
              nside = nside, pixelfile = pixelfile, pixelvar = pixelvar, hpx_dftsetup_savefile = hpx_dftsetup_savefile, $
              savefilebase = savefilebase_in, weight_savefilebase = weight_savefilebase_in, refresh = refresh, $
              dft_refresh_data = dft_refresh_data, dft_refresh_weight = dft_refresh_weight, dft_fchunk = dft_fchunk, $
              no_weighting = no_weighting, std_power = std_power, no_kzero = no_kzero, log_kpar = log_kpar, $
              log_kperp = log_kperp, kperp_bin = kperp_bin, kpar_bin = kpar_bin, log_k1d = log_k1d, k1d_bin = k1d_bin, $
              no_weighted_averaging = no_weighted_averaging, input_units = input_units, fill_holes = fill_holes, quiet = quiet
;, clean_type = clean_type

  if n_params() ne 6 then message, 'Wrong number of parameters passed. Required parameters are: datafile (string), ' + $
                                 'datavar (string), weightfile (string), weightvar (string), ' + $
                                 'freqencies (float or double array), max_baseline (float or double)'

  if keyword_set(healpix) and (keyword_set(dft_refresh_data) or keyword_set(dft_refresh_weight)) then refresh=1

  datafile_test = file_test(datafile)
  if datafile_test eq 0 then message, 'datafile not found'

  if not keyword_set(no_weighting) then begin
     weightfile_test = file_test(weightfile)
     if weightfile_test eq 0 then message, 'weightfile not found'
  endif     

  if n_elements(input_units) eq 0 then input_units = 'jansky'
  units_enum = ['jansky', 'mk']
  wh = where(units_enum eq input_units, count)
  if count eq 0 then message, 'input units not recognized, options are: ' + units_enum

  if n_elements(fill_holes) eq 0 then fill_holes = 0

  ;; clean_type_enum = ['hmf', 'iterate', 'fit']
  ;; if n_elements(clean_type) ne 0 then begin
  ;;    wh = where(clean_type_enum eq clean_type, count)
  ;;    if count eq 0 then message, 'Clean type not recognized'
  ;; endif

  if n_elements(savefilebase_in) eq 0 then begin
     froot = file_dirname(datafile, /mark_directory)
     infilebase = file_basename(datafile)
     temp2 = strpos(infilebase, '.', /reverse_search)
     savefilebase = strmid(infilebase, 0, temp2)
  endif else begin
     temp = file_dirname(savefilebase_in, /mark_directory)
     if temp ne '.' then froot = temp else froot = file_dirname(datafile, /mark_directory)
     savefilebase = file_basename(savefilebase_in)
  endelse

  save_file = froot + savefilebase + '_power.idlsave'
  fits_savefile = froot + savefilebase + '_power.fits'

  ;;if n_elements(clean_type) ne 0 then if clean_type ne 'fit' then model_save = froot + savefilebase + '_modeluv.idlsave'

  test_save = file_test(save_file) *  (1 - file_test(save_file, /zero_length))
  if test_save eq 0 or keyword_set(refresh) then begin
     
     datafile_obj = obj_new('IDL_Savefile', datafile)
     datafile_names = datafile_obj->names()
     datavar = strupcase(datavar)
     wh = where(datafile_names eq datavar, count)
     if count eq 0 then message, 'specified datavar is not present in datafile (datafile=' + datafile + ', datavar=' + datavar + ')'
     data_dims = datafile_obj->size(datavar, /dimensions)
     obj_destroy, datafile_obj
     
     if not keyword_set(no_weighting) then begin
        weightfile_obj = obj_new('IDL_Savefile', weightfile)
        weightfile_names = weightfile_obj->names()
        weightvar = strupcase(weightvar)
        wh = where(weightfile_names eq weightvar, count)
        if count eq 0 then message, 'specified weightvar is not present in weightfile (weightfile=' + datafile + ', weightvar=' + $
                                    weightvar + ')'
        weight_dims = weightfile_obj->size(weightvar, /dimensions)
        obj_destroy, weightfile_obj
        
        if total(abs(data_dims - weight_dims)) ne 0 then message, 'data and weight dimensions do not match'
        undefine, weight_dims
     endif
     dims = data_dims
     undefine, data_dims
     
     if keyword_set(healpix) then n_freq = dims[1] else n_freq = dims[2]
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
        nominal_freqs = dindgen(floor(((max(frequencies)-min(frequencies))/freq_resolution))+1)*freq_resolution + min(frequencies)
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
     delay_max = delay_delta * n_freq/2d ;; factor of 2 b/c of neg/positive
     delay_params = [delay_delta, delay_max]     

     z_mpc_length = max(comov_dist_los) - min(comov_dist_los) + z_mpc_delta
     kz_mpc_range =  (2d*!dpi) / (z_mpc_delta)
     kz_mpc_delta = (2d*!dpi) / z_mpc_length
     kz_mpc = dindgen(round(kz_mpc_range / kz_mpc_delta)) * kz_mpc_delta - kz_mpc_range/2d
     if n_elements(n_kz) ne 0 then begin
        if n_elements(kz_mpc) ne n_kz then stop
     endif else n_kz = n_elements(kz_mpc)
    
     print, 'z delta: ', z_mpc_delta
     print,  'kz delta: ', kz_mpc_delta

     if input_units eq 'jansky' then begin
        ;; beam_diameter_rad = (3d * 10^8d) / (frequencies * 10^6d * max_baseline)
        ;; beam_area_str = !pi * beam_diameter_rad^2d /4d
        
        ;; conv_factor = (10^(double(-26+16+3-12+23)) * 9d) / (beam_area_str * 2d * frequencies^2d * 1.38)
        ;; if max(conv_factor-conv_factor[0]) gt 1e-8 then stop else conv_factor = conv_factor[0]
        conv_factor = 2d * max_baseline^2d / (!dpi * 1.38)
     endif else conv_factor = 1d

     ;; data and weights are floats, avoid increasing to doubles
     data_cube = getvar_savefile(datafile, datavar) * float(conv_factor)
     if not keyword_set(no_weighting) then begin
        weights_cube = getvar_savefile(weightfile, weightvar)
        if max(abs(imaginary(weights_cube))) eq 0 then weights_cube = real_part(weights_cube) $
        else stop
     endif else weights_cube = real_part(data_cube)*0. + 1.
        
     if keyword_set(healpix) then begin
        if n_elements(pixelfile) eq 0 or n_elements(pixelvar) eq 0 then message, $
           'pixelfile and pixelvar keywords must be specified for healpix'

        if n_elements(nside) eq 0 then message, 'nside keyword must be specified for healpix'

        pixelfile_test = file_test(pixelfile)
        if pixelfile_test eq 0 then message, 'pixelfile not found'

        uvf_data_savefile = froot + savefilebase + '_uvf.idlsave'
        test_uvf = file_test(uvf_data_savefile) *  (1 - file_test(uvf_data_savefile, /zero_length))
      
        if not keyword_set(no_weighting) then begin
           if n_elements(weight_savefilebase_in) eq 0 then begin
              weight_savefilebase = savefilebase + '_weights'
              wt_froot = froot
           endif else begin
              temp = file_dirname(weight_savefilebase_in, /mark_directory)
              if temp ne '.' then wt_froot = froot else wt_froot = temp
              weight_savefilebase = file_basename(weight_savefilebase_in)
           endelse

           uvf_weight_savefile = wt_froot + weight_savefilebase + '_uvf.idlsave'
           test_wt_uvf = file_test(uvf_weight_savefile) *  (1 - file_test(uvf_weight_savefile, /zero_length))
        endif        

        if test_uvf eq 0 or test_wt_uvf eq 0 or keyword_set(dft_refresh_data) or keyword_set(dft_refresh_weight) then begin 
           pixelfile_obj = obj_new('IDL_Savefile', pixelfile)
           pixelfile_names = pixelfile_obj->names()
           pixelvar = strupcase(pixelvar)
           wh = where(pixelfile_names eq pixelvar, count)
           if count eq 0 then message, 'specified pixelvar is not present in pixelfile (pixelfile=' + pixelfile + ', pixelvar=' + $
                                       pixelvar + ')'
           pixel_dims = pixelfile_obj->size(pixelvar, /dimensions)
           if total(abs(dims - pixel_dims)) ne 0 then message, 'pixel and data dimensions do not match'
           undefine, pixel_dims
           obj_destroy, pixelfile_obj

           pixel_nums = getvar_savefile(pixelfile, pixelvar)

           if n_elements(hpx_dftsetup_savefile) eq 0 then hpx_dftsetup_savefile = froot + savefilebase + '_dftsetup.idlsave'
           test_setup = file_test(hpx_dftsetup_savefile) *  (1 - file_test(hpx_dftsetup_savefile, /zero_length))
           if test_setup eq 0 then begin
              ;; figure out k values to calculate dft
              healpix_setup_ft, pixel_nums, nside, new_pix_vec, limits, kx_rad_vals, ky_rad_vals, /quiet
              save, file = hpx_dftsetup_savefile, new_pix_vec, limits, kx_rad_vals, ky_rad_vals
           endif else restore, hpx_dftsetup_savefile

           ;; do FT.
           if test_uvf eq 0 or keyword_set(dft_refresh_data) then begin
              arr = temporary(data_cube)
              transform = discrete_ft_2D_fast(new_pix_vec[*,0], new_pix_vec[*,1], arr, kx_rad_vals, ky_rad_vals, timing = ft_time, $
                                             fchunk = dft_fchunk)
              data_cube = transform

              save, file=uvf_data_savefile, kx_rad_vals, ky_rad_vals, data_cube
           endif else restore, uvf_data_savefile

           if test_wt_uvf eq 0 or keyword_set(dft_refresh_weight) then begin
              arr = temporary(weights_cube)
              transform = discrete_ft_2D_fast(new_pix_vec[*,0], new_pix_vec[*,1], arr, kx_rad_vals, ky_rad_vals, timing = ft_time, $
                                             fchunk = dft_fchunk)

              weights_cube = abs(transform) ;; make weights real (amplitude)
 
              save, file=uvf_weight_savefile, kx_rad_vals, ky_rad_vals, weights_cube
           endif else restore, uvf_weight_savefile
           undefine, transform, new_pix_vec, arr
         
        endif else begin
           restore, uvf_weight_savefile
           restore, uvf_data_savefile
        endelse

        dims = size(data_cube, /dimensions)
        n_kx = dims[0]
        kx_mpc = temporary(kx_rad_vals) / z_mpc_mean
        
        n_ky = dims[1]
        ky_mpc = temporary(ky_rad_vals) / z_mpc_mean


     endif else begin
        if n_elements(degpix) eq 0 then message, 'degpix keyword must be specified unless in healpix'

        x_rad_delta = abs(degpix) * !pi / 180d
        n_kx = dims[0]
        x_rad_length = dims[0] * x_rad_delta
        x_mpc_delta = x_rad_delta * z_mpc_mean
        x_mpc_length = x_rad_length * z_mpc_mean
        kx_mpc_range = (2d*!pi) / x_mpc_delta
        kx_mpc_delta = (2d*!pi) / x_mpc_length
        kx_mpc = (dindgen(n_kx)-n_kx/2) * kx_mpc_delta
        
        y_rad_delta = abs(degpix) * !pi / 180d
        n_ky = dims[1]
        y_rad_length = dims[1] * y_rad_delta
        y_mpc_delta = y_rad_delta * z_mpc_mean
        y_mpc_length = y_rad_length * z_mpc_mean
        ky_mpc_range = (2d*!pi) / y_mpc_delta
        ky_mpc_delta = (2d*!pi) / y_mpc_length
        ky_mpc = (dindgen(n_ky)-n_ky/2) * ky_mpc_delta
        
     endelse

     sigma2_cube = 1d/(weights_cube)
     wh_wt0 = where(weights_cube eq 0, count_wt0)
     ;; wh = where(weights_cube le 1e-10, count)
     if count_wt0 ne 0 then sigma2_cube[wh_wt0] = 0
     undefine, weights_cube

     mask = intarr(n_kx, n_ky, n_kz) + 1
     if count_wt0 gt 0 then mask[wh_wt0] = 0
     ;; n_pix_contrib = total(total(mask, 2), 1)
     n_freq_contrib = total(mask, 3)
     wh_nofreq = where(n_freq_contrib eq 0, count_nofreq)
     undefine, mask

      print, 'pre-weighting sum(data_cube^2)*z_delta:', total(abs(data_cube)^2d)*z_mpc_delta

     ;; divide by weights (~array beam) to estimate true sky
     data_cube = data_cube * sigma2_cube
     if count_wt0 ne 0 then data_cube[wh_wt0] = 0
     
     print, 'sum(data_cube^2)*z_delta (after weighting):', total(abs(data_cube)^2d)*z_mpc_delta

     ;; save some slices of the data cube
     uf_savefile = froot + savefilebase + '_uf_plane.idlsave'
     uf_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 1, $
                          slice_inds = n_ky/2, slice_savefile = uf_savefile)

     vf_savefile = froot + savefilebase + '_vf_plane.idlsave'
     vf_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 0, $
                          slice_inds = n_kx/2, slice_savefile = vf_savefile)

     if max(abs(vf_slice)) eq 0 then begin
        nloop = 0
        while max(abs(vf_slice)) eq 0 do begin
           nloop = nloop+1
           vf_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, $
                                slice_axis = 0, slice_inds = n_kx/2+nloop, slice_savefile = vf_savefile)
        endwhile
     endif

     uv_savefile = froot + savefilebase + '_uv_plane.idlsave'
     uv_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 2, $
                          slice_inds = 0, slice_savefile = uv_savefile)

     ;; need to cut uvf cubes in half because image is real -- we'll cut in v
     data_cube = data_cube[*, n_ky/2:n_ky-1,*]
     sigma2_cube = sigma2_cube[*, n_ky/2:n_ky-1,*]
     n_freq_contrib = n_freq_contrib[*, n_ky/2:n_ky-1]
 
     ky_mpc = ky_mpc[n_ky/2:n_ky-1]
     n_ky = n_elements(ky_mpc)

     print, 'sum(data_cube^2)*z_delta (after cut):', total(abs(data_cube)^2d)*z_mpc_delta

     ;; now take FFT
     data_ft = fft(data_cube, dimension=3) * n_freq * z_mpc_delta / (2d*!dpi)
     ;; put k0 in middle of cube
     data_ft = shift(data_ft, [0,0,n_kz/2])
     undefine, data_cube

     print, 'full_ft^2 integral:', total(abs(data_ft)^2d)
     print, 'full_ft^2 integral * 2pi*delta_k^2:', total(abs(data_ft)^2d) * kz_mpc_delta * 2d * !dpi
 
     ;; factor to go to eor theory FT convention
     ;; Only 1 factor of 2pi b/c only transforming along z
     ;; factor = (2d*!pi)
     ;; if not keyword_set(eor_test) then data_ft = factor * temporary(data_ft)

     ;; print, 'full_ft^2d integral (after theory factor):', total(abs(data_ft)^2d)

     n_val = kz_mpc / kz_mpc_delta
     a_0 = 2d * data_ft[*,*,where(n_val eq 0)]
     a_n = data_ft[*,*, where(n_val gt 0)] + data_ft[*,*, reverse(where(n_val lt 0))]
     b_n = complex(0,1) * (data_ft[*,*, where(n_val gt 0)] - data_ft[*,*, reverse(where(n_val lt 0))])
     undefine, data_ft
   
     kz_mpc = kz_mpc[where(n_val ge 0)]
     n_kz = n_elements(kz_mpc)

     if keyword_set(std_power) then begin
        ;; for standard power calc. just need ft of sigma2
        sigma2_ft = fft(sigma2_cube, dimension=3) * n_freq * z_mpc_delta / (2d*!dpi)
        sigma2_ft = shift(sigma2_ft, [0,0,n_kz/2])

        undefine, sigma2_cube

        ;; factor to go to eor theory FT convention
        ;;sigma2_ft = factor * temporary(sigma2_ft)

        power_3d = dblarr(n_kx, n_ky, n_kz)
        power_3d[*,*,0] = (a_0 * conj(a_0))/4d
        power_3d[*,*,1:n_kz-1] = ((a_n * conj(a_n)) + (b_n * conj(b_n)))/2d
    
        ;; power_3d[*,*,0] = real_part(data_ft[*,*,where(n_val eq 0)] * conj(data_ft[*,*,where(n_val eq 0)]))
        ;; power_3d[*,*,where(n_val gt 0)] = real_part(data_ft[*,*,where(n_val gt 0)] * conj(data_ft[*,*,where(n_val gt 0)])) + $
        ;;                   reverse(real_part(data_ft[*,*,where(n_val lt 0)] * conj(data_ft[*,*,where(n_val lt 0)])), 3)
        ;;undefine, data_ft

        sigma_a0 = 2d * abs(sigma2_ft[*,*,where(n_val eq 0)])
        sigma_an_bn = sqrt(abs(sigma2_ft[*,*, where(n_val gt 0)])^2d + abs(sigma2_ft[*,*, reverse(where(n_val lt 0))])^2d)

        sigma2_3d = dblarr(n_kx, n_ky, n_kz)
        sigma2_3d[*,*,0] = abs(a_0) * sigma_a0 / 2d
        sigma2_3d[*,*,1:n_kz-1] = sqrt(abs(a_n)^2d + abs(b_n)^2d) * sigma_an_bn
        undefine, sigma2_ft

        weights_3d = 1d/sigma2_3d
        wh_sig0 = where(sigma2_3d eq 0, count_sig0)
        if count_sig0 gt 0 then weights_3d[wh_sig0] = 0
        sigma2_3d=0

     endif else begin  
 
        ;; drop pixels with less than 1/3 of the frequencies (set weights to 0)
        wh_fewfreq = where(n_freq_contrib lt ceil(n_freq/3d), count_fewfreq)
        if count_fewfreq gt 0 then begin
           mask_fewfreq = n_freq_contrib * 0 + 1
           mask_fewfreq[wh_fewfreq] = 0
           mask_fewfreq = rebin(temporary(mask_fewfreq), n_kx, n_ky, n_kz)

           a_0 = temporary(a_0) * mask_fewfreq[*,*,0]
           a_n = temporary(a_n) * mask_fewfreq[*,*,1:*]
           b_n = temporary(b_n) * mask_fewfreq[*,*,1:*]
        endif

        data_cos = complex(dblarr(n_kx, n_ky, n_kz))
        data_sin = complex(dblarr(n_kx, n_ky, n_kz))
        data_cos[*, *, 0] = a_0 /2d
        data_cos[*, *, 1:n_kz-1] = a_n
        data_sin[*, *, 1:n_kz-1] = b_n
 
        ;; for new power calc, need cos2, sin2, cos*sin transforms
        ;; have to do this in a for loop for memory's sake
        covar_cos = dblarr(n_kx, n_ky, n_kz)
        covar_sin = dblarr(n_kx, n_ky, n_kz)
        covar_cross = dblarr(n_kx, n_ky, n_kz)
        
        ;; comov_dist_los goes from large to small z
        z_relative = dindgen(n_freq)*z_mpc_delta
        freq_kz_arr = rebin(reform(rebin(reform(kz_mpc, 1, n_kz), n_freq, n_kz) * $
                                   rebin(z_relative, n_freq, n_kz), 1, n_freq, n_kz), n_ky, n_freq, n_kz)
        
        cos_arr = cos(freq_kz_arr)
        sin_arr = sin(freq_kz_arr)
        ;;z_exp_arr = exp(-1*dcomplex(0,1)*freq_kz_arr)
        
        for i=0, n_kx-1 do begin
           if max(n_freq_contrib[i,*]) eq 0 then continue
           wh_divide = where(n_freq_contrib[i,*] gt 0, count_divide)
           n_divide = double(rebin(reform(n_freq_contrib[i,wh_divide]),count_divide, n_kz))

           sigma2_arr = rebin(reform(sigma2_cube[i,wh_divide,*]), count_divide, n_freq, n_kz) 

           covar_cos[i,wh_divide,*] = total(sigma2_arr*cos_arr[wh_divide, *, *]^2d, 2)/n_divide
           covar_sin[i,wh_divide,*] = total(sigma2_arr*sin_arr[wh_divide, *, *]^2d, 2)/n_divide
           covar_cross[i,wh_divide,*] = total(sigma2_arr*cos_arr[wh_divide, *, *]*sin_arr[wh_divide, *, *], 2)/n_divide
        endfor
        
        ;; drop pixels with less than 1/3 of the frequencies
        if count_fewfreq gt 0 then begin
           covar_cos = temporary(covar_cos) * mask_fewfreq
           covar_sin = temporary(covar_sin) * mask_fewfreq
           covar_cross = temporary(covar_cross) * mask_fewfreq
           undefine, mask_fewfreq
        endif
        
        ;; cos 0 term has different normalization
        covar_cos[*,*,0] = covar_cos[*,*,0]/4d
      
        undefine, sigma2_cube
        undefine, freq_kz_arr
        undefine, cos_arr
        undefine, sin_arr
        undefine, sigma2_arr

        ;; factor to go to eor theory FT convention
        ;; I don't think I need these factors in the covariance
        ;; matrix because I've use the FT & inv FT -- should cancel
        ;; covar_cos = factor * temporary(covar_cos2)
        ;; covar_sin = factor * temporary(covar_sin2)
        ;; covar_cross = factor * temporary(covar_cross)

        ;; get rotation angle to diagonalize covariance block
        theta = atan(2d*covar_cross, covar_cos - covar_sin)/2d
        cos_theta = cos(theta)
        sin_theta = sin(theta)
        undefine, theta

        sigma1_2 = covar_cos*cos_theta^2 + 2d*covar_cross*cos_theta*sin_theta + covar_sin*sin_theta^2d
        sigma2_2 = covar_cos*sin_theta^2d - 2d*covar_cross*cos_theta*sin_theta + covar_sin*cos_theta^2d

        undefine, covar_cos
        undefine, covar_sin
        undefine, covar_cross
       
        data1 = data_cos*cos_theta + data_sin*sin_theta
        data2 = (-1d)*data_cos*sin_theta + data_sin*cos_theta
        undefine, data_cos
        undefine, data_sin

        ;;weights_1 = (1/sigma1_2)^2d
        ;;weights_1 = 1d/(4d*real_part(data1 * conj(data1))*sigma1_2)
        weights_1 = 1d/(4*(sigma1_2)^2d)
        term1 = real_part(data1 * conj(data1))*weights_1
        wh_sig1_0 = where(sigma1_2^2d eq 0, count_sig1_0)
        if count_sig1_0 ne 0 then begin
           weights_1[wh_sig1_0] = 0
           term1[wh_sig1_0] = 0
        endif
        
        ;;weights_2 = (1/sigma2_2)^2d
        ;;weights_2 = 1d/(4d*real_part(data2 * conj(data2))*sigma2_2)
        weights_2 = 1d/(4*(sigma2_2)^2d)
        term2 = real_part(data2 * conj(data2))*weights_2
        wh_sig2_0 = where(sigma2_2^2d eq 0, count_sig2_0)
        if count_sig2_0 ne 0 then begin
           weights_2[wh_sig2_0] = 0
           term2[wh_sig2_0] = 0
        endif
        undefine, data1
        undefine, data2
        undefine, sigma1_2
        undefine, sigma2_2   

        weights_3d = weights_1 + weights_2
        ;;weights_3d = 1d/(1d/weights_1 + 1d/weights_2)
        if count_sig1_0 gt 0 then weights_3d[wh_sig1_0] = weights_2[wh_sig1_0]
        if count_sig2_0 gt 0 then weights_3d[wh_sig2_0] = weights_1[wh_sig2_0]
        power_3d = (term1 + term2) / weights_3d
        ;;power_error = 1/weights
        wh_wt0 = where(weights_3d eq 0, count_wt0)
        if count_wt0 ne 0 then begin
           power_3d[wh_wt0] = 0
        endif

        undefine, term1
        undefine, term2
        undefine, weights1
        undefine, weights2

     endelse

     save, file = save_file, power_3d, weights_3d, kx_mpc, ky_mpc, kz_mpc, kperp_lambda_conv, delay_params, hubble_param, $
           n_freq_contrib

     write_ps_fits, fits_savefile, power_3d, weights_3d, kx_mpc, ky_mpc, kz_mpc, kperp_lambda_conv, delay_params, hubble_param

  endif else restore, save_file

  print, 'power integral:', total(power_3d)

  n_kx = n_elements(kx_mpc)
  n_ky = n_elements(ky_mpc)
  n_kz = n_elements(kz_mpc)

  if keyword_set (no_kzero) then begin 
     ;; leave out kz=0 -- full of foregrounds
     kz_mpc = kz_mpc[1:*]
     power_3d = temporary(power_3d[*, *, 1:*])
     weights_3d = temporary(weights_3d[*,*,1:*])
     n_kz = n_elements(kz_mpc)
  endif

  fadd = ''
  if keyword_set(no_weighted_averaging) then fadd = fadd + '_nowtave'
  if keyword_set(no_kzero) then fadd = fadd + '_nok0'

  fadd_2d = ''
  if keyword_set(fill_holes) then fadd_2d = fadd_2d + '_nohole'
  if keyword_set(log_kpar) then fadd_2d = fadd_2d + '_logkpar'
  if keyword_set(log_kperp) then fadd_2d = fadd_2d + '_logkperp'

  savefile = froot + savefilebase + fadd + fadd_2d + '_2dkpower.idlsave'

  print, 'Binning to 2D power spectrum'

  if keyword_set(no_weighting) then $
     power_rebin = kspace_rebinning_2d(power_3D, kx_mpc, ky_mpc, kz_mpc, kperp_edges_mpc, kpar_edges_mpc, log_kpar = log_kpar, $
                                       log_kperp = log_kperp, kperp_bin = kperp_bin, kpar_bin = kpar_bin, $
                                       binned_weights = binned_weights, fill_holes = fill_holes) $
  else begin
     power_rebin = kspace_rebinning_2d(power_3D, kx_mpc, ky_mpc, kz_mpc, kperp_edges_mpc, kpar_edges_mpc, log_kpar = log_kpar, $
                                       log_kperp = log_kperp, kperp_bin = kperp_bin, kpar_bin = kpar_bin, $
                                       weights = weights_3d, binned_weights = binned_weights, fill_holes = fill_holes)
     if keyword_set(no_weighted_averaging) then $
        power_rebin = kspace_rebinning_2d(power_3D, kx_mpc, ky_mpc, kz_mpc, kperp_edges_mpc, kpar_edges_mpc, log_kpar = log_kpar, $
                                       log_kperp = log_kperp, kperp_bin = kperp_bin, kpar_bin = kpar_bin, fill_holes = fill_holes)
  endelse

  power = power_rebin
  kperp_edges = kperp_edges_mpc
  kpar_edges = kpar_edges_mpc
  weights = binned_weights

  wh_good_kperp = where(total(weights, 2) gt 0, count)
  if count eq 0 then stop
  kperp_plot_range = [min(kperp_edges[wh_good_kperp]), max(kperp_edges[wh_good_kperp+1])]
  
  save, file = savefile, power, weights, kperp_edges, kpar_edges, kperp_bin, kpar_bin, kperp_lambda_conv, delay_params, hubble_param

  if not keyword_set(quiet) then begin
     kpower_2d_plots, savefile, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
                      data_range = data_range
     kpower_2d_plots, savefile, /plot_weights, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
                      window_num = 2, title = 'Weights'
  endif

  ;; now do slices    
  yslice_savefile = froot + savefilebase + '_xz_plane.idlsave'
  yslice_power = kpower_slice(power_3d, kx_mpc, ky_mpc, kz_mpc, weights_3d, kperp_lambda_conv, delay_params, hubble_param, $
                              slice_axis = 1, slice_inds = 0, slice_weights = yslice_weights, slice_savefile = yslice_savefile)

  xslice_savefile = froot + savefilebase + '_yz_plane.idlsave'
  xslice_power = kpower_slice(power_3d, kx_mpc, ky_mpc, kz_mpc, weights_3d, kperp_lambda_conv, delay_params, hubble_param, $
                              slice_axis = 0, slice_inds = n_kx/2, slice_weights = xslice_weights, slice_savefile = xslice_savefile)
  if max(xslice_power) eq 0 then begin
     nloop = 0
     while max(xslice_power) eq 0 do begin
         nloop = nloop+1
         xslice_power = kpower_slice(power_3d, kx_mpc, ky_mpc, kz_mpc, weights_3d, kperp_lambda_conv, delay_params, hubble_param, $
                                     slice_axis = 0, slice_inds = n_kx/2+nloop, slice_weights = xslice_weights, $
                                     slice_savefile = xslice_savefile)
     endwhile
  endif

  zslice_savefile = froot + savefilebase + '_xy_plane.idlsave'
  zslice_power = kpower_slice(power_3d, kx_mpc, ky_mpc, kz_mpc, weights_3d, kperp_lambda_conv, delay_params, hubble_param, $
                              slice_axis = 2, slice_inds = 1, slice_weights = zslice_weights, slice_savefile = zslice_savefile)



  print, 'Binning to 1D power spectrum'
 
  if keyword_set(no_weighting) then $
     power_1d = kspace_rebinning_1d(power_3d, kx_mpc, ky_mpc, kz_mpc, k_edges_mpc, k_bin = k1d_bin, log_k = log_k1d, $
                                    binned_weights = weights_1d, mask = mask, pixelwise_mask = pixelwise_mask, k1_mask = k1_mask, $
                                    k2_mask = k2_mask,  k3_mask = k3_mask, edge_on_grid = edge_on_grid, match_datta = match_datta)$
  else begin
     power_1d = kspace_rebinning_1d(power_3d, kx_mpc, ky_mpc, kz_mpc, k_edges_mpc, k_bin = k1d_bin, log_k = log_k1d, $
                                    weights = weights_3d, binned_weights = weights_1d, mask = mask, pixelwise_mask = pixelwise_mask, $
                                    k1_mask = k1_mask, k2_mask = k2_mask,  k3_mask = k3_mask, edge_on_grid = edge_on_grid, $
                                    match_datta = match_datta)
     if keyword_set(no_weighted_averaging) then $
        power_1d = kspace_rebinning_1d(power_3d, kx_mpc, ky_mpc, kz_mpc, k_edges_mpc, k_bin = k1d_bin, log_k = log_k1d, mask = mask, $
                                       pixelwise_mask = pixelwise_mask, k1_mask = k1_mask, k2_mask = k2_mask,  k3_mask = k3_mask, $
                                       edge_on_grid = edge_on_grid, match_datta = match_datta)
  endelse

  power = power_1d
  weights = weights_1d
  k_edges = k_edges_mpc
  k_bin = k1d_bin

  fadd_1d = ''
  if keyword_set(log_k) then fadd_1d = fadd_1d + '_logk'

  savefile = froot + savefilebase + fadd + fadd_1d + '_1dkpower.idlsave'
  save, file = savefile, power, weights, k_edges, k_bin

  if not keyword_set(quiet) then begin
     kpower_1d_plots, savefile, window_num = 5
  endif

  ;; eor_file_1d = base_path() + 'power_spectrum/eor_data/eor_power_1d.idlsave'
  ;; file_arr = [savefile, eor_file_1d]
  ;; if keyword_set(eor_only) then begin
  ;;    if keyword_set(eor_test) then names_arr = 'Input EoR' else names_arr = 'Simulated EoR'
  ;; endif else names_arr = 'Simulation PS'
  ;; names_arr = [names_arr, 'EoR signal']
  ;; colors_arr = [0, 254]

  ;;   if not keyword_set(quiet) then begin
  ;;      kpower_1d_plots, file_arr, window_num = 5, names = names_arr, colors = colors_arr
  ;;   endif
end
