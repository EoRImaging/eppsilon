pro fhd_kcube, file_struct, dft_refresh_data = dft_refresh_data, dft_refresh_weight = dft_refresh_weight, $
               dft_fchunk = dft_fchunk, freq_ch_range = freq_ch_range, spec_window_type = spec_window_type, $
               std_power = std_power, input_units = input_units, quiet = quiet

  if n_elements(file_struct.nside) ne 0 then healpix = 1 else healpix = 0
  nfiles = n_elements(file_struct.datafile)

  if n_elements(input_units) eq 0 then input_units = 'jansky'
  units_enum = ['jansky', 'mk']
  wh = where(units_enum eq input_units, count)
  if count eq 0 then message, 'input units not recognized, options are: ' + units_enum
  
  datavar = strupcase(file_struct.datavar)
  if datavar ne '' then begin
     for i=0, nfiles-1 do begin
        datafile_obj = obj_new('IDL_Savefile', file_struct.datafile[i])
        datafile_names = datafile_obj->names()
        wh = where(datafile_names eq datavar, count)
        if count eq 0 then message, 'specified datavar is not present in datafile (datafile=' + file_struct.datafile[i] + $
                                    ', datavar=' + file_struct.datavar + ')'
        data_dims = datafile_obj->size(file_struct.datavar, /dimensions)
        obj_destroy, datafile_obj
  
        if i gt 0 then if total(abs(data_dims - dims)) ne 0 then message, 'data dimensions in files do not match'
        
        weightfile_obj = obj_new('IDL_Savefile', file_struct.weightfile[i])
        weightfile_names = weightfile_obj->names()
        weightvar = strupcase(file_struct.weightvar)
        wh = where(weightfile_names eq weightvar, count)
        if count eq 0 then message, 'specified weightvar is not present in weightfile (weightfile=' + file_struct.weightfile[i] + $
                                    ', weightvar=' + file_struct.weightvar + ')'
        weight_dims = weightfile_obj->size(weightvar, /dimensions)
        obj_destroy, weightfile_obj
        
        if total(abs(data_dims - weight_dims)) ne 0 then message, 'data and weight dimensions do not match'
        undefine, weight_dims
        
        
        variancefile_obj = obj_new('IDL_Savefile', file_struct.variancefile[i])
        variancefile_names = variancefile_obj->names()
        variancevar = strupcase(file_struct.variancevar)
        wh = where(variancefile_names eq variancevar, count)
        if count eq 0 then begin
           print, 'specified variancevar is not present in variancefile (variancefile=' + file_struct.variancefile[i] $
                  +  ', variancevar=' + file_struct.variancevar + '). Weights will be used instead' 
           no_var = 1
        endif else begin
           if n_elements(no_var) eq 0 then no_var = 0

           variance_dims = variancefile_obj->size(variancevar, /dimensions)
           if total(abs(data_dims - variance_dims)) ne 0 then message, 'data and variance dimensions do not match'
           undefine, variance_dims
        endelse
        obj_destroy, variancefile_obj
     
        dims = data_dims
        undefine, data_dims
     endfor


     if healpix then n_freq = dims[1] else n_freq = dims[2]
  endif else begin
     ;; working with a 'derived' cube (ie residual cube) that is constructed from uvf_savefiles
     input_uvf_files = reform(file_struct.res_uvf_inputfiles, nfiles, 2)
     for i=0, n_elements(input_uvf_files)-1 do begin
        datafile_obj = obj_new('IDL_Savefile', input_uvf_files[i])
        datafile_names = datafile_obj->names()
        wh = where(datafile_names eq 'data_cube', count)
        if count eq 0 then message, 'specified res_uvf_inputfile does not contain a data cube (res_uvf_inputfile=' + $
                                    input_uvf_files[i] + ')'
        data_dims = datafile_obj->size('data_cube', /dimensions)
        obj_destroy, datafile_obj

        dims = data_dims
        if i gt 0 then if total(abs(data_dims - dims)) ne 0 then message, 'data dimensions in files do not match'
        undefine, data_dims
     endfor

     input_uvf_wtfiles = file_struct.uvf_weight_savefile
     for i=0, n_elements(input_uvf_wtfiles)-1 do begin
        weightfile_obj = obj_new('IDL_Savefile', input_uvf_wtfiles[i])
        weightfile_names = datafile_obj->names()
        wh = where(weightfile_names eq 'weights_cube', count)
        if count eq 0 then message, 'specified uvf_weight_savefile does not contain a weights cube (res_uvf_inputfile=' + $
                                    input_uvf_wtfiles[i] + ')'
        weights_dims = datafile_obj->size('weights_cube', /dimensions)
        wh = where(weightfile_names eq 'variance_cube', count)
        if count eq 0 then message, 'specified uvf_weight_savefile does not contain a variance cube (res_uvf_inputfile=' + $
                                    input_uvf_wtfiles[i] + ')'
        variance_dims = datafile_obj->size('variance_cube', /dimensions)
        obj_destroy, datafile_obj
     
        if total(abs(dims - weight_dims)) ne 0 then message, 'data and weight dimensions do not match'
        if total(abs(dims - variance_dims)) ne 0 then message, 'data and variance dimensions do not match'
     endfor
     undefine, weights_dims, variance_dims


     n_freq = dims[2]
  endelse     

  frequencies = file_struct.frequencies
  if n_elements(frequencies) ne n_freq then message, 'number of frequencies does not match frequency dimension of data'

  if n_elements(freq_ch_range) ne 0 then frequencies = frequencies[min(freq_ch_range):max(freq_ch_range)]
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

     ;; converting from Jy (in u,v,f) to mK*str
     conv_factor = float((3e8)^2 / (2. * (frequencies*1e6)^2. * 1.38065))

     ;; mK str -> mK Mpc^2
     conv_factor = conv_factor * z_mpc_mean^2.

  endif else conv_factor = 1. + fltarr(n_freq)
  
  ;;t_sys = 440. ; K
  t_sys = 280. * sqrt(2.) * ((1+redshifts)/7.5)^2.3 ;; from skew w/ stu + srt(2) for single pol
  ;;eff_area = 16. ; m^2
  eff_area = 21. ; m^2 -- from Aaron's memo
  df = file_struct.freq_resolution ; Hz -- native visibility resolution NOT cube resolution
  tau = file_struct.time_resolution ; seconds
  vis_sigma = (2. * (1.38065e-23) * 1e26) * t_sys / (eff_area * sqrt(df * tau)) ;; in Jy
  vis_sigma = float(vis_sigma)


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
           if datavar eq '' then begin
              ;; working with a 'derived' cube (ie residual cube) that is constructed from uvf_savefiles
              restore, input_uvf_files[i,0]
              kx_dirty = temporary(kx_rad_vals)
              ky_dirty = temporary(ky_rad_vals)
              dirty_cube = temporary(data_cube)

              restore, input_uvf_files[i,1]
              model_cube = temporary(data_cube)
              
              if total(abs(kx_rad_vals - kx_dirty)) ne 0 then message, 'kx_rad_vals for dirty and model cubes must match'
              if total(abs(ky_rad_vals - ky_dirty)) ne 0 then message, 'kx_rad_vals for dirty and model cubes must match'
              undefine, kx_dirty, ky_dirty
stop
              data_cube = temporary(dirty_cube) - temporary(model_cube)
              save, file = file_struct.uvf_savefile[i], kx_rad_vals, ky_rad_vals, data_cube
              undefine, data_cube
           endif else begin        
              test_setup = file_test(file_struct.hpx_dftsetup_savefile) * $
                           (1 - file_test(file_struct.hpx_dftsetup_savefile, /zero_length))
              if test_setup eq 0 then begin
                 ;; figure out k values to calculate dft
                 healpix_setup_ft, pixel_nums1, file_struct.nside, new_pix_vec, limits, kx_rad_vals, ky_rad_vals, /quiet
                 save, file = file_struct.hpx_dftsetup_savefile, new_pix_vec, limits, kx_rad_vals, ky_rad_vals
              endif else restore, file_struct.hpx_dftsetup_savefile
              
              ;; drop kperp >> max_baseline to save time on DFT
              ;; go a little beyond max_baseline to account for expansion due to w projection
              max_kperp_rad = (file_struct.max_baseline_lambda/kperp_lambda_conv) * z_mpc_mean * 1.1
              
              wh_kx_good = where(abs(kx_rad_vals) le max_kperp_rad, count_kx)
              wh_ky_good = where(abs(ky_rad_vals) le max_kperp_rad, count_ky)
              
              if count_kx gt 0 then kx_rad_vals = kx_rad_vals[wh_kx_good] else stop
              if count_ky gt 0 then ky_rad_vals = ky_rad_vals[wh_ky_good] else stop
              
              ;; need to cut uvf cubes in half because image is real -- we'll cut in v
              ;; drop the unused half before the DFT to save time
              n_ky = n_elements(ky_rad_vals)
              ky_rad_vals = ky_rad_vals[n_ky/2:n_ky-1]
              
              ;; do DFT.
              if test_uvf eq 0 or keyword_set(dft_refresh_data) then begin
                 arr = getvar_savefile(file_struct.datafile[i], file_struct.datavar)
                 if n_elements(freq_ch_range) ne 0 then arr = arr[*, min(freq_ch_range):max(freq_ch_range)]
                 
                 transform = discrete_ft_2D_fast(new_pix_vec[*,0], new_pix_vec[*,1], arr, kx_rad_vals, ky_rad_vals, $
                                                 max_k_mag = max_kperp_rad, timing = ft_time, fchunk = dft_fchunk)
                 data_cube = temporary(transform)
                 undefine, arr
                 
                 save, file = file_struct.uvf_savefile[i], kx_rad_vals, ky_rad_vals, data_cube
                 undefine, data_cube
              endif
              
              if test_wt_uvf eq 0 or keyword_set(dft_refresh_weight) then begin
                 arr = getvar_savefile(file_struct.weightfile[i], file_struct.weightvar)
                 if n_elements(freq_ch_range) ne 0 then arr = arr[*, min(freq_ch_range):max(freq_ch_range)]
                 
                 transform = discrete_ft_2D_fast(new_pix_vec[*,0], new_pix_vec[*,1], arr, kx_rad_vals, ky_rad_vals, timing = ft_time, $
                                                 fchunk = dft_fchunk)            
                 weights_cube = temporary(transform)
                 
                 if not no_var then begin
                    arr = getvar_savefile(file_struct.variancefile[i], file_struct.variancevar)
                    if n_elements(freq_ch_range) ne 0 then arr = arr[*, min(freq_ch_range):max(freq_ch_range)]
                    
                    transform = discrete_ft_2D_fast(new_pix_vec[*,0], new_pix_vec[*,1], arr, kx_rad_vals, ky_rad_vals, $
                                                    timing = ft_time, fchunk = dft_fchunk)            
                    variance_cube = abs(temporary(transform)) ;; make variances real, positive definite (amplitude)
                    undefine, arr
                 endif
 
                 save, file = file_struct.uvf_weight_savefile[i], kx_rad_vals, ky_rad_vals, weights_cube, variance_cube
                 undefine, new_pix_vec, weights_cube, variance_cube
              endif
           endelse
        endif else begin
           kx_rad_vals = getvar_savefile(file_struct.uvf_savefile[0], 'kx_rad_vals')
           ky_rad_vals = getvar_savefile(file_struct.uvf_savefile[0], 'ky_rad_vals')
        endelse
     endfor
  endif

  if healpix then begin
     n_kx = n_elements(kx_rad_vals)
     kx_rad_delta = kx_rad_vals[1] - kx_rad_vals[0]
     kx_mpc = temporary(kx_rad_vals) / z_mpc_mean
     kx_mpc_delta = kx_mpc[1] - kx_mpc[0]     

     n_ky = n_elements(ky_rad_vals)
     ky_rad_delta = ky_rad_vals[1] - ky_rad_vals[0]
     ky_mpc = temporary(ky_rad_vals) / z_mpc_mean
     ky_mpc_delta = ky_mpc[1] - ky_mpc[0]     
   
     ;; Angular resolution is given in Healpix paper in units of arcminutes, need to convert to radians
     ang_resolution = sqrt(3./!pi) * 3600./file_struct.nside * (1./60.) * (!pi/180.)
     pix_area_rad = ang_resolution^2. ;; by definition of ang. resolution in Healpix paper
     
     pix_area_mpc = pix_area_rad * z_mpc_mean^2.
     
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
     
     pix_area_mpc = x_mpc_delta * y_mpc_delta

  endelse

  if nfiles eq 2 then begin
     if no_var then begin
        ;; use 1/abs(weights) instead
        if healpix then begin
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

        sigma2_cube1 = 1./abs(weights_cube1)
        sigma2_cube2 = 1./abs(weights_cube2)

     endif else begin
        if healpix then begin
           variance_cube1 = getvar_savefile(file_struct.uvf_weight_savefile[0], 'variance_cube')
           variance_cube2 = getvar_savefile(file_struct.uvf_weight_savefile[1], 'variance_cube')
           weights_cube1 = getvar_savefile(file_struct.uvf_weight_savefile[0], 'weights_cube')
           weights_cube2 = getvar_savefile(file_struct.uvf_weight_savefile[1], 'weights_cube')

           if min(ky_mpc) lt 0 then begin
              ;; calculate integral of window function before cut for comparison
              window_int_orig = [total(variance_cube1)*pix_area_rad/file_struct.n_vis[0], $
                                 total(variance_cube2)*pix_area_rad/file_struct.n_vis[1]]

              ;; negative ky values haven't been cut yet
              ;; need to cut uvf cubes in half because image is real -- we'll cut negative ky
              variance_cube1 = variance_cube1[*, n_ky/2:n_ky-1,*]
              variance_cube2 = variance_cube2[*, n_ky/2:n_ky-1,*]
              weights_cube1 = weights_cube1[*, n_ky/2:n_ky-1,*]
              weights_cube2 = weights_cube2[*, n_ky/2:n_ky-1,*]
           endif

        endif else begin
           variance_cube1 = getvar_savefile(file_struct.variancefile[0], file_struct.variancevar)
           variance_cube2 = getvar_savefile(file_struct.variancefile[1], file_struct.variancevar)
           weights_cube1 = getvar_savefile(file_struct.weightfile[0], file_struct.weightvar)
           weights_cube2 = getvar_savefile(file_struct.weightfile[1], file_struct.weightvar)

           ;; need to cut uvf cubes in half because image is real -- we'll cut negative ky
           variance_cube1 = variance_cube1[*, n_ky/2:n_ky-1,*]
           variance_cube2 = variance_cube2[*, n_ky/2:n_ky-1,*]
           weights_cube1 = weights_cube1[*, n_ky/2:n_ky-1,*]
           weights_cube2 = weights_cube2[*, n_ky/2:n_ky-1,*]
       endelse

        ;; Also need to drop 1/2 of ky=0 line -- drop ky=0, kx<0
        weights_cube1[0:n_kx/2-1, 0, *] = 0
        weights_cube2[0:n_kx/2-1, 0, *] = 0
        variance_cube1[0:n_kx/2-1, 0, *] = 0
        variance_cube2[0:n_kx/2-1, 0, *] = 0
        
        ;; calculate integral of window function (use pix_area_rad for FT normalization)
        ;; already cut out negative ky, so multiply by 2
        window_int = 2*[total(variance_cube1)*pix_area_rad/file_struct.n_vis[0], $
                        total(variance_cube2)*pix_area_rad/file_struct.n_vis[1]]

        sigma2_cube1 = temporary(variance_cube1) / (abs(weights_cube1)^2.*pix_area_rad)
        sigma2_cube2 = temporary(variance_cube2) / (abs(weights_cube2)^2.*pix_area_rad)
     endelse

     wh_wt1_0 = where(abs(weights_cube1) eq 0, count_wt1_0)
     if count_wt1_0 ne 0 then sigma2_cube1[wh_wt1_0] = 0
     wh_wt2_0 = where(abs(weights_cube2) eq 0, count_wt2_0)
     if count_wt2_0 ne 0 then sigma2_cube2[wh_wt2_0] = 0
 
     if min(sigma2_cube1) lt 0 or min(sigma2_cube2) lt 0 then message, 'sigma2 should be positive definite.'
     if total(abs(sigma2_cube1)) le 0 or total(abs(sigma2_cube2)) le 0 then message, 'one or both sigma2 cubes is all zero'

     ;; now get data cubes
     if healpix then begin
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
        data_cube1 = getvar_savefile(file_struct.datafile[0], file_struct.datavar)
        data_cube2 = getvar_savefile(file_struct.datafile[1], file_struct.datavar)

        ;; need to cut uvf cubes in half because image is real -- we'll cut negative ky
        data_cube1 = data_cube1[*, n_ky/2:n_ky-1,*]
        data_cube2 = data_cube2[*, n_ky/2:n_ky-1,*]
        
        ky_mpc = ky_mpc[n_ky/2:n_ky-1]
        n_ky = n_elements(ky_mpc)
     endelse
   
    ;; Also need to drop 1/2 of ky=0 line -- drop ky=0, kx<0
    data_cube1[0:n_kx/2-1, 0, *] = 0
    data_cube2[0:n_kx/2-1, 0, *] = 0

  
  endif else begin
     ;; single file mode
     if healpix then begin
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
           window_int = 2*total(variance_cube1)*pix_area_rad/file_struct.n_vis

           sigma2_cube1 = temporary(variance_cube1) / abs(weights_cube1)^2.
        endif else sigma2_cube1 = 1./abs(weights_cube1)

        if min(ky_mpc) lt 0 then begin
           ;; negative ky values haven't been cut yet
           ;; need to cut uvf cubes in half because image is real -- we'll cut negative ky
           ky_mpc = ky_mpc[n_ky/2:n_ky-1]
           n_ky = n_elements(ky_mpc)
        endif

        wh_wt1_0 = where(abs(weights_cube1) eq 0, count_wt1_0)
        if count_wt1_0 ne 0 then sigma2_cube1[wh_wt1_0] = 0
     endif else begin
        data_cube1 = getvar_savefile(file_struct.datafile, file_struct.datavar)
        weights_cube1 = getvar_savefile(file_struct.weightfile, file_struct.weightvar)

        ;; need to cut uvf cubes in half because image is real -- we'll cut negative ky
        weights_cube1 = weights_cube1[*, n_ky/2:n_ky-1,*]
        data_cube1 = data_cube1[*, n_ky/2:n_ky-1,*]        

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
           window_int = 2*total(variance_cube1)*pix_area_rad/file_struct.n_vis

           sigma2_cube1 = temporary(variance_cube1) / abs(weights_cube1)^2.
        endif else sigma2_cube1 = 1./abs(weights_cube1)
  
        ;; need to cut uvf cubes in half because image is real -- we'll cut negative ky
        ky_mpc = ky_mpc[n_ky/2:n_ky-1]
        n_ky = n_elements(ky_mpc)

        wh_wt1_0 = where(abs(weights_cube1) eq 0, count_wt1_0)
        if count_wt1_0 ne 0 then sigma2_cube1[wh_wt1_0] = 0
     endelse
  endelse
  
 
  ;; multiply data and sigma cubes by pixel_area_mpc to get proper units from DFT
  ;; divide by (2*!pi)^2d to get right FT convention
  ;; (not squared for sigmas because they weren't treated as units squared in FHD code)
  data_cube1 = data_cube1 * pix_area_mpc / (2.*!pi)^2.
  sigma2_cube1 = sigma2_cube1 * pix_area_mpc / (2.*!pi)^2.
  if nfiles eq 2 then begin
     data_cube2 = data_cube2 * pix_area_mpc / (2.*!pi)^2.
     sigma2_cube2 = sigma2_cube2 * pix_area_mpc / (2.*!pi)^2.
  endif

  ;; get sigma into Jy
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

  ;; divide by weights
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
                          slice_inds = 0, slice_savefile = file_struct.uf_savefile[i])
     
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

  ;; now take FFT
  data_sum_ft = fft(data_sum, dimension=3) * n_freq * z_mpc_delta / (2.*!pi)
  ;; put k0 in middle of cube
  data_sum_ft = shift(data_sum_ft, [0,0,n_kz/2])
  undefine, data_sum
  if nfiles eq 2 then begin
     data_diff_ft = fft(data_diff, dimension=3) * n_freq * z_mpc_delta / (2.*!pi)
     ;; put k0 in middle of cube
     data_diff_ft = shift(data_diff_ft, [0,0,n_kz/2])
     undefine, data_diff
  endif

  ;; factor to go to eor theory FT convention
  ;; Only 1 factor of 2pi b/c only transforming along z
  ;; factor = (2d*!pi)
  ;; if not keyword_set(eor_test) then data_ft = factor * temporary(data_ft)
  
  ;; print, 'full_ft^2d integral (after theory factor):', total(abs(data_ft)^2d)
  
  n_val = round(kz_mpc_orig / kz_mpc_delta)
  kz_mpc_orig[where(n_val eq 0)] = 0
  a1_0 = 2. * data_sum_ft[*,*,where(n_val eq 0)]
  a1_n = data_sum_ft[*,*, where(n_val gt 0)] + data_sum_ft[*,*, reverse(where(n_val lt 0))]
  b1_n = complex(0,1) * (data_sum_ft[*,*, where(n_val gt 0)] - data_sum_ft[*,*, reverse(where(n_val lt 0))])
  undefine, data_sum_ft

  if nfiles gt 1 then begin
     a2_0 = 2. * data_diff_ft[*,*,where(n_val eq 0)]
     a2_n = data_diff_ft[*,*, where(n_val gt 0)] + data_diff_ft[*,*, reverse(where(n_val lt 0))]
     b2_n = complex(0,1) * (data_diff_ft[*,*, where(n_val gt 0)] - data_diff_ft[*,*, reverse(where(n_val lt 0))])
     undefine, data_diff_ft
  endif

  kz_mpc = kz_mpc_orig[where(n_val ge 0)]
  n_kz = n_elements(kz_mpc)
  
  if keyword_set(std_power) then begin
     ;; for standard power calc. just need ft of sigma2 (sigma has squared units relative to data, so use z_mpc_delta^2d)
     sigma2_ft = fft(sum_sigma2, dimension=3) * n_freq * z_mpc_delta^2. / (2.*!pi)
     sigma2_ft = shift(sigma2_ft, [0,0,n_kz/2])
     undefine, sum_sigma2
     
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
