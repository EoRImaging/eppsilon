pro ps_kcube, file_struct, dft_refresh_data = dft_refresh_data, dft_refresh_weight = dft_refresh_weight, refresh_beam = refresh_beam, dft_ian = dft_ian, $
    dft_fchunk = dft_fchunk, freq_ch_range = freq_ch_range, freq_flags = freq_flags, $
    spec_window_type = spec_window_type, cut_image = cut_image, delta_uv_lambda = delta_uv_lambda, max_uv_lambda = max_uv_lambda, sim=sim, $
    std_power = std_power, input_units = input_units, uvf_input = uvf_input, uv_avg = uv_avg, uv_img_clip = uv_img_clip, no_dft_progress = no_dft_progress
    
  if tag_exist(file_struct, 'nside') ne 0 then healpix = 1 else healpix = 0
  if keyword_set(uvf_input) or tag_exist(file_struct, 'uvf_savefile') eq 0 then uvf_input = 1 else uvf_input = 0
  nfiles = n_elements(file_struct.datafile)
  if tag_exist(file_struct, 'no_var') ne 0 then no_var = 1 else no_var = 0
  
  if tag_exist(file_struct, 'beam_savefile') eq 0 then refresh_beam = 0
  
  if n_elements(input_units) eq 0 then input_units = 'jansky'
  units_enum = ['jansky', 'mk']
  wh = where(units_enum eq input_units, count)
  if count eq 0 then message, 'input units not recognized, options are: ' + units_enum
  
  
  datavar = strupcase(file_struct.datavar)
  if datavar eq '' then begin
    ;; working with a 'derived' cube (ie residual cube) that is constructed from uvf_savefiles
    input_uvf_files = reform(file_struct.res_uvf_inputfiles, nfiles, 2)
    input_uvf_varname = reform(file_struct.res_uvf_varname, nfiles, 2)
    
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
    if max(freq_diff_hist)/float(n_freq) lt .5 then stop else begin
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
      if n_elements(vis_sigma_adam) ne n_freq_orig then stop
      vis_sigma_adam = vis_sigma_adam[min(freq_ch_range):max(freq_ch_range)]
    endif else if n_elements(vis_sigma_adam) ne n_freq then stop
    
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
  
  t_sys_used = (eff_area * sqrt(df * tau) * vis_sigma) / ((2. * (1.38065e-23) * 1e26))  ;; in K
  
  n_vis = reform(file_struct.n_vis)
  n_vis_freq = reform(file_struct.n_vis_freq)
  if n_elements(freq_ch_range) ne 0 then n_vis = total(n_vis_freq[*, min(freq_ch_range):max(freq_ch_range)], 2)
  
  if healpix or not uvf_input then begin
  
    for i=0, nfiles-1 do begin
      test_uvf = file_test(file_struct.uvf_savefile[i]) *  (1 - file_test(file_struct.uvf_savefile[i], /zero_length))
      
      test_wt_uvf = file_test(file_struct.uvf_weight_savefile[i]) * (1 - file_test(file_struct.uvf_weight_savefile[i], /zero_length))
      
      if tag_exist(file_struct, 'beam_savefile') then $
        test_beam = file_test(file_struct.beam_savefile[i]) * ( 1- file_test(file_struct.beam_savefile[i], /zero_length)) $
      else test_beam = 1
      
      if test_uvf eq 1 and n_elements(freq_mask) ne 0 then begin
        old_freq_mask = getvar_savefile(file_struct.uvf_savefile[i], 'freq_mask')
        if total(abs(old_freq_mask - freq_mask)) ne 0 then test_uvf = 0
      endif
      
      if test_wt_uvf eq 1 and n_elements(freq_mask) ne 0 then begin
        old_freq_mask = getvar_savefile(file_struct.uvf_savefile[i], 'freq_mask')
        if total(abs(old_freq_mask - freq_mask)) ne 0 then test_uvf = 0
      endif
      
      
      if test_uvf eq 0 and not keyword_set(dft_refresh_data) and (n_elements(freq_ch_range) ne 0 or n_elements(freq_flags) ne 0) then begin
        ;; if this is a limited freq. range cube, check for the full cube to avoid redoing the DFT
        full_uvf_file = file_struct.uvf_savefile[i]
        if n_elements(freq_ch_range) ne 0 then full_uvf_file = strjoin(strsplit(full_uvf_file, '_ch[0-9]+-[0-9]+', /regex, /extract))
        if n_elements(freq_flags) ne 0 then full_uvf_file = strjoin(strsplit(full_uvf_file, '_flag[a-z0-9]+', /regex, /extract))
        test_full_uvf = file_test(full_uvf_file) *  (1 - file_test(full_uvf_file, /zero_length))
        
        if test_full_uvf eq 1 then begin
          restore, full_uvf_file
          
          if n_elements(freq_flags) ne 0 then data_cube = data_cube * rebin(reform(freq_mask, 1, 1, n_elements(file_struct.frequencies)), size(data_cube, /dimension), /sample)
          if n_elements(freq_ch_range) ne 0 then data_cube = data_cube[*, *, min(freq_ch_range):max(freq_ch_range)]
          
          if keyword_set(dft_ian) then save, file = file_struct.uvf_savefile[i], u_lambda_vals, v_lambda_vals, data_cube, freq_mask, uvf_git_hash $
          else save, file = file_struct.uvf_savefile[i], kx_rad_vals, ky_rad_vals, data_cube, freq_mask, uvf_git_hash
          undefine, data_cube
          
          test_uvf = 1
        endif
        
      endif
      
      if test_wt_uvf eq 0 and not keyword_set(dft_refresh_weight) and (n_elements(freq_ch_range) ne 0 or n_elements(freq_flags) ne 0) then begin
        ;; if this is a limited freq. range cube, check for the full cube to avoid redoing the DFT
        full_uvf_wt_file = file_struct.uvf_weight_savefile[i]
        if n_elements(freq_ch_range) ne 0 then full_uvf_wt_file = strjoin(strsplit(full_uvf_wt_file, '_ch[0-9]+-[0-9]+', /regex, /extract))
        if n_elements(freq_flags) ne 0 then full_uvf_wt_file = strjoin(strsplit(full_uvf_wt_file, '_flag[a-z0-9]+', /regex, /extract))
        test_full_wt_uvf = file_test(full_uvf_wt_file) *  (1 - file_test(full_uvf_wt_file, /zero_length))
        if test_full_wt_uvf eq 1 then begin
          restore, full_uvf_wt_file
          
          if n_elements(freq_flags) ne 0 then begin
            flag_arr = rebin(reform(freq_mask, 1, 1, n_elements(file_struct.frequencies)), size(weights_cube, /dimension), /sample)
            weights_cube = weights_cube * flag_arr
            variance_cube = variance_cube * flag_arr
          endif
          
          if n_elements(freq_ch_range) ne 0 then begin
            weights_cube = weights_cube[*, *, min(freq_ch_range):max(freq_ch_range)]
            variance_cube = variance_cube[*, *, min(freq_ch_range):max(freq_ch_range)]
          endif
          
          if keyword_set(dft_ian) then $
            save, file = file_struct.uvf_weight_savefile[i], u_lambda_vals, v_lambda_vals, weights_cube, variance_cube, freq_mask, uvf_wt_git_hash else $
            save, file = file_struct.uvf_weight_savefile[i], kx_rad_vals, ky_rad_vals, weights_cube, variance_cube, freq_mask, uvf_wt_git_hash
          undefine, weights_cube, variance_cube
          
          test_wt_uvf = 1
        endif
        
      endif
      
      if test_uvf eq 0 or test_wt_uvf eq 0 or test_beam eq 0 or keyword_set(dft_refresh_data) or keyword_set(dft_refresh_weight) or keyword_set(refresh_beam) then begin
        if datavar eq '' then begin
          ;; working with a 'derived' cube (ie residual cube) that is constructed from uvf_savefiles
        
          dirty_cube = getvar_savefile(input_uvf_files[i,0], input_uvf_varname[i,0])
          kx_dirty = getvar_savefile(input_uvf_files[i,0], 'kx_rad_vals')
          ky_dirty = getvar_savefile(input_uvf_files[i,0], 'ky_rad_vals')
          
          model_cube = getvar_savefile(input_uvf_files[i,1], input_uvf_varname[i,1])
          kx_rad_vals = getvar_savefile(input_uvf_files[i,1], 'kx_rad_vals')
          ky_rad_vals = getvar_savefile(input_uvf_files[i,1], 'ky_rad_vals')
          
          if n_elements(freq_mask) ne 0 then begin
            dirty_freq_mask = getvar_savefile(input_uvf_files[i,0], 'freq_mask')
            model_freq_mask = getvar_savefile(input_uvf_files[i,1], 'freq_mask')
            if total(abs(dirty_freq_mask - freq_mask)) ne 0 then message, 'freq_mask of dirty file does not match current freq_mask'
            if total(abs(model_freq_mask - freq_mask)) ne 0 then message, 'freq_mask of model file does not match current freq_mask'
          endif
          
          if total(abs(kx_rad_vals - kx_dirty)) ne 0 then message, 'kx_rad_vals for dirty and model cubes must match'
          if total(abs(ky_rad_vals - ky_dirty)) ne 0 then message, 'kx_rad_vals for dirty and model cubes must match'
          undefine, kx_dirty, ky_dirty
          
          uvf_git_hash_dirty = getvar_savefile(input_uvf_files[i,0], 'uvf_git_hash')
          uvf_git_hash = getvar_savefile(input_uvf_files[i,1], 'uvf_git_hash')
          if uvf_git_hash_dirty ne uvf_git_hash then print, 'git hashes for dirty and model cubes does not match'
          
          data_cube = temporary(dirty_cube) - temporary(model_cube)
          if keyword_set(dft_ian) then save, file = file_struct.uvf_savefile[i], u_lambda_vals, v_lambda_vals, data_cube, freq_mask, uvf_git_hash $
          else save, file = file_struct.uvf_savefile[i], kx_rad_vals, ky_rad_vals, data_cube, freq_mask, uvf_git_hash
          undefine, data_cube
        endif else begin
        
          if healpix then begin
            pixel_nums1 = getvar_savefile(file_struct.pixelfile[0], file_struct.pixelvar[0])
            
            if nfiles eq 2 then begin
              ;; check that they have the same set of healpix pixels
              pixel_nums2 = getvar_savefile(file_struct.pixelfile[1], file_struct.pixelvar[0])
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
            
            ;; then rotate to make as rectangular as possible
            pred_angle = healpix_rot(new_pix_vec[*,0], new_pix_vec[*,1])
            
            x_rot = new_pix_vec[*,0] * cos(pred_angle) - new_pix_vec[*,1] * sin(pred_angle)
            y_rot = new_pix_vec[*,0] * sin(pred_angle) + new_pix_vec[*,1] * cos(pred_angle)
            
            
          endif else begin
            ;; gridded image to dft to parallel Healpix computation
            pix_size_rad = abs(file_struct.degpix) * !pi / 180d
            
            data_size = getvar_savefile(file_struct.datafile[i], file_struct.datavar, /return_size)
            
            dims = data_size[1:data_size(0)]
            
            x_vec = (findgen(dims[0]) - dims[0]/2.) * pix_size_rad
            y_vec = (findgen(dims[1]) - dims[1]/2.) * pix_size_rad
            
            x_rot = fltarr(dims[0]*dims[1])
            y_rot = fltarr(dims[0]*dims[1])
            x_rot = reform(rebin(x_vec, dims[0], dims[1], /sample), dims[0]*dims[1])
            y_rot = reform(rebin(reform(y_vec, 1, dims[1]), dims[0], dims[1], /sample), dims[0]*dims[1])
            
          endelse
          
          ;; figure out k values to calculate dft
          uv_cellsize_m = 5 ;; based on calculations of beam FWHM by Aaron
          if keyword_set(dft_ian) then begin
          
            ;;delta_u_lambda = file_struct.kpix
            if n_elements(delta_uv_lambda) gt 0 then delta_u_lambda = delta_uv_lambda $
            else delta_u_lambda = uv_cellsize_m * mean(frequencies*1e6) / (3e8)
            
            ;; go a little beyond max_baseline to account for expansion due to w projection
            ;;max_u_lambda = (file_struct.max_baseline_lambda) * 1.1
            ;; use kspan of Ian's cubes
            if tag_exist(file_struct, 'kspan') then begin
              max_kperp_rad = min([file_struct.kspan/2.,file_struct.max_baseline_lambda])* (2.*!pi)
            endif else max_kperp_rad = min([file_struct.max_baseline_lambda])* (2.*!pi)
            
            if n_elements(max_uv_lambda) gt 0 then max_kperp_rad = min([max_kperp_rad, max_uv_lambda])
            
            if keyword_set(cut_image) then begin
              ;; limit field of view to match calculated k-modes
              xy_len = 1/delta_u_lambda
              
              ;; image may be smaller than expected, may need to adjust delta_kperp_rad
              
              if healpix then begin
                ;; if file_struct.kpix = 1./sqrt(2.) or less (ie touching horizon) then use lims = +/- 1/(sqrt(2))
                if file_struct.kpix lt 1./sqrt(2.) then begin
                  limits = [-1., -1., 1. ,1.] / sqrt(2.)
                endif else begin
                  ;; get surrounding pixels
                  dists = sqrt((pix_center_vec[*,0]-vec_mid[0])^2d + (pix_center_vec[*,1]-vec_mid[1])^2d + (pix_center_vec[*,2]-vec_mid[2])^2d)
                  radius = max(dists)*1.1
                  query_disc, file_struct.nside, vec_mid, radius, listpix, nlist, /inc
                  min_pix = min([pixel_nums1, listpix])
                  wh2 = where(histogram(listpix, min=min_pix) eq 0 and histogram(pixel_nums1, min=min_pix) gt 0, count2)
                  while count2 gt 0 do begin
                    radius = radius*1.1
                    query_disc, nside, vec_mid, radius, listpix, nlist, /inc
                    min_pix = min([pixel_nums1, listpix])
                    wh2 = where(histogram(listpix, min=min_pix) eq 0 and histogram(pixel_nums1, min=min_pix) gt 0, count2)
                  endwhile
                  
                  ;; remove pixels from listpix that are in my image -- only want nearby pixels not in my image
                  min_pix = min([pixel_nums1, listpix])
                  max_pix = max([pixel_nums1, listpix])
                  wh = where(histogram(listpix, min = min_pix, max = max_pix) gt 0 and histogram(pixel_nums1, min = min_pix, max = max_pix) eq 0, count_out)
                  if count_out gt 0 then outside_pix = wh + min_pix else stop
                  pix2vec_ring, file_struct.nside, outside_pix, out_center_vec
                  new_out_vec = rot_matrix ## out_center_vec
                  x_out_rot = new_out_vec[*,0] * cos(pred_angle) - new_out_vec[*,1] * sin(pred_angle)
                  y_out_rot = new_out_vec[*,0] * sin(pred_angle) + new_out_vec[*,1] * cos(pred_angle)
                  limits = healpix_limits(x_rot, y_rot, x_out_rot, y_out_rot)
                endelse
                
              endif else begin
                limits = [min(x_vec), min(y_vec), max(x_vec), max(y_vec)]
              endelse
              
              image_len = min([limits[2]-limits[0],limits[3]-limits[1]])
              if image_len lt xy_len then begin
                print, 'Image FoV is smaller than expected, increasing delta kperp to match image FoV'
                delta_kperp_rad = 1./image_len
                
                wh_close = where(x_rot le image_len/2. and x_rot ge -1.*image_len/2. and $
                  y_rot le image_len/2. and y_rot ge -1.*image_len/2., count_close, $
                  ncomplement = count_far)
              endif else wh_close = where(x_rot le xy_len/2. and x_rot ge -1.*xy_len/2. and $
                y_rot le xy_len/2. and y_rot ge -1.*xy_len/2., count_close, $
                ncomplement = count_far)
                
            endif else count_far = 0
            
            n_u = round(max_u_lambda / delta_u_lambda) * 2 + 1
            u_lambda_vals = (findgen(n_u) - (n_u-1)/2) * delta_u_lambda
            
            ;; need to cut uvf cubes in half because image is real -- we'll cut in v
            ;; drop the unused half before the DFT to save time
            v_lambda_vals = u_lambda_vals[n_u/2:n_u-1]
            
          endif else begin
            ;;delta_kperp_rad = file_struct.kpix * z_mpc_mean / kperp_lambda_conv
          
            if n_elements(delta_uv_lambda) gt 0 then delta_kperp_rad = delta_uv_lambda * (2.*!pi) else begin
              if n_elements(freq_ch_range) gt 0 then delta_kperp_rad = uv_cellsize_m * mean(file_struct.frequencies*1e6) * (2.*!pi) / (3e8) $
              else delta_kperp_rad = uv_cellsize_m * mean(frequencies*1e6) * (2.*!pi) / (3e8)
            endelse
            ;; go a little beyond max_baseline to account for expansion due to w projection
            ;; max_kperp_rad = (file_struct.max_baseline_lambda/kperp_lambda_conv) * z_mpc_mean * 1.1
            ;; use kspan of Ian's cubes
            if tag_exist(file_struct, 'kspan') then begin
              max_kperp_rad = min([file_struct.kspan/2.,file_struct.max_baseline_lambda])* (2.*!pi)
            endif else max_kperp_rad = min([file_struct.max_baseline_lambda])* (2.*!pi)
            
            if n_elements(max_uv_lambda) gt 0 then max_kperp_rad = min([max_kperp_rad, max_uv_lambda * (2.*!pi)])
            
            if keyword_set(cut_image) then begin
              ;; limit field of view to match calculated k-modes
              xy_len = 2*!pi/delta_kperp_rad
              
              ;; image may be smaller than expected, may need to adjust delta_kperp_rad
              
              if healpix then begin
                if file_struct.kpix lt 1./sqrt(2.) then begin
                  limits = [-1., -1., 1. ,1.] / sqrt(2.)
                endif else begin
                  ;; get surrounding pixels
                  dists = sqrt((pix_center_vec[*,0]-vec_mid[0])^2d + (pix_center_vec[*,1]-vec_mid[1])^2d + (pix_center_vec[*,2]-vec_mid[2])^2d)
                  radius = max(dists)*1.1
                  query_disc, file_struct.nside, vec_mid, radius, listpix, nlist, /inc
                  min_pix = min([pixel_nums1, listpix])
                  wh2 = where(histogram(listpix, min=min_pix) eq 0 and histogram(pixel_nums1, min=min_pix) gt 0, count2)
                  while count2 gt 0 do begin
                    radius = radius*1.1
                    query_disc, file_struct.nside, vec_mid, radius, listpix, nlist, /inc
                    min_pix = min([pixel_nums1, listpix])
                    wh2 = where(histogram(listpix, min=min_pix) eq 0 and histogram(pixel_nums1, min=min_pix) gt 0, count2)
                  endwhile
                  
                  ;; remove pixels from listpix that are in my image -- only want nearby pixels not in my image
                  min_pix = min([pixel_nums1, listpix])
                  max_pix = max([pixel_nums1, listpix])
                  wh = where(histogram(listpix, min = min_pix, max = max_pix) gt 0 and histogram(pixel_nums1, min = min_pix, max = max_pix) eq 0, count_out)
                  if count_out gt 0 then outside_pix = wh + min_pix else stop
                  pix2vec_ring, file_struct.nside, outside_pix, out_center_vec
                  new_out_vec = rot_matrix ## out_center_vec
                  x_out_rot = new_out_vec[*,0] * cos(pred_angle) - new_out_vec[*,1] * sin(pred_angle)
                  y_out_rot = new_out_vec[*,0] * sin(pred_angle) + new_out_vec[*,1] * cos(pred_angle)
                  limits = healpix_limits(x_rot, y_rot, x_out_rot, y_out_rot)
                endelse
              endif else begin
                limits = [min(x_vec), min(y_vec), max(x_vec), max(y_vec)]
              endelse
              
              image_len = min([limits[2]-limits[0],limits[3]-limits[1]])
              if image_len lt xy_len then begin
                print, 'Image FoV is smaller than expected, increasing delta kperp to match image FoV'
                delta_kperp_rad = 2*!pi/image_len
                
                wh_close = where(x_rot le image_len/2. and x_rot ge -1.*image_len/2. and $
                  y_rot le image_len/2. and y_rot ge -1.*image_len/2., count_close, $
                  ncomplement = count_far)
              endif else wh_close = where(x_rot le xy_len/2. and x_rot ge -1.*xy_len/2. and $
                y_rot le xy_len/2. and y_rot ge -1.*xy_len/2., count_close, $
                ncomplement = count_far)
                
            endif else count_far = 0
            
            if n_elements(wh_close) eq 0 and count_far eq 0 then wh_close = lindgen(n_elements(x_rot))
            
            ;; for deciding on pixel sets:
            ;            consv_delta_kperp_rad = 4.5* mean(frequencies*1e6) * z_mpc_mean / (3e8 * kperp_lambda_conv) ;use 4.5m to be conservative
            ;            consv_xy_len = 2*!pi/consv_delta_kperp_rad
            ;            radius = consv_xy_len/2.*sqrt(2)*1.1
            ;            query_disc, file_struct.nside, vec_mid, radius, listpix, nlist, /inc
            ;            pix2vec_ring, file_struct.nside, listpix, list_center_vec
            ;            new_list_vec = rot_matrix ## list_center_vec
            ;            x_list_rot = new_list_vec[*,0] * cos(pred_angle) - new_list_vec[*,1] * sin(pred_angle)
            ;            y_list_rot = new_list_vec[*,0] * sin(pred_angle) + new_list_vec[*,1] * cos(pred_angle)
            ;            cgplot, x_list_rot, y_list_rot, psym=3
            ;            consv_lims = [-1*consv_xy_len/2., -1*consv_xy_len/2., consv_xy_len/2., consv_xy_len/2.]
            ;            cgpolygon, reform(rebin(consv_lims[[0,2]], 2,2),4), reform(rebin(reform(consv_lims[[1,3]],1,2), 2,2),4), color='aqua'
            ;            wh_listpix_close = where(x_list_rot ge consv_lims[0] and x_list_rot le consv_lims[2] and $
            ;              y_list_rot ge consv_lims[1] and y_list_rot le consv_lims[3], count_list_close)
            ;            hpx_inds = listpix[wh_listpix_close]
            ;            nside = file_struct.nside
            ;            stop
            ;            save, file='/Users/bryna/Documents/Physics/FHD/Observations/EoR1_low_healpix_inds.idlsave', nside, hpx_inds
            
            
            n_kperp = round(max_kperp_rad / delta_kperp_rad) * 2 + 1
            kx_rad_vals = (findgen(n_kperp) - (n_kperp-1)/2) * delta_kperp_rad
            
            ;; need to cut uvf cubes in half because image is real -- we'll cut in v
            ;; drop the unused half before the DFT to save time
            ky_rad_vals = kx_rad_vals[n_kperp/2:n_kperp-1]
            
          endelse
          
          ;; get beam if needed
          if (test_beam eq 0 or keyword_set(refresh_beam)) and tag_exist(file_struct, 'beam_savefile') then begin
            arr = getvar_savefile(file_struct.beamfile[i], file_struct.beamvar)
            if count_far ne 0 then arr = arr[wh_close, *]
            if n_elements(freq_ch_range) ne 0 then arr = arr[*, min(freq_ch_range):max(freq_ch_range)]
            if n_elements(freq_flags) ne 0 then arr = arr * rebin(reform(freq_mask, 1, n_elements(file_struct.frequencies)), size(arr, /dimension), /sample)
            
            pixels = pixel_nums1[wh_close]
            
            if max(arr) le 1.1 then begin
              ;; beam is peak normalized to 1
              temp = arr * rebin(reform(n_vis_freq[i, *], 1, n_freq), count_close, n_freq, /sample)
            endif else if max(arr) le file_struct.n_obs[i]*1.1 then begin
              ;; beam is peak normalized to 1 for each obs, then summed over obs so peak is ~ n_obs
              temp = (arr/file_struct.n_obs[i]) * rebin(reform(n_vis_freq[i, *], 1, n_freq), count_close, n_freq, /sample)
            endif else begin
              ;; beam is peak normalized to 1 then multiplied by n_vis_freq for each obs & summed
              temp = arr
            endelse
            
            avg_beam = total(temp, 2) / total(n_vis_freq[i, *])
            
            nside = file_struct.nside
            
            git, repo_path = ps_repository_dir(), result=beam_git_hash
            
            save, file=file_struct.beam_savefile[i], avg_beam, pixels, nside, beam_git_hash
            
          endif
          
          ;; do DFT.
          if test_uvf eq 0 or keyword_set(dft_refresh_data) then begin
            print, 'calculating DFT for ' + file_struct.datavar + ' in ' + file_struct.datafile[i]
            
            time0 = systime(1)
            arr = getvar_savefile(file_struct.datafile[i], file_struct.datavar)
            time1 = systime(1)
            
            if time1 - time0 gt 60 then print, 'Time to restore cube was ' + number_formatter((time1 - time0)/60.) + ' minutes'
            
            if size(arr,/type) eq 6 or size(arr,/type) eq 9 then $
              if max(abs(imaginary(arr))) eq 0 then arr = real_part(arr) else stop
            if not healpix then arr = reform(arr, dims[0]*dims[1], n_freq)
            
            if count_far ne 0 then arr = arr[wh_close, *]
            if n_elements(freq_ch_range) ne 0 then arr = arr[*, min(freq_ch_range):max(freq_ch_range)]
            if n_elements(freq_flags) ne 0 then arr = arr * rebin(reform(freq_mask, 1, n_elements(file_struct.frequencies)), size(arr, /dimension), /sample)
            
            if keyword_set(dft_ian) then begin
              transform = discrete_ft_2D_fast(x_rot[wh_close], y_rot[wh_close], arr, u_lambda_vals, v_lambda_vals, /exp2pi, $
                timing = ft_time, fchunk = dft_fchunk, no_progress = no_dft_progress)
            endif else begin
              transform = discrete_ft_2D_fast(x_rot[wh_close], y_rot[wh_close], arr, kx_rad_vals, ky_rad_vals, $
                timing = ft_time, fchunk = dft_fchunk, no_progress = no_dft_progress)
            endelse
            data_cube = temporary(transform)
            undefine, arr
            
            git, repo_path = ps_repository_dir(), result=uvf_git_hash
            
            if keyword_set(dft_ian) then save, file = file_struct.uvf_savefile[i], u_lambda_vals, v_lambda_vals, data_cube, freq_mask, uvf_git_hash $
            else save, file = file_struct.uvf_savefile[i], kx_rad_vals, ky_rad_vals, data_cube, freq_mask, uvf_git_hash
            undefine, data_cube
          endif
          
          if test_wt_uvf eq 0 or keyword_set(dft_refresh_weight) then begin
            print, 'calculating DFT for ' + file_struct.weightvar + ' in ' + file_struct.weightfile[i]
            
            time0 = systime(1)
            arr = getvar_savefile(file_struct.weightfile[i], file_struct.weightvar)
            time1 = systime(1)
            
            if time1 - time0 gt 60 then print, 'Time to restore cube was ' + number_formatter((time1 - time0)/60.) + ' minutes'
            
            if size(arr,/type) eq 6 or size(arr,/type) eq 9 then $
              if max(abs(imaginary(arr))) eq 0 then arr = real_part(arr) else stop
            if not healpix then arr = reform(arr, dims[0]*dims[1], n_freq)
            
            if count_far ne 0 then arr = arr[wh_close, *]
            if n_elements(freq_ch_range) ne 0 then arr = arr[*, min(freq_ch_range):max(freq_ch_range)]
            if n_elements(freq_flags) ne 0 then arr = arr * rebin(reform(freq_mask, 1, n_elements(file_struct.frequencies)), size(arr, /dimension), /sample)
            
            if keyword_set(dft_ian) then begin
              transform = discrete_ft_2D_fast(x_rot[wh_close], y_rot[wh_close], arr, u_lambda_vals, v_lambda_vals, /exp2pi, $
                timing = ft_time, fchunk = dft_fchunk, no_progress = no_dft_progress)
            endif else begin
              transform = discrete_ft_2D_fast(x_rot[wh_close], y_rot[wh_close], arr, kx_rad_vals, ky_rad_vals, $
                timing = ft_time, fchunk = dft_fchunk, no_progress = no_dft_progress)
            endelse
            
            weights_cube = temporary(transform)
            
            if not no_var then begin
              print, 'calculating DFT for ' + file_struct.variancevar + ' in ' + file_struct.variancefile[i]
              
              time0 = systime(1)
              arr = getvar_savefile(file_struct.variancefile[i], file_struct.variancevar)
              time1 = systime(1)
              
              if time1 - time0 gt 60 then print, 'Time to restore cube was ' + number_formatter((time1 - time0)/60.) + ' minutes'
              
              if size(arr,/type) eq 6 or size(arr,/type) eq 9 then $
                if max(abs(imaginary(arr))) eq 0 then arr = real_part(arr) else stop
              if not healpix then arr = reform(arr, dims[0]*dims[1], n_freq)
              
              if count_far ne 0 then arr = arr[wh_close, *]
              if n_elements(freq_ch_range) ne 0 then arr = arr[*, min(freq_ch_range):max(freq_ch_range)]
              if n_elements(freq_flags) ne 0 then arr = arr * rebin(reform(freq_mask, 1, n_elements(file_struct.frequencies)), size(arr, /dimension), /sample)
              
              if keyword_set(dft_ian) then begin
                transform = discrete_ft_2D_fast(x_rot[wh_close], y_rot[wh_close], arr, u_lambda_vals, v_lambda_vals, /exp2pi, $
                  timing = ft_time, fchunk = dft_fchunk, no_progress = no_dft_progress)
              endif else begin
                transform = discrete_ft_2D_fast(x_rot[wh_close], y_rot[wh_close], arr, kx_rad_vals, ky_rad_vals, $
                  timing = ft_time, fchunk = dft_fchunk, no_progress = no_dft_progress)
              endelse
              variance_cube = abs(temporary(transform)) ;; make variances real, positive definite (amplitude)
              undefine, arr
            endif
            
            
            git, repo_path = ps_repository_dir(), result=uvf_wt_git_hash
            
            if keyword_set(dft_ian) then $
              save, file = file_struct.uvf_weight_savefile[i], u_lambda_vals, v_lambda_vals, weights_cube, variance_cube, freq_mask, uvf_wt_git_hash else $
              save, file = file_struct.uvf_weight_savefile[i], kx_rad_vals, ky_rad_vals, weights_cube, variance_cube, freq_mask, uvf_wt_git_hash
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

  if healpix or not keyword_set(uvf_input) then begin
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
    
  endif
  
  if healpix or not keyword_set(uvf_input) then begin
    weights_cube1 = getvar_savefile(file_struct.uvf_weight_savefile[0], 'weights_cube')
    if nfiles eq 2 then weights_cube2 = getvar_savefile(file_struct.uvf_weight_savefile[1], 'weights_cube')
    
    ave_weights = total(total(abs(weights_cube1),2),1)/(n_kx*n_ky)
    if nfiles eq 2 then ave_weights = transpose([[ave_weights], [total(total(abs(weights_cube2),2),1)/(n_kx*n_ky)]])
    
    void = getvar_savefile(file_struct.uvf_weight_savefile[0], names = uvf_varnames)
    wh_hash = where(uvf_varnames eq 'uvf_wt_git_hash', count_hash)
    if count_hash gt 0 then begin
      uvf_wt_git_hashes = getvar_savefile(file_struct.uvf_weight_savefile[0], 'uvf_wt_git_hash')
      if nfiles eq 2 then uvf_wt_git_hashes = [uvf_wt_git_hashes, getvar_savefile(file_struct.uvf_weight_savefile[1], 'uvf_wt_git_hash')]
    endif else uvf_wt_git_hashes = strarr(nfiles)
    
    if min(ky_mpc) lt 0 then begin
      ;; negative ky values haven't been cut yet
      ;; need to cut uvf cubes in half because image is real -- we'll cut negative ky
      weights_cube1 = weights_cube1[*, n_ky/2:n_ky-1,*]
      if nfiles eq 2 then weights_cube2 = weights_cube2[*, n_ky/2:n_ky-1,*]
    endif
    
  endif else begin
    ;; uvf_input
  
    weights_cube1 = getvar_savefile(file_struct.weightfile[0], file_struct.weightvar)
    if nfiles eq 2 then weights_cube2 = getvar_savefile(file_struct.weightfile[1], file_struct.weightvar)
    
    uvf_wt_git_hashes = strarr(nfiles)
    
    wt_size = size(weights_cube1)
    if wt_size[n_elements(wt_size)-2] eq 10 then begin
      ;; weights cube is a pointer
      dims2 = size(*weights_cube1[0], /dimension)
      temp = complex(fltarr([dims2, n_freq]))
      if nfiles eq 2 then temp2 = complex(fltarr([dims2, n_freq]))
      for i = 0, n_freq-1 do begin
        temp[*,*,i] = *weights_cube1[file_struct.pol_index, i]
        if nfiles eq 2 then temp2[*,*,i] = *weights_cube2[file_struct.pol_index, i]
      endfor
      undefine_fhd, weights_cube1, weights_cube2
      
      weights_cube1 = temporary(temp)
      if nfiles eq 2 then weights_cube2 = temporary(temp2)
    endif else dims2 = size(weights_cube1, /dimension)
    
    n_kx = dims2[0]
    if abs(file_struct.kpix-1/(n_kx[0] * (abs(file_struct.degpix) * !pi / 180d)))/file_struct.kpix gt 1e-4 then stop
    kx_mpc_delta = (2.*!pi)*file_struct.kpix / z_mpc_mean
    kx_mpc = (dindgen(n_kx)-n_kx/2) * kx_mpc_delta
    
    n_ky = dims2[1]
    ky_mpc_delta = (2.*!pi)*file_struct.kpix / z_mpc_mean
    ky_mpc = (dindgen(n_ky)-n_ky/2) * kx_mpc_delta
    
    ave_weights = total(total(abs(weights_cube1),2),1)/(n_kx*n_ky)
    if nfiles eq 2 then ave_weights = transpose([[ave_weights], [total(total(abs(weights_cube2),2),1)/(n_kx*n_ky)]])
    
    pix_area_rad = (!dtor*file_struct.degpix)^2.
    pix_area_mpc = pix_area_rad * z_mpc_mean^2.
    
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
      nkx_new = floor(n_kx / uv_avg)
      temp = complex(fltarr(nkx_new, n_ky, n_freq))
      if nfiles eq 2 then temp2 = complex(fltarr(nkx_new, n_ky, n_freq))
      temp_kx = fltarr(nkx_new)
      for i=0, nkx_new-1 do begin
        temp[i,*,*] = total(weights_cube1[i*uv_avg:(i+1)*uv_avg-1,*,*], 1) / uv_avg
        if nfiles eq 2 then temp2[i,*,*] = total(weights_cube2[i*uv_avg:(i+1)*uv_avg-1,*,*], 1) / uv_avg
        temp_kx[i] = total(kx_mpc[i*uv_avg:(i+1)*uv_avg-1]) / uv_avg
      endfor
      
      nky_new = floor(n_ky / uv_avg)
      temp3 = complex(fltarr(nkx_new, nky_new, n_freq))
      if nfiles eq 2 then temp4 = complex(fltarr(nkx_new, nky_new, n_freq))
      temp_ky = fltarr(nky_new)
      for i=0, nky_new-1 do begin
        temp3[*,i,*] = total(temp[*,i*uv_avg:(i+1)*uv_avg-1,*], 2) / uv_avg
        if nfiles eq 2 then temp4[*,i,*] = total(temp2[*,i*uv_avg:(i+1)*uv_avg-1,*], 2) / uv_avg
        temp_ky[i] = total(ky_mpc[i*uv_avg:(i+1)*uv_avg-1]) / uv_avg
      endfor
      undefine, temp, temp2
      
      ;; averging reduces the value of total(weights) ~ n_vis needed for the window int calculation
      n_vis = n_vis/(uv_avg)^2.
      
      weights_cube1 = temporary(temp3)
      if nfiles eq 2 then weights_cube2 = temporary(temp4)
      kx_mpc = temporary(temp_kx)
      ky_mpc = temporary(temp_ky)
      n_kx = nkx_new
      n_ky = nky_new
    endif
    
    if keyword_set(uv_img_clip) then begin
      kx_mpc_delta_old = kx_mpc_delta
      ky_mpc_delta_old = ky_mpc_delta
      temp = shift(fft(fft(shift(weights_cube1,dims2[0]/2,dims2[1]/2,0), dimension=1),dimension=2),dims2[0]/2,dims2[1]/2,0)
      temp = temp[(dims2[0]/2)-(dims2[0]/uv_img_clip)/2:(dims2[0]/2)+(dims2[0]/uv_img_clip)/2-1, *, *]
      temp = temp[*, (dims2[1]/2)-(dims2[1]/uv_img_clip)/2:(dims2[1]/2)+(dims2[1]/uv_img_clip)/2-1, *]
      if nfiles eq 2 then begin
        temp2 = shift(fft(fft(shift(weights_cube2,dims2[0]/2,dims2[1]/2,0), dimension=1), dimension=2),dims2[0]/2,dims2[1]/2,0)
        temp2 = temp2[(dims2[0]/2)-(dims2[0]/uv_img_clip)/2:(dims2[0]/2)+(dims2[0]/uv_img_clip)/2-1, *, *]
        temp2 = temp2[*, (dims2[1]/2)-(dims2[1]/uv_img_clip)/2:(dims2[1]/2)+(dims2[1]/uv_img_clip)/2-1, *]
      endif
      temp_dims = size(temp, /dimension)
      n_kx = temp_dims[0]
      n_ky = temp_dims[1]
      kx_mpc_delta = kx_mpc_delta * uv_img_clip
      kx_mpc = (dindgen(n_kx)-n_kx/2) * kx_mpc_delta
      ky_mpc_delta = ky_mpc_delta * uv_img_clip
      ky_mpc = (dindgen(n_kx)-n_kx/2) * ky_mpc_delta
      
      temp = shift(fft(fft(shift(temp,n_kx/2,n_ky/2,0), dimension=1, /inverse), dimension=2, /inverse),n_kx/2,n_ky/2,0)
      if nfiles eq 2 then temp2 = shift(fft(fft(shift(temp2,n_kx/2,n_ky/2,0), dimension=1, /inverse), dimension=2, /inverse),n_kx/2,n_ky/2,0)
      
      weights_cube1 = temp
      if nfiles eq 2 then weights_cube2 = temp2
    endif
    
    if n_elements(freq_ch_range) ne 0 then begin
      weights_cube1 = weights_cube1[*, *, min(freq_ch_range):max(freq_ch_range)]
      if nfiles eq 2 then weights_cube2 = weights_cube2[*, *, min(freq_ch_range):max(freq_ch_range)]
    endif
    if n_elements(freq_flags) ne 0 then begin
      flag_arr = rebin(reform(freq_mask, 1, 1, n_elements(file_struct.frequencies)), size(weights_cube1,/dimension), /sample)
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
  
  
  if not no_var then begin
    if healpix or not keyword_set(uvf_input) then begin
      variance_cube1 = getvar_savefile(file_struct.uvf_weight_savefile[0], 'variance_cube')
      if nfiles eq 2 then variance_cube2 = getvar_savefile(file_struct.uvf_weight_savefile[1], 'variance_cube')
      
      if min(ky_mpc) lt 0 then begin
        ;; calculate integral of window function before cut for comparison
        if nfiles eq 2 then window_int_orig = [total(variance_cube1)*pix_area_rad/n_vis[0], $
          total(variance_cube2)*pix_area_rad/n_vis[1]] $
        else window_int_orig = total(variance_cube1)*pix_area_rad/n_vis[0]
        
        ;; negative ky values haven't been cut yet
        ;; need to cut uvf cubes in half because image is real -- we'll cut negative ky
        variance_cube1 = variance_cube1[*, n_ky/2:n_ky-1,*]
        if nfiles eq 2 then variance_cube2 = variance_cube2[*, n_ky/2:n_ky-1,*]
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
          if nfiles eq 2 then beam_git_hashes = [beam_git_hashes, getvar_savefile(file_struct.beam_savefile[1], 'beam_git_hash')]
        endif else beam_git_hashes = strarr(nfiles)
        
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
      
      
    endif else begin
      ;; uvf_input
      variance_cube1 = getvar_savefile(file_struct.variancefile[0], file_struct.variancevar)
      if nfiles eq 2 then variance_cube2 = getvar_savefile(file_struct.variancefile[1], file_struct.variancevar)
      
      var_size = size(variance_cube1)
      if var_size[n_elements(var_size)-2] eq 10 then begin
        ;; variance cube is a pointer
        dims2 = size(*variance_cube1[0], /dimension)
        temp = complex(fltarr([dims2, n_freq]))
        if nfiles eq 2 then temp2 = complex(fltarr([dims2, n_freq]))
        for i = 0, n_freq-1 do begin
          temp[*,*,i] = *variance_cube1[file_struct.pol_index, i]
          if nfiles eq 2 then temp2[*,*,i] = *variance_cube2[file_struct.pol_index, i]
        endfor
        undefine_fhd, variance_cube1, variance_cube2
        
        variance_cube1 = temporary(temp)
        if nfiles eq 2 then variance_cube2 = temporary(temp2)
      endif
      
      if keyword_set(uv_avg) then begin
        temp = complex(fltarr(nkx_new, dims2[1], n_freq))
        if nfiles eq 2 then temp2 = complex(fltarr(nkx_new, dims2[1], n_freq))
        for i=0, nkx_new-1 do begin
          temp[i,*,*] = total(variance_cube1[i*uv_avg:(i+1)*uv_avg-1,*,*], 1) / uv_avg
          if nfiles eq 2 then temp2[i,*,*] = total(variance_cube2[i*uv_avg:(i+1)*uv_avg-1,*,*], 1) / uv_avg
        endfor
        
        temp3 = complex(fltarr(nkx_new, nky_new, n_freq))
        if nfiles eq 2 then temp4 = complex(fltarr(nkx_new, nky_new, n_freq))
        for i=0, nky_new-1 do begin
          temp3[*,i,*] = total(temp[*,i*uv_avg:(i+1)*uv_avg-1,*], 2) / uv_avg
          if nfiles eq 2 then temp4[*,i,*] = total(temp2[*,i*uv_avg:(i+1)*uv_avg-1,*], 2) / uv_avg
        endfor
        undefine, temp, temp2
        
        variance_cube1 = temporary(temp3)
        if nfiles eq 2 then variance_cube2 = temporary(temp4)
      endif
      
      if keyword_set(uv_img_clip) then begin
        temp = shift(fft(fft(shift(variance_cube1,dims2[0]/2,dims2[1]/2,0), dimension=1), dimension=2),dims2[0]/2,dims2[1]/2,0)
        temp = temp[(dims2[0]/2)-(dims2[0]/uv_img_clip)/2:(dims2[0]/2)+(dims2[0]/uv_img_clip)/2-1, *, *]
        temp = temp[*, (dims2[1]/2)-(dims2[1]/uv_img_clip)/2:(dims2[1]/2)+(dims2[1]/uv_img_clip)/2-1, *]
        temp = shift(fft(fft(shift(temp,n_kx/2,n_ky/2,0), dimension=1, /inverse), dimension=2, /inverse),n_kx/2,n_ky/2,0)
        if nfiles eq 2 then begin
          temp2 = shift(fft(fft(shift(variance_cube2,dims2[0]/2,dims2[1]/2,0), dimension=1), dimension=2),dims2[0]/2,dims2[1]/2,0)
          temp2 = temp2[(dims2[0]/2)-(dims2[0]/uv_img_clip)/2:(dims2[0]/2)+(dims2[0]/uv_img_clip)/2-1, *, *]
          temp2 = temp2[*, (dims2[1]/2)-(dims2[1]/uv_img_clip)/2:(dims2[1]/2)+(dims2[1]/uv_img_clip)/2-1, *]
          temp2 = shift(fft(fft(shift(temp2,n_kx/2,n_ky/2,0), dimension=1, /inverse), dimension=2, /inverse),n_kx/2,n_ky/2,0)
        endif
        
        variance_cube1 = temp
        if nfiles eq 2 then variance_cube2 = temp2
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
        flag_arr = rebin(reform(freq_mask, 1, 1, n_elements(file_struct.frequencies)), size(variance_cube1,/dimension), /sample)
        variance_cube1 = variance_cube1 * flag_arr
        if nfiles eq 2 then variance_cube2 = variance_cube2 * flag_arr
      endif
      
      ;; Also need to drop 1/2 of ky=0 line -- drop ky=0, kx<0
      variance_cube1[0:n_kx/2-1, 0, *] = 0
      if nfiles eq 2 then variance_cube2[0:n_kx/2-1, 0, *] = 0
      
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
      if nfiles eq 2 then uvf_git_hashes = [uvf_git_hashes, getvar_savefile(file_struct.uvf_savefile[1], 'uvf_git_hash')]
    endif else uvf_git_hashes = strarr(nfiles)
    
    if min(ky_mpc) lt 0 then begin
      ;; negative ky values haven't been cut yet
      ;; need to cut uvf cubes in half because image is real -- we'll cut negative ky
      data_cube1 = data_cube1[*, n_ky/2:n_ky-1,*]
      if nfiles eq 2 then data_cube2 = data_cube2[*, n_ky/2:n_ky-1,*]
      
      ky_mpc = ky_mpc[n_ky/2:n_ky-1]
      n_ky = n_elements(ky_mpc)
    endif
    
  endif else begin
    if datavar eq '' then begin
      ;; working with a 'derived' cube (ie residual cube) that is constructed from other cubes
      ;; need to cut uvf cubes in half because image is real -- we'll cut negative ky
      dirty_cube1 = getvar_savefile(input_uvf_files[0,0], input_uvf_varname[0,0])
      dirty_cube1 = dirty_cube1[*, n_ky/2:n_ky-1,*]
      model_cube1 = getvar_savefile(input_uvf_files[0,1], input_uvf_varname[0,1])
      model_cube1 = model_cube1[*, n_ky/2:n_ky-1,*]
      data_cube1 = temporary(dirty_cube1) - temporary(model_cube1)
      
      if nfiles eq 2 then begin
        dirty_cube2 = getvar_savefile(input_uvf_files[1,0], input_uvf_varname[1,0])
        dirty_cube2 = dirty_cube2[*, n_ky/2:n_ky-1,*]
        model_cube2 = getvar_savefile(input_uvf_files[1,1], input_uvf_varname[1,1])
        model_cube2 = model_cube2[*, n_ky/2:n_ky-1,*]
        data_cube2 = temporary(dirty_cube2) - temporary(model_cube2)
      endif
      
      void = getvar_savefile(file_struct.input_uvf_files[0,0], names = uvf_varnames)
      wh_hash = where(uvf_varnames eq 'uvf_git_hash', count_hash)
      if count_hash gt 0 then begin
        uvf_git_hashes = getvar_savefile(file_struct.input_uvf_files[0,0], 'uvf_git_hash')
        if nfiles eq 2 then uvf_git_hashes = [uvf_git_hashes, getvar_savefile(file_struct.input_uvf_files[0,1], 'uvf_git_hash')]
      endif else uvf_git_hashes = strarr(nfiles)
      
      
      if n_elements(freq_ch_range) ne 0 then begin
        data_cube1 = data_cube1[*, *, min(freq_ch_range):max(freq_ch_range)]
        if nfiles eq 2 then data_cube2 = data_cube2[*, *, min(freq_ch_range):max(freq_ch_range)]
      endif
      if n_elements(freq_flags) ne 0 then begin
        flag_arr = rebin(reform(freq_mask, 1, 1, n_elements(file_struct.frequencies)), size(data_cube1,/dimension), /sample)
        data_cube1 = data_cube1 * flag_arr
        if nfiles eq 2 then data_cube2 = data_cube2 * flag_arr
      endif
      
    endif else begin
      ;; uvf_input
      data_cube1 = getvar_savefile(file_struct.datafile[0], file_struct.datavar)
      if nfiles eq 2 then data_cube2 = getvar_savefile(file_struct.datafile[1], file_struct.datavar)
      
      uvf_git_hashes = strarr(nfiles)
      
      data_size = size(data_cube1)
      if data_size[n_elements(data_size)-2] eq 10 then begin
        ;; weights cube is a pointer
        dims2 = size(*data_cube1[0], /dimension)
        temp = complex(fltarr([dims2, n_freq]))
        if nfiles eq 2 then temp2 = complex(fltarr([dims2, n_freq]))
        for i = 0, n_freq-1 do begin
          temp[*,*,i] = *data_cube1[file_struct.pol_index, i]
          if nfiles eq 2 then temp2[*,*,i] = *data_cube2[file_struct.pol_index, i]
        endfor
        undefine_fhd, data_cube1, data_cube2
        
        data_cube1 = temporary(temp)
        if nfiles eq 2 then data_cube2 = temporary(temp2)
      endif
      
      if keyword_set(uv_avg) then begin
        temp = complex(fltarr(nkx_new, dims2[1], n_freq))
        if nfiles eq 2 then temp2 = complex(fltarr(nkx_new, dims2[1], n_freq))
        for i=0, nkx_new-1 do begin
          temp[i,*,*] = total(data_cube1[i*uv_avg:(i+1)*uv_avg-1,*,*], 1) / uv_avg
          if nfiles eq 2 then temp2[i,*,*] = total(data_cube2[i*uv_avg:(i+1)*uv_avg-1,*,*], 1) / uv_avg
        endfor
        
        temp3 = complex(fltarr(nkx_new, nky_new, n_freq))
        if nfiles eq 2 then temp4 = complex(fltarr(nkx_new, nky_new, n_freq))
        for i=0, nky_new-1 do begin
          temp3[*,i,*] = total(temp[*,i*uv_avg:(i+1)*uv_avg-1,*], 2) / uv_avg
          if nfiles eq 2 then temp4[*,i,*] = total(temp2[*,i*uv_avg:(i+1)*uv_avg-1,*], 2) / uv_avg
        endfor
        undefine, temp, temp2
        
        data_cube1 = temporary(temp3)
        if nfiles eq 2 then data_cube2 = temporary(temp4)
      endif
      
      if keyword_set(uv_img_clip) then begin
        temp = shift(fft(fft(shift(data_cube1,dims2[0]/2,dims2[1]/2,0), dimension=1), dimension=2),dims2[0]/2,dims2[1]/2,0)
        temp = temp[(dims2[0]/2)-(dims2[0]/uv_img_clip)/2:(dims2[0]/2)+(dims2[0]/uv_img_clip)/2-1, *, *]
        temp = temp[*, (dims2[1]/2)-(dims2[1]/uv_img_clip)/2:(dims2[1]/2)+(dims2[1]/uv_img_clip)/2-1, *]
        temp = shift(fft(fft(shift(temp,n_kx/2,n_ky/2,0), dimension=1, /inverse), dimension=2, /inverse),n_kx/2,n_ky/2,0)
        if nfiles eq 2 then begin
          temp2 = shift(fft(fft(shift(data_cube2,dims2[0]/2,dims2[1]/2,0), dimension=1), dimension=2),dims2[0]/2,dims2[1]/2,0)
          temp2 = temp2[(dims2[0]/2)-(dims2[0]/uv_img_clip)/2:(dims2[0]/2)+(dims2[0]/uv_img_clip)/2-1, *, *]
          temp2 = temp2[*, (dims2[1]/2)-(dims2[1]/uv_img_clip)/2:(dims2[1]/2)+(dims2[1]/uv_img_clip)/2-1, *]
          temp2 = shift(fft(fft(shift(temp2,n_kx/2,n_ky/2,0), dimension=1, /inverse), dimension=2, /inverse),n_kx/2,n_ky/2,0)
        endif
        
        data_cube1 = temp
        if nfiles eq 2 then data_cube2 = temp2
      endif
      
      ;; need to cut uvf cubes in half because image is real -- we'll cut negative ky
      data_cube1 = data_cube1[*, n_ky/2:n_ky-1,*]
      if nfiles eq 2 then data_cube2 = data_cube2[*, n_ky/2:n_ky-1,*]
      
      if n_elements(freq_ch_range) ne 0 then begin
        data_cube1 = data_cube1[*, *, min(freq_ch_range):max(freq_ch_range)]
        if nfiles eq 2 then data_cube2 = data_cube2[*, *, min(freq_ch_range):max(freq_ch_range)]
      endif
      if n_elements(freq_flags) ne 0 then begin
        flag_arr = rebin(reform(freq_mask, 1, 1, n_elements(file_struct.frequencies)), size(data_cube1,/dimension), /sample)
        data_cube1 = data_cube1 * flag_arr
        if nfiles eq 2 then data_cube2 = data_cube2 * flag_arr
      endif
      
    endelse
    
    ky_mpc = ky_mpc[n_ky/2:n_ky-1]
    n_ky = n_elements(ky_mpc)
  endelse
  
  ;; Also need to drop 1/2 of ky=0 line -- drop ky=0, kx<0
  data_cube1[0:n_kx/2-1, 0, *] = 0
  if nfiles eq 2 then data_cube2[0:n_kx/2-1, 0, *] = 0
  
  ;; save some slices of the raw data cube (before dividing by weights) & weights
  for i=0, nfiles-1 do begin
    if i eq 0 then data_cube = data_cube1 else data_cube = data_cube2
    if i eq 0 then weights_cube = weights_cube1 else weights_cube = weights_cube2
    
    uf_tot = total(total(abs(weights_cube),3),1)
    wh_uf_n0 = where(uf_tot gt 0, count_uf_n0)
    if count_uf_n0 eq 0 then stop
    min_dist_uf_n0 = min(wh_uf_n0, min_loc)
    uf_slice_ind = wh_uf_n0[min_loc]
    
    uf_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 1, $
      slice_inds = uf_slice_ind, slice_savefile = file_struct.uf_raw_savefile[i])
      
    uf_weight_slice = uvf_slice(weights_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 1, $
      slice_inds = uf_slice_ind, slice_savefile = file_struct.uf_weight_savefile[i])
      
      
    vf_tot = total(total(abs(weights_cube),3),2)
    wh_vf_n0 = where(vf_tot gt 0, count_vf_n0)
    if count_vf_n0 eq 0 then stop
    min_dist_vf_n0 = min(abs(n_kx/2-wh_vf_n0), min_loc)
    vf_slice_ind = wh_vf_n0[min_loc]
    
    vf_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 0, $
      slice_inds = vf_slice_ind, slice_savefile = file_struct.vf_raw_savefile[i])
      
    vf_weight_slice = uvf_slice(weights_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 0, $
      slice_inds = vf_slice_ind, slice_savefile = file_struct.vf_weight_savefile[i])
      
    if max(abs(vf_slice)) eq 0 then stop
    
    uv_tot = total(total(abs(weights_cube),2),1)
    wh_uv_n0 = where(uv_tot gt 0, count_uv_n0)
    if count_uv_n0 eq 0 then stop
    min_dist_uv_n0 = min(wh_uv_n0, min_loc)
    uv_slice_ind = wh_uv_n0[min_loc]
    
    uv_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 2, $
      slice_inds = uv_slice_ind, slice_savefile = file_struct.uv_raw_savefile[i])
      
    uv_weight_slice = uvf_slice(weights_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 2, $
      slice_inds = uv_slice_ind, slice_savefile = file_struct.uv_weight_savefile[i])
      
    if max(abs(uv_slice)) eq 0 then stop
    
  endfor
  
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
      
    vf_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 0, $
      slice_inds = vf_slice_ind, slice_savefile = file_struct.vf_savefile[i])
      
    uv_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 2, $
      slice_inds = uv_slice_ind, slice_savefile = file_struct.uv_savefile[i])
      
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
    if tag_exist(file_struct, 'beam_savefile') then print, 'window integral from beam: ' + number_formatter(window_int_beam[0], format='(e10.4)')
    if tag_exist(file_struct, 'beam_int') then print, 'window integral from obs.beam_integral: ' + number_formatter(window_int_beam_obs[0], format='(e10.4)')
    
    if tag_exist(file_struct, 'beam_int') then window_int = window_int_beam_obs $
    else if tag_exist(file_struct, 'beam_savefile') then window_int = window_int_beam else window_int = window_int_k
  ;if keyword_set(sim) then window_int = 2.39e9 + fltarr(nfiles)
  endif else begin
    window_int_k = window_int * (z_mpc_delta * n_freq) * (2.*!pi)^2. / (kx_mpc_delta * ky_mpc_delta)
    print, 'var_cube multiplier: ', (z_mpc_delta * n_freq) * (2.*!pi)^2. / (kx_mpc_delta * ky_mpc_delta)
    print, 'window integral from variances: ' + number_formatter(window_int_k[0], format='(e10.4)')
    if tag_exist(file_struct, 'beam_savefile') then print, 'window integral from beam: ' + number_formatter(window_int_beam[0], format='(e10.4)')
    if tag_exist(file_struct, 'beam_int') then print, 'window integral from obs.beam_integral: ' + number_formatter(window_int_beam_obs[0], format='(e10.4)')
    
    if tag_exist(file_struct, 'beam_int') then window_int = window_int_beam_obs $
    else if tag_exist(file_struct, 'beam_savefile') then window_int = window_int_beam else window_int = window_int_k
    ;if keyword_set(sim) then window_int = 2.39e9 + fltarr(nfiles)
    
    if keyword_set(sim) then if stregex(file_struct.savefilebase, 'yy', /boolean) then begin
    
      ;volume_factor_2 = ((1./file_struct.kpix)^2. * z_mpc_mean^2.)*(z_mpc_delta * n_freq)
      ;window_int = volume_factor*16
    
      ;conv_factor = conv_factor_adrian/(2.*!pi)
      conv_factor = conv_factor / z_mpc_mean
      window_int = bandwidth_factor ;bandwidth_factor = z_mpc_delta * n_freq
      
    endif
    
  endelse
  
  
  ;; old convention
  ;; take care of FT convention for EoR (uv -> kx,ky)
  ;  data_cube1 = data_cube1 / (2.*!pi)^2.
  ;  sigma2_cube1 = sigma2_cube1 / (2.*!pi)^4.
  ;  if nfiles eq 2 then begin
  ;    data_cube2 = data_cube2 / (2.*!pi)^2.
  ;    sigma2_cube2 = sigma2_cube2 / (2.*!pi)^4.
  ;  endif
  
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
  
  ;; temp = data_cube1[*,*,10]
  ;; temp_exp = complex(dblarr(165, 83+82))
  ;; temp_exp[*, 82:*] = temp
  ;; temp_exp[*, 0:81] = conj(reverse(reverse(temp[*,1:*]),2))
  ;; temp_exp[0:81, 82] =  reverse(conj(temp[83:*,0]))
  ;; temp_exp = temp_exp[0:163, 0:163]
  ;; temp_img = fft_shift(fft(fft_shift(temp_exp)))
  ;; temp_2 = fft_shift(fft(fft_shift(temp_img), /inverse))
  ;; quick_image, abs(temp_exp), title='107',/log
  
  
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
    
    wt_ave_power_freq = fltarr(2, n_freq)
    wt_ave_power_freq[0,*] = total(total(sum_weights1 * abs(data_cube1)^2., 2), 1)/total(total(sum_weights1, 2), 1)
    wt_ave_power_freq[1,*] = total(total(sum_weights2 * abs(data_cube2)^2., 2), 1)/total(total(sum_weights2, 2), 1)
    ave_power_freq = fltarr(2, n_freq)
    for i=0, n_freq-1 do ave_power_freq[*, i] = [mean(abs((data_cube1[*,*,i])[where(sum_weights1[*,*,i] ne 0),*])^2.), $
      mean(abs((data_cube2[*,*,i])[where(sum_weights2[*,*,i] ne 0),*])^2.)]
      
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
    sum_weights1 = 1./sigma2_cube1
    wh_sig1_0 = where(sigma2_cube1 eq 0, count_sig1_0)
    if count_sig1_0 ne 0 then sum_weights1[wh_sig1_0] = 0
    
    wt_ave_power_freq = total(total(sum_weights1 * abs(data_cube1)^2., 2), 1)/total(total(sum_weights1, 2), 1)
    ave_power_freq = fltarr(n_freq)
    for i=0, n_freq-1 do ave_power_freq[i] = mean(abs((data_cube1[*,*,i])[where(sum_weights1[*,*,i] ne 0),*])^2.)
    undefine, sum_weights1
    
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
    slice_inds = uf_slice_ind, slice_savefile = file_struct.uf_sum_savefile)
  if nfiles eq 2 then $
    uf_slice = uvf_slice(data_diff, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 1, $
    slice_inds = uf_slice_ind, slice_savefile = file_struct.uf_diff_savefile)
    
  vf_slice = uvf_slice(data_sum, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 0, $
    slice_inds = vf_slice_ind, slice_savefile = file_struct.vf_sum_savefile)
  if nfiles eq 2 then $
    vf_slice2 = uvf_slice(data_diff, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 0, $
    slice_inds = vf_slice_ind, slice_savefile = file_struct.vf_diff_savefile) $
  else vf_slice2 = vf_slice
  
  uv_slice = uvf_slice(data_sum, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 2, $
    slice_inds = uv_slice_ind, slice_savefile = file_struct.uv_sum_savefile)
  if nfiles eq 2 then $
    uv_slice2 = uvf_slice(data_diff, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, hubble_param, slice_axis = 2, $
    slice_inds = uv_slice_ind, slice_savefile = file_struct.uv_diff_savefile) $
  else uv_slice2 = uv_slice
  
  
  ;; apply spectral windowing function if desired
  if n_elements(spec_window_type) ne 0 then begin
    window = spectral_window(n_freq, type = spec_window_type, /periodic)
    
    norm_factor = sqrt(n_freq/total(window^2.))
    
    window = window * norm_factor
    
    window_expand = rebin(reform(window, 1, 1, n_freq), n_kx, n_ky, n_freq, /sample)
    
    data_sum = data_sum * window_expand
    if nfiles eq 2 then data_diff = data_diff * window_expand
    
    sum_sigma2 = sum_sigma2 * temporary(window_expand^2.)
  endif
  
  ;; now take frequency FT
  if even_freq then begin
    ;; evenly spaced, just use fft
    ;; old ft convention
    ; data_sum_ft = fft(data_sum, dimension=3) * n_freq * z_mpc_delta / (2.*!pi)
    data_sum_ft = fft(data_sum, dimension=3) * n_freq * z_mpc_delta
    
    ;; put k0 in middle of cube
    data_sum_ft = shift(data_sum_ft, [0,0,n_kz/2])
    
    undefine, data_sum
    if nfiles eq 2 then begin
      ;; old ft convention
      ; data_diff_ft = fft(data_diff, dimension=3) * n_freq * z_mpc_delta / (2.*!pi)
      data_diff_ft = fft(data_diff, dimension=3) * n_freq * z_mpc_delta
      ;; put k0 in middle of cube
      data_diff_ft = shift(data_diff_ft, [0,0,n_kz/2])
      undefine, data_diff
    endif
  endif else begin
    ;; Not evenly spaced. Do a dft
    z_exp =  exp(-1.*complex(0,1)*matrix_multiply(comov_dist_los, kz_mpc_orig, /btranspose))
    
    ;; old ft convention
    ;    data_sum_ft = z_mpc_delta/(2.*!pi) * $
    ;      reform(matrix_multiply(reform(temporary(data_sum), n_kx*n_ky, n_freq), z_exp), n_kx, n_ky, n_kz)
    ;    if nfiles eq 2 then $
    ;      data_diff_ft = z_mpc_delta/(2.*!pi) * $
    ;      reform(matrix_multiply(reform(temporary(data_diff), n_kx*n_ky, n_freq), z_exp), n_kx, n_ky, n_kz)
    data_sum_ft = z_mpc_delta * $
      reform(matrix_multiply(reform(temporary(data_sum), n_kx*n_ky, n_freq), z_exp), n_kx, n_ky, n_kz)
    if nfiles eq 2 then $
      data_diff_ft = z_mpc_delta * $
      reform(matrix_multiply(reform(temporary(data_diff), n_kx*n_ky, n_freq), z_exp), n_kx, n_ky, n_kz)
  endelse
  
  n_val = round(kz_mpc_orig / kz_mpc_delta)
  kz_mpc_orig[where(n_val eq 0)] = 0
  
  kz_mpc = kz_mpc_orig[where(n_val ge 0)]
  n_kz = n_elements(kz_mpc)
  
  if keyword_set(std_power) then begin
    ;; comov_dist_los goes from large to small z
    z_relative = dindgen(n_freq)*z_mpc_delta
    freq_kz_arr = rebin(reform(kz_mpc_orig, 1, n_elements(kz_mpc_orig)), n_freq, n_elements(kz_mpc_orig)) * rebin(z_relative, n_freq, n_elements(kz_mpc_orig))
    
    ;kz_arr = rebin(reform(kz_mpc_orig, 1, n_elements(kz_mpc_orig)), n_elements(kz_mpc_orig), n_elements(kz_mpc_orig)) - $
    ;  rebin(kz_mpc_orig, n_elements(kz_mpc_orig), n_elements(kz_mpc_orig))
    ;full_freq_kz_arr = rebin(reform(kz_arr, 1, n_elements(kz_mpc_orig)^2), n_freq, n_elements(kz_mpc_orig)^2) * rebin(z_relative, n_freq, n_elements(kz_mpc_orig)^2)
    
    ;; for standard power calc. just need ft of sigma2 (sigma has squared units relative to data, so use z_mpc_delta^2d)
    ;; old ft convention
    ;     sigma2_ft = matrix_multiply(reform(sum_sigma2, n_kx*n_ky, n_freq), exp(-1.*complex(0,1)*freq_kz_arr)*exp(1.*complex(0,1)*freq_kz_arr)) * (z_mpc_delta^2. / (2.*!pi))^2d
    sigma2_ft = matrix_multiply(reform(sum_sigma2, n_kx*n_ky, n_freq), exp(-1.*complex(0,1)*freq_kz_arr)*exp(1.*complex(0,1)*freq_kz_arr)) * (z_mpc_delta^2.)^2d
    sigma2_ft = reform(sigma2_ft, n_kx, n_ky, n_elements(kz_mpc_orig))
    
    ;; old ft convention
    ;    sigma2_ft_cov = matrix_multiply(reform(sum_sigma2, n_kx*n_ky, n_freq), exp(-1.*complex(0,1)*freq_kz_arr)*exp(-1.*complex(0,1)*freq_kz_arr)) * (z_mpc_delta^2. / (2.*!pi))^2d
    sigma2_ft_cov = matrix_multiply(reform(sum_sigma2, n_kx*n_ky, n_freq), exp(-1.*complex(0,1)*freq_kz_arr)*exp(-1.*complex(0,1)*freq_kz_arr)) * (z_mpc_delta^2.)^2d
    sigma2_ft_cov = reform(sigma2_ft_cov, n_kx, n_ky, n_elements(kz_mpc_orig))
    sigma2_ft_cov[*,*, where(n_val eq 0)] = 0
    
    ;sigma2_ft_full = matrix_multiply(reform(sum_sigma2, n_kx*n_ky, n_freq), exp(1.*complex(0,1)*full_freq_kz_arr)) * (z_mpc_delta^2. / (2.*!pi))^2d
    ;sigma2_ft_full = reform(sigma2_ft_full, n_kx, n_ky, n_elements(kz_mpc_orig), n_elements(kz_mpc_orig))
    undefine, sum_sigma2
    
    ;; diagonalize to get rid of covariance
    theta = atan(2*sigma2_ft_cov[*,*,where(n_val ge 0)], sigma2_ft[*,*,where(n_val ge 0)]-sigma2_ft[*,*,where(n_val le 0)])/2.
    
    theta[where(sigma2_ft_cov[*,*,where(n_val ge 0)] eq 0)] = 0
    
    cos_theta = cos(theta)
    sin_theta = sin(theta)
    undefine, theta
    
    sigma2_1 = sigma2_ft[*,*,where(n_val ge 0)]*cos_theta^2. + 2.*sigma2_ft_cov[*,*,where(n_val ge 0)]*cos_theta*sin_theta + sigma2_ft[*,*,where(n_val le 0)]*sin_theta^2.
    sigma2_2 = sigma2_ft[*,*,where(n_val ge 0)]*sin_theta^2. - 2.*sigma2_ft_cov[*,*,where(n_val ge 0)]*cos_theta*sin_theta + sigma2_ft[*,*,where(n_val le 0)]*cos_theta^2.
    sigma2_2[*,*,0] = 0.
    undefine, sigma2_ft, sigma2_ft_cov
    
    data_sum_1 = data_sum_ft[*,*,where(n_val ge 0)]*cos_theta + data_sum_ft[*,*,where(n_val le 0)]*sin_theta
    data_sum_2 = (-1d)*data_sum_ft[*,*,where(n_val ge 0)]*sin_theta + data_sum_ft[*,*,where(n_val le 0)]*cos_theta
    data_sum_2[*,*,0] = 0.
    undefine, data_sum_ft
    if nfiles eq 2 then begin
      data_diff_1 = data_diff_ft[*,*,where(n_val ge 0)]*cos_theta + data_diff_ft[*,*,where(n_val le 0)]*sin_theta
      data_diff_2 = (-1d)*data_diff_ft[*,*,where(n_val ge 0)]*sin_theta + data_diff_ft[*,*,where(n_val le 0)]*cos_theta
      data_diff_2[*,*,0] = 0.
      undefine, data_diff_ft
    endif
    
    git, repo_path = ps_repository_dir(), result=kcube_git_hash
    git_hashes = {uvf:uvf_git_hashes, uvf_wt:uvf_wt_git_hashes, beam:beam_git_hashes, kcube:kcube_git_hash}
    
    save, file = file_struct.kcube_savefile, data_sum_1, data_sum_2, data_diff_1, data_diff_2, sigma2_1, sigma2_2, n_val, $
      kx_mpc, ky_mpc, kz_mpc, kperp_lambda_conv, delay_params, hubble_param, n_freq_contrib, freq_mask, $
      vs_name, vs_mean, window_int, wt_meas_ave, wt_meas_min, ave_weights, wt_ave_power_freq, ave_power_freq, git_hashes
      
  endif else begin
    ;; these an and bn calculations don't match the standard
    ;; convention (they down by a factor of 2) but they make more sense
    ;; and remove factors of 2 we'd otherwise have in the power
    ;; and variance calculations
    ;; note that the 0th mode will have higher noise because there's half as many measurements going into it
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
    data_sum_cos[*, *, 0] = a1_0 ;/2. changed 3/12/14.
    data_sum_cos[*, *, 1:n_kz-1] = a1_n
    data_sum_sin[*, *, 1:n_kz-1] = b1_n
    if nfiles gt 1 then begin
      data_diff_cos = complex(fltarr(n_kx, n_ky, n_kz))
      data_diff_sin = complex(fltarr(n_kx, n_ky, n_kz))
      data_diff_cos[*, *, 0] = a2_0 ;/2. changed 3/12/14.
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
    ;; old ft convention
    ;    covar_cos = matrix_multiply(sum_sigma2, cos_arr^2d) * (z_mpc_delta / (2.*!pi))^2.
    ;    covar_sin = matrix_multiply(sum_sigma2, sin_arr^2d) * (z_mpc_delta / (2.*!pi))^2.
    ;    covar_cross = matrix_multiply(sum_sigma2, cos_arr*sin_arr) * (z_mpc_delta / (2.*!pi))^2.
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
    
    ;; cos 0 term does NOT have a different normalization
    ;; changed 3/12/14 -- the noise is higher in 0th bin because there are only half as many measurements
    ;;covar_cos[*,*,0] = covar_cos[*,*,0]/4d
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
    
    ;    sigma2_1 = covar_cos
    ;    sigma2_2 = covar_sin
    
    undefine, covar_cos, covar_sin, covar_cross
    
    data_sum_1 = data_sum_cos*cos_theta + data_sum_sin*sin_theta
    data_sum_2 = (-1d)*data_sum_cos*sin_theta + data_sum_sin*cos_theta
    undefine, data_sum_cos, data_sum_sin
    if nfiles eq 2 then begin
      data_diff_1 = data_diff_cos*cos_theta + data_diff_sin*sin_theta
      data_diff_2 = (-1d)*data_diff_cos*sin_theta + data_diff_sin*cos_theta
      undefine, data_diff_cos, data_diff_sin
    endif
    
    ;    data_sum_1 = data_sum_cos
    ;    data_sum_2 = data_sum_sin
    ;    undefine, data_sum_cos, data_sum_sin
    ;    if nfiles eq 2 then begin
    ;      data_diff_1 = data_diff_cos
    ;      data_diff_2 = data_diff_sin
    ;      undefine, data_diff_cos, data_diff_sin
    ;    endif
    
    git, repo_path = ps_repository_dir(), result=kcube_git_hash
    git_hashes = {uvf:uvf_git_hashes, uvf_wt:uvf_wt_git_hashes, beam:beam_git_hashes, kcube:kcube_git_hash}
    
    save, file = file_struct.kcube_savefile, data_sum_1, data_sum_2, data_diff_1, data_diff_2, sigma2_1, sigma2_2, $
      kx_mpc, ky_mpc, kz_mpc, kperp_lambda_conv, delay_params, hubble_param, n_freq_contrib, freq_mask, $
      vs_name, vs_mean, window_int, wt_meas_ave, wt_meas_min, ave_weights, wt_ave_power_freq, ave_power_freq, git_hashes
  endelse
  
end
