

function kspace_rebinning_1d, power, k1_mpc, k2_mpc, k3_mpc, k_edges_mpc, k_bin = k_bin, log_k = log_k, $
    noise_expval = noise_expval, binned_noise_expval = noise_expval_1d, weights = weights, $
    binned_weights = weights_1d, kperp_range = kperp_range, kx_range = kx_range, ky_range = ky_range, $
    kpar_range = kpar_range, wedge_amp = wedge_amp, $
    nbins_1d = bin_hist, var_power_1d = var_power_1d, mean_var_1d = mean_var_1d, $
    coarse_harm0 = coarse_harm0, coarse_width = coarse_width, bin_arr_3d = bin_arr_3d, noise_frac_3d = noise_frac_3d, $
    edge_on_grid = edge_on_grid, match_datta = match_datta, kpar_power = kpar_power, kperp_power = kperp_power, $
    kperp_density_measure = kperp_density_measure, kperp_density_cutoff = kperp_density_cutoff
    
  power_size = size(power, /dimensions)
  power_dim = n_elements(power_size)
  
  if keyword_set(kpar_power) and keyword_set(kperp_power) then message, 'Only one of kpar_power and kperp_power can be set'
  
  if power_dim lt 2 or power_dim gt 3 then $
    message, 'power array must be 2 or 3 dimensional and ordered (kperpendicular, kparallel) or (kx,ky,kz)'
    
    
  if n_elements(k_bin) gt 1 then message, 'k_bin must be a scalar because it specifies the linear bin size ' + $
    'or the log bin size (bins per decade = 1/log bin size) if log_k keyword is set. '
    
  ;; for log binning, k_bins define log binsize (1/bins per decade). Default to 0.1 (10 bins per decade)
  if keyword_set(log_k) then begin
    if n_elements(k_bin) eq 0 then k_bin = 0.1d else if k_bin gt 1 then $
      message, 'k_bin must be less than one if log_k keyword is set (so there will be > 1 bin per decade)'
  endif else begin
    ;; for linear binning default to min of (kx, ky, kz) binsize
    if n_elements(k_bin) eq 0 then begin
      if power_dim eq 3 then k_bin = max([k1_mpc[1]-k1_mpc[0], k2_mpc[1]-k2_mpc[0], k3_mpc[1]-k3_mpc[0]]) $
      else k_bin = max([k1_mpc[1]-k1_mpc[0], k2_mpc[1]-k2_mpc[0]])
    endif
  endelse
  
  if n_elements(kperp_density_measure) gt 0 and n_elements(kperp_density_cutoff) gt 0 then begin
    if total(abs(size(kperp_density_measure, /dimensions) - power_size[0:power_dim-2])) ne 0 then $
      message, 'If kperp_density_measure is provided, it must have the same kperp dimensions as the power'
    kperp_density_measure_use = kperp_density_measure
  endif
  
  if n_elements(weights) ne 0 then begin
    if total(abs(size(weights, /dimensions) - power_size)) ne 0 then $
      message, 'If weights array is provided, it must have the same dimensionality as the power'
    weights_use = weights
    power_use = power
  endif else begin
    weights_use = dblarr(power_size) + 1d
    power_use = power
  endelse
  
  if n_elements(noise_expval) ne 0 then begin
    if total(abs(size(noise_expval, /dimensions) - power_size)) ne 0 then $
      message, 'If noise_expval array is provided, it must have the same dimensionality as the power'
    noise_expval_use = noise_expval
  endif else begin
    noise_expval_use = 1/sqrt(weights_use)
    wh_wt0 = where(weights_use eq 0, count_wt0)
    if count_wt0 gt 0 then noise_expval_use[wh_wt0] = 0
  endelse
  
  if power_dim eq 3 then begin
  
    kx_mpc = k1_mpc
    ky_mpc = k2_mpc
    kz_mpc = k3_mpc
    kpar_ind_map = indgen(n_elements(k3_mpc))
    
    n_kx = power_size[0]
    n_ky = power_size[1]
    n_kz = power_size[2]
    
    if n_elements(kx_mpc) ne n_kx then message, 'Length of k1_mpc must be the same as the first dimension of the power array'
    if n_elements(ky_mpc) ne n_ky then message, 'Length of k2_mpc must be the same as the second dimension of the power array'
    if n_elements(kz_mpc) ne n_kz then message, 'Length of k3_mpc must be the same as the third dimension of the power array'
    
    if n_elements(bin_arr_3d) eq 0 then begin
      bin_arr_3d = lonarr(n_kx, n_ky, n_kz) + 1
      
      if n_elements(kpar_range) gt 0 then begin
        if n_elements(kpar_range) ne 2 then message, 'kpar_range must be a 2 element vector'
        
        wh_good = where(abs(kz_mpc) ge kpar_range[0] and abs(kz_mpc) le kpar_range[1], count_good, ncomplement = count_bad, complement = wh_bad)
        if count_good eq 0 then message, 'No kz values within kpar_range'
        
        if count_bad gt 0 then bin_arr_3d[*,*,wh_bad] = 0
      endif else kpar_range = minmax(kz_mpc)
      
      if n_elements(coarse_harm0) gt 0 then begin
        n_harm = floor(n_kz/float(coarse_harm0))+1
        n_width = coarse_width*2-1
        kpar_ch_bad = rebin(coarse_harm0 * (findgen(n_harm)+1), n_harm, n_width) + $
          rebin(reform(findgen(n_width)-(coarse_width-1), 1, n_width), n_harm, n_width)
          
        match,kpar_ind_map,kpar_ch_bad[*],suba,subb
        bin_arr_3d[*,*,suba] = 0
      endif
      
      if n_elements(kx_range) gt 0 then begin
        if n_elements(kx_range) ne 2 then message, 'kx_range must be a 2 element vector'
        
        wh_good = where(abs(kx_mpc) ge kx_range[0] and abs(kx_mpc) le kx_range[1], count_good, ncomplement = count_bad, complement = wh_bad)
        if count_good eq 0 then message, 'No kx values within kx_range'
        
        if count_bad gt 0 then bin_arr_3d[wh_bad,*,*] = 0
      endif else kx_range = minmax(kx_mpc)
      
      if n_elements(ky_range) gt 0 then begin
        if n_elements(ky_range) ne 2 then message, 'ky_range must be a 2 element vector'
        
        wh_good = where(abs(ky_mpc) ge ky_range[0] and abs(ky_mpc) le ky_range[1], count_good, ncomplement = count_bad, complement = wh_bad)
        if count_good eq 0 then message, 'No ky values within ky_range'
        
        if count_bad gt 0 then bin_arr_3d[*,wh_bad,*] = 0
      endif else ky_range = minmax(ky_mpc)
      
      kperp_vals = sqrt(rebin(kx_mpc^2d, n_kx, n_ky) + rebin(reform(ky_mpc^2d, 1, n_ky), n_kx, n_ky))
      if n_elements(kperp_range) gt 0 then begin
        if n_elements(kperp_range) ne 2 then message, 'kperp_range must be a 2 element vector'
        
        wh_good = where(kperp_vals ge kperp_range[0] and kperp_vals le kperp_range[1], count_good, ncomplement = count_bad, complement = wh_bad)
        if count_good eq 0 then message, 'No kperp values within kperp_range'
        
        if count_bad gt 0 then begin
          bin_arr_3d = reform(bin_arr_3d, n_kx*n_ky, n_kz)
          bin_arr_3d[wh_bad, *] = 0
          bin_arr_3d = reform(bin_arr_3d, n_kx, n_ky, n_kz)
        endif
      endif else kperp_range = minmax(kperp_vals)
      
      if n_elements(wedge_amp) gt 0 then begin
      
        temp_kpar = rebin(reform(kz_mpc, 1, 1, n_kz), n_kx, n_ky, n_kz)
        temp_kperp = rebin(kperp_vals, n_kx, n_ky, n_kz)
        wh_above_wedge = where(temp_kpar gt temp_kperp*max(wedge_amp), count_above_wedge, ncomplement = count_below_wedge, complement = wh_below_wedge)
        undefine, temp_kperp, temp_kpar
        
        if count_above_wedge eq 0 then begin
          message, 'no pixels above wedge, using full volume'
        endif else begin
          bin_arr_3d[wh_below_wedge] = 0
        endelse
      endif
      
    endif
    
    
    ;; apply density correction
    if n_elements(kperp_density_measure_use) gt 0 and n_elements(kperp_density_cutoff) gt 0 then begin
      wh_kperp_dense = where(kperp_density_measure_use gt kperp_density_cutoff, count_kperp_dense, ncomplement = count_sparse, complement = wh_sparse)
      if count_kperp_dense eq 0 then begin
        print, 'no kperp values exceed kperp_density_cutoff, using full volume'
      endif else begin
        bin_arr_3d = reform(bin_arr_3d, n_kx*n_ky, n_kz)
        bin_arr_3d[wh_sparse, *] = 0
        bin_arr_3d = reform(bin_arr_3d, n_kx, n_ky, n_kz)
        
        power_use = power_use/0.5
        if n_elements(noise_expval_use) gt 0 then noise_expval_use = noise_expval_use/0.5
        if n_elements(weights_use) gt 0 then weights_use = weights_use*0.5^2.
      endelse
    endif
    if max(bin_arr_3d) eq 0 then message, 'All voxels have been excluded (max of bin_arr_3d is zero)'
    
    ;; if bin_arr_3d is just a mask (only 0s & 1s) then need to get 1d bin values
    if max(bin_arr_3d) eq 1 then begin
      mask_arr_3d = bin_arr_3d
      bin_arr_3d = lonarr(n_kx, n_ky, n_kz)
      
      if keyword_set(kpar_power) then begin
        ;; generate kpar values for histogram
        temp = rebin(reform(kz_mpc, 1, 1, n_kz), n_kx, n_ky, n_kz)
      endif else if keyword_set(kperp_power) then begin
        ;; generate kperp values for histogram
        temp = sqrt(rebin(kx_mpc^2d, n_kx, n_ky, n_kz) + rebin(reform(ky_mpc^2d, 1, n_ky), n_kx, n_ky, n_kz))
      endif else begin
        ;; generate k values for histogram
        temp = sqrt(rebin(kx_mpc^2d, n_kx, n_ky, n_kz) + rebin(reform(ky_mpc^2d, 1, n_ky), n_kx, n_ky, n_kz) + $
          rebin(reform(kz_mpc^2d, 1, 1, n_kz), n_kx, n_ky, n_kz))
      endelse
      
      if not keyword_set(log_k) then begin
        ;; linear binning
        if keyword_set(edge_on_grid) then k_min = floor(min(temp) / k_bin) * k_bin $
        else k_min = floor((min(temp) + k_bin/2d) / k_bin) * k_bin - k_bin/2d
        
        ;; Use histogram with reverse indicies to bin in k
        k_hist = histogram(temp, binsize = k_bin, min = k_min, omax = k_max, locations = lin_k_locs, $
          reverse_indices = k_ri)
          
      endif else begin
        ;; log binning
        ;; get min value, adjust so centers of bins will be on major grid.
        min_val = sqrt(min(kx_mpc^2d) + min(ky_mpc^2d) + min(kz_mpc^2d))
        if min_val gt 0 then begin
          min_non_zero = min_val
          if keyword_set(edge_on_grid) then $
            k_min = floor(alog10(min_non_zero) / k_bin) * k_bin $
          else $
            k_min = floor((alog10(min_non_zero) + k_bin/2d) / k_bin) * k_bin - k_bin/2d
        endif else begin
          min_candidates = sqrt([min(kx_mpc[where(kx_mpc^2d ne 0)]^2d) + min(ky_mpc^2d) + min(kz_mpc^2d), $
            min(kx_mpc^2d) + min(ky_mpc[where(ky_mpc^2d ne 0)]^2d) + min(kz_mpc^2d), $
            min(kx_mpc^2d) + min(ky_mpc^2d) + min(kz_mpc[where(kz_mpc^2d ne 0)]^2d)])
          min_non_zero = min(min_candidates)
          
          ;; Add an extra bin to put the k=0 mode into.
          if keyword_set(edge_on_grid) then $
            k_min = floor(alog10(min_non_zero) / k_bin - 1) * k_bin $
          else $
            k_min = floor((alog10(min_non_zero) + k_bin/2d) / k_bin - 1) * $
            k_bin - k_bin/2d
          ;; if keyword_set(match_datta) then begin
          ;;    if keyword_set(edge_on_grid) then k_min = -3d $
          ;;    else k_min = -3d - k_log_binsize/2d
          ;; endif
            
          ;;Set 0 mode to bottom value so it doesn't get left out.
          wh_k0 = where(temp eq 0, n_k0)
          if n_k0 gt 1 then stop
          if n_k0 gt 0 then temp[wh_k0] = 10^k_min
          
        endelse
        
        ;; calculate log k at each location in 1/Mpc.
        k_array = alog10(temporary(temp))
        
        print, 'Histogramming array of log k values'
        ;; Use histogram with reverse indicies to bin in log k (want even spacing in log space.)
        k_hist = histogram((k_array), binsize = k_bin, min = k_min, omax = k_max, locations = log_k_locs, $
          reverse_indices = k_ri)
        undefine, k_array
      endelse
      
      mask_hist = lonarr(n_elements(k_hist))
      for i=0, n_elements(k_hist)-1 do begin
        if k_hist[i] ne 0 then begin
          inds =  k_ri[k_ri[i] : k_ri[i+1]-1]
          bin_arr_3d[inds] = i+1
          mask_hist[i] = total(mask_arr_3d[inds])
        endif
      endfor
      bin_arr_3d = bin_arr_3d * mask_arr_3d
      undefine, mask_arr_3d
      undefine, k_ri, k_hist
      
      if keyword_set(log_k) then begin
        k_centers_mpc = 10^(log_k_locs + k_bin/2d)
        k_edges_mpc = [10^(log_k_locs), 10^(max(log_k_locs) + k_bin)]
      endif else begin
        k_centers_mpc = lin_k_locs + k_bin/2d
        k_edges_mpc = [lin_k_locs, max(lin_k_locs) + k_bin]
      endelse
      
      ;; strip off any bins with no unmasked pixels
      while mask_hist[-1] eq 0 do begin
        k_centers_mpc = k_centers_mpc[0:-2]
        k_edges_mpc = k_edges_mpc[0:-2]
        mask_hist = mask_hist[0:-2]
      endwhile

      
    endif
    
    ;; now bin_arr_3d gives us the 1d bins to put each voxel into. If this is passed in, k_edges_mpc need to also be passed in.
    if n_elements(k_edges_mpc) eq 0 then message, 'if a bin array is passed in with values greater than one (indicating which bin to put each voxel in), ' + $
      'k_edges_mpc also needs to be passed in because it is not calculated.'
    ;; we need reverse indices, so use histogram on bin_arr_3d
    bin_hist = histogram(bin_arr_3d, binsize = 1, min = 1, omax = bin_max, reverse_indices = bin_ri)
    
    full_ind_arr = lindgen(n_kx, n_ky, n_kz)
    noise_frac_3d = fltarr(n_kx, n_ky, n_kz)
    
    n_bins = n_elements(bin_hist)
    power_1d = dblarr(n_bins)
    if n_elements(noise_expval_use) gt 0 then nev_1d = dblarr(n_bins)
    weights_1d = dblarr(n_bins)
    var_power_1d = dblarr(n_bins)
    mean_var_1d = dblarr(n_bins)
    
    for i=0, n_bins-1 do begin
      if bin_hist[i] ne 0 then begin
        inds =  bin_ri[bin_ri[i] : bin_ri[i+1]-1]
        power_1d[i] = total(power_use[inds] * weights_use[inds])
        if n_elements(noise_expval_use) gt 0 then nev_1d[i] = total(noise_expval_use[inds] * weights_use[inds])
        if n_elements(weights_use) gt 0 then weights_1d[i] = total(weights_use[inds]) $
        else weights_1d[i] = bin_hist[i]
        
        temp = weights_use[inds]
        temp_sigma = 1./temp
        zero_i=where(temp eq 0,n_zero)
        IF n_zero GT 0 THEN temp_sigma[zero_i] = 0
        
        if min(weights_use[inds]) eq 0 then begin
          wt_n0 = where(weights_use[inds] gt 0, count_n0)
          if count_n0 gt 0 then begin
            var_power_1d[i] = total(power_use[inds[wt_n0]]^2.)/count_n0
            mean_var_1d[i] = mean(temp_sigma[wt_n0])
          endif
        endif else begin
          var_power_1d[i] = total(power_use[inds]^2.)/bin_hist[i]
          mean_var_1d[i] = mean(temp_sigma)
        endelse
        
        full_inds = full_ind_arr[inds]
        
        where_noise = where(power_use[inds]/weights_use[inds] le 1/sqrt(weights_use[inds]) and weights_use[inds] gt 0, count_noise)
        noise_frac_3d[full_inds] = count_noise / float(bin_hist[i])
      endif
    endfor
    undefine, bin_ri
    
    
  endif else begin
  
    kperp_mpc = k1_mpc
    kpar_mpc = k2_mpc
    k3_mpc = 0
    
    n_kperp = power_size[0]
    n_kpar = power_size[1]
    
    if n_elements(kperp_mpc) ne n_kperp then message, 'Length of k1_mpc must be the same as the first dimension of the power array'
    if n_elements(kpar_mpc) ne n_kpar then message, 'Length of k2_mpc must be the same as the second dimension of the power array'
    
    if n_elements(kperp_range) gt 0 then begin
      if  n_elements(kperp_range) ne 2 then message, 'kperp_range must be a 2 element vector'
      
      wh_good = where(kperp_mpc ge kperp_range[0] and kperp_mpc le kperp_range[1], count_good, ncomplement = count_bad)
      if count_good eq 0 then message, 'No kperp values within kperp_range'
      
      if count_bad gt 0 then begin
        power_use = power_use[wh_good,*]
        if n_elements(weights_use) gt 0 then weights_use = weights_use[wh_good,*]
        if n_elements(noise_expval_use) gt 0 then noise_expval_use = noise_expval_use[wh_good,*]
        kperp_mpc = kperp_mpc[wh_good]
        n_kperp = n_elements(kperp_mpc)
      endif
    endif else kperp_range = minmax(kperp_mpc)
    
    if n_elements(kpar_range) gt 0 then begin
      if  n_elements(kpar_range) ne 2 then message, 'kpar_range must be a 2 element vector'
      
      wh_good = where(kpar_mpc ge kpar_range[0] and kpar_mpc le kpar_range[1], count_good, ncomplement = count_bad)
      if count_good eq 0 then message, 'No kpar values within kpar_range'
      
      if count_bad gt 0 then begin
        power_use = power_use[*,wh_good]
        if n_elements(weights_use) gt 0 then weights_use = weights_use[*,wh_good]
        if n_elements(noise_expval_use) gt 0 then noise_expval_use = noise_expval_use[*,wh_good]
        kpar_mpc = kpar_mpc[wh_good]
        n_kz = n_elements(kpar_mpc)
      endif
    endif else kpar_range = minmax(kpar_mpc)
    
    if keyword_set(kpar_power) then begin
      ;; generate kpar values for histogram
      temp = rebin(reform(kpar_mpc, 1, n_kpar), n_kperp, n_kpar)
    endif else if keyword_set(kperp_power) then begin
      ;; generate kperp values for histogram
      temp =  rebin(kperp_mpc, n_kperp, n_kpar)
    endif else begin
      temp = sqrt(rebin(kperp_mpc, n_kperp, n_kpar)^2 + rebin(reform(kpar_mpc, 1, n_kpar), n_kperp, n_kpar)^2)
    endelse
    temp_kperp = rebin(kperp_mpc, n_kperp, n_kpar)
    temp_kpar = rebin(reform(kpar_mpc, 1, n_kpar), n_kperp, n_kpar)
    
    if n_elements(kperp_density_measure_use) gt 0 and n_elements(kperp_density_cutoff) gt 0 then begin
      wh_kperp_good = where(kperp_density_measure_use gt kperp_density_cutoff, count_kperp_good)
      
      if count_kperp_good gt 0 then begin
        temp = temp[wh_kperp_good, *]
        temp_kperp = temp_kperp[wh_kperp_good, *]
        temp_kpar = temp_kpar[wh_kperp_good, *]
        power_use = power_use[wh_good,*]/0.5
        if n_elements(weights_use) gt 0 then weights_use = weights_use[wh_good,*]*0.5^2.
        if n_elements(noise_expval_use) gt 0 then noise_expval_use = noise_expval_use[wh_good,*]/0.5
      endif else print, 'no kperp values exceed kperp_density_cutoff, using full volume'
      
    endif
    
    if n_elements(wedge_amp) gt 0 then begin
      wh_above_wedge = where(temp_kpar gt temp_kperp*max(wedge_amp), count_above_wedge)
      if count_above_wedge gt 0 then begin
        temp = temp[wh_above_wedge]
        power_use = power_use[wh_above_wedge]
        if n_elements(weights_use) gt 0 then weights_use = weights_use[wh_above_wedge]
        if n_elements(noise_expval_use) gt 0 then noise_expval_use = noise_expval_use[wh_above_wedge]
      endif else print, 'no pixels above wedge, using full volume'
    endif
    
    if not keyword_set(log_k) then begin
    
      if keyword_set(edge_on_grid) then k_min = floor(min(temp) / k_bin) * k_bin $
      else kperp_min = floor((min(temp) + k_bin/2d) / k_bin) * k_bin - k_bin/2d
      
      ;; Use histogram with reverse indicies to bin in kperp
      k_hist = histogram(temp, binsize = k_bin, min = k_min, omax = k_max, locations = lin_k_locs, $
        reverse_indices = k_ri)
        
    endif else begin
      ;; get min value, adjust so centers of bins will be on major grid.
      if min(temp) gt 0 then begin
        min_non_zero = min(temp)
        if keyword_set(edge_on_grid) then $
          k_min = floor(alog10(min_non_zero) / k_bin) * k_bin $
        else $
          k_min = floor((alog10(min_non_zero) + k_bin/2d) / k_bin) * $
          k_bin - k_bin/2d
      endif else begin
        wh0 = where(temp eq 0, n_k0, complement = wh_non0)
        min_non_zero = min(temp[wh_non0])
        ;; Add an extra bin to put the k=0 mode into.
        if keyword_set(edge_on_grid) then $
          k_min = floor(alog10(min_non_zero) / k_bin - 1) * k_bin $
        else $
          k_min = floor((alog10(min_non_zero) + k_bin/2d) / k_bin - 1) * $
          k_bin - k_bin/2d
          
        ;;Set 0 mode to bottom value so it doesn't get left out.
        if n_k0 gt 0 then temp[wh0] = 10^k_min
      endelse
      
      ;; calculate log k at each location in 1/Mpc.
      k_array = alog10(temporary(temp))
      
      ;; Use histogram with reverse indicies to bin in log k (want even spacing in log space.)
      k_hist = histogram(temporary(k_array), binsize = k_bin, min = k_min, omax = k_max, locations = log_k_locs, $
        reverse_indices = k_ri)
    endelse
    
    if n_elements(wedge_amp) gt 0 then if count_above_wedge gt 0 then begin
      power_use = power_use[wh_above_wedge]
      if n_elements(noise_expval_use) gt 0 then noise_expval_use = noise_expval_use[wh_above_wedge]
      if n_elements(weights_use) gt 0 then weights_use = weights_use[wh_above_wedge]
    endif
    bin_hist = k_hist
    
    n_k = n_elements(k_hist)
    power_1d = dblarr(n_k)
    if n_elements(noise_expval_use) gt 0 then nev_1d = dblarr(n_k)
    weights_1d = dblarr(n_k)
    var_power_1d = dblarr(n_k)
    mean_var_1d = dblarr(n_k)
    
    for i=0, n_k-1 do begin
      if k_hist[i] ne 0 then begin
        inds = k_ri[k_ri[i] : k_ri[i+1]-1]
        power_1d[i] = total(power_use[inds] * weights_use[inds])
        if n_elements(noise_expval_use) gt 0 then nev_1d[i] = total(noise_expval_use[inds] * weights_use[inds])
        if n_elements(weights_use) gt 0 then weights_1d[i] = total(weights_use[inds]) $
        else weights_1d[i] = k_hist[i]
        
        temp = weights_use[inds]
        temp_sigma = 1./temp
        temp_sigma[where(temp eq 0)] = 0
        
        if min(weights_use[inds]) eq 0 then begin
          wt_n0 = where(weights_use[inds] gt 0, count_n0)
          if count_n0 gt 0 then begin
            var_power_1d[i] = total(power_use[inds[wt_n0]]^2.)/count_n0
            mean_var_1d[i] = mean(temp_sigma[wh_n0])
          endif
        endif else begin
          var_power_1d[i] = total(power_use[inds]^2.)/k_hist[i]
          mean_var_1d[i] = mean(temp_sigma)
        endelse
        
      endif
    endfor
    k_ri=0
    
    if keyword_set(log_k) then begin
      k_centers_mpc = 10^(log_k_locs + k_bin/2d)
      k_edges_mpc = [10^(log_k_locs), 10^(max(log_k_locs) + k_bin)]
    endif else begin
      k_centers_mpc = lin_k_locs + k_bin/2d
      k_edges_mpc = [lin_k_locs, max(lin_k_locs) + k_bin]
    endelse
  endelse
  
  power_ave = power_1d/weights_1d
  if n_elements(nev_1d) gt 0 then noise_expval_1d = nev_1d/weights_1d
  wh = where(weights_1d eq 0, count)
  if count ne 0 then begin
    power_ave[wh] = 0d
    if n_elements(nev_1d) gt 0 then noise_expval_1d[wh] = 0d
  endif
  
  return, power_ave
  
end
