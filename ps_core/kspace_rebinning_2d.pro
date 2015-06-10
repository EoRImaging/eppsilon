;; Function to rebin a 3D power to 2D kperp,kpar space, returns 2D power
;;
;; Inputs:
;;   power_3d: 3D power array
;;   kx_mpc: list of kx values, same size as 1st dimesion of power_3d
;;   ky_mpc: list of ky values, same size as 2nd dimesion of power_3d
;;   kz_mpc: list of kz values, same size as 3rd dimesion of power_3d
;;
;; Outputs:
;;   kperp_edges_mpc: edges of kperp bins, 1 element longer than 1st dimension of power_2d
;;   kpar_edges_mpc: edges of kpar bins, 1 element longer than 2nd dimension of power_2d
;;
;; Keywords:
;;   weights: array of inverse variances, same size as power_3d
;;   binned_weights: array of binned inverse variances, same size as power_2d
;;   noise_expval: array of expected noise values, same size as power_3d
;;   binned_weights: array of binned expected noise values, same size as power_2d
;;   nbins_2d: 2 element output vector returning [n_kperp, n_kpar]
;;   fill_holes: flag to fill in locations with no contributing data by copying data from lower kperp or kpar bin (default to set)
;;   kx_lims: limits to apply to kx values before binning -- only data between kx_lims will be binned
;;   ky_lims: limits to apply to ky values before binning -- only data between ky_lims will be binned
;;   log_kpar: bin with logrithmic kpar bins (default to linear bins)
;;   log_kperp: bin with logrithmic kperp bins (default to linear bins)
;;   kperp_bin: kperp bin size (for log bins per decade = 1/kperp_bin). If not set, the values used will be returned in this keyword
;;   kpar_bin: kpar bin size (for log bins per decade = 1/kpar_bin). If not set, the values used will be returned in this keyword
;;

function kspace_rebinning_2d, power_3d, kx_mpc, ky_mpc, kz_mpc, kperp_edges_mpc, kpar_edges_mpc, weights = weights, $
    binned_weights = weights_2d, noise_expval = noise_expval, binned_noise_expval = noise_expval_2d, $
    nbins_2d = nvox_2d, edge_on_grid = edge_on_grid, match_datta = match_datta, $
    fill_holes = fill_holes, kx_lims = kx_lims, ky_lims = ky_lims, $
    log_kpar = log_kpar, log_kperp = log_kperp, kperp_bin = kperp_bin, kpar_bin = kpar_bin, $
    kperp_density_measure = kperp_density_measure, kperp_density_cutoff = kperp_density_cutoff
    
  if n_elements(fill_holes) eq 0 then fill_holes = 1
  
  if n_elements(kpar_bin) gt 1 then message, 'kpar_bin must be a scalar because it specifies the linear bin size ' + $
    'or the log bin size (bins per decade = 1/log bin size) if log_kpar keyword is set. '
  if n_elements(kperp_bin) gt 1 then message, 'kperp_bin must be a scalar because it specifies the linear bin size ' + $
    'or the log bin size (bins per decade = 1/log bin size) if log_kperp keyword is set. '
    
  ;; for log binning, kpar/kperp_bins define log binsize (1/bins per decade). Default to 0.1 (10 bins per decade)
  if keyword_set(log_kpar) then begin
    if n_elements(kpar_bin) eq 0 then kpar_bin = 0.1d else if kpar_bin gt 1 then $
      message, 'kpar_bin must be less than one if log_kpar keyword is set (so there will be > 1 bin per decade)'
  endif else begin
    ;; for linear binning default to kz binsize
    if n_elements(kpar_bin) eq 0 then kpar_bin = kz_mpc[1]-kz_mpc[0]
  endelse
  
  if keyword_set(log_kperp) then begin
    if n_elements(kperp_bin) eq 0 then kperp_bin = 0.1d else if kperp_bin gt 1 then $
      message, 'kperp_bin must be less than one if log_kperp keyword is set (so there will be > 1 bin per decade)'
  endif else begin
    ;; for linear binning default to binsize = min(kx, ky binsize)
    if n_elements(kperp_bin) eq 0 then kperp_bin = max([kx_mpc[1]-kx_mpc[0], ky_mpc[1]-ky_mpc[0]])
  endelse
  
  dims = size(power_3d, /dimensions)
  input_type = size(power_3d, /type)
  
  if n_elements(dims) ne 3 then message, 'input power array must be 3 dimensional and ordered (x,y,z)'
  
  n_kx = dims[0]
  n_ky = dims[1]
  n_kz = dims[2]
  
  if n_elements(kx_mpc) ne n_kx then message, 'Length of kx_mpc must be the same as the first dimension of the power array'
  if n_elements(ky_mpc) ne n_ky then message, 'Length of ky_mpc must be the same as the second dimension of the power array'
  if n_elements(kz_mpc) ne n_kz then message, 'Length of kz_mpc must be the same as the third dimension of the power array'
  
  if n_elements(kperp_density_measure) gt 0 and n_elements(kperp_density_cutoff) gt 0 then begin
    if total(abs(size(kperp_density_measure, /dimensions) - dims[0:1])) ne 0 then $
      message, 'If kperp_density_measure is provided, it must have the same kperp dimensions as the power'
  endif
  
  if n_elements(weights) ne 0 then begin
    if total(abs(size(weights, /dimensions) - dims)) ne 0 then $
      message, 'If weights array is provided, it must have the same dimensionality as the power' $
    else weighted_power = weights * power_3d
  endif else begin
    weights = dblarr(dims) + 1d
    weighted_power = power_3d
  endelse
  
  if n_elements(noise_expval) ne 0 then begin
    if total(abs(size(noise_expval, /dimensions) - dims)) ne 0 then $
      message, 'If noise_expval array is provided, it must have the same dimensionality as the power' $
    else weighted_noise_expval = weights * noise_expval
  endif else begin
    noise_expval = 1/sqrt(weights)
    wh_wt0 = where(weights eq 0, count_wt0)
    if count_wt0 gt 0 then noise_expval[wh_wt0] = 0
    weighted_noise_expval = sqrt(weights)
  endelse
  
  ;; apply any pre-binning cuts
  if n_elements(kx_lims) gt 0 then begin
    if n_elements(kx_lims) gt 2 then message, 'kx_lims has too many elements' $
    else if n_elements(kx_lims) eq 1 then kx_lims = [-1, 1] * abs(kx_lims)
    
    wh = where(kx_mpc ge min(kx_lims) and kx_mpc le max(kx_lims), count)
    
    if count eq 0 then message, 'No data between kx_lims' else begin
      weighted_power = weighted_power[wh, *, *]
      weighted_noise_expval = weighted_noise_expval[wh, *, *]
      weights_use = weights[wh, *, *]
      kx_mpc_use = kx_mpc[wh]
      n_kx = count
      if n_elements(kperp_density_measure) gt 0 and n_elements(kperp_density_cutoff) gt 0 $
        then kperp_density_measure_use = kperp_density_measure[wh, *]
    endelse
  endif else begin
    kx_mpc_use = kx_mpc
    weights_use = weights
    if n_elements(kperp_density_measure) gt 0 and n_elements(kperp_density_cutoff) gt 0 $
      then kperp_density_measure_use = kperp_density_measure
  endelse
  
  if n_elements(ky_lims) gt 0 then begin
    if n_elements(ky_lims) gt 2 then message, 'ky_lims has too many elements' $
    else if n_elements(ky_lims) eq 1 then ky_lims = [-1, 1] * abs(ky_lims)
    
    wh = where(ky_mpc ge min(ky_lims) and ky_mpc le max(ky_lims), count)
    if count eq 0 then message, 'No data between ky_lims' else begin
      weighted_power = weighted_power
      weighted_noise_expval = weighted_noise_expval[*, wh , *]
      weights_use = weights_use[*, wh, *]
      ky_mpc_use = ky_mpc[wh]
      n_ky = count
      if n_elements(kperp_density_measure_use) gt 0 and n_elements(kperp_density_cutoff) gt 0 $
          then kperp_density_measure_use = kperp_density_measure_use[*, wh]
    endelse
  endif else ky_mpc_use = ky_mpc
  
  ;; apply density corrections if desired
  if n_elements(kperp_density_measure) gt 0 and n_elements(kperp_density_cutoff) gt 0 then begin
    wh_kperp_dense = where(kperp_density_measure_use gt kperp_density_cutoff, count_kperp_dense)
    
    if count_kperp_dense gt 0 then begin
      temp = reform(weighted_power, n_kx*n_ky, n_kz)
      temp[wh_kperp_dense, *] /= 0.5
      weighted_power = reform(temp, n_kx, n_ky, n_kz)
    endif else print, 'no kperp values exceed kperp_density_cutoff'
        
  endif
  
  if not keyword_set(log_kperp) then begin
  
    temp = sqrt(rebin(kx_mpc_use, n_kx, n_ky)^2 + rebin(reform(ky_mpc_use, 1, n_ky), n_kx, n_ky)^2)
    if keyword_set(edge_on_grid) then kperp_min = floor(min(temp) / kperp_bin) * kperp_bin $
    else kperp_min = floor((min(temp) + kperp_bin/2d) / kperp_bin) * kperp_bin - kperp_bin/2d
    
    kperp_array = temp
    
    ;; Use histogram with reverse indicies to bin in kperp
    kperp_hist = histogram(kperp_array, binsize = kperp_bin, min = kperp_min, omax = kperp_max, locations = lin_kperp_locs, $
      reverse_indices = kperp_ri)
      
    kperp_centers_mpc = lin_kperp_locs + kperp_bin/2d
    kperp_edges_mpc = [lin_kperp_locs, max(lin_kperp_locs) + kperp_bin]
    
  endif else begin
    ;; bin in logarithmic k_perp space
    kperp_log_binsize = kperp_bin
    temp = sqrt(rebin(kx_mpc_use, n_kx, n_ky)^2 + rebin(reform(ky_mpc_use, 1, n_ky), n_kx, n_ky)^2)
    
    ;; get min value, adjust so centers of bins will be on major grid.
    if min(temp) gt 0 then begin
      if keyword_set(edge_on_grid) then $
        kperp_min = floor((alog10(min(temp))) / kperp_log_binsize) * kperp_log_binsize $
      else $
        kperp_min = floor((alog10(min(temp)) + kperp_log_binsize/2d) / kperp_log_binsize) * kperp_log_binsize - $
        kperp_log_binsize/2d
      kperp_array = alog10(temp)
    endif else begin
      ;; Add an extra bin to put the k=0 mode into.
      if keyword_set(edge_on_grid) then $
        kperp_min = floor((alog10(min(temp[where(temp gt 0)]))) / kperp_log_binsize - 1) * kperp_log_binsize $
      else kperp_min = alog10(min(temp[where(temp gt 0)]))-3d*kperp_log_binsize/2d
      if keyword_set(match_datta) then begin
        if keyword_set(edge_on_grid) then kperp_min = -3d $
        else kperp_min = -3d - kperp_log_binsize/2d
      endif
      
      ;; calculate log kperp at each location in 1/Mpc. Set 0 mode to bottom value so it doesn't get left out.
      kperp_array = alog10(temp)
      kperp_array[where(temp le 0)] = kperp_min
    endelse
    
    ;; Use histogram with reverse indicies to bin in log kperp (want even spacing in log space.)
    kperp_hist = histogram(kperp_array, binsize = kperp_log_binsize, min = kperp_min, omax = kperp_max, locations = log_kperp_locs, $
      reverse_indices = kperp_ri)
      
    kperp_centers_mpc = 10^(log_kperp_locs + kperp_log_binsize/2d)
    kperp_edges_mpc = [10^(log_kperp_locs), 10^(max(log_kperp_locs) + kperp_log_binsize)]
  endelse
  
  n_kperp = n_elements(kperp_hist)
  weighted_power_mid = make_array(n_kperp, n_kz, type=input_type)
  weighted_nev_mid = make_array(n_kperp, n_kz, type=input_type)
  weights_mid = dblarr(n_kperp, n_kz)
  nvox_mid = lonarr(n_kperp, n_kz)
  
  reformed_wtpower = reform(temporary(weighted_power), n_kx*n_ky, n_kz)
  reformed_wtnev = reform(temporary(weighted_noise_expval), n_kx*n_ky, n_kz)
  reformed_weights = reform(temporary(weights_use), n_kx*n_ky, n_kz)
  for i=0, n_kperp-1 do begin
    if kperp_hist[i] gt 0 then begin
      weighted_power_mid[i, *] = total(reform(reformed_wtpower[kperp_ri[kperp_ri[i] : kperp_ri[i+1]-1], *], kperp_hist[i], n_kz), 1)
      weighted_nev_mid[i, *] = total(reform(reformed_wtnev[kperp_ri[kperp_ri[i] : kperp_ri[i+1]-1], *], kperp_hist[i], n_kz), 1)
      
      temp = reformed_weights[kperp_ri[kperp_ri[i] : kperp_ri[i+1]-1], *]
      weights_mid[i, *] = total(reform(temp, kperp_hist[i], n_kz), 1)
      
      if min(temp) le 0 then begin
        for j=0, n_kz -1 do begin
          wh_n0 = where(temp[*,j] gt 0, count_n0)
          nvox_mid[i, j] = count_n0
        endfor
      endif else nvox_mid[i, *] = kperp_hist[i]
      
    endif
  endfor
  
  undefine, reformed_wtpower
  undefine, reformed_wtnev
  undefine, reformed_weights
  
  if keyword_set(fill_holes) then begin
    wh = where(weights_mid eq 0, n0)
    for i=0, n0-1 do begin
      perp_ind = wh[i] mod n_kperp
      z_ind = wh[i] / n_kperp
      if perp_ind ne 0 then begin
        weighted_power_mid[perp_ind, z_ind] = weighted_power_mid[perp_ind-1, z_ind]
        weighted_nev_mid[perp_ind, z_ind] = weighted_nev_mid[perp_ind-1, z_ind]
        weights_mid[perp_ind, z_ind] = weights_mid[perp_ind-1, z_ind]
      endif
    endfor
  endif
  
  kpar_centers_mpc = kz_mpc
  kz_delta = kz_mpc[1] - kz_mpc[0]
  if not keyword_set(log_kpar) then begin
  
    ;; kpar_edges_mpc = [kz_mpc[0] - kz_delta/2d, kz_mpc + kz_delta/2d]
  
    ;; weighted_power_2d = temporary(weighted_power_mid)
    ;; weights_2d = temporary(weights_mid)
  
    if keyword_set(edge_on_grid) then kpar_min = floor(min(kz_mpc) / kpar_bin) * kpar_bin $
    else kpar_min = floor((min(kz_mpc) + kpar_bin/2d) / kpar_bin) * kpar_bin - kpar_bin/2d
    
    ;; Use histogram with reverse indicies to bin in kpar
    kpar_hist = histogram(kz_mpc, binsize = kpar_bin, min = kpar_min, omax = kpar_max, locations = lin_kpar_locs, $
      reverse_indices = kpar_ri)
      
    kpar_centers_mpc = lin_kpar_locs + kpar_bin/2d
    kpar_edges_mpc = [lin_kpar_locs, max(lin_kpar_locs) + kpar_bin]
    
  endif else begin
    ;; now make kz bins logrithmic by resampling
    kpar_log_binsize = kpar_bin
    
    ;; get min value, adjust so centers of bins will be on major grid
    if min(kz_mpc) gt 0 then begin
      if keyword_set(edge_on_grid) then $
        kpar_min = floor((alog10(min(kz_mpc))) / kpar_log_binsize) * kpar_log_binsize $
      else $
        kpar_min = floor((alog10(min(kz_mpc)) + kpar_log_binsize/2d) / kpar_log_binsize) * $
        kpar_log_binsize - kpar_log_binsize/2d
        
      kpar_arr = alog10(kz_mpc)
    endif else begin
      if keyword_set(edge_on_grid) then $
        kpar_min = floor((alog10(min(kz_mpc[where(kz_mpc gt 0)]))) / kpar_log_binsize - 1) * kpar_log_binsize $
      else kpar_min = alog10(min(kz_mpc[where(kz_mpc gt 0)]))-3d*kpar_log_binsize/2d
      kpar_arr = alog10(kz_mpc)
      kpar_arr[where(kz_mpc le 0)] = kpar_min
    endelse
    
    kpar_hist = histogram(kpar_arr, binsize = kpar_log_binsize, min = kpar_min, omax = kpar_max, locations = log_kpar_locs, $
      reverse_indices = kpar_ri)
    kpar_centers_mpc = 10^(log_kpar_locs + kpar_log_binsize/2d)
    kpar_edges_mpc = [10^(log_kpar_locs), 10^(max(log_kpar_locs) + kpar_log_binsize)]
    kz_mpc_edges = [kz_mpc, max(kz_mpc) + kz_mpc[1] - kz_mpc[0]]
  endelse
  
  n_kpar = n_elements(kpar_hist)
  weighted_power_2d = make_array(n_kperp, n_kpar, type=input_type)
  weighted_nev_2d = make_array(n_kperp, n_kpar, type=input_type)
  weights_2d = dblarr(n_kperp, n_kpar)
  nvox_2d = lonarr(n_kperp, n_kpar)
  
  for j=0, n_kpar-1 do begin
    if kpar_hist[j] ne 0 then begin
      weighted_power_2d[*,j] = total(reform(weighted_power_mid[*, kpar_ri[kpar_ri[j] : kpar_ri[j+1]-1]], n_kperp, kpar_hist[j]),2)
      weighted_nev_2d[*,j] = total(reform(weighted_nev_mid[*, kpar_ri[kpar_ri[j] : kpar_ri[j+1]-1]], n_kperp, kpar_hist[j]),2)
      
      temp = reform(weights_mid[*, kpar_ri[kpar_ri[j] : kpar_ri[j+1]-1]], n_kperp, kpar_hist[j])
      weights_2d[*, j] = total(temp, 2)
      
      if min(temp) le 0 then begin
        wh_n0 = where(temp gt 0, count_n0)
        mask = intarr(n_kperp, kpar_hist[j])
        if count_n0 gt 0 then mask[wh_n0] = 1
        
        nvox_2d[*, j] = total(reform(nvox_mid[*, kpar_ri[kpar_ri[j] : kpar_ri[j+1]-1]]*mask, n_kperp, kpar_hist[j]),2)
      endif else nvox_2d[*, j] = total(reform(nvox_mid[*, kpar_ri[kpar_ri[j] : kpar_ri[j+1]-1]], n_kperp, kpar_hist[j]),2)
    endif
  endfor
  
  undefine, weighted_power_mid
  undefine, weighted_nev_mid
  undefine, weights_mid
  undefine, nvox_mid
  
  if keyword_set(fill_holes) then begin
    wh = where(weights_2d eq 0, n0)
    for j=0, n0-1 do begin
      perp_ind = wh[j] mod n_kperp
      par_ind = wh[j] / n_kperp
      if par_ind ne 0 then begin
        weighted_power_2d[perp_ind, par_ind] = weighted_power_2d[perp_ind, par_ind-1]
        weighted_nev_2d[perp_ind, par_ind] = weighted_nev_2d[perp_ind, par_ind-1]
        weights_2d[perp_ind, par_ind] = weights_2d[perp_ind, par_ind-1]
      endif
    endfor
  endif
  
  power_ave = weighted_power_2d/weights_2d
  noise_expval_2d = weighted_nev_2d/weights_2d
  wh = where(weights_2d eq 0, count)
  if count ne 0 then begin
    power_ave[wh] = 0d
    noise_expval_2d[wh] = 0d
  endif
  
  return, power_ave
  
end
