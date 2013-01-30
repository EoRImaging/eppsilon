

function kspace_rebinning_2d, power_3d, kx_mpc, ky_mpc, kz_mpc, kperp_edges_mpc, kpar_edges_mpc, weights = weights, $
                              binned_weights = weights_2d, bins_per_decade = bins_per_decade, edge_on_grid = edge_on_grid, $
                              match_datta = match_datta, fill_holes = fill_holes, kx_lims = kx_lims, ky_lims = ky_lims, $
                              linear_kpar = linear_kpar, linear_kperp = linear_kperp, linkperp_bin = linkperp_bin

  if n_elements(fill_holes) eq 0 then fill_holes = 1

  if n_elements(bins_per_decade) eq 0 then bins_per_decade = 10d
  if bins_per_decade lt 1 or n_elements(bins_per_decade) ne 1 then begin
     print, 'Number of kspace bins per decade must be a scalar greater than 1. Defaulting to 10.'
     bins_per_decade = 10
  endif

  dims = size(power_3d, /dimensions)
  input_type = size(power_3d, /type)

  if n_elements(dims) ne 3 then message, 'input power array must be 3 dimensional and ordered (x,y,z)'

  n_kx = dims[0]
  n_ky = dims[1]
  n_kz = dims[2]

  if n_elements(kx_mpc) ne n_kx then message, 'Length of kx_mpc must be the same as the first dimension of the power array'
  if n_elements(ky_mpc) ne n_ky then message, 'Length of ky_mpc must be the same as the second dimension of the power array'
  if n_elements(kz_mpc) ne n_kz then message, 'Length of kz_mpc must be the same as the third dimension of the power array'

  if n_elements(weights) ne 0 then begin
     if total(abs(size(weights, /dimensions) - dims)) ne 0 then $
        message, 'If weights array is provided, it must have the same dimensionality as the power' $
     else weighted_power = weights * power_3d
  endif else begin
     weights = dblarr(dims) + 1d
     weighted_power = power_3d
  endelse

  ;; apply any pre-binning cuts
  if n_elements(kx_lims) gt 0 then begin
     if n_elements(kx_lims) gt 2 then message, 'kx_lims has too many elements' $
     else if n_elements(kx_lims) eq 1 then kx_lims = [-1, 1] * abs(kx_lims)

     wh = where(kx_mpc ge min(kx_lims) and kx_mpc le max(kx_lims), count)

     if count eq 0 then message, 'No data between kx_lims' $
     else begin
        weighted_power = weighted_power[wh, *, *]
        weights_use = weights[wh, *, *]
        kx_mpc_use = kx_mpc[wh]
        n_kx = count
     endelse
  endif else begin
     kx_mpc_use = kx_mpc
     weights_use = weights
  endelse

  if n_elements(ky_lims) gt 0 then begin
     if n_elements(ky_lims) gt 2 then message, 'ky_lims has too many elements' $
     else if n_elements(ky_lims) eq 1 then ky_lims = [-1, 1] * abs(ky_lims)

     wh = where(ky_mpc ge min(ky_lims) and ky_mpc le max(ky_lims), count)
     if count eq 0 then message, 'No data between ky_lims' $
     else begin
        weighted_power = weighted_power[*, wh , *]
        weights_use = weights_use[*, wh, *]
        ky_mpc_use = ky_mpc[wh]
        n_ky = count
     endelse
  endif else ky_mpc_use = ky_mpc

  if keyword_set(linear_kperp) then begin
     ;; bin in linear kperp. default binsize = min(kx, ky binsize)
     if n_elements(kperp_lin_binsize) eq 0 then kperp_lin_binsize = min([kx_mpc_use[1]-kx_mpc_use[0], ky_mpc_use[1]-ky_mpc_use[0]])

     temp = sqrt(rebin(kx_mpc_use, n_kx, n_ky)^2 + rebin(reform(ky_mpc_use, 1, n_ky), n_kx, n_ky)^2)
     if keyword_set(edge_on_grid) then kperp_min = floor(min(temp) / kperp_lin_binsize) * kperp_lin_binsize $
     else kperp_min = floor((min(temp) + kperp_lin_binsize/2d) / kperp_lin_binsize) * kperp_lin_binsize - kperp_lin_binsize/2d

     kperp_array = temp

     ;; Use histogram with reverse indicies to bin in kperp
     kperp_hist = histogram(kperp_array, binsize = kperp_lin_binsize, min = kperp_min, omax = kperp_max, locations = lin_kperp_locs, $
                            reverse_indices = kperp_ri)

     kperp_centers_mpc = lin_kperp_locs + kperp_lin_binsize/2d
     kperp_edges_mpc = [lin_kperp_locs, max(lin_kperp_locs) + kperp_lin_binsize]
 
  endif else begin
     ;; bin in logarithmic k_perp space
     kperp_log_binsize = 1d/bins_per_decade
     temp = sqrt(rebin(kx_mpc_use, n_kx, n_ky)^2 + rebin(reform(ky_mpc_use, 1, n_ky), n_kx, n_ky)^2)
     
     ;; get min value, adjust so centers of bins will be on major grid. 
     if min(temp) gt 0 then begin
        if keyword_set(edge_on_grid) then $
           kperp_min = floor((alog10(min(temp))) / kperp_log_binsize) * kperp_log_binsize $ 
        else kperp_min = floor((alog10(min(temp)) + kperp_log_binsize/2d) / kperp_log_binsize) * $
                         kperp_log_binsize - kperp_log_binsize/2d
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
  weights_mid = dblarr(n_kperp, n_kz)

  reformed_wtpower = reform(temporary(weighted_power), n_kx*n_ky, n_kz)
  reformed_weights = reform(temporary(weights_use), n_kx*n_ky, n_kz)
  for i=0, n_kperp-1 do begin
     if kperp_hist[i] gt 0 then begin
        weighted_power_mid[i, *] = total(reform(reformed_wtpower[kperp_ri[kperp_ri[i] : kperp_ri[i+1]-1], *], kperp_hist[i], n_kz), 1)
        weights_mid[i, *] = total(reform(reformed_weights[kperp_ri[kperp_ri[i] : kperp_ri[i+1]-1], *], kperp_hist[i], n_kz), 1)
     endif 
  endfor
  undefine, reformed_wtpower
  undefine, reformed_weights

  if keyword_set(fill_holes) then begin
     wh = where(weights_mid eq 0, n0)
     for i=0, n0-1 do begin
        perp_ind = wh[i] mod n_kperp
        z_ind = wh[i] / n_kperp
        if perp_ind ne 0 then begin
              weighted_power_mid[perp_ind, z_ind] = weighted_power_mid[perp_ind-1, z_ind]
              weights_mid[perp_ind, z_ind] = weights_mid[perp_ind-1, z_ind]
           endif
     endfor
  endif

  kpar_centers_mpc = kz_mpc
  kz_delta = kz_mpc[1] - kz_mpc[0]
  if keyword_set(linear_kpar) then begin
     kpar_edges_mpc = [kz_mpc[0] - kz_delta/2d, kz_mpc + kz_delta/2d]

     weighted_power_2d = temporary(weighted_power_mid)
     weights_2d = temporary(weights_mid)
  endif else begin
     ;; now make kz bins logrithmic by resampling
     kpar_log_binsize = 1d/bins_per_decade
     
     ;; get min value, adjust so centers of bins will be on major grid
     if min(kz_mpc) gt 0 then begin
        if keyword_set(edge_on_grid) then $
           kpar_min = floor((alog10(min(kz_mpc))) / kpar_log_binsize) * kpar_log_binsize $
        else kpar_min = floor((alog10(min(kz_mpc)) + kpar_log_binsize/2d) / kpar_log_binsize) * $
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
     n_kpar = n_elements(kpar_hist)
     kpar_centers_mpc = 10^(log_kpar_locs + kpar_log_binsize/2d)
     kpar_edges_mpc = [10^(log_kpar_locs), 10^(max(log_kpar_locs) + kpar_log_binsize)]
     kz_mpc_edges = [kz_mpc, max(kz_mpc) + kz_mpc[1] - kz_mpc[0]] 
     
     weighted_power_2d = make_array(n_kperp, n_kpar, type=input_type)
     weights_2d = dblarr(n_kperp, n_kpar)
     
     

     for j=0, n_kpar-1 do begin
        if kpar_hist[j] ne 0 then begin
           weighted_power_2d[*,j] = total(reform(weighted_power_mid[*, kpar_ri[kpar_ri[j] : kpar_ri[j+1]-1]], n_kperp, kpar_hist[j]),2)
           weights_2d[*, j] = total(reform(weights_mid[*, kpar_ri[kpar_ri[j] : kpar_ri[j+1]-1]], n_kperp, kpar_hist[j]), 2)
        endif
     endfor
     undefine, weighted_power_mid
     undefine, weights_mid

     if keyword_set(fill_holes) then begin
           wh = where(weights_2d eq 0, n0)
           for j=0, n0-1 do begin
              perp_ind = wh[j] mod n_kperp
              par_ind = wh[j] / n_kperp
              if par_ind ne 0 then begin
                 weighted_power_2d[perp_ind, par_ind] = weighted_power_2d[perp_ind, par_ind-1]
                 weights_2d[perp_ind, par_ind] = weights_2d[perp_ind, par_ind-1]
              endif
           endfor
        endif

  endelse     
     
  power_ave = weighted_power_2d/weights_2d
  wh = where(weights_2d eq 0, count)
  if count ne 0 then power_ave[wh] = 0d 
 
  return, power_ave

end 
