

function kspace_rebinning_1d, power, k1_mpc, k2_mpc, k3_mpc, k_edges_mpc, k_bin = k_bin, log_k = log_k, $
                              noise_expval = noise_expval, binned_noise_expval = noise_expval_1d, weights = weights, $
                              binned_weights = weights_1d, mask = mask, pixelwise_mask = pixelwise_mask, k1_mask = k1_mask, $
                              k2_mask = k2_mask, k3_mask = k3_mask, edge_on_grid = edge_on_grid, match_datta = match_datta

  power_size = size(power, /dimensions)
  power_dim = n_elements(power_size)
  
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


  if n_elements(weights) ne 0 then begin
     if total(abs(size(weights, /dimensions) - power_size)) ne 0 then $
        message, 'If weights array is provided, it must have the same dimensionality as the power' $
     else weighted_power = weights * power
  endif else begin
     weights = dblarr(power_size) + 1d
     weighted_power = power
  endelse
 
 if n_elements(noise_expval) ne 0 then begin
     if total(abs(size(noise_expval, /dimensions) - power_size)) ne 0 then $
        message, 'If noise_expval array is provided, it must have the same dimensionality as the power' $
     else weighted_noise_expval = weights * noise_expval
  endif else begin
     noise_expval = 1/sqrt(weights)
     weighted_noise_expval = sqrt(weights)
  endelse

  if n_elements(mask) ne 0 then begin
     mask_size = size(mask, /dimensions)
     mask_dim = n_elements(mask_size)

     if mask_dim gt power_dim then message, 'Mask array cannot have more dimensions than power array'
     if mask_dim lt 2 then message, 'Mask array must be 2 or 3 dimensional and ordered (kperpendicular, kparallel) or (kx,ky,kz)'

     wh = where(mask ne 0 and mask ne 1, count)
     if size(mask, /type) ne 2 or count ne 0 then message, "Mask must be an integer array with binary values (only 0's and 1's)"

     if keyword_set(pixelwise_mask) then begin
        if total(abs(mask_size - power_size)) gt 0 then $
           message, 'If pixelwise_mask keyword is set, the mask and power arrays must have the same size ' + $
                    '(so the mask can be applied pixel-by-pixel to the power array)'
     endif else begin
        ;; if not pixelwise, then kvalues should be edges (not centers) so there should be n+1 k elements
        if mask_dim eq 3 then begin
           mask_kx_edges = k1_mask
           mask_ky_edges = k2_mask
           mask_kz_edges = k3_mask
           
           n_mask_kx = mask_size[0]
           n_mask_ky = mask_size[1]
           n_mask_kz = mask_size[2]
           
           if n_elements(mask_kx_edges)-1 ne n_mask_kx then $
              message, 'Length of mask_k1 must be one greater than the first dimension of the mask array (giving k bin edges)'
           if n_elements(mask_ky_edges)-1 ne n_mask_ky then $
              message, 'Length of mask_k2 must be one greater than the second dimension of the mask array (giving k bin edges)'
           if n_elements(mask_kz_edges)-1 ne n_mask_kz then $
              message, 'Length of mask_k3 must be one greater than the third dimension of the mask array (giving k bin edges)'
        endif else begin
           mask_kperp_edges = k1_mask
           mask_kpar_edges = k2_mask
           
           n_mask_kperp = mask_size[0]
           n_mask_kpar = mask_size[1]
           
           if n_elements(mask_kperp_edges)-1 ne n_mask_kperp then $
              message, 'Length of mask_k1 must be one greater than the first dimension of the mask array (giving k bin edges)'
           if n_elements(mask_kpar_edges)-1 ne n_mask_kpar then $
              message, 'Length of mask_k2 must be one greater than the second dimension of the mask array (giving k bin edges)'
           
        endelse
     endelse
  endif


  if power_dim eq 3 then begin

     kx_mpc = k1_mpc
     ky_mpc = k2_mpc
     kz_mpc = k3_mpc

     n_kx = power_size[0]
     n_ky = power_size[1]
     n_kz = power_size[2]
     
     if n_elements(kx_mpc) ne n_kx then message, 'Length of k1_mpc must be the same as the first dimension of the power array'
     if n_elements(ky_mpc) ne n_ky then message, 'Length of k2_mpc must be the same as the second dimension of the power array'
     if n_elements(kz_mpc) ne n_kz then message, 'Length of k3_mpc must be the same as the third dimension of the power array'
     
     print, 'Generating array of k values'
     temp = sqrt(rebin(kx_mpc^2d, n_kx, n_ky, n_kz) + rebin(reform(ky_mpc^2d, 1, n_ky), n_kx, n_ky, n_kz) + $
                 rebin(reform(kz_mpc^2d, 1, 1, n_kz), n_kx, n_ky, n_kz))

     if not keyword_set(log_k) then begin

        if keyword_set(edge_on_grid) then k_min = floor(min(temp) / k_bin) * k_bin $
        else kperp_min = floor((min(temp) + k_bin/2d) / k_bin) * k_bin - k_bin/2d

        ;; Use histogram with reverse indicies to bin in kperp
        k_hist = histogram(temp, binsize = k_bin, min = k_min, omax = k_max, locations = lin_k_locs, $
                               reverse_indices = k_ri)

     endif else begin
        ;; get min value, adjust so centers of bins will be on major grid.
        min_val = sqrt(min(kx_mpc^2d) + min(ky_mpc^2d) + min(kz_mpc^2d))
        if min_val gt 0 then begin
           min_non_zero = min_val
           if keyword_set(edge_on_grid) then $
              k_min = floor(alog10(min_non_zero) / k_bin) * k_bin $ 
           else k_min = floor((alog10(min_non_zero) + k_bin/2d) / k_bin) * $
                        k_bin - k_bin/2d
        endif else begin
           min_candidates = sqrt([min(kx_mpc[where(kx_mpc^2d ne 0)]^2d) + min(ky_mpc^2d) + min(kz_mpc^2d), $
                                  min(kx_mpc^2d) + min(ky_mpc[where(ky_mpc^2d ne 0)]^2d) + min(kz_mpc^2d), $
                                  min(kx_mpc^2d) + min(ky_mpc^2d) + min(kz_mpc[where(kz_mpc^2d ne 0)]^2d)])
           min_non_zero = min(min_candidates)
           
           ;; Add an extra bin to put the k=0 mode into.
           if keyword_set(edge_on_grid) then $
              k_min = floor(alog10(min_non_zero) / k_bin - 1) * k_bin $
           else k_min = floor((alog10(min_non_zero) + k_bin/2d) / k_bin - 1) * $
                        k_bin - k_bin/2d
           ;; if keyword_set(match_datta) then begin
           ;;    if keyword_set(edge_on_grid) then k_min = -3d $
           ;;    else k_min = -3d - k_log_binsize/2d
           ;; endif
           
           ;;Set 0 mode to bottom value so it doesn't get left out.
           wh_kx_0 = where(kx_mpc^2d eq 0, nkx_0)
           wh_ky_0 = where(ky_mpc^2d eq 0, nky_0)
           wh_kz_0 = where(kz_mpc^2d eq 0, nkz_0)

           if nkx_0*nky_0*nkz_0 gt 1 then stop
        
           temp[wh_kx_0, wh_ky_0, wh_kz_0] = 10^k_min
        endelse

        ;; calculate log kperp at each location in 1/Mpc.
        k_array = alog10(temporary(temp))
        
        print, 'Histogramming array of log k values'
        ;; Use histogram with reverse indicies to bin in log k (want even spacing in log space.)
        k_hist = histogram(temporary(k_array), binsize = k_bin, min = k_min, omax = k_max, locations = log_k_locs, $
                           reverse_indices = k_ri)
     endelse
     if n_elements(mask) ne 0 then begin
        if keyword_set(pixelwise_mask) then pixel_mask = mask $
        else begin
           pixel_mask = intarr(size(power, /dimensions))
           if mask_dim eq 2 then perp_arr = sqrt(rebin(kx_mpc, n_kx, n_ky)^2d + rebin(reform(ky_mpc, 1, n_ky), n_kx, n_ky)^2)
           wh_mask_gt0 = where(mask gt 0, n_gt0)
           for i=0, n_gt0 -1 do begin
              mask_inds = array_indices(mask, wh_mask_gt0[i])
              
              if mask_dim eq 3 then begin
                 x_inds = where(kx_mpc gt mask_kx_edges[mask_inds[0]] and kx_mpc lt mask_kx_edges[mask_inds[0]+1], nx_inds)
                 y_inds = where(ky_mpc gt mask_ky_edges[mask_inds[1]] and ky_mpc lt mask_ky_edges[mask_inds[1]+1], ny_inds)
                 z_inds = where(kz_mpc gt mask_kz_edges[mask_inds[2]] and kz_mpc lt mask_kz_edges[mask_inds[2]+1], nz_inds)
                 
                 if (nx_inds*ny_inds*nz_inds) ne 0 then begin
                    inds = rebin(x_inds, nx_inds, ny_inds, nz_inds) + $
                           rebin(reform(y_inds, 1, ny_inds), nx_inds, ny_inds, nz_inds)*n_kx + $
                           rebin(reform(z_inds, 1, 1, nz_inds), nx_inds, ny_inds, nz_inds)*n_kx*n_ky
                    inds = reform(inds, n_elements(inds))
                    
                    pixel_mask[inds] = mask[wh_mask_gt0[i]]
                 endif 
              endif else begin
                 perp_inds = where(perp_arr gt mask_kperp_edges[mask_inds[0]] and perp_arr lt mask_kperp_edges[mask_inds[0]+1], $
                                   nperp_inds)
                 if nperp_inds ne 0 then begin
                    temp = array_indices(perp_arr, perp_inds)
                    x_inds = reform(temp[0,*])
                    y_inds = reform(temp[1,*])
                    
                    z_inds = where(kz_mpc gt mask_kpar_edges[mask_inds[1]] and kz_mpc lt mask_kpar_edges[mask_inds[1]+1], nz_inds)
                    
                    if nz_inds ne 0 then begin
                       for j=0, nperp_inds-1 do begin
                          if j eq 0 then inds = rebin([x_inds[j]], nz_inds) + rebin([y_inds[j]], nz_inds)*n_kx + z_inds*n_kx*n_ky $
                          else inds = [inds, rebin([x_inds[j]], nz_inds) + rebin([y_inds[j]], nz_inds)*n_kx + z_inds*n_kx*n_ky]
                       endfor
                       inds = reform(inds, n_elements(inds))
                       
                       pixel_mask[inds] = mask[wh_mask_gt0[i]]
                    endif
                 endif
              endelse               
           endfor
        endelse
     endif
     
     n_k = n_elements(k_hist)
     power_1d = dblarr(n_k)
     nev_1d = dblarr(n_k)
     weights_1d = dblarr(n_k)
     norm = dblarr(n_k)

     for i=0, n_k-1 do begin
        if k_hist[i] ne 0 then begin
           inds =  k_ri[k_ri[i] : k_ri[i+1]-1]
           if n_elements(mask) ne 0 then begin
              power_1d[i] = total(weighted_power[inds] * double(pixel_mask[inds]))
              nev_1d[i] = total(weighted_noise_expval[inds] * double(pixel_mask[inds]))
              weights_1d[i] = total(weights[inds] * double(pixel_mask[inds]))
              norm[i] = total(pixel_mask[inds])
           endif else begin
              power_1d[i] = total(weighted_power[inds])
              nev_1d[i] = total(weighted_noise_expval[inds])
              weights_1d[i] = total(weights[inds])
              norm[i] = k_hist[i]
           endelse
        endif 
     endfor
     k_ri=0

  endif else begin
     
     kperp_mpc = k1_mpc
     kpar_mpc = k2_mpc
     k3_mpc = 0
  
     n_kperp = power_size[0]
     n_kpar = power_size[1]
     
     if n_elements(kperp_mpc) ne n_kperp then message, 'Length of k1_mpc must be the same as the first dimension of the power array'
     if n_elements(kpar_mpc) ne n_kpar then message, 'Length of k2_mpc must be the same as the second dimension of the power array'
     
     temp = rebin(kperp_mpc, n_kperp, n_kpar)^2 + rebin(reform(kpar_mpc, 1, n_kpar), n_kperp, n_kpar)^2

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
           else k_min = floor((alog10(min_non_zero) + k_bin/2d) / k_bin) * $
                        k_bin - k_bin/2d
        endif else begin
           wh0 = where(temp eq 0, complement = wh_non0)
           min_non_zero = min(temp[wh_non0])
           ;; Add an extra bin to put the k=0 mode into.
           if keyword_set(edge_on_grid) then $
              k_min = floor(alog10(min_non_zero) / k_bin - 1) * k_bin $
           else k_min = floor((alog10(min_non_zero) + k_bin/2d) / k_bin - 1) * $
                        k_bin - k_bin/2d
           
           ;;Set 0 mode to bottom value so it doesn't get left out.
           temp[wh0] = 10^k_min
        endelse
        
        ;; calculate log kperp at each location in 1/Mpc.
        k_array = alog10(temporary(temp))
        
        ;; Use histogram with reverse indicies to bin in log k (want even spacing in log space.)
        k_hist = histogram(temporary(k_array), binsize = k_bin, min = k_min, omax = k_max, locations = log_k_locs, $
                           reverse_indices = k_ri)
     endelse

     if n_elements(mask) ne 0 then begin
        if keyword_set(pixelwise_mask) then pixel_mask = mask $
        else begin
           pixel_mask = intarr(size(power, /dimensions))
           wh_mask_gt0 = where(mask gt 0, n_gt0)
           for i=0, n_gt0 -1 do begin
              mask_inds = array_indices(mask, wh_mask_gt0[i])
              
              x_inds = where(kx_mpc gt mask_kx_edges[mask_inds[0]] and kx_mpc lt mask_kx_edges[mask_inds[0]+1], nx_inds)
              y_inds = where(ky_mpc gt mask_ky_edges[mask_inds[1]] and ky_mpc lt mask_ky_edges[mask_inds[1]+1], ny_inds)
              
              if (nx_inds*ny_inds) ne 0 then begin
                 inds = rebin(x_inds, nx_inds, ny_inds) + rebin(reform(y_inds, 1, ny_inds), nx_inds, ny_inds)*n_kx 
                 inds = reform(inds, n_elements(inds))
                 
                 pixel_mask[inds] = mask[wh_mask_gt0[i]]
              endif
           endfor
        endelse
     endif


     n_k = n_elements(k_hist)
     power_1d = dblarr(n_k)
     nev_1d = dblarr(n_k)
     weights_1d = dblarr(n_k)
     norm = dblarr(n_k)

     for i=0, n_k-1 do begin
        if k_hist[i] ne 0 then begin
           inds = k_ri[k_ri[i] : k_ri[i+1]-1]
           if n_elements(mask) ne 0 then begin
              power_1d[i] = total(weighted_power[inds] * double(pixel_mask[inds]))
              nev_1d[i] = total(weighted_noise_expval[inds] * double(pixel_mask[inds]))
              weights_1d[i] = total(weights[inds] * double(pixel_mask[inds]))
              norm[i] = total(pixel_mask[inds])
           endif else begin
              power_1d[i] = total(weighted_power[inds])
              nev_1d[i] = total(weighted_noise_expval[inds])
              weights_1d[i] = total(weights[inds])
              norm[i] = k_hist[i]
           endelse
        endif 
     endfor
     k_ri=0
          
  endelse

  if keyword_set(log_k) then begin
     k_centers_mpc = 10^(log_k_locs + k_bin/2d)
     k_edges_mpc = [10^(log_k_locs), 10^(max(log_k_locs) + k_bin)]
  endif else begin
     k_centers_mpc = lin_k_locs + k_bin/2d
     k_edges_mpc = [lin_k_locs, max(lin_k_locs) + k_bin]
  endelse


  power_ave = power_1d/weights_1d
  noise_expval_1d = nev_1d/weights_1d
  wh = where(weights_1d eq 0, count)
  if count ne 0 then begin
     power_ave[wh] = 0d
     noise_expval_1d[wh] = 0d
  endif

  return, power_ave
  
end 
