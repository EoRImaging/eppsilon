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
;;   kx_lims: limits to apply to kx values before binning -- only data between kx_lims will be binned
;;   ky_lims: limits to apply to ky values before binning -- only data between ky_lims will be binned
;;   binning_2d_options: keywords that control binning, defined in create_binning_2d_options

function kspace_rebinning_2d, power_3d, kx_mpc, ky_mpc, kz_mpc, kperp_edges_mpc, $
    kpar_edges_mpc, weights = weights, binned_weights = weights_2d, $
    noise_expval = noise_expval, binned_noise_expval = noise_expval_2d, $
    var_power_2d = var_power_2d, mean_var_2d = mean_var_2d, $
    nbins_2d = nvox_2d, edge_on_grid = edge_on_grid, match_datta = match_datta, $
    kx_lims = kx_lims, ky_lims = ky_lims, $
    kperp_density_measure = kperp_density_measure, kperp_density_cutoff = kperp_density_cutoff, $
    binning_2d_options = binning_2d_options

  if tag_exist(binning_2d_options, 'kpar_bin') then begin
    if n_elements(binning_2d_options.kpar_bin) gt 1 then begin
      message, 'kpar_bin must be a scalar because it specifies the linear bin size ' + $
        'or the log bin size (bins per decade = 1/log bin size) if log_kpar keyword is set. '
    endif
  endif
  if tag_exist(binning_2d_options, 'kperp_bin') then begin
    if n_elements(binning_2d_options.kperp_bin) gt 1 then begin
      message, 'kperp_bin must be a scalar because it specifies the linear bin size ' + $
        'or the log bin size (bins per decade = 1/log bin size) if log_kperp keyword is set. '
    endif
  endif
  ;; for log binning, kpar/kperp_bins define log binsize (1/bins per decade).
  ;; Default to 0.1 (10 bins per decade)
  if binning_2d_options.log_kpar then begin
    if not tag_exist(binning_2d_options, 'kpar_bin') then begin
      binning_2d_options = create_struct(binning_2d_options, 'kpar_bin', 0.1d)
    endif else if binning_2d_options.kpar_bin gt 1 then begin
      message, 'kpar_bin must be less than one if log_kpar keyword is set ' + $
        '(so there will be > 1 bin per decade)'
    endif
  endif else begin
    ;; for linear binning default to kz binsize
    if not tag_exist(binning_2d_options, 'kpar_bin') then begin
      binning_2d_options = create_struct(binning_2d_options, 'kpar_bin', kz_mpc[1]-kz_mpc[0])
    endif
  endelse

  if binning_2d_options.log_kperp then begin
    if not tag_exist(binning_2d_options, 'kperp_bin') then begin
      binning_2d_options = create_struct(binning_2d_options, 'kperp_bin', 0.1d)
    endif else if binning_2d_options.kperp_bin gt 1 then begin
      message, 'kperp_bin must be less than one if log_kperp keyword is set ' + $
        '(so there will be > 1 bin per decade)'
    endif
  endif else begin
    ;; for linear binning default to binsize = min(kx, ky binsize)
    if not tag_exist(binning_2d_options, 'kperp_bin') then begin
      binning_2d_options = create_struct(binning_2d_options, 'kperp_bin', $
        max([kx_mpc[1]-kx_mpc[0], ky_mpc[1]-ky_mpc[0]]))
    endif
  endelse

  dims = size(power_3d, /dimensions)
  input_type = size(power_3d, /type)

  if n_elements(dims) ne 3 then begin
    message, 'input power array must be 3 dimensional and ordered (x,y,z)'
  endif

  n_kx = dims[0]
  n_ky = dims[1]
  n_kz = dims[2]

  if n_elements(kx_mpc) ne n_kx then begin
    message, 'Length of kx_mpc must be the same as the first dimension of the power array'
  endif
  if n_elements(ky_mpc) ne n_ky then begin
    message, 'Length of ky_mpc must be the same as the second dimension of the power array'
  endif
  if n_elements(kz_mpc) ne n_kz then begin
    message, 'Length of kz_mpc must be the same as the third dimension of the power array'
  endif

  if n_elements(kperp_density_measure) gt 0 and n_elements(kperp_density_cutoff) gt 0 then begin
    if total(abs(size(kperp_density_measure, /dimensions) - dims[0:1])) ne 0 then $
      message, 'If kperp_density_measure is provided, it must have the same kperp dimensions as the power'
  endif

  if n_elements(weights) ne 0 then begin
    if total(abs(size(weights, /dimensions) - dims)) ne 0 then begin
      message, 'If weights array is provided, it must have the same dimensionality as the power'
    endif

    power_use = power_3d
    weights_use = weights
  endif else begin
    weights_use = dblarr(dims) + 1d
    power_use = power_3d
  endelse

  if n_elements(noise_expval) ne 0 then begin
    if total(abs(size(noise_expval, /dimensions) - dims)) ne 0 then begin
      message, 'If noise_expval array is provided, it must have the same ' + $
        'dimensionality as the power'
    endif else begin
      noise_expval_use = noise_expval
    endelse
  endif else begin
    noise_expval_use = 1/sqrt(weights_use)
    wh_wt0 = where(weights_use eq 0, count_wt0)
    if count_wt0 gt 0 then noise_expval_use[wh_wt0] = 0
  endelse

  ;; apply any pre-binning cuts
  if n_elements(kx_lims) gt 0 then begin
    if n_elements(kx_lims) gt 2 then begin
      message, 'kx_lims has too many elements'
    endif else begin
      if n_elements(kx_lims) eq 1 then kx_lims = [-1, 1] * abs(kx_lims)
    endelse

    wh = where(kx_mpc ge min(kx_lims) and kx_mpc le max(kx_lims), count)

    if count eq 0 then message, 'No data between kx_lims' else begin
      power_use = power_use[wh, *, *]
      noise_expval_use = noise_expval_use[wh, *, *]
      weights_use = weights_use[wh, *, *]
      kx_mpc_use = kx_mpc[wh]
      n_kx = count
      if n_elements(kperp_density_measure) gt 0 and n_elements(kperp_density_cutoff) gt 0 $
        then kperp_density_measure_use = kperp_density_measure[wh, *]
    endelse
  endif else begin
    kx_mpc_use = kx_mpc

    if n_elements(kperp_density_measure) gt 0 and n_elements(kperp_density_cutoff) gt 0 $
      then kperp_density_measure_use = kperp_density_measure
  endelse

  if n_elements(ky_lims) gt 0 then begin
    if n_elements(ky_lims) gt 2 then message, 'ky_lims has too many elements' $
    else if n_elements(ky_lims) eq 1 then ky_lims = [-1, 1] * abs(ky_lims)

    wh = where(ky_mpc ge min(ky_lims) and ky_mpc le max(ky_lims), count)
    if count eq 0 then message, 'No data between ky_lims' else begin
      power_use = power_use[*, wh , *]
      noise_expval_use = noise_expval_use[*, wh , *]
      weights_use = weights_use[*, wh, *]
      ky_mpc_use = ky_mpc[wh]
      n_ky = count
      if n_elements(kperp_density_measure_use) gt 0 and n_elements(kperp_density_cutoff) gt 0 $
        then kperp_density_measure_use = kperp_density_measure_use[*, wh]
    endelse
  endif else ky_mpc_use = ky_mpc

  ;; apply density corrections if desired
  if n_elements(kperp_density_measure) gt 0 and n_elements(kperp_density_cutoff) gt 0 then begin
    wh_kperp_dense = where(kperp_density_measure_use gt kperp_density_cutoff, count_kperp_dense, $
      complement = wh_kperp_nd, ncomplement=count_kperp_nd)

    if count_kperp_dense gt 0 then begin
      temp = reform(power_use, n_kx*n_ky, n_kz)
      temp[wh_kperp_dense, *] /= 0.5
      power_use = reform(temp, n_kx, n_ky, n_kz)

      temp = reform(weights_use, n_kx*n_ky, n_kz)
      temp[wh_kperp_dense, *] *= 0.5^2.
      weights_use = reform(temp, n_kx, n_ky, n_kz)

      temp = reform(noise_expval_use, n_kx*n_ky, n_kz)
      temp[wh_kperp_dense, *] /= 0.5
      noise_expval_use = reform(temp, n_kx, n_ky, n_kz)
    endif else print, 'no kperp values exceed kperp_density_cutoff'

  endif

  if not binning_2d_options.log_kperp then begin

    temp = sqrt(rebin(kx_mpc_use, n_kx, n_ky)^2 + rebin(reform(ky_mpc_use, 1, n_ky), n_kx, n_ky)^2)
    if keyword_set(edge_on_grid) then begin
      kperp_min = floor(min(temp) / binning_2d_options.kperp_bin) * binning_2d_options.kperp_bin
    endif else begin
      kperp_min = floor((min(temp) + binning_2d_options.kperp_bin/2d) / binning_2d_options.kperp_bin) * binning_2d_options.kperp_bin - binning_2d_options.kperp_bin/2d
    endelse

    kperp_array = temp

    ;; Use histogram with reverse indicies to bin in kperp
    kperp_hist = histogram(kperp_array, binsize = binning_2d_options.kperp_bin, min = kperp_min, $
      omax = kperp_max, locations = lin_kperp_locs, reverse_indices = kperp_ri)

    kperp_centers_mpc = lin_kperp_locs + binning_2d_options.kperp_bin/2d
    kperp_edges_mpc = [lin_kperp_locs, max(lin_kperp_locs) + binning_2d_options.kperp_bin]

  endif else begin
    ;; bin in logarithmic k_perp space
    kperp_log_binsize = binning_2d_options.kperp_bin
    temp = sqrt(rebin(kx_mpc_use, n_kx, n_ky)^2 + rebin(reform(ky_mpc_use, 1, n_ky), n_kx, n_ky)^2)

    ;; get min value, adjust so centers of bins will be on major grid.
    if min(temp) gt 0 then begin
      if keyword_set(edge_on_grid) then begin
        kperp_min = floor((alog10(min(temp))) / kperp_log_binsize) * kperp_log_binsize
      endif else begin
        kperp_min = floor((alog10(min(temp)) + kperp_log_binsize/2d) / kperp_log_binsize) * $
          kperp_log_binsize - kperp_log_binsize/2d
      endelse
      kperp_array = alog10(temp)
    endif else begin
      ;; Add an extra bin to put the k=0 mode into.
      if max(temp) eq 0 then message, 'kperps are all zero'
      if keyword_set(edge_on_grid) then begin
        kperp_min = floor((alog10(min(temp[where(temp gt 0)]))) / kperp_log_binsize - 1) * kperp_log_binsize
      endif else begin
        kperp_min = alog10(min(temp[where(temp gt 0)]))-3d*kperp_log_binsize/2d
      endelse
      if keyword_set(match_datta) then begin
        if keyword_set(edge_on_grid) then kperp_min = -3d $
        else kperp_min = -3d - kperp_log_binsize/2d
      endif

      ;; calculate log kperp at each location in 1/Mpc. Set 0 mode to bottom
      ;; value so it doesn't get left out.
      kperp_array = alog10(temp)
      kperp_array[where(temp le 0)] = kperp_min
    endelse

    ;; Use histogram with reverse indicies to bin in log kperp (want even spacing in log space.)
    kperp_hist = histogram(kperp_array, binsize = kperp_log_binsize, min = kperp_min, $
      omax = kperp_max, locations = log_kperp_locs, reverse_indices = kperp_ri)

    kperp_centers_mpc = 10^(log_kperp_locs + kperp_log_binsize/2d)
    kperp_edges_mpc = [10^(log_kperp_locs), 10^(max(log_kperp_locs) + kperp_log_binsize)]
  endelse

  n_kperp = n_elements(kperp_hist)
  weighted_power_mid = make_array(n_kperp, n_kz, type=input_type)
  weighted_nev_mid = make_array(n_kperp, n_kz, type=input_type)
  weights_mid = dblarr(n_kperp, n_kz)
  nvox_mid = lonarr(n_kperp, n_kz)
  var_power_mid = make_array(n_kperp, n_kz, type=input_type)
  mean_var_mid = make_array(n_kperp, n_kz, type=input_type)

  reformed_power = reform(temporary(power_use), n_kx*n_ky, n_kz)
  reformed_nev = reform(temporary(noise_expval_use), n_kx*n_ky, n_kz)
  reformed_weights = reform(temporary(weights_use), n_kx*n_ky, n_kz)
  reformed_wtpower = reformed_power * reformed_weights
  reformed_wtnev = reformed_nev * reformed_weights
  for i=0, n_kperp-1 do begin
    if kperp_hist[i] gt 0 then begin
      weighted_power_mid[i, *] = total(reform(reformed_wtpower[kperp_ri[kperp_ri[i] : kperp_ri[i+1]-1], *], kperp_hist[i], n_kz), 1)
      weighted_nev_mid[i, *] = total(reform(reformed_wtnev[kperp_ri[kperp_ri[i] : kperp_ri[i+1]-1], *], kperp_hist[i], n_kz), 1)

      temp = reformed_weights[kperp_ri[kperp_ri[i] : kperp_ri[i+1]-1], *]
      weights_mid[i, *] = total(reform(temp, kperp_hist[i], n_kz), 1)

      temp_sigma = 1./temp
      wh_temp0 = where(temp eq 0, count_temp0)
      if count_temp0 gt 0 then temp_sigma[wh_temp0] = 0
      temp_power = reformed_power[kperp_ri[kperp_ri[i] : kperp_ri[i+1]-1], *]

      if min(temp) le 0 then begin
        for j=0, n_kz -1 do begin
          wh_n0 = where(temp[*,j] gt 0, count_n0)
          nvox_mid[i, j] = count_n0

          if count_n0 gt 0 then begin
            ;; calculate variance for an (assumed) zero mean signal
            ;; IDL variance function is unbiased sample variance for the case
            ;; where the mean is unknown
            var_power_mid[i, j] = total(temp_power[wh_n0, j]^2.)/count_n0
            ;if count_n0 gt 1 then var_power_mid[i,j] = variance(temp_power[wh_n0, j])
            mean_var_mid[i, j] = mean(temp_sigma[wh_n0, j])
          endif
        endfor
      endif else begin
        nvox_mid[i, *] = kperp_hist[i]

        for j=0, n_kz -1 do begin
          ;; calculate variance for an (assumed) zero mean signal
          ;; IDL variance function is unbiased sample variance for the case where
          ;; the mean is unknown
          var_power_mid[i, j] = total(temp_power[*, j]^2.)/kperp_hist[i]
          ;if kperp_hist[i] gt 1 then var_power_mid[i, j] = variance(temp_power[*, j])
          mean_var_mid[i, j] = mean(temp_sigma[*, j])
        endfor
      endelse
    endif
  endfor

  undefine, reformed_power, reformed_nev, reformed_wtpower, reformed_wtnev, reformed_weights

  if not binning_2d_options.log_kpar then begin

    if keyword_set(edge_on_grid) then begin
      kpar_min = floor(min(kz_mpc) / binning_2d_options.kpar_bin) * binning_2d_options.kpar_bin
    endif else begin
      kpar_min = floor((min(kz_mpc) + binning_2d_options.kpar_bin/2d) / binning_2d_options.kpar_bin) * binning_2d_options.kpar_bin - binning_2d_options.kpar_bin/2d
    endelse

    ;; Use histogram with reverse indicies to bin in kpar
    kpar_hist = histogram(kz_mpc, binsize = binning_2d_options.kpar_bin, min = kpar_min, $
      omax = kpar_max, locations = lin_kpar_locs, reverse_indices = kpar_ri)

    kpar_centers_mpc = lin_kpar_locs + binning_2d_options.kpar_bin/2d
    kpar_edges_mpc = [lin_kpar_locs, max(lin_kpar_locs) + binning_2d_options.kpar_bin]

  endif else begin
    ;; now make kz bins logrithmic by resampling
    kpar_log_binsize = binning_2d_options.kpar_bin

    ;; get min value, adjust so centers of bins will be on major grid
    if min(kz_mpc) gt 0 then begin
      if keyword_set(edge_on_grid) then begin
        kpar_min = floor((alog10(min(kz_mpc))) / kpar_log_binsize) * kpar_log_binsize
      endif else begin
        kpar_min = floor((alog10(min(kz_mpc)) + kpar_log_binsize/2d) / kpar_log_binsize) * $
          kpar_log_binsize - kpar_log_binsize/2d
      endelse

      kpar_arr = alog10(kz_mpc)
    endif else begin
      if max(kz_mpc) eq 0 then message, 'kz_mpc is all zero'
      if keyword_set(edge_on_grid) then begin
        kpar_min = floor((alog10(min(kz_mpc[where(kz_mpc gt 0)]))) / kpar_log_binsize - 1) * kpar_log_binsize
      endif else begin
        kpar_min = alog10(min(kz_mpc[where(kz_mpc gt 0)]))-3d*kpar_log_binsize/2d
      endelse
      kpar_arr = alog10(kz_mpc)
      kpar_arr[where(kz_mpc le 0)] = kpar_min
    endelse

    kpar_hist = histogram(kpar_arr, binsize = kpar_log_binsize, min = kpar_min, $
      omax = kpar_max, locations = log_kpar_locs, reverse_indices = kpar_ri)
    kpar_centers_mpc = 10^(log_kpar_locs + kpar_log_binsize/2d)
    kpar_edges_mpc = [10^(log_kpar_locs), 10^(max(log_kpar_locs) + kpar_log_binsize)]
    kz_mpc_edges = [kz_mpc, max(kz_mpc) + kz_mpc[1] - kz_mpc[0]]
  endelse

  n_kpar = n_elements(kpar_hist)
  weighted_power_2d = make_array(n_kperp, n_kpar, type=input_type)
  weighted_nev_2d = make_array(n_kperp, n_kpar, type=input_type)
  weights_2d = dblarr(n_kperp, n_kpar)
  nvox_2d = lonarr(n_kperp, n_kpar)
  var_power_2d = make_array(n_kperp, n_kpar, type=input_type)
  mean_var_2d = make_array(n_kperp, n_kpar, type=input_type)

  for j=0, n_kpar-1 do begin
    if kpar_hist[j] ne 0 then begin
      weighted_power_2d[*,j] = total(reform(weighted_power_mid[*, kpar_ri[kpar_ri[j] : kpar_ri[j+1]-1]], $
        n_kperp, kpar_hist[j]),2)
      weighted_nev_2d[*,j] = total(reform(weighted_nev_mid[*, kpar_ri[kpar_ri[j] : kpar_ri[j+1]-1]], $
        n_kperp, kpar_hist[j]),2)

      temp = reform(weights_mid[*, kpar_ri[kpar_ri[j] : kpar_ri[j+1]-1]], n_kperp, kpar_hist[j])
      weights_2d[*, j] = total(temp, 2)

      ;; var_power_2d is a weighted average of var_power_mid with nvox_mid as weights
      var_power_2d[*, j] = total(reform(var_power_mid[*, kpar_ri[kpar_ri[j] : kpar_ri[j+1]-1]] * $
        nvox_mid[*, kpar_ri[kpar_ri[j] : kpar_ri[j+1]-1]], n_kperp, kpar_hist[j]),2) $
        / total(reform(nvox_mid[*, kpar_ri[kpar_ri[j] : kpar_ri[j+1]-1]], n_kperp, kpar_hist[j]),2)

      ;; mean_var_2d is a weighted average of mean_var_mid with nvox_mid as weights
      mean_var_2d[*,j] = total(reform(mean_var_mid[*, kpar_ri[kpar_ri[j] : kpar_ri[j+1]-1]] * $
        nvox_mid[*, kpar_ri[kpar_ri[j] : kpar_ri[j+1]-1]], n_kperp, kpar_hist[j]),2) $
        / total(reform(nvox_mid[*, kpar_ri[kpar_ri[j] : kpar_ri[j+1]-1]], n_kperp, kpar_hist[j]),2)

      if min(temp) le 0 then begin
        wh_n0 = where(temp gt 0, count_n0)
        mask = intarr(n_kperp, kpar_hist[j])
        if count_n0 gt 0 then mask[wh_n0] = 1

        nvox_2d[*, j] = total(reform(nvox_mid[*, kpar_ri[kpar_ri[j] : kpar_ri[j+1]-1]] * mask, $
          n_kperp, kpar_hist[j]),2)
      endif else begin
        nvox_2d[*, j] = total(reform(nvox_mid[*, kpar_ri[kpar_ri[j] : kpar_ri[j+1]-1]], $
          n_kperp, kpar_hist[j]),2)
      endelse
    endif
  endfor

  undefine, weighted_power_mid, weighted_nev_mid, weights_mid, nvox_mid, var_power_mid, mean_var_mid

  power_ave = weighted_power_2d/weights_2d
  noise_expval_2d = weighted_nev_2d/weights_2d
  wh = where(weights_2d eq 0, count)
  if count ne 0 then begin
    power_ave[wh] = 0d
    noise_expval_2d[wh] = 0d
  endif

  return, power_ave

end
