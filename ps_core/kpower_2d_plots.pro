pro kpower_2d_plots, power_savefile, multi_pos = multi_pos, start_multi_params = start_multi_params, window_num = window_num, $
    plot_weights = plot_weights, plot_noise = plot_noise, plot_sigma = plot_sigma, plot_exp_noise = plot_exp_noise, $
    ratio = ratio, diff = diff, snr = snr, nnr = nnr, $
    kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, data_range = data_range, $
    color_profile = color_profile, log_cut_val = log_cut_val, plotfile = plotfile, png = png, eps = eps, $
    no_title = no_title, title = title, title_prefix = title_prefix, norm_2d = norm_2d, $
    norm_factor = norm_factor, grey_scale = grey_scale, wedge_amp = wedge_amp, plot_wedge_line = plot_wedge_line, $
    baseline_axis = baseline_axis, delay_axis = delay_axis, kperp_linear_axis = kperp_linear_axis, $
    kpar_linear_axis = kpar_linear_axis, no_units = no_units, hinv = hinv, charsize = charsize_in, $
    cb_size = cb_size_in, margin=margin_in, cb_margin = cb_margin_in
      
  if n_elements(plotfile) gt 0 or keyword_set(png) or keyword_set(eps) then pub = 1 else pub = 0
  if pub eq 1 then begin
    if not (keyword_set(png) or keyword_set(eps)) then begin
      basename = cgRootName(plotfile, directory=directory, extension=extension)
      
      case extension of
        'eps': eps=1
        'png': png=1
        '': png = 1
        else: begin
          print, 'Unrecognized extension, using png'
          png = 1
        end
      endcase
      
    endif
    if n_elements(plotfile) eq 0 and n_elements(multi_pos) eq 0 then begin
      if keyword_set(eps) then plotfile = 'idl_kpower_2d_plots.eps' else plotfile = 'idl_kpower_2d_plots'
      cd, current = current_dir
      print, 'no filename specified for kpower_2d_plots output. Using ' + current_dir + path_sep() + plotfile
    endif
    
    if keyword_set(png) and keyword_set(eps) then begin
      print, 'both eps and png cannot be set, using png'
      eps = 0
    endif
    
    case 1 of
      png: begin
        plot_exten = '.png'
        delete_ps = 1
      end
      eps: begin
        plot_exten = '.eps'
        delete_ps = 0
      end
    endcase
  endif
  
  
  if n_elements(window_num) eq 0 then window_num = 1
  
  if n_elements(start_multi_params) gt 0 and n_elements(multi_pos) gt 0 then message, 'If start_multi_params are passed, ' + $
    'multi_pos cannot be passed because then it is used as an output to pass back the positions for future plots.'
    
  if n_elements(multi_pos) gt 0 then begin
    if n_elements(multi_pos) ne 4 then message, 'multi_pos must be a 4 element plot position vector'
    if max(multi_pos) gt 1 or min(multi_pos) lt 0 then message, 'multi_pos must be in normalized coordinates (between 0 & 1)'
    if multi_pos[2] le multi_pos[0] or multi_pos[3] le multi_pos[1] then $
      message, 'In multi_pos, x1 must be greater than x0 and y1 must be greater than y0 '
  endif
  
  color_profile_enum = ['log_cut', 'sym_log', 'abs']
  if n_elements(color_profile) eq 0 then color_profile = 'log_cut'
  
  wh_prof = where(color_profile_enum eq color_profile, count)
  if count eq 0 then message, 'Color profile must be one of: ' + strjoin(color_profile_enum, ', ')
  
  if total([keyword_set(plot_weights), keyword_set(plot_sigma), keyword_set(plot_exp_noise), $
    keyword_set(plot_noise), keyword_set(snr), keyword_set(nnr)]) gt 1 then $
    message, 'only one of [plot_noise, plot_sigma, plot_exp_noise, plot_weights, snr, nnr] keywords can be set'
    
    
  if keyword_set(ratio) or keyword_set(diff) then begin
    if keyword_set(ratio) and keyword_set(diff) then message, 'Ratio and diff keywords cannot be set simultaneously'
    if n_elements(power_savefile) gt 2 then message, 'Only 2 files can be specified in power_savefile if ratio or diff keyword is set'
    restore, power_savefile[0]
    power1 = power
    weights1 = weights
    weights1 = weights
    noise_expval1 = noise_expval
    if n_elements(noise) ne 0 then noise1 = noise
    kperp1 = kperp_edges
    kpar1 = kpar_edges
    kpar_bin1 = kpar_bin
    kperp_bin1 = kperp_bin
    
    restore, power_savefile[1]
    power2 = power
    if n_elements(noise) ne 0 then noise2 = noise
    weights2 = weights
    noise_expval2 = noise_expval
    
    if total(abs(kperp1 - kperp_edges)) ne 0 or total(abs(kpar1 - kpar_edges)) ne 0 or $
      total(abs(size(power1,/dimension) - size(power2, /dimension))) ne 0 then $
      message, 'dimensions and kperp/kpar edges must be the same in both files'
      
    if keyword_set(snr) then message, 'snr keyword cannot be used with ratio keyword'
    
    if (keyword_set(plot_noise) or keyword_set(nnr)) and (n_elements(noise1) eq 0 or n_elements(noise2) eq 0) then $
      message, 'Noise is not included in one or both files'
      
    if keyword_set(diff) then begin
      power = power1 - power2
      plot_type = 'diff'
    endif
    
    if keyword_set(plot_weights) then begin
      power = weights1 / weights2
      plot_type = 'weight_ratio'
    endif
    
    if keyword_set(plot_sigma) then begin
      power = sqrt(weights2 / weights1)
      plot_type = 'sigma_ratio'
    endif
    
    if keyword_set(plot_noise_exp) then begin
      power = noise_expval1 / noise_expval2
      plot_type = 'exp_noise_ratio'
    endif
    
    if keyword_set(plot_noise) then begin
      power = noise1 / noise2
      plot_type = 'noise_ratio'
    endif
    
    if n_elements(plot_type) eq 0 then begin
      power = power1 / power2
      plot_type = 'ratio'
    endif
    
  endif else begin
    if n_elements(power_savefile) gt 1 then message, $
      'Only 1 file can be specified in power_savefile unless ratio or diff keyword is set'
      
    restore, power_savefile
    
    if keyword_set(snr) then begin
      power = power / noise_expval
      wh_err0 = where(noise_expval eq 0, count_err0)
      if count_err0 gt 0 then power[wh_err0] = 0
      
      plot_type = 'snr'
    endif
    
    if keyword_set(plot_weights) then begin
      power = weights
      plot_type = 'weight'
    endif
    
    if keyword_set(plot_sigma) then begin
      power = 1/sqrt(weights)
      wh_wt0 = where(weights eq 0, count_wt0)
      if count_wt0 gt 0 then power[wh_wt0 ] = 0
      plot_type = 'sigma'
    endif
    
    if keyword_set(plot_exp_noise) then begin
      power = noise_expval
      plot_type = 'exp_noise'
    endif
    
    if keyword_set(plot_noise) then begin
      if n_elements(noise) eq 0 then message, 'noise is undefined in this file'
      power = noise
      plot_type = 'noise'
    endif
    
    if keyword_set(nnr) then begin
      power = noise / noise_expval
      wh_err0 = where(noise_expval eq 0, count_err0)
      if count_err0 gt 0 then power[wh_err0] = 0
      
      plot_type = 'nnr'
    endif
    
    if n_elements(plot_type) eq 0 then plot_type = 'power'
    
  endelse
  
  
  dims = size(power, /dimension)
  n_kperp = dims[0]
  n_kpar = dims[1]
  
  if keyword_set(hinv) then begin
    kperp_edges = kperp_edges / hubble_param
    kpar_edges = kpar_edges / hubble_param
    if plot_type eq 'power' or plot_type eq 'noise' or plot_type eq 'sigma' or plot_type eq 'exp_noise' then $
      power = power * (hubble_param)^3d
  endif
  
  if n_elements(kperp_plot_range) eq 0 then kperp_plot_range = minmax(kperp_edges)
  if n_elements(kpar_plot_range) eq 0 then kpar_plot_range = minmax(kpar_edges)
  
  units_str = ''
  case plot_type of
    'power': begin
      if keyword_set(hinv) then units_str = textoidl(' (mK^2 h^{-3} Mpc^3)', font = font) $
      else units_str = textoidl(' (mK^2 Mpc^3)', font = font)
      
      plot_title = textoidl('P_k', font = font)
      if pub then plotfile_add = '_2dkpower' + plot_exten
    end
    'ratio': begin
      units_str = ''
      plot_title = 'Power Ratio'
      if pub then plotfile_add = '_ratio' + plot_exten
    end
    'diff': begin
      if keyword_set(hinv) then units_str = textoidl(' (mK^2 h^{-3} Mpc^3)', font = font) $
      else units_str = textoidl(' (mK^2 Mpc^3)', font = font)
      
      plot_title = textoidl('P_k difference', font = font)
      if pub then plotfile_add = '_diff' + plot_exten
    end
    'weight': begin
      units_str = ''
      plot_title = 'Weights'
      if pub then plotfile_add = '_2dweight' + plot_exten
    end
    'weight_ratio': begin
      units_str = ''
      plot_title = 'Weights Ratio'
      if pub then plotfile_add = '_weight_ratio' + plot_exten
    end
    'sigma': begin
      if keyword_set(hinv) then units_str = textoidl(' (mK^2 h^{-3} Mpc^3)', font = font) $
      else units_str = textoidl(' (mK^2 Mpc^3)', font = font)
      
      plot_title = 'Error (sigma)'
      if pub then plotfile_add = '_2dsigma' + plot_exten
    end
    'sigma_ratio': begin
      units_str = ''
      
      if keyword_set(pub) then sigma_char = textoidl('\sigma', font = 1) else sigma_char = textoidl('\sigma', font = -1)
      plot_title = 'Error (sigma) Ratio'
      if pub then plotfile_add = '_sigma_ratio' + plot_exten
    end
    'exp_noise': begin
      if keyword_set(hinv) then units_str = textoidl(' (mK^2 h^{-3} Mpc^3)', font = font) $
      else units_str = textoidl(' (mK^2 Mpc^3)', font = font)
      
      plot_title = 'Expected Noise'
      if pub then plotfile_add = '_2dnoise_expval' + plot_exten
    end
    'exp_noise_ratio': begin
      units_str = ''
      plot_title = 'Expected Noise ratio'
      if pub then plotfile_add = '_noise_expval_ratio' + plot_exten
    end
    'snr': begin
      units_str = ''
      plot_title = 'SNR (' + textoidl('P_k/N_E', font = font) + ')'
      if pub then plotfile_add = '_2dsnr' + plot_exten
    end
    'noise': begin
      if keyword_set(hinv) then units_str = textoidl(' (mK^2 h^{-3} Mpc^3)', font = font) $
      else units_str = textoidl(' (mK^2 Mpc^3)', font = font)
      
      plot_title = 'Observed Noise'
      if pub then plotfile_add = '_2dnoise' + plot_exten
    end
    'nnr': begin
      units_str = ''
      plot_title = 'Noise Ratio (' + textoidl('N_O/N_E', font = font) + ')'
      if pub then plotfile_add = '_2dnnr' + plot_exten
    end
    'noise_ratio': begin
      units_str = ''
      plot_title = 'Observed Noise Ratio'
      if pub then plotfile_add = '_noise_ratio' + plot_exten
    end
  endcase
  
  if pub then begin
    if n_elements(plotfile) eq 0 then plotfile = strsplit(power_savefile[0], '.idlsave', /regex, /extract) + plotfile_add $
    else if strcmp(strmid(plotfile, strlen(plotfile)-4), plot_exten, /fold_case) eq 0 then plotfile = plotfile + plot_exten
  endif
  
  if keyword_set(no_units) then units_str = ''
  
  wh_kperp_inrange = where(kperp_edges ge kperp_plot_range[0] and kperp_edges[1:*] le kperp_plot_range[1], n_kperp_plot)
  wh_kpar_inrange = where(kpar_edges ge kpar_plot_range[0] and kpar_edges[1:*] le kpar_plot_range[1], n_kpar_plot)
  
  if n_kperp_plot eq 0 or n_kpar_plot eq 0 then message, 'No data in plot k range'
  
  if n_kperp_plot ne n_kperp then begin
    power = power[wh_kperp_inrange, *]
    temp = [wh_kperp_inrange, wh_kperp_inrange[n_kperp_plot-1]+1]
    kperp_edges =kperp_edges[temp]
    n_kperp = n_kperp_plot
  endif
  if n_kpar_plot ne n_kpar then begin
    power = power[*, wh_kpar_inrange]
    temp = [wh_kpar_inrange, wh_kpar_inrange[n_kpar_plot-1]+1]
    kpar_edges = kpar_edges[temp]
    n_kpar = n_kpar_plot
  endif
  
  ;; Check whether binning is log or not
  log_bins = [1, 1]
  kperp_log_diffs = (alog10(kperp_edges) - shift(alog10(kperp_edges), 1))[2:*]
  if total(abs(kperp_log_diffs - kperp_log_diffs[0])) gt n_kperp*1e-15 then log_bins[0] = 0
  kpar_log_diffs = (alog10(kpar_edges) - shift(alog10(kpar_edges), 1))[2:*]
  if total(abs(kpar_log_diffs - kpar_log_diffs[0])) gt n_kpar*1e-15 then log_bins[1] = 0
  
  log_bins = [0, 0]
  kperp_diffs = (kperp_edges - shift(kperp_edges, 1))[1:*]
  if total(abs(kperp_diffs - kperp_diffs[0])) gt n_kperp*1e-7 then log_bins[0] = 1
  kpar_diffs = (kpar_edges - shift(kpar_edges, 1))[1:*]
  if total(abs(kpar_diffs - kpar_diffs[0])) gt n_kpar*1e-7 then log_bins[1] = 1
  
  
  
  tvlct, r, g, b, /get
  
  if keyword_set(grey_scale) then begin
    cgloadct, 0, /reverse
    color_range = [0, 255]
    background_color = 'white'
    annotate_color = 'black'
  endif else begin
    cgloadct, 25, /brewer, /reverse
    color_range = [0, 255]
    background_color = 'white'
    annotate_color = 'black'
  endelse
  n_colors = color_range[1] - color_range[0]
  
  
  if keyword_set(norm_2d) then begin
    if n_elements(norm_factor) eq 0 then norm_factor = 1/max(power)
    power = power * norm_factor
  endif
  
  if n_elements(data_range) eq 0 then begin
    no_input_data_range = 1
    data_range = minmax(power)
  endif else begin
    if n_elements(data_range) ne 2 then message, 'data_range must be a 2 element vector'
    no_input_data_range = 0
  endelse
  
  if data_range[1] lt data_range[0] then message, 'data_range[0] must be less than data_range[1]'
  wh = where(power gt 0d, count)
  if count gt 0 then min_pos = min(power[wh]) else if data_range[0] gt 0 then min_pos = data_range[0] else $
    if data_range[1] gt 0 then min_pos = data_range[1]/10d else min_pos = 0.01d
    
  ;; check whether we need to have varying bin sizes (log/linear options don't match)
  log_axes = [1,1]
  if keyword_set(kperp_linear_axis) then log_axes[0] = 0
  if keyword_set(kpar_linear_axis) then log_axes[1] = 0
  
  if total(abs(log_bins-log_axes)) ne 0 then begin
    ;; need to make a new image array with varying bin sizes
    if log_bins[0] ne log_axes[0] then begin
      if log_bins[0] eq 0 then begin
        ;; linear binning, log axes
        wh_kperp0 = where(kperp_edges lt 0, count_kperp0, complement = wh_kperp_good)
        if count_kperp0 gt 1 then stop
        
        kperp_log_edges = alog10(kperp_edges)
        if count_kperp0 eq 1 then begin
          kperp_log_diffs = (kperp_log_edges[1:*] - shift(kperp_log_edges[1:*], 1))[1:*]
          kperp_log_diffs = [kperp_log_diffs[0], kperp_log_diffs]
          kperp_log_edges[wh_kperp0] = kperp_log_edges[wh_kperp0+1] - kperp_log_diffs[wh_kperp0]
        endif
        
        image_kperp_delta = min(kperp_log_diffs)
        kperp_bin_widths = round(kperp_log_diffs / image_kperp_delta)
      endif else begin
        ;; log binning, linear axes
        kperp_diffs = (kperp_edges[1:*] - shift(kperp_edges[1:*], 1))[1:*]
        image_kperp_delta = min(kperp_diffs)/2d
        kperp_bin_widths = round(kperp_diffs / image_kperp_delta)
      endelse
      rebin_x = 1
    endif else begin
      ;; axes and binning agree
      if log_axes[0] eq 1 then kperp_log_edges = alog10(kperp_edges)
      rebin_x = 0
    endelse
    
    if log_bins[1] ne log_axes[1] then begin
      if log_bins[1] eq 0 then begin
        ;; linear binning, log axes
        wh_kpar0 = where(kpar_edges lt 0, count_kpar0, complement = wh_kpar_good)
        if count_kpar0 gt 1 then stop
        
        kpar_log_edges = alog10(kpar_edges)
        if count_kpar0 eq 1 then begin
          kpar_log_diffs = (kpar_log_edges[1:*] - shift(kpar_log_edges[1:*], 1))[1:*]
          kpar_log_diffs = [kpar_log_diffs[0], kpar_log_diffs]
          kpar_log_edges[wh_kpar0] = kpar_log_edges[wh_kpar0+1] - kpar_log_diffs[wh_kpar0]
        endif
        
        image_kpar_delta = min(kpar_log_diffs)/2d
        kpar_bin_widths = round(kpar_log_diffs / image_kpar_delta)
      endif else begin
        ;; log binning, linear axes
        kpar_diffs = (kpar_edges[1:*] - shift(kpar_edges[1:*], 1))[1:*]
        image_kpar_delta = min(kpar_diffs)/2d
        kpar_bin_widths = round(kpar_diffs / image_kpar_delta)
      endelse
      rebin_y = 1
    endif else begin
      if log_bins[1] eq 1 then kpar_log_edges = alog10(kpar_edges)
      rebin_y = 0
    endelse
    
    if rebin_x eq 1 then begin
      ;; now get width for each input bin in image array
      nkperp_image = total(kperp_bin_widths)
      
      h_kperp = histogram(total(kperp_bin_widths,/cumulative)-1,binsize=1, min=0, reverse_indices=ri_kperp)
      undefine, h_kperp
      kperp_inds = ri_kperp[0:nkperp_image-1]-ri_kperp[0]
    endif else begin
      nkperp_image = n_kperp
      kperp_inds = indgen(nkperp_image)
    endelse
    
    if rebin_y eq 1 then begin
      nkpar_image = total(kpar_bin_widths)
      
      h_kpar = histogram(total(kpar_bin_widths,/cumulative)-1,binsize=1, min=0, reverse_indices=ri_kpar)
      undefine, h_kpar
      kpar_inds = rebin(reform(ri_kpar[0:nkpar_image-1]-ri_kpar[0], 1, nkpar_image), nkperp_image, nkpar_image)
    endif else begin
      nkpar_image = n_kpar
      kpar_inds = rebin(reform(indgen(nkpar_image), 1, nkpar_image), nkperp_image, nkpar_image)
    endelse
    
    kperp_inds = rebin(kperp_inds, nkperp_image, nkpar_image)
    power_plot = power[kperp_inds, kpar_inds]
    
    ;; now expand array in any non-rebinned direction to prevent interpolation
    if rebin_x eq 0 then power_plot = congrid(power_plot, nkperp_image*10, nkpar_image)
    if rebin_y eq 0 then power_plot = congrid(power_plot, nkperp_image, nkpar_image*10)
    
  endif else begin
    ;; axes & binning agree for both directions
    ;; expand image array to prevent interpolation in postscript
    power_plot = congrid(power, n_kperp*10, n_kpar*10)
    if log_axes[0] eq 1 then kperp_log_edges = alog10(kperp_edges)
    if log_axes[1] eq 1 then kpar_log_edges = alog10(kpar_edges)
  endelse
  
  
  case color_profile of
    'log_cut': begin
    
      if data_range[1] lt 0 then message, 'log_cut color profile will not work for entirely negative arrays.'
      
      if n_elements(log_cut_val) eq 0 then begin
        if data_range[0] gt 0 then log_cut_val = alog10(data_range[0]) else $
          log_cut_val = alog10(min_pos)
      endif
      
      log_data_range = [log_cut_val, alog10(data_range[1])]
      
      ;; if log_cut_val is less than min_pos then indicate it in the color bar
      if min(power_plot) le 0 then begin
      
      
        wh_zero = where(power_plot eq 0, count_zero)
        if count_zero gt 0 then begin
          min_pos_color = 2
          zero_color=1
          zero_val = n_colors/(log_data_range[1]-log_data_range[0])
          
          wh_zero = where(power_plot eq 0, count)
          if count ne 0 then power_log[wh_zero] = zero_val
          
          if log_cut_val gt min_pos then log_data_range = [log_cut_val-2*n_colors/(log_data_range[1]-log_data_range[0]), alog10(data_range[1])]
          
        endif else begin
          min_pos_color = 1
          if log_cut_val gt min_pos then log_data_range = [log_cut_val-n_colors/(log_data_range[1]-log_data_range[0]), alog10(data_range[1])]
          
        endelse
        neg_color=0
        neg_val = log_data_range[0]
        
        power_log = alog10(power_plot)
        wh_under = where(power_plot lt 10^double(log_cut_val), count)
        if count ne 0 then power_log[wh_under] = log_data_range[0]
        wh_neg = where(power_plot lt 0, count)
        if count ne 0 then power_log[wh_neg] = neg_val
        wh_over = where(power_log gt log_data_range[1], count)
        if count ne 0 then power_log[wh_over] = log_data_range[1]
        
        oob_low = 0
        
        cgloadct, 25, /brewer, /reverse, BOTTOM = min_pos_color, NCOLORS = 256-min_pos_color
        if count_zero gt 0 then begin
          cgLoadCT, 0, CLIP=[50, 50], BOTTOM=0, NCOLORS=1
          cgLoadCT, 0, CLIP=[250, 250], BOTTOM=zero_color, NCOLORS=1
        endif else cgLoadCT, 0, CLIP=[50, 50], BOTTOM=0, NCOLORS=1
        
      endif else begin
      
        min_pos_color = 0
        
        power_log = alog10(power_plot)
        wh_under = where(power_plot lt 10^double(log_cut_val), count)
        if count ne 0 then power_log[wh_under] = log_data_range[0]
        wh_over = where(power_log gt log_data_range[1], count)
        if count ne 0 then power_log[wh_over] = log_data_range[1]
        
      endelse
      
      if no_input_data_range eq 1 then data_range = 10^log_data_range
    end
    'sym_log': begin
    
      ;; want log-like behavior on each side of 0: log(pos. vals +1) & log((-1)*neg. vals+1)
      neg_inds = where(power_plot lt 0, n_neg, complement = pos_inds, ncomplement = n_pos)
      
      power_log = alog10(power_plot+1)
      if n_neg gt 0 then power_log[neg_inds] = (-1) * alog10((-1)*power_plot[neg_inds] + 1)
      
      log_data_range = alog10(data_range + 1)
      if data_range[0] lt 0 then begin
        temp = alog10((-1)*data_range[0]+1)
        if temp gt log_data_range[1] then log_data_range = [(-1)*temp, temp] $
        else log_data_range[0] = (-1)*log_data_range[1]
      endif
    end
    'abs': begin
    
      abs_power_plot = abs(power_plot)
      log_data_range = dblarr(2)
      if data_range[0] lt 0 then log_data_range[0] = alog10(min(abs_power_plot[where(abs_power_plot gt 0)])) $
      else log_data_range[0] = alog10(data_range[0])
      log_data_range[1] = alog10(max(abs(data_range)))
      
      power_log = alog10(abs_power_plot)
      wh_zero = where(power_plot eq 0, count)
      if count ne 0 then power_log[wh_zero] = log_data_range[0]
      
      abs_power_plot = 0
      
    end
  endcase
  
  power_log_norm = (power_log-log_data_range[0])*n_colors/(log_data_range[1]-log_data_range[0]) + color_range[0]
  
  max_ysize = 1000
  max_xsize = 1200
  base_size = 600
  if n_elements(multi_pos) eq 4 then begin
    ;; work out positions scaled to the area allowed in multi_pos with proper aspect ratio
    multi_xlen = (multi_pos[2]-multi_pos[0])
    multi_ylen = (multi_pos[3]-multi_pos[1])
    multi_center = [multi_pos[0] + multi_xlen/2d, multi_pos[1] + multi_ylen/2d]
    
    multi_size = [!d.x_vsize*multi_xlen, !d.y_vsize*multi_ylen]
  endif
  
  ;; Work out plot & colorbar positions
  ;; in units of plot area (incl. margins)
  if n_elements(cb_size_in) eq 0 then cb_size = 0.025 else cb_size = cb_size_in
  if n_elements(margin_in) lt 4 then begin
    margin = [0.2, 0.17, 0.02, 0.1]
    if keyword_set(baseline_axis) and not keyword_set(no_title) then margin[3] = 0.15
    if keyword_set(delay_axis) then margin[2] = 0.07
  endif else margin = margin_in
  
  if n_elements(cb_margin_in) lt 2 then begin
    cb_margin = [0.2, 0.02]
  ;; if units_str ne '' then begin
  ;;    cb_margin = [0.2, 0.02]
  ;; endif else begin
  ;;    ;; no label on colorbar in this case
  ;;    cb_margin = [0.1, 0.02]
  ;; endelse
  endif else cb_margin = cb_margin_in
  
  plot_pos = [margin[0], margin[1], (1-cb_margin[1]-cb_size-cb_margin[0]-margin[2]), (1-margin[3])]
  cb_pos = [(1-cb_margin[1]-cb_size), margin[1], (1-cb_margin[1]), (1-margin[3])]
  
  plot_len = [plot_pos[2]-plot_pos[0], plot_pos[3] - plot_pos[1]]
  if min(plot_len) le 0 then stop
  
  plot_aspect = (plot_pos[3] - plot_pos[1]) / (plot_pos[2] - plot_pos[0])
  
  if log_axes[0] eq 0 then begin
    kperp_length = max(kperp_edges) - min(kperp_edges)
    kpar_length = max(kpar_edges) - min(kpar_edges)
  endif else begin
    if log_axes[1] eq 0 then kpar_length = alog10(max(kpar_edges)) - alog10(min(kpar_edges[where(kpar_edges gt 0)])) $
    else kpar_length = max(kpar_log_edges) - min(kpar_log_edges)
    if min(kperp_edges) le 0 then begin
      wh_zero = where(kperp_edges eq 0, n_zero)
      if n_zero ne 0 then stop
      
      wh_pos = where(kperp_edges ge 0, n_pos, complement = wh_neg, ncomplement = n_neg)
      if n_neg gt 0 then neg_leng = max(alog10((-1)*kperp_edges[wh_neg])) - min(alog10((-1)*kperp_edges[wh_neg]))
      if n_pos gt 0 then pos_leng = max(alog10(kperp_edges[wh_pos])) - min(alog10(kperp_edges[wh_pos]))
      
      kperp_length = neg_leng + pos_leng + pos_leng/(n_pos-1)
      
    endif else kperp_length = max(kperp_log_edges) - min(kperp_log_edges)
  endelse
  data_aspect = (kpar_length / kperp_length)
  
  aspect_ratio =  data_aspect /plot_aspect
  if aspect_ratio gt 1 then begin
    y_factor = aspect_ratio
    x_factor = 1.
  endif else begin
    y_factor = 1.
    x_factor = 1./aspect_ratio
  endelse
  
  if n_elements(multi_pos) eq 4 or n_elements(start_multi_params) gt 0 then begin
    if n_elements(start_multi_params) gt 0 then begin
      ;; calculate desired window size and positions for all plots
      ncol = start_multi_params.ncol
      nrow = start_multi_params.nrow
      
      multi_pos = fltarr(4, ncol*nrow)
      
      if tag_exist(start_multi_params, 'ordering') eq 0 then ordering = 'row' $
      else ordering = start_multi_params.ordering
      
      case ordering of
        'col': begin
          ;; col-major values
          col_val = reform(rebin(reform(indgen(ncol), 1, ncol), nrow, ncol), ncol*nrow)
          row_val = reverse(reform(rebin(indgen(nrow), nrow, ncol), ncol*nrow))
        end
        'row': begin
          ;; row-major values
          col_val = reform(rebin(indgen(ncol), ncol, nrow), ncol*nrow)
          row_val = reverse(reform(rebin(reform(indgen(nrow), 1, nrow), ncol, nrow), ncol*nrow))
        end
        else: message, 'unrecognized ordering value in start_multi_params, use "col" or "row" '
      endcase
      
      multi_pos[0,*] = col_val/double(ncol)
      multi_pos[1,*] = row_val/double(nrow)
      multi_pos[2,*] = (col_val+1)/double(ncol)
      multi_pos[3,*] = (row_val+1)/double(nrow)
      
      ;; define window size based on aspect ratio
      base_size_use = base_size
      xsize = round(base_size * x_factor * ncol)
      ysize = round(base_size * y_factor * nrow)
      while (ysize gt max_ysize) or (xsize gt max_xsize) do begin
        base_size_use = base_size_use - 100
        xsize = round(base_size_use * x_factor * ncol)
        ysize = round(base_size_use * y_factor * nrow)
      endwhile
      
      ;; make or set window
      if windowavailable(window_num) then begin
        wset, window_num
        if !d.x_size ne xsize or !d.y_size ne ysize then make_win = 1 else make_win = 0
      endif else make_win = 1
      if make_win eq 1 then window, window_num, xsize = xsize, ysize = ysize
      cgerase, background_color
      
      ;; if pub is set, start ps output
      if keyword_set(pub) then cgps_open, plotfile, /font, encapsulated=eps
      
      ;; calculate multi_size & multi x/ylen not calculated earlier
      multi_xlen = (multi_pos[2,0]-multi_pos[0,0])
      multi_ylen = (multi_pos[3,0]-multi_pos[1,0])
      multi_center = [multi_pos[0,0] + multi_xlen/2d, multi_pos[1,0] + multi_ylen/2d]
      
      multi_size = [!d.x_vsize*multi_xlen, !d.y_vsize*multi_ylen]
      
      multi_pos_use = multi_pos[*,0]
    endif else multi_pos_use = multi_pos
    
    multi_aspect = multi_size[1]/float(multi_size[0])
    
    new_aspect = aspect_ratio/multi_aspect
    if new_aspect gt 1 then begin
      y_factor = 1.
      x_factor = 1/new_aspect
    endif else begin
      y_factor = new_aspect
      x_factor = 1.
    endelse
    
    new_xlen = multi_xlen*x_factor
    new_ylen = multi_ylen*y_factor
    new_multi = [multi_center[0] - new_xlen/2d, multi_center[1] - new_ylen*y_factor/2d, $
      multi_center[0] + new_xlen/2d, multi_center[1] + new_ylen*y_factor/2d]
      
    new_pos = [new_xlen * plot_pos[0] + new_multi[0], new_ylen * plot_pos[1] + new_multi[1], $
      new_xlen * plot_pos[2] + new_multi[0], new_ylen * plot_pos[3] + new_multi[1]]
      
    new_cb_pos = [new_xlen * cb_pos[0] + new_multi[0], new_ylen * cb_pos[1] + new_multi[1], $
      new_xlen * cb_pos[2] + new_multi[0], new_ylen * cb_pos[3] + new_multi[1]]
      
    plot_pos = new_pos
    cb_pos = new_cb_pos
    
    no_erase = 1
  endif else begin
    xsize = round(base_size * x_factor)
    ysize = round(base_size * y_factor)
    
    if not keyword_set(pub) then begin
      while (ysize gt max_ysize) or (xsize gt max_xsize) do begin
        base_size = base_size - 100
        xsize = round(base_size * x_factor)
        ysize = round(base_size * y_factor)
      endwhile
    endif
    
    no_erase = 0
  endelse
  
  if not keyword_set(no_title) then begin
    xloc_title = (plot_pos[2] - plot_pos[0])/2. + plot_pos[0]
    if n_elements(multi_pos) gt 0 then yloc_title = plot_pos[3] + 0.6* (multi_pos_use[3]-plot_pos[3]) $
    else yloc_title = plot_pos[3] + 0.6* (1-plot_pos[3])
  endif
  
  if n_elements(multi_pos) gt 0 then begin
    xloc_lambda = plot_pos[0] - 0.15* (plot_pos[0]-multi_pos_use[0])
    yloc_lambda = plot_pos[3] + 0.15* (multi_pos_use[3]-plot_pos[3])
    
    ;; xloc_delay = plot_pos[2] + 0.2* (multi_pos[2]-plot_pos[2])
    ;; yloc_delay = (plot_pos[3] - plot_pos[1])/2d + plot_pos[1]
    xloc_delay = plot_pos[2] + 0.15 * (multi_pos_use[2]-plot_pos[2])
    yloc_delay = plot_pos[1] + 0.1* (plot_pos[1]-multi_pos_use[1])
  endif else begin
    xloc_lambda = plot_pos[0] - 0.1* (plot_pos[0]-0)
    yloc_lambda = plot_pos[3] + 0.1* (1-plot_pos[3])
    
    xloc_delay = plot_pos[2] + 0.1*(1-plot_pos[2])
    yloc_delay = plot_pos[1] +0.1*(plot_pos[1]-0)
  endelse
  
  if keyword_set(pub) then begin
    charthick = 3
    thick = 3
    xthick = 3
    ythick = 3
    if n_elements(charsize_in) eq 0 then begin
      if n_elements(multi_pos) gt 0 then begin
        charsize = 1.5d * (multi_size[0]/10000.)
      endif else charsize = 2
    endif else charsize = charsize_in
    
    font = 1
    
    if n_elements(multi_pos) eq 0 then begin
      window, window_num, xsize = xsize, ysize = ysize
      cgps_open, plotfile, /font, encapsulated=eps
    endif
    
    DEVICE, /ISOLATIN1
    perp_char = '!9' + String("136B) + '!X' ;"
    
    
  endif else begin
    charthick = 1
    thick = 1
    xthick = 1
    ythick = 1
    font = -1
    if n_elements(charsize_in) eq 0 then begin
      if n_elements(multi_pos) gt 0 then begin
        charsize = 2d * (multi_size[0]/float(base_size))
      endif else charsize = 2
    endif else charsize = charsize_in
    
    if n_elements(multi_pos) eq 0 then begin
      perp_char = '!9' + string(120B) + '!X'
      
      if windowavailable(window_num) then begin
        wset, window_num
        if !d.x_size ne xsize or !d.y_size ne ysize then make_win = 1 else make_win = 0
      endif else make_win = 1
      
      if make_win eq 1 then window, window_num, xsize = xsize, ysize = ysize
    endif
  endelse
  if log_axes[0] eq 1 then plot_kperp = 10^kperp_log_edges else plot_kperp = kperp_edges
  if log_axes[1] eq 1 then plot_kpar = 10^kpar_log_edges else plot_kpar = kpar_edges
  
  ;; if plot title includes sigma need to replace 'sigma' with appropriate character
  ;; (textoidl has to be called after cgps_open)
  if keyword_set(plot_sigma) then plot_title = repstr(plot_title, 'sigma', textoidl('\sigma', font=font))
  
  if keyword_set(title_prefix) then plot_title = title_prefix + ' ' + plot_title
  if n_elements(title) ne 0 then plot_title = title
  if keyword_set(no_title) then undefine, plot_title
  
  if keyword_set (hinv) then xtitle = textoidl('k_{perp} (h Mpc^{-1})', font = font) $
  else xtitle = textoidl('k_{perp} (Mpc^{-1})', font = font)
  xtitle = repstr(xtitle, 'perp', perp_char)
  if keyword_set (hinv) then ytitle = textoidl('k_{||} (h Mpc^{-1})', font = font) $
  else ytitle = textoidl('k_{||} (Mpc^{-1})', font = font)
  
  if keyword_set(no_title) or keyword_set(baseline_axis) then initial_title = '' else initial_title = plot_title
  
  ;;cgplot, plot_kperp, plot_kpar, /xlog, /ylog, /nodata, xstyle=5, ystyle=5, title = initial_title, $
  ;;        position = plot_pos, thick = thick, charthick = charthick, xthick = xthick, ythick = ythick, charsize = charsize, $
  ;;        font = font, noerase = no_erase, axiscolor = annotate_color, background = background_color
  
  ;;cgimage, power_log_norm, /nointerp,/overplot,/noerase
  
  if log_axes[0] eq 1 then xtickformat = 'exponent' else begin
    nticks = 5
    
    log_size = round(min(alog10(plot_kperp)))
    xtick_width = round((max(plot_kperp) - min(plot_kperp))/(nticks-1.)/(10.^log_size))*(10.^log_size)
    xtick_start = round(min(plot_kperp)/xtick_width)*(10.^log_size)
    xticks_in = round((dindgen(nticks)*xtick_width + xtick_start)/(10.^log_size))*(10.^log_size)
    
    x_nticks = nticks-1
    n_minor = 4
    
    if keyword_set(baseline_axis) then begin
      if keyword_set(hinv) then baseline_range = minmax(plot_kperp * hubble_param * kperp_lambda_conv) $
      else baseline_range = minmax(plot_kperp* kperp_lambda_conv)
      
      log_size2 = round(min(alog10(baseline_range)))
      xtick_width2 = round((max(baseline_range) - min(baseline_range))/(nticks-1.)/(10.^log_size2))*(10.^log_size2)
      xtick_start2 = round(min(baseline_range)/xtick_width2)*(10.^log_size2)
      xticks_in2 = round((dindgen(nticks)*xtick_width2 + xtick_start2)/(10.^log_size2))*(10.^log_size2)
    endif
    
  endelse
  if log_axes[1] eq 1 then ytickformat = 'exponent'
  
  
  axkeywords = {xlog: log_axes[0], ylog: log_axes[1], xstyle: 5, ystyle: 5, thick: thick, charthick: charthick, xthick: xthick, $
    ythick: ythick, charsize: charsize, font: font}
  cgimage, power_log_norm, /nointerp, xrange = minmax(plot_kperp), yrange = minmax(plot_kpar), $
    title=initial_title, position = plot_pos, noerase = no_erase, color = annotate_color, background = background_color, $
    axkeywords = axkeywords, /axes
    
  if keyword_set(plot_wedge_line) then begin
    n_lines = n_elements(wedge_amp)
    sorted_amp = reverse(wedge_amp[sort(wedge_amp)])
    if n_lines gt 1 then linestyles = [0, 2, 1] else linestyles=2
    
    for i=0, n_lines-1 do cgplot, /overplot, plot_kperp, plot_kperp * sorted_amp[i], color = annotate_color, thick = thick+1, $
      psym=-0, linestyle = linestyles[i]
  endif
  
  cgaxis, xaxis=0, xtick_get = xticks, xtickv = xticks_in, xticks = x_nticks, xminor=n_minor, xrange = minmax(plot_kperp), $
    xtitle = xtitle, charthick = charthick, xthick = xthick, ythick = ythick, charsize = charsize, font = font, $
    xtickformat = xtickformat, xstyle = 1, color = annotate_color
    
  cgaxis, yaxis=0, ytick_get = yticks, ytitle = ytitle, yrange = minmax(plot_kpar), $
    charthick = charthick, xthick = xthick, ythick = ythick, charsize = charsize, font = font, $
    ytickformat = ytickformat, ystyle = 1, color = annotate_color
  if keyword_set(baseline_axis) then begin
    ;; baselines don't care about hinv -- take it back out.
    if keyword_set(hinv) then baseline_range = minmax(plot_kperp * hubble_param * kperp_lambda_conv) $
    else baseline_range = minmax(plot_kperp* kperp_lambda_conv)
    
    if keyword_set(no_title) then xtitle = textoidl('(\lambda)', font = font) else undefine, xtitle
    
    cgaxis, xaxis=1, xtickv = xticks_in2, xticks = x_nticks, xminor=n_minor, xrange = baseline_range, xtickformat = xtickformat, $
      xthick = xthick, xtitle = xtitle, $
      charthick = charthick, ythick = ythick, charsize = charsize, font = font, xstyle = 1, color = annotate_color
      
    if not keyword_set(no_title) then begin
      cgtext, xloc_title, yloc_title, plot_title, /normal, alignment=0.5, charsize=1.2 * charsize, $
        color = annotate_color, font = font
      cgtext, xloc_lambda, yloc_lambda, textoidl('(\lambda)', font = font), /normal, alignment=0.5, charsize=charsize, $
        color = annotate_color, font = font
    endif
  endif else $
    cgaxis, xaxis=1, xrange = minmax(plot_kperp), xtickv = xticks, xtickname = replicate(' ', n_elements(xticks)), $
    charthick = charthick, xthick = xthick, ythick = ythick, charsize = charsize, font = font, xstyle = 1, $
    color = annotate_color
  if keyword_set(delay_axis) then begin
    min_delay_plot = 2d*alog10(delay_params[0]) - alog10(delay_params[0]*2) ;; in analogy with min kperp/par with 0 bins
    cgaxis, yaxis=1, yrange = [min_delay_plot, delay_params[1]], ytickformat = ytickformat, charthick = charthick, xthick = xthick, $
      ythick = ythick, charsize = charsize, font = font, ystyle = 1, color = annotate_color
      
    cgtext, xloc_delay, yloc_delay, '(ns)', /normal, alignment=0.5, charsize=charsize*0.9, $
      color = annotate_color, font = font
      
  endif else $
    cgaxis, yaxis=1, yrange = minmax(plot_kpar), ytickv = yticks, ytickname = replicate(' ', n_elements(yticks)), $
    charthick = charthick, xthick = xthick, ythick = ythick, charsize = charsize, font = font, ystyle = 1, $
    color = annotate_color
    
  tick_vals = loglevels(10d^[floor(log_data_range[0])-.1, ceil(log_data_range[1])+.1], coarse=0)
  
  wh_keep = where(tick_vals gt 10^(log_data_range[0]-0.001) and tick_vals lt 10^(log_data_range[1]+0.001), count_keep)
  if count_keep gt 0 then tick_vals = tick_vals[wh_keep] else undefine, tick_vals
  
  ;; want minor tick marks if there aren't very many loglevels.
  ;; unfortunately cgcolorbar can't do log minor tick marks
  ;; with specified tick locations (which I have to do b/c
  ;; can't use divisions=0 with formatting keyword)
  ;; solution: add regular tickmarks without labels for the minor ones.
  
  if n_elements(tick_vals) lt 4 then begin
    tick_vals_use = loglevels(10d^[floor(log_data_range[0])-.1, ceil(log_data_range[1])+.1], coarse=0)
    
    if n_elements(tick_vals) lt 2 then minor_multipliers = dindgen(8)+2 else minor_multipliers = (dindgen(4)+1)*2d
    n_minor_mult = n_elements(minor_multipliers)
    n_major = n_elements(tick_vals_use)
    
    minor_tick_vals = reform(rebin(minor_multipliers, n_minor_mult, n_major), n_minor_mult*n_major) $
      * reform(rebin(reform(tick_vals_use, 1, n_major), n_minor_mult, n_major), n_minor_mult*n_major)
    wh_keep = where(minor_tick_vals gt 10^log_data_range[0] and minor_tick_vals lt 10^log_data_range[1], count_keep)
    
    if count_keep gt 0 then begin
      minor_tick_vals = minor_tick_vals[wh_keep]
      n_minor = n_elements(minor_tick_vals)
      minor_tick_names = strarr(n_minor) + ' '
    endif else n_minor = 0
    
    if n_elements(tick_vals) eq 0 then begin
      ;; range doesn't include any decade levels. Use minor ticks as major ticks
      if n_minor gt 0 then begin
        tick_vals = minor_tick_vals
        n_minor = 0
      endif else stop
    endif
    
  endif else n_minor = 0
  
  nloop = 0
  while(n_elements(tick_vals) gt 8) do begin
    nloop = nloop + 1
    factor = double(nloop+1)
    if color_profile eq 'sym_log' then begin
      pos_exp_vals = dindgen(ceil((alog10(max(tick_vals))-alog10(min(tick_vals)) + 1)/(2d*factor)) + 1)*factor
      if max(pos_exp_vals) gt temp[1] then pos_exp_vals = pos_exp_vals[0:n_elements(pos_exp_vals)-2]
      
      exp_vals = [(-1)*reverse(pos_exp_vals[1:*]), pos_exp_vals]
    endif else begin
      exp_vals = (dindgen(ceil((alog10(max(tick_vals))-alog10(min(tick_vals)) + 1)/factor) + 1)*factor + alog10(min(tick_vals)))
      if max(exp_vals) gt alog10(max(tick_vals)) then exp_vals = exp_vals[0:n_elements(exp_vals)-2]
    endelse
    tick_vals = 10^exp_vals
  endwhile
  
  if color_profile eq 'sym_log' then begin
    names = strarr(n_elements(tick_vals))
    
    wh_neg = where(tick_vals lt 0, count_neg)
    wh_pos = where(tick_vals gt 0, count_pos)
    wh_zero = where(tick_vals eq 0, count_zero)
    if count_pos gt 0 then names[wh_pos] = number_formatter(tick_vals[wh_pos], format = '(e0)',/print_exp)
    if count_neg gt 0 then names[wh_neg] = '-' + number_formatter(abs(tick_vals[wh_neg]), format = '(e0)',/print_exp)
    if count_zero gt 0 then names[wh_zero] = '0'
  endif else begin
  
    if min_pos_color gt 0 then begin
      new_tick_vals = [10^log_data_range[0], tick_vals]
      
      min_tick_color = (alog10(min(tick_vals))-log_data_range[0])*n_colors/(log_data_range[1]-log_data_range[0]) + color_range[0]
      if min_tick_color lt 10 then names = ['<0', ' ', number_formatter(tick_vals[1:*], format = '(e0)',/print_exp)] $
      else names = ['<0',  number_formatter(tick_vals, format = '(e0)',/print_exp)]
      
      tick_vals = new_tick_vals
    endif else names = number_formatter(tick_vals, format = '(e0)',/print_exp)
    
  endelse
  
  if n_minor gt 0 then begin
    temp_ticks = [tick_vals, minor_tick_vals]
    temp_names = [names, minor_tick_names]
    order = sort(temp_ticks)
    
    tick_vals = temp_ticks[order]
    names = temp_names[order]
  endif
  
  
  if (alog10(tick_vals[0]) - log_data_range[0]) lt 10^(-3d) then begin
    cb_ticknames = [' ', names]
    cb_ticks = [color_range[0]-1, (alog10(tick_vals) - log_data_range[0]) * n_colors / $
      (log_data_range[1] - log_data_range[0]) + color_range[0]] - color_range[0]
      
  endif else begin
    cb_ticknames = names
    cb_ticks = ((alog10(tick_vals) - log_data_range[0]) * (n_colors+1) / $
      (log_data_range[1] - log_data_range[0]) + color_range[0]) - color_range[0]
  endelse
  
  if (log_data_range[1] - alog10(max(tick_vals))) gt 10^(-3d) then begin
    cb_ticknames = [cb_ticknames, ' ']
    cb_ticks = [cb_ticks, color_range[1]-color_range[0]+1]
  endif
  
  min_pos_cb_val = ((alog10(min_pos) - log_data_range[0]) * n_colors / $
    (log_data_range[1] - log_data_range[0]) + color_range[0]) - color_range[0]
    
  cgcolorbar, color = annotate_color, /vertical, position = cb_pos, bottom = color_range[0], ncolors = n_colors, minor = 0, $
    ticknames = cb_ticknames, ytickv = cb_ticks, yticks = n_elements(cb_ticks) -1, title = units_str, $
    charsize = charsize, font = font, oob_low = oob_low
  ;; cgcolorbar, color = annotate_color, /vertical, position = cb_pos, bottom = color_range[0]+1, ncolors = n_colors, /ylog, $
  ;;             range = 10d^log_data_range, format = 'exponent', minor=minor, charsize = charsize, font = font, $
  ;;             divisions = 0
    
  if keyword_set(pub) and n_elements(multi_pos) eq 0 then begin
    cgps_close, png = png, delete_ps = delete_ps
    wdelete, window_num
  endif
  
  tvlct, r, g, b
  
end
