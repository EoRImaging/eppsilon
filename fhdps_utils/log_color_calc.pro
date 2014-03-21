pro log_color_calc, data, data_log_norm, cb_ticks, cb_ticknames, color_range, n_colors, data_range = data_range, $
    color_profile = color_profile, log_cut_val = log_cut_val, grey_scale = grey_scale, oob_low = oob_low, $
    missing_value = missing_value, missing_color = missing_color
    
    
  color_profile_enum = ['log_cut', 'sym_log', 'abs']
  if n_elements(color_profile) eq 0 then color_profile = 'log_cut'
  
  wh_prof = where(color_profile_enum eq color_profile, count)
  if count eq 0 then message, 'Color profile must be one of: ' + strjoin(color_profile_enum, ', ')
  
  if n_elements(data_range) eq 0 then begin
    no_input_data_range = 1
    data_range = minmax(data)
  endif else begin
    if n_elements(data_range) ne 2 then message, 'data_range must be a 2 element vector'
    no_input_data_range = 0
  endelse
  
  if data_range[1] lt data_range[0] then message, 'data_range[0] must be less than data_range[1]'
  
  if color_profile eq 'sym_log' and data_range[0] gt 0 then begin
    print, 'sym_log profile cannot be selected with an entirely positive data range. Switching to log_cut'
    color_profile = 'log_cut'
  endif
  
  color_range = [0, 255]
  if n_elements(missing_value) ne 0 then begin
    wh_missing = where(data eq missing_value, count_missing)
    if count_missing gt 0 then begin
      missing_color = 255
      data_color_range = [0, 254]
    endif else data_color_range = color_range
  endif else data_color_range = [0, 255]
  n_colors = color_range[1] - color_range[0] + 1
  data_n_colors = data_color_range[1] - data_color_range[0] + 1
  
  wh = where(data gt 0d, count)
  if count gt 0 then min_pos = min(data[wh]) else if data_range[0] gt 0 then min_pos = data_range[0] else $
    if data_range[1] gt 0 then min_pos = data_range[1]/10d else min_pos = 0.01d
  wh = where(data lt 0d, count)
  if count gt 0 then max_neg = max(data[wh]) else if data_range[1] lt 0 then max_neg = data_range[1] else $
    if data_range[0] lt 0 then max_neg = data_range[0]/10d
    
  case color_profile of
    'log_cut': begin
    
      if data_range[1] lt 0 then message, 'log_cut color profile will not work for entirely negative arrays.'
      
      if n_elements(log_cut_val) eq 0 then begin
        if data_range[0] gt 0 then log_cut_val = alog10(data_range[0]) else $
          log_cut_val = alog10(min_pos)
      endif
      
      log_data_range = [log_cut_val, alog10(data_range[1])]
      
      ;; if log_cut_val is less than min_pos then indicate it in the color bar
      if min(data) lt 0 then begin
      
      
        wh_zero = where(data eq 0, count_zero)
        if count_zero gt 0 then begin
          min_pos_color = 2
          zero_color=1
          zero_val = data_n_colors/(log_data_range[1]-log_data_range[0])
          
          if log_cut_val gt min_pos then log_data_range = [log_cut_val-2*data_n_colors/(log_data_range[1]-log_data_range[0]), alog10(data_range[1])]
          
        endif else begin
          min_pos_color = 1
          if log_cut_val gt min_pos then log_data_range = [log_cut_val-data_n_colors/(log_data_range[1]-log_data_range[0]), alog10(data_range[1])]
          
        endelse
        neg_color=0
        neg_val = log_data_range[0]
        
        data_log = alog10(data)
        wh_under = where(data lt 10^double(log_cut_val), count)
        if count ne 0 then data_log[wh_under] = log_data_range[0]
        wh_neg = where(data lt 0, count)
        if count ne 0 then data_log[wh_neg] = neg_val
        wh_over = where(data_log gt log_data_range[1], count)
        if count ne 0 then data_log[wh_over] = log_data_range[1]
        if count_zero ne 0 then data_log[wh_zero] = zero_val
        
        oob_low = 0
        
        cgloadct, 25, /brewer, /reverse, BOTTOM = min_pos_color, NCOLORS = 256-min_pos_color, clip = [0, 235]
        if count_zero gt 0 then begin
          cgLoadCT, 25, /brewer, /reverse, CLIP=[255, 255], BOTTOM=0, NCOLORS=1
          cgLoadCT, 25, /brewer, /reverse, CLIP=[245, 245], BOTTOM=zero_color, NCOLORS=1
        endif else cgLoadCT, 25, /brewer, /reverse, CLIP=[255, 255], BOTTOM=0, NCOLORS=1
        
      ;        cgloadct, 27, /brewer, BOTTOM = min_pos_color, NCOLORS = 256-min_pos_color, clip = [40, 255]
      ;        if count_zero gt 0 then begin
      ;          cgLoadCT, 27, /brewer, CLIP=[0, 0], BOTTOM=0, NCOLORS=1
      ;          cgLoadCT, 27, /brewer, CLIP=[1, 1], BOTTOM=zero_color, NCOLORS=1
      ;        endif else cgLoadCT, 27, /brewer, CLIP=[0, 0], BOTTOM=0, NCOLORS=1
        
      endif else begin
      
        min_pos_color = 0
        
        data_log = alog10(data)
        wh_under = where(data lt 10^double(log_cut_val), count)
        if count ne 0 then data_log[wh_under] = log_data_range[0]
        wh_over = where(data_log gt log_data_range[1], count)
        if count ne 0 then data_log[wh_over] = log_data_range[1]
        
        cgloadct, 25, /brewer, /reverse
        
      endelse
      
      if no_input_data_range eq 1 then data_range = 10^log_data_range
      
      data_log_norm = (data_log-log_data_range[0])*data_n_colors/(log_data_range[1]-log_data_range[0]) + data_color_range[0]
      
    end
    'sym_log': begin
      ;; find the middle of our color range
      if n_elements(data_range) gt 0 then max_abs = max(abs(data_range)) else max_abs = max(abs(data[where(abs(data) ne 0)]))
      min_abs = min(abs(data[where(abs(data) gt 0)]))
      
      log_data_range = alog10([min_abs, max_abs])
      
      neg_color_range = [0, floor(data_n_colors/2)-1]
      zero_color = floor(data_n_colors/2)
      pos_color_range = [floor(data_n_colors/2),data_n_colors]
      n_pos_neg_colors = floor(data_n_colors/2)
      
      cgLoadCT, 16, /brewer, /reverse, clip=[20, 220], bottom=0, ncolors=floor(data_n_colors/2)
      cgloadct, 0, clip = [255, 255], bottom = zero_color, ncolors = 1
      cgLoadCT, 13, /brewer, clip=[20, 220], bottom=zero_color+1, ncolors=data_n_colors-(zero_color+1)
      
      
      ;; construct 2 separate data logs (for negative & positive)
      pos_data_log = alog10(data)
      neg_data_log = alog10(-1*data)
      
      data_log_norm = data*0.
      wh_pos = where(data gt 0, count_pos)
      if count_pos gt 0 then data_log_norm[wh_pos] = (pos_data_log[wh_pos]-log_data_range[0])*n_pos_neg_colors/(log_data_range[1]-log_data_range[0]) + pos_color_range[0]
      wh_neg = where(data lt 0, count_neg)
      if count_neg gt 0 then data_log_norm[wh_neg] = (log_data_range[1] - neg_data_log[wh_neg])*n_pos_neg_colors/(log_data_range[1]-log_data_range[0]) + neg_color_range[0]
      wh_0 = where(data eq 0, count_0)
      if count_0 gt 0 then data_log_norm[wh_0] = zero_color
      
    end
    'abs': begin
    
      abs_data = abs(data)
      log_data_range = dblarr(2)
      if data_range[0] lt 0 then log_data_range[0] = alog10(min(abs_data[where(abs_data gt 0)])) $
      else log_data_range[0] = alog10(data_range[0])
      log_data_range[1] = alog10(max(abs(data_range)))
      
      data_log = alog10(abs_data)
      wh_zero = where(data eq 0, count)
      if count ne 0 then data_log[wh_zero] = log_data_range[0]
      
      abs_data = 0
      
      data_log_norm = (data_log-log_data_range[0])*data_n_colors/(log_data_range[1]-log_data_range[0]) + data_color_range[0]
      cgloadct, 25, /brewer, /reverse
      
    end
  endcase
  
  if n_elements(missing_value) gt 0 then begin
    if count_missing gt 0 then data_log_norm[wh_missing] = missing_color
  endif
  
  
  if color_profile eq 'sym_log' then begin
  
    pos_tick_vals = loglevels(10d^[floor(log_data_range[0])-.1, ceil(log_data_range[1])+.1], coarse=0)
    wh_keep = where(pos_tick_vals gt 10^(log_data_range[0]-0.001) and pos_tick_vals lt 10^(log_data_range[1]+0.001), count_keep)
    if count_keep gt 0 then pos_tick_vals = pos_tick_vals[wh_keep] else stop
    
    neg_tick_vals = reverse(pos_tick_vals)
    
    nloop = 0
    while(n_elements(pos_tick_vals) gt 4) do begin
      nloop = nloop + 1
      factor = double(nloop+1)
      
      pos_exp_vals = (dindgen(ceil((alog10(max(pos_tick_vals))-alog10(min(pos_tick_vals)) + 1)/factor) + 1)*factor + alog10(min(pos_tick_vals)))
      if max(pos_exp_vals) gt alog10(max(pos_tick_vals)) then pos_exp_vals = pos_exp_vals[0:n_elements(pos_exp_vals)-2]
      pos_tick_vals = 10^pos_exp_vals
      
      neg_tick_vals = reverse(pos_tick_vals)
    endwhile
    
    top_neg_color = max((log_data_range[1] - alog10(neg_tick_vals)) * n_pos_neg_colors / (log_data_range[1] - log_data_range[0]) + neg_color_range[0])
    bot_pos_color = min((alog10(pos_tick_vals) - log_data_range[0]) * n_pos_neg_colors / (log_data_range[1] - log_data_range[0]) + pos_color_range[0])
    if abs(zero_color - top_neg_color) lt 15 then neg_tick_vals = neg_tick_vals[0:n_elements(neg_tick_vals)-2]
    if abs(zero_color - bot_pos_color) lt 15 then pos_tick_vals = pos_tick_vals[1:*]
    
    neg_names = '-' + number_formatter(neg_tick_vals, format = '(e0)',/print_exp)
    pos_names = number_formatter(pos_tick_vals, format = '(e0)',/print_exp)
    zero_name = number_formatter(0) + '  '
    
    if (log_data_range[1] - alog10(max(neg_tick_vals))) gt 10^(-3d) then begin
      cb_ticknames_neg = [' ', neg_names]
      cb_ticks_neg = [neg_color_range[0]-1, (log_data_range[1] - alog10(neg_tick_vals)) * n_pos_neg_colors / $
        (log_data_range[1] - log_data_range[0]) + neg_color_range[0]]
        
    endif else begin
      cb_ticknames_neg = neg_names
      cb_ticks_neg = (log_data_range[1] - alog10(neg_tick_vals)) * (n_pos_neg_colors+1) / $
        (log_data_range[1] - log_data_range[0]) + neg_color_range[0]
    endelse
    
    if (log_data_range[1] - alog10(max(pos_tick_vals))) gt 10^(-3d) then begin
      cb_ticknames_pos = [pos_names, ' ']
      cb_ticks_pos = [(alog10(pos_tick_vals) - log_data_range[0]) * n_pos_neg_colors / $
        (log_data_range[1] - log_data_range[0]) + pos_color_range[0], pos_color_range[1]+1]
    endif else begin
      cb_ticknames_pos = pos_names
      cb_ticks_pos = (alog10(pos_tick_vals) - log_data_range[0]) * (n_pos_neg_colors+1) / $
        (log_data_range[1] - log_data_range[0]) + pos_color_range[0]
    endelse
    
    cb_ticks = [cb_ticks_neg, zero_color, cb_ticks_pos]
    cb_ticknames = [cb_ticknames_neg, zero_name, cb_ticknames_pos]
    
  endif else begin
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
      exp_vals = (dindgen(ceil((alog10(max(tick_vals))-alog10(min(tick_vals)) + 1)/factor) + 1)*factor + alog10(min(tick_vals)))
      if max(exp_vals) gt alog10(max(tick_vals)) then exp_vals = exp_vals[0:n_elements(exp_vals)-2]
      tick_vals = 10^exp_vals
    endwhile
    
    if min_pos_color gt 0 then begin
      new_tick_vals = [10^log_data_range[0], tick_vals]
      
      min_tick_color = (alog10(min(tick_vals))-log_data_range[0])*n_colors/(log_data_range[1]-log_data_range[0]) + color_range[0]
      if min_tick_color lt 10 then names = ['<0', ' ', number_formatter(tick_vals[1:*], format = '(e0)',/print_exp)] $
      else names = ['<0',  number_formatter(tick_vals, format = '(e0)',/print_exp)]
      
      tick_vals = new_tick_vals
    endif else names = number_formatter(tick_vals, format = '(e0)',/print_exp)
    
    
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
  endelse
  
end