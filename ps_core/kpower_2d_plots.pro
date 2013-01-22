pro kpower_2d_plots, power_savefile, multi_pos = multi_pos, multi_aspect = multi_aspect, plot_weights = plot_weights, $
                     kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, data_range = data_range, $
                     color_profile = color_profile, log_cut_val = log_cut_val, pub = pub, plotfile = plotfile, no_title = no_title, $
                     window_num = window_num, title = title, norm_2d = norm_2d, norm_factor = norm_factor, grey_scale = grey_scale, $
                     wedge_amp = wedge_amp, plot_wedge_line = plot_wedge_line, baseline_axis = baseline_axis, ratio = ratio, $
                     no_units = no_units, charsize = charsize_in, cb_size = cb_size_in, margin=margin_in, cb_margin = cb_margin_in

  if n_elements(window_num) eq 0 then window_num = 1
 
  if n_elements(multi_pos) gt 0 then begin
     if n_elements(multi_pos) ne 4 then message, 'multi_pos must be a 4 element plot position vector'
     if max(multi_pos) gt 1 or min(multi_pos) lt 0 then message, 'multi_pos must be in normalized coordinates (between 0 & 1)'
     if multi_pos[2] le multi_pos[0] or multi_pos[3] le multi_pos[1] then $
        message, 'In multi_pos, x1 must be greater than x0 and y1 must be greater than y0 '

     if n_elements(multi_aspect) eq 0 then begin
        print, 'No aspect ratio for multi_pos supplied. Assuming aspect = 1'
        multi_aspect = 1
     endif else if n_elements(multi_aspect) gt 1 then message, 'too many elements in multi_aspect'
  endif

  color_profile_enum = ['log_cut', 'sym_log', 'abs']
  if n_elements(color_profile) eq 0 then color_profile = 'log_cut'

  wh_prof = where(color_profile_enum eq color_profile, count)
  if count eq 0 then message, 'Color profile must be one of: ' + strjoin(color_profile_enum, ', ')

  if keyword_set(ratio) then begin
     if n_elements(power_savefile) gt 2 then message, 'Only 2 files can be specified in power_savefile if ratio keyword is set'
     restore, power_savefile[0]
     power1 = power
     weights1 = weights
     kperp1 = kperp_edges
     kpar1 = kpar_edges
     bins1 = bins_per_decade
     
     restore, power_savefile[1]
     power2 = power
     weights2 = weights
     if total(abs(kperp1 - kperp_edges)) ne 0 or total(abs(kpar1 - kpar_edges)) ne 0 or $
        total(abs(size(power1,/dimension) - size(power2, /dimension))) ne 0 then $
           message, 'dimensions and kperp/kpar edges must be the same in both files'

     power = power1 / power2

  endif else begin
     if n_elements(power_savefile) gt 1 then message, 'Only 1 file can be specified in power_savefile unless ratio keyword is set'

     restore, power_savefile
  endelse

  if keyword_set(plot_weights) then begin
     if n_elements(weights) ne 0 then if keyword_set(ratio) then power = weights1 / weights2 else power = weights $
     else message, 'No weights array included in this file'
  endif

  dims = size(power, /dimension)
  n_kperp = dims[0]
  n_kpar = dims[1]

  if n_elements(kperp_plot_range) eq 0 then kperp_plot_range = minmax(kperp_edges)
  if n_elements(kpar_plot_range) eq 0 then kpar_plot_range = minmax(kpar_edges)

  if n_elements(plotfile) eq 0 then begin
     if keyword_set(ratio) then plotfile = strsplit(power_savefile[0], '.idlsave', /regex, /extract) + '_ratio.eps' $
     else plotfile = strsplit(power_savefile, '.idlsave', /regex, /extract) + '_2dkplot.eps'
  endif else if strcmp(strmid(plotfile, strlen(plotfile)-4), '.eps', /fold_case) eq 0 then plotfile = plotfile + '.eps'
  
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
  kperp_log_diffs = (alog10(kperp_edges) - shift(alog10(kperp_edges), 1))[1:*]
  if total(abs(kperp_log_diffs - kperp_log_diffs[0])) gt n_kperp*1e-15 then log_bins[0] = 0
  kpar_log_diffs = (alog10(kpar_edges) - shift(alog10(kpar_edges), 1))[2:*]
  if total(abs(kpar_log_diffs - kpar_log_diffs[0])) gt n_kpar*1e-15 then log_bins[1] = 0

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

  if n_elements(data_range) eq 0 then data_range = minmax(power)
  wh = where(power gt 0d, count)
  if count gt 0 then min_pos = min(power[wh]) else if data_range[0] gt 0 then min_pos = data_range[0] else $
     if data_range[1] gt 0 then min_pos = data_range[1]/10d else min_pos = 0.01d

  if total(log_bins) eq 2 then begin
     ;; expand image array to prevent interpolation in postscript, normalize
     power_plot = congrid(power, n_kperp*10, n_kpar*10)
     kperp_log_edges = alog10(kperp_edges)
     kpar_log_edges = alog10(kpar_edges)
  endif else begin
     ;; make a new image array to allow log axes
     
     if log_bins[0] eq 0 then begin
        wh_kperp0 = where(kperp_edges lt 0, count_kperp0, complement = wh_kperp_good)
        if count_kperp0 gt 1 then stop
        
        kperp_log_edges = alog10(kperp_edges)
        if count_kperp0 eq 1 then begin
           kperp_log_diffs = (kperp_log_edges[1:*] - shift(kperp_log_edges[1:*], 1))[1:*]
           kperp_log_diffs = [kperp_log_diffs[0], kperp_log_diffs]
           kperp_log_edges[wh_kperp0] = kperp_log_edges[wh_kperp0+1] - kperp_log_diffs[wh_kperp0]
        endif

        image_kperp_delta = min(kperp_log_diffs)
     endif else begin
        kperp_log_edges = alog10(kperp_edges)
        image_kperp_delta = min(kperp_log_diffs)/10d
     endelse

     if log_bins[1] eq 0 then begin
        wh_kpar0 = where(kpar_edges lt 0, count_kpar0, complement = wh_kpar_good)
        if count_kpar0 gt 1 then stop
        
        kpar_log_edges = alog10(kpar_edges)
        if count_kpar0 eq 1 then begin
           kpar_log_diffs = (kpar_log_edges[1:*] - shift(kpar_log_edges[1:*], 1))[1:*]
           kpar_log_diffs = [kpar_log_diffs[0], kpar_log_diffs]
           kpar_log_edges[wh_kpar0] = kpar_log_edges[wh_kpar0+1] - kpar_log_diffs[wh_kpar0]
        endif

        image_kpar_delta = min(kpar_log_diffs)/2d
     endif else begin
        kpar_log_edges = alog10(kpar_edges)
        image_kpar_delta = min(kpar_log_diffs)/10d
     endelse

     ;; now get width for each input bin in image array
     kperp_bin_widths = round(kperp_log_diffs / image_kperp_delta)
     nkperp_image = total(kperp_bin_widths)
     kpar_bin_widths = round(kpar_log_diffs / image_kpar_delta) 
     nkpar_image = total(kpar_bin_widths)
    
     h_kperp = histogram(total(kperp_bin_widths,/cumulative)-1,binsize=1, min=0, reverse_indices=ri_kperp)
     undefine, h_kperp
     kperp_inds = rebin(ri_kperp[0:nkperp_image-1]-ri_kperp[0], nkperp_image, nkpar_image)
     
     h_kpar = histogram(total(kpar_bin_widths,/cumulative)-1,binsize=1, min=0, reverse_indices=ri_kpar)
     undefine, h_kpar
     kpar_inds = rebin(reform(ri_kpar[0:nkpar_image-1]-ri_kpar[0], 1, nkpar_image), nkperp_image, nkpar_image)

     power_plot = power[kperp_inds, kpar_inds]
  endelse

  case color_profile of
     'log_cut': begin

        if data_range[1] lt 0 then message, 'log_cut color profile will not work for entirely negative arrays.'
        
        if n_elements(log_cut_val) eq 0 then begin 
           if data_range[0] gt 0 then log_cut_val = alog10(data_range[0]) else $
              log_cut_val = alog10(min_pos)
        endif

        log_data_range = [log_cut_val, alog10(data_range[1])]

        power_log = alog10(power_plot)
        wh_under = where(power_plot lt 10^double(log_cut_val), count)
        if count ne 0 then power_log[wh_under] = log_data_range[0]
        wh_over = where(power_log gt log_data_range[1], count)
        if count ne 0 then power_log[wh_over] = log_data_range[1]
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


  if n_elements(multi_pos) eq 4 then begin
     ;; work out positions scaled to the area allowed in multi_pos with proper aspect ratio
     multi_xlen = (multi_pos[2]-multi_pos[0])
     multi_ylen = (multi_pos[3]-multi_pos[1])
     multi_center = [multi_pos[0] + multi_xlen/2d, multi_pos[1] + multi_ylen/2d]

     min_len = min([multi_xlen, multi_ylen])
  endif

  ;; Work out plot & colorbar positions
  ;; in units of plot area (incl. margins)
  if n_elements(cb_size_in) eq 0 then cb_size = 0.025 else cb_size = cb_size_in
  if n_elements(margin_in) lt 4 then begin
     margin = [0.15, 0.15, 0.02, 0.1] 
     if keyword_set(baseline_axis) and not keyword_set(no_title) then $
        if n_elements(multi_pos) gt 0 then margin[3] = 0.15 else margin[3] = 0.15
  endif else margin = margin_in

  if n_elements(cb_margin_in) lt 2 then begin
     cb_margin = [0.08, 0.02] 
     if n_elements(multi_pos) gt 0 then cb_margin[0] = 0.08/multi_xlen
  endif else cb_margin = cb_margin_in 
  
  plot_pos = [margin[0], margin[1], (1-cb_margin[1]-cb_size-cb_margin[0]-margin[2]), (1-margin[3])]
  cb_pos = [(1-cb_margin[1]-cb_size), margin[1], (1-cb_margin[1]), (1-margin[3])]
  
  plot_aspect = (plot_pos[3] - plot_pos[1]) / (plot_pos[2] - plot_pos[0])
  
  kpar_length_log = max(kpar_log_edges) - min(kpar_log_edges)
  if min(kperp_edges) le 0 then begin
     wh_zero = where(kperp_edges eq 0, n_zero)
     if n_zero ne 0 then stop
     
     wh_pos = where(kperp_edges ge 0, n_pos, complement = wh_neg, ncomplement = n_neg)
     if n_neg gt 0 then neg_leng = max(alog10((-1)*kperp_edges[wh_neg])) - min(alog10((-1)*kperp_edges[wh_neg]))
     if n_pos gt 0 then pos_leng = max(alog10(kperp_edges[wh_pos])) - min(alog10(kperp_edges[wh_pos]))
     
     kperp_length_log = neg_leng + pos_leng + pos_leng/(n_pos-1)
     
  endif else kperp_length_log = max(kperp_log_edges) - min(kperp_log_edges)
  data_aspect = (kpar_length_log / kperp_length_log)

  aspect_ratio =  data_aspect /plot_aspect
  if aspect_ratio gt 1 then begin
     y_factor = aspect_ratio
     x_factor = 1.
  endif else begin   
     y_factor = 1.
     x_factor = 1./aspect_ratio
  endelse
  
  if n_elements(multi_pos) eq 4 then begin
     
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
     base_size = 600
     xsize = round(base_size * x_factor)
     ysize = round(base_size * y_factor)
     
     if not keyword_set(pub) then begin
        while ysize gt 1100 do begin
           base_size = base_size - 100
           xsize = round(base_size * x_factor)
           ysize = round(base_size * y_factor)
        endwhile
     endif

     no_erase = 0
  endelse
  
  if not keyword_set(no_title) then begin
     xloc_title = (plot_pos[2] - plot_pos[0])/2. + plot_pos[0]
     if n_elements(multi_pos) gt 0 then yloc_title = plot_pos[3] + 0.6* (multi_pos[3]-plot_pos[3]) $
     else yloc_title = plot_pos[3] + 0.6* (1-plot_pos[3])
  endif

  if n_elements(multi_pos) gt 0 then begin
     xloc_lambda = plot_pos[0] - 0.3* (plot_pos[0]-multi_pos[0])
     yloc_lambda = plot_pos[3] + 0.15* (multi_pos[3]-plot_pos[3])
  endif else begin
     xloc_lambda = plot_pos[0] - 0.1* (plot_pos[0]-0)
     yloc_lambda = plot_pos[3] + 0.1* (1-plot_pos[3])
  endelse
  if keyword_set(no_title) then xloc_lambda = (plot_pos[2] - plot_pos[0])/2d + plot_pos[0]

  if keyword_set(pub) then begin
     charthick = 3
     thick = 3
     xthick = 3
     ythick = 3
     if n_elements(charsize_in) eq 0 then begin
        if n_elements(multi_pos) gt 0 then begin
           charsize = 5d * min_len
        endif else charsize = 2
     endif else charsize = charsize_in
     font = 1
     ;;perp_char = '!Z(22A5)'
     ;;perp_char = '!Z(27C2)'
     perp_char = 'perp'

     if n_elements(multi_pos) eq 0 then begin
        window, window_num, xsize = xsize, ysize = ysize
        pson, file = plotfile, /eps 
     endif
     ;;Device, /ISOLATIN1, Set_Font='Symbol', /TT_Font

  endif else begin
     charthick = 1
     thick = 1
     xthick = 1
     ythick = 1
     font = -1
     if n_elements(charsize_in) eq 0 then charsize=1 else charsize = charsize_in
     if n_elements(multi_pos) eq 0 then begin
        perp_char = '!9' + string(120B) + '!X'

        if windowavailable(window_num) then begin 
           wset, window_num
           if !d.x_size ne xsize or !d.y_size ne ysize then make_win = 1 else make_win = 0
        endif else make_win = 1
        
        if make_win eq 1 then window, window_num, xsize = xsize, ysize = ysize
     endif
  endelse
  plot_kperp = 10^kperp_log_edges
  plot_kpar = 10^kpar_log_edges
 
  units_str = ''
  if not keyword_set(ratio) and not keyword_set(no_units) then units_str = textoidl(' (mK^2 Mpc^3)', font = font)
 
  if not keyword_set(no_title) then begin
     if keyword_set(ratio) then begin
        if n_elements(title) ne 0 then plot_title = title else begin
           if keyword_set(plot_weights) then plot_title = 'Weights Ratio' else plot_title = 'Power Ratio'
        endelse
     endif else begin
        if n_elements(title) ne 0 then plot_title = title else begin
           if keyword_set(plot_weights) then plot_title = 'Weights' else plot_title = textoidl('P_k', font = font)
        endelse
     endelse
  endif

  xtitle = textoidl('k_{perp} (Mpc^{-1})', font = font)
  xtitle = repstr(xtitle, 'perp', perp_char)
  ytitle = textoidl('k_{||} (Mpc^{-1})', font = font)

  if keyword_set(no_title) or keyword_set(baseline_axis) then initial_title = '' else initial_title = plot_title

  ;;cgplot, plot_kperp, plot_kpar, /xlog, /ylog, /nodata, xstyle=5, ystyle=5, title = initial_title, $
  ;;        position = plot_pos, thick = thick, charthick = charthick, xthick = xthick, ythick = ythick, charsize = charsize, $
  ;;        font = font, noerase = no_erase, axiscolor = annotate_color, background = background_color

  ;;cgimage, power_log_norm, /nointerp,/overplot,/noerase

  axkeywords = {xlog: 1, ylog: 1, xstyle: 5, ystyle: 5, thick: thick, charthick: charthick, xthick: xthick, ythick: ythick, $
                charsize: charsize, font: font} 
  cgimage, power_log_norm, /nointerp, xrange = minmax(plot_kperp), yrange = minmax(plot_kpar), $
           title=initial_title, position = plot_pos, noerase = no_erase, color = annotate_color, background = background_color, $
           axkeywords = axkeywords, /axes
           


  if keyword_set(plot_wedge_line) then $
     cgplot, /overplot, plot_kperp, plot_kperp * wedge_amp, color = annotate_color, thick = thick+1, psym=-0, linestyle = 2

  cgaxis, xaxis=0, xtick_get = xticks, xtitle = xtitle, xrange = minmax(plot_kperp), $
        charthick = charthick, xthick = xthick, ythick = ythick, charsize = charsize, font = font, $
        xtickformat = 'exponent', xstyle = 1, color = annotate_color

  cgaxis, yaxis=0, ytick_get = yticks, ytitle = ytitle, yrange = minmax(plot_kpar), $
        charthick = charthick, xthick = xthick, ythick = ythick, charsize = charsize, font = font, $
        ytickformat = 'exponent', ystyle = 1, color = annotate_color
  if keyword_set(baseline_axis) then begin
     cgaxis, xaxis=1, xrange = minmax(plot_kperp* kperp_lambda_conv), xtickformat = 'exponent', xthick = xthick, $
             charthick = charthick, ythick = ythick, charsize = charsize, font = font, xstyle = 1, color = annotate_color

     if not keyword_set(no_title) then cgtext, xloc_title, yloc_title, plot_title, /normal, alignment=0.5, charsize=1.2 * charsize, $
                                               color = annotate_color, font = font
     cgtext, xloc_lambda, yloc_lambda, textoidl('(\lambda)', font = font), /normal, alignment=0.5, charsize=charsize, $
             color = annotate_color, font = font
  endif else $
     cgaxis, xaxis=1, xrange = minmax(plot_kperp), xtickv = xticks, xtickname = replicate(' ', n_elements(xticks)), $
             charthick = charthick, xthick = xthick, ythick = ythick, charsize = charsize, font = font, xstyle = 1, $
             color = annotate_color

  cgaxis, yaxis=1, yrange = minmax(plot_kpar), ytickv = yticks, ytickname = replicate(' ', n_elements(yticks)), $
          charthick = charthick, xthick = xthick, ythick = ythick, charsize = charsize, font = font, ystyle = 1, color = annotate_color

  tick_vals = loglevels(10d^log_data_range, coarse=2)

  ;; want minor tick marks if there aren't very many loglevels.
  ;; unfortunately cgcolorbar can't do log minor tick marks
  ;; with specified tick locations (which I have to do b/c
  ;; can't use divisions=0 with formatting keyword)
  ;; solution: add regular tickmarks without labels for the minor ones.

  if n_elements(tick_vals) lt 5 then begin
     if n_elements(tick_vals) lt 2 then minor_multipliers = dindgen(8)+2 else minor_multipliers = (dindgen(4)+1)*2d
     n_minor_mult = n_elements(minor_multipliers)
     n_major = n_elements(tick_vals)

     minor_tick_vals = reform(rebin(minor_multipliers, n_minor_mult, n_major), n_minor_mult*n_major) $
                       * reform(rebin(reform(tick_vals, 1, n_major), n_minor_mult, n_major), n_minor_mult*n_major)
     wh_keep = where(minor_tick_vals gt 10^log_data_range[0] and minor_tick_vals lt 10^log_data_range[1], count_keep)

     if count_keep gt 0 then begin
        minor_tick_vals = minor_tick_vals[wh_keep]
        n_minor = n_elements(minor_tick_vals)
        minor_tick_names = strarr(n_minor) + ' '
     endif else n_minor = 0

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

     ;; if log_data_range[0] lt min_pos then begin
     ;;    wh = where(tick_vals lt min_pos, count, complement = wh_keep, ncomplement = count_keep) 
     ;;    if count_keep gt 0 then tick_vals = [min_pos, tick_vals[wh_keep]]

     ;;    names = ['<0', '10!U' + strtrim(string(round(alog10(tick_vals))), 2) + '!N']

     ;; endif else names = '10!U' + strtrim(string(round(alog10(tick_vals))), 2) + '!N'
     if min(power) lt 0 and min(tick_vals) lt min_pos then begin
        wh = where(tick_vals lt min_pos, count, complement = wh_keep, ncomplement = count_keep) 

        if count lt 1 then stop $
        else if count eq 1 then names = ['<0', number_formatter(tick_vals[wh_keep], format = '(e0)',/print_exp)] $
        else names = [strarr(count-1), '<0', number_formatter(tick_vals[wh_keep], format = '(e0)',/print_exp)]

     endif else names = number_formatter(tick_vals, format = '(e0)',/print_exp)
  endelse

  if n_minor gt 0 then begin
     temp_ticks = [tick_vals, minor_tick_vals]
     temp_names = [names, minor_tick_names]
     order = sort(temp_ticks)

     tick_vals = temp_ticks[order]
     names = temp_names[order]
  endif

  if (alog10(tick_vals[0]) - log_data_range[0]) gt 10^(-3d) then begin
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
              charsize = charsize, font = font
   ;; cgcolorbar, color = annotate_color, /vertical, position = cb_pos, bottom = color_range[0]+1, ncolors = n_colors, /ylog, $
   ;;             range = 10d^log_data_range, format = 'exponent', minor=minor, charsize = charsize, font = font, $
   ;;             divisions = 0

  if keyword_set(pub) and n_elements(multi_pos) eq 0 then begin
     psoff
     wdelete, window_num
  endif
  
  tvlct, r, g, b
  
 end
