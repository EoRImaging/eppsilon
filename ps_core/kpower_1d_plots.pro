

pro kpower_1d_plots, power_savefile, plot_weights = plot_weights, multi_pos = multi_pos, data_range = data_range, k_range = k_range, $
                     pub = pub, plotfile = plotfile, window_num = window_num, colors = colors, names = names, save_text = save_text, $
                     delta = delta

  if n_elements(window_num) eq 0 then window_num = 2

  nfiles = n_elements(power_savefile)
  if n_elements(names) gt 0 and n_elements(names) ne nfiles then message, 'Number of names does not match number of files'

  if n_elements(multi_pos) gt 0 then begin
     if n_elements(multi_pos) ne 4 then message, 'multi_pos must be a 4 element plot position vector'
     if max(multi_pos) gt 1 or min(multi_pos) lt 0 then message, 'multi_pos must be in normalized coordinates (between 0 & 1)'
     if multi_pos[2] le multi_pos[0] or multi_pos[3] le multi_pos[1] then $
        message, 'In multi_pos, x1 must be greater than x0 and y1 must be greater than y0 '

     no_erase = 1

     xlen = (multi_pos[2]-multi_pos[0])
     ylen = (multi_pos[3]-multi_pos[1])

     big_margin = [0.15, 0.2]      ;; in units of plot area (incl. margins)
     small_margin = [0.05, 0.1]    ;; in units of plot area (incl. margins)

     plot_pos = [big_margin[0], big_margin[1], (1-small_margin[0]), (1-small_margin[1])]

     new_pos = [xlen * plot_pos[0] + multi_pos[0], ylen * plot_pos[1] + multi_pos[1], $
                xlen * plot_pos[2] + multi_pos[0], ylen * plot_pos[3] + multi_pos[1]]

  endif else no_erase = 0

  if n_elements(plotfile) eq 0 then plotfile = strsplit(power_savefile[0], '.idlsave', /regex, /extract) + '_1dkplot.eps' $
  else if strcmp(strmid(plotfile, strlen(plotfile)-4), '.eps', /fold_case) eq 0 then plotfile = plotfile + '.eps'

  
  if n_elements(colors) eq 0 then colors = indgen(nfiles)*254/(nfiles-1)

  if keyword_set(save_text) then begin
     text_filename = strsplit(plotfile, '.eps', /regex, /extract) + '.txt'
     if nfiles gt 1 then if n_elements(names) ne 0 then text_labels = names else text_labels = strarr(nfiles)
        
     openw, lun, text_filename, /get_lun
  endif

  for i=0, nfiles-1 do begin
     restore, power_savefile[i]
     n_k = n_elements(power)

     if keyword_set(plot_weights) then begin
        if n_elements(weights) ne 0 then power = weights $
        else message, 'No weights array included in this file'
     endif
     
     if n_elements(k_centers) ne 0 then k_mid = k_centers $
     else k_mid = 10^(alog10(k_edges[1:*]) - (1d/bins_per_decade)/2.)
     theory_delta = (power * k_mid^3d / (2d*!pi^2d)) ^(1/2d)

     if keyword_set(save_text) then begin
        printf, lun,  text_labels[i]+ ' k'
        printf, lun, transpose(k_mid)
        printf, lun, ''
        if keyword_set(delta) then begin 
           printf, lun,  text_labels[i]+ ' delta (sqrt(k^3 Pk/(2pi^2)) -- mk)'
           printf, lun, transpose(theory_delta)
        endif else begin
           printf, lun,  text_labels[i]+ ' power (mk^2 Mpc^3)'
           printf, lun, transpose(power)
        endelse
        printf, lun, ''
     endif

     if keyword_set(delta) then power = theory_delta

     wh = where(power eq 0, count, complement = wh_good, ncomplement = count_good)
     if count_good eq 0 then message, 'No non-zero power'
     if count gt 0 then begin
        power = power[wh_good]
        k_mid = k_mid[wh_good]
     endif

     tag = 'f' + strsplit(string(i),/extract)
     if i eq 0 then begin
        power_plot = create_struct(tag, power)
        k_plot = create_struct(tag, k_mid)
  
        if n_elements(data_range) eq 0 then yrange = minmax(power) else yrange = data_range
        if n_elements(k_range) eq 0 then xrange = minmax(k_edges) else xrange = k_range

     endif else begin
        if n_elements(data_range) eq 0 then yrange = minmax([yrange, power])
        if n_elements(k_range) eq 0 then if n_elements(k_edges) gt 0 then xrange = minmax([xrange, k_edges])

        power_plot = create_struct(tag, power, power_plot)
        k_plot = create_struct(tag, k_mid, k_plot)
     endelse

     undefine, power
     if n_elements(k_edges) ne 0 then undefine, k_edges
     if n_elements(k_centers) ne 0 then undefine, k_centers
  endfor

  if keyword_set(save_text) then free_lun, lun

  tvlct, r, g, b, /get

  if keyword_set(pub) then begin
     charthick = 3
     thick = 3
     xthick = 3
     ythick = 3
     charsize = 2.5
     font = 1
     if nfiles gt 3 then legend_charsize = charsize / (nfiles/4.5d)  else legend_charsize = 2

    if n_elements(multi_pos) eq 0 then begin
       window, window_num
       pson, file = plotfile, /eps 
    endif
  endif else if n_elements(multi_pos) eq 0 then begin
     if windowavailable(window_num) then wset, window_num else window, window_num
  endif
  

  ;;plot, k_plot, power_plot, /ylog, /xlog, xrange = xrange, xstyle=1
  plot_order = sort(tag_names(power_plot))
  if keyword_set(delta) then ytitle = textoidl('(k^3 P_k /(2\pi^2))^{1/2} (mK)', font = font) $
  else ytitle = textoidl('P_k (mK^2 Mpc^3)', font = font)
  xtitle = textoidl('k (Mpc^{-1})', font = font)

  plot, k_plot.(plot_order[0]), power_plot.(plot_order[0]), position = new_pos, /ylog, /xlog, xrange = xrange, yrange = yrange, $
        xstyle=1, ystyle=1, xtitle = xtitle, ytitle = ytitle, psym=-3, xtickformat = 'exponent', ytickformat = 'exponent', $
        thick = thick, charthick = charthick, xthick = xthick, ythick = ythick, charsize = charsize, font = font, noerase = no_erase
  for i=0, nfiles - 1 do oplot, k_plot.(plot_order[i]), power_plot.(plot_order[i]), psym=-3, color = colors[i], thick = thick
  if n_elements(names) ne 0 then $
     al_legend, names, textcolor = colors, box = 0, /right,/bottom, charsize = legend_charsize, charthick = charthick


  if keyword_set(pub) and n_elements(multi_pos) eq 0 then begin
     psoff
     wdelete, window_num
  endif
  
  tvlct, r, g, b
  

end
