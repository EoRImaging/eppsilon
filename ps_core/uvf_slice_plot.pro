pro uvf_slice_plot, slice_savefile, multi_pos = multi_pos, multi_aspect = multi_aspect, plot_xrange = plot_xrange, $
                    plot_yrange = plot_yrange, data_range = data_range, type = type, pub = pub, plotfile = plotfile, $
                    window_num = window_num, title = title, grey_scale = grey_scale, baseline_axis = baseline_axis, mark_0 = mark_0, $
                    image_space = image_space, color_0amp = color_0amp, $
                    charsize = charsize_in, cb_size = cb_size_in, margin=margin_in, cb_margin = cb_margin_in

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

  type_enum = ['abs', 'phase', 'real', 'imaginary']
  if n_elements(type) eq 0 then if keyword_set(image_space) then type = 'real' else type = 'phase'
  wh = where(type_enum eq type, count)
  if count eq 0 then message, 'unknown type. Use one of: ' + type_enum

  if keyword_set(image_space) and keyword_set(baseline_axis) then begin
     print, 'baseline_axis keyword cannot be used with image_space keyword'
     baseline_axis = 0
  endif

  restore, slice_savefile
  dims = size(uvf_slice, /dimension)
  if n_elements(dims) gt 2 then uvf_slice = total(uvf_slice, slice_axis+1) / n_elements(slice_inds)

  if max(abs(uvf_slice)) eq 0 then all_zero = 1

  xdelta = xarr[1] - xarr[0]
  ydelta = yarr[1] - yarr[0]

  if keyword_set(image_space) then begin
     old_xarr = xarr
     old_xdelta = xdelta

     xdelta = 1d/(kperp_lambda_conv * dims[0] * old_xdelta)
     xarr = (dindgen(dims[0])-dims[0]/2) * xdelta

     if slice_axis eq 2 then begin
        old_yarr = yarr
        old_ydelta = ydelta
        
        ydelta = 1d/(kperp_lambda_conv * dims[1] * old_ydelta)
        yarr = (dindgen(dims[1])-dims[1]/2) * ydelta

        image = shift(fft(shift(uvf_slice,dims[0]/2,dims[1]/2)),dims[0]/2,dims[1]/2)
     endif else image = shift(fft(shift(uvf_slice,dims[0]/2,0)),dims[0]/2,0)
     old_slice = uvf_slice
     uvf_slice = image

  endif

  xarr_edges = [xarr - xdelta/2, max(xarr) + xdelta/2]
  yarr_edges = [yarr - ydelta/2, max(yarr) + ydelta/2]

  if n_elements(plot_xrange) eq 0 then begin
     temp = where(total(abs(uvf_slice),2) gt 0, count)
     if count gt 1 then plot_xrange =  minmax(xarr_edges[[temp, max(temp)+1]]) $
     else plot_xrange = minmax(xarr_edges)
  endif

  if n_elements(plot_yrange) eq 0 then begin
     temp = where(total(abs(uvf_slice),1) gt 0, count)
     
     if count gt 1 then plot_yrange = minmax(yarr_edges[[temp, max(temp)+1]]) $
     else plot_yrange = minmax(yarr_edges)
  endif

  ;; if slice_axis eq 2 then begin
  ;;    plot_xrange = minmax([plot_xrange, plot_yrange])
  ;;    plot_yrange = plot_xrange
  ;; endif

  wh_x_inrange = where(xarr ge plot_xrange[0] and xarr + xdelta le plot_xrange[1], n_x_plot)
  wh_y_inrange = where(yarr ge plot_yrange[0] and yarr + ydelta le plot_yrange[1], n_y_plot)

  if n_elements(plotfile) eq 0 then $
      plotfile = base_path() + 'power_spectrum/plots/' + plane_name + ' plane.eps' $
   else if strcmp(strmid(plotfile, strlen(plotfile)-4), '.eps', /fold_case) eq 0 then plotfile = plotfile + '.eps'
    
  if n_x_plot eq 0 or n_y_plot eq 0 then message, 'No data in plot k range'

  if n_x_plot ne n_elements(xarr) then begin
     uvf_slice = uvf_slice[wh_x_inrange, *]
     xarr = xarr[wh_x_inrange]
  endif
  if n_y_plot ne n_elements(yarr) then begin
     uvf_slice = uvf_slice[*, wh_y_inrange]
     yarr = yarr[wh_y_inrange]
  endif

  case type of
     'abs': begin
        plot_slice = abs(uvf_slice)
        cb_title = 'Magnitude'
     end
     'phase': begin
        plot_slice = atan(uvf_slice, /phase)
        cb_title = 'Phase (degrees)'
     end
     'real': begin
        plot_slice = real_part(uvf_slice)
        cb_title = 'Real Part'
     end
     'imaginary':  begin
        plot_slice = imaginary(uvf_slice)
        cb_title = 'Imaginary Part'       
     end
 endcase

  xarr_edges = [xarr - xdelta/2, max(xarr) + xdelta/2]
  yarr_edges = [yarr - ydelta/2, max(yarr) + ydelta/2]

  tvlct, r, g, b, /get

  if keyword_set(grey_scale) then begin
     cgloadct, 0, /reverse
     background_color = 'white'
     annotate_color = 'black'
  endif else begin
     cgloadct, 25, /brewer, /reverse
     background_color = 'white'
     annotate_color = 'black'
  endelse
  n_colors = 256
  
  if n_elements(data_range) eq 0 then if not keyword_set(all_zero) then data_range = minmax(plot_slice) else data_range = [-1*!pi, !pi]

  slice_plot = congrid(plot_slice, n_x_plot*10, n_y_plot * 10)
     
  slice_plot_norm = (slice_plot-data_range[0])*n_colors/(data_range[1]-data_range[0])
 
  tvlct, r2, g2, b2, /get
  plot_dims = size(slice_plot_norm, /dimension)
  slice_plot_norm_rgb = bytarr(3, plot_dims[0], plot_dims[1])
  slice_plot_norm_rgb[0,*,*] = r2[slice_plot_norm]
  slice_plot_norm_rgb[1,*,*] = g2[slice_plot_norm]
  slice_plot_norm_rgb[2,*,*] = b2[slice_plot_norm]

  if keyword_set(color_0amp) then begin
     wh_0amp = where(abs(slice_plot) eq 0, count_0amp)
     if count_0amp ne 0 then begin

        for i=0, 2 do begin
           temp = slice_plot_norm_rgb[i,*,*]
           temp[wh_0amp] = 255
           slice_plot_norm_rgb[i,*,*] = temp
        endfor

     endif
  endif

  ;; Work out plot & colorbar positions
  ;; in units of plot area (incl. margins)
  if n_elements(cb_size_in) eq 0 then cb_size = 0.025 else cb_size = cb_size_in
  if n_elements(margin_in) lt 4 then begin
     margin = [0.15, 0.15, 0.02, 0.1]
     if keyword_set(baseline_axis) then if (n_elements(multi_pos) gt 0 and keyword_set(pub)) then $
        margin[3] = 0.3 else margin[3] = 0.15
     if keyword_set(baseline_axis) and slice_axis eq 2 then margin[2] = 0.1
  endif else margin = margin_in

  if n_elements(cb_margin_in) lt 2 then begin
     cb_margin = [0.08, 0.02] 
     if n_elements(multi_pos) gt 0 then cb_margin[0] = 0.18
  endif else cb_margin = cb_margin_in 

  plot_pos = [margin[0], margin[1], (1-cb_margin[1]-cb_size-cb_margin[0]-margin[2]), (1-margin[3])]
  cb_pos = [(1-cb_margin[1]-cb_size), margin[1], (1-cb_margin[1]), (1-margin[3])]

  plot_aspect = (plot_pos[3] - plot_pos[1]) / (plot_pos[2] - plot_pos[0])

  plot_xlength = max(xarr_edges) - min(xarr_edges)
  if slice_axis eq 2 then plot_ylength = (max(yarr_edges) - min(yarr_edges)) $
  else plot_ylength = plot_xlength


  data_aspect = (plot_ylength / plot_xlength)
  aspect_ratio =  data_aspect /plot_aspect
  
  if aspect_ratio gt 1 then begin
     y_factor = 1.
     x_factor = 1/aspect_ratio
  endif else begin

     y_factor = aspect_ratio
     x_factor = 1.
  endelse
  
  if n_elements(multi_pos) eq 4 then begin
     ;; work out positions scaled to the area allowed in multi_pos with proper aspect ratio
     multi_xlen = (multi_pos[2]-multi_pos[0])
     multi_ylen = (multi_pos[3]-multi_pos[1])
     multi_center = [multi_pos[0] + multi_xlen/2d, multi_pos[1] + multi_ylen/2d]
     
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
     if n_elements(multi_pos) gt 0 then yloc_title = plot_pos[3] + 0.8* (multi_pos[3]-plot_pos[3]) $
     else yloc_title = plot_pos[3] + 0.6* (1-plot_pos[3])
  endif

  ;; if n_elements(multi_pos) gt 0 then begin
  ;;    xloc_lambda = plot_pos[0] - 0.1* (plot_pos[0]-multi_pos[0])
  ;;    yloc_lambda = plot_pos[3] + 0.1* (multi_pos[3]-plot_pos[3])

  ;;    xloc2_lambda = plot_pos[2] + 0.1* (multi_pos[2]-plot_pos[2])
  ;;    yloc2_lambda = plot_pos[1] - 0.1* (plot_pos[1]-multi_pos[1])
  ;; endif else begin
  ;;    xloc_lambda = plot_pos[0] - 0.1* (plot_pos[0]-0)
  ;;    yloc_lambda = plot_pos[3] + 0.1* (1-plot_pos[3])

  ;;    xloc2_lambda = plot_pos[2] + 0.1* (1-plot_pos[2])
  ;;    yloc2_lambda = plot_pos[1] - 0.1* (plot_pos[1]-0)
  ;; endelse

  if n_elements(multi_pos) gt 0 then begin
     xloc_lambda = plot_pos[0] - 0.1* (plot_pos[0]-multi_pos[0])
     yloc_lambda = plot_pos[3] + 0.1* (multi_pos[3]-plot_pos[3])

     xloc2_lambda = plot_pos[2] + 0.1* (multi_pos[2]-plot_pos[2])
     yloc2_lambda = plot_pos[1] - 0.1* (plot_pos[1]-multi_pos[1])
  endif else begin
     xloc_lambda = (plot_pos[2] - plot_pos[0])/2d + plot_pos[0]
     yloc_lambda = plot_pos[3] + 0.3* (1-plot_pos[3])

     xloc2_lambda = plot_pos[2] + 0.3* (1-plot_pos[2])
     yloc2_lambda =  (plot_pos[3] - plot_pos[1])/2d + plot_pos[1]
  endelse

  if keyword_set(pub) then begin
     charthick = 3
     thick = 3
     xthick = 3
     ythick = 3
      if n_elements(charsize_in) eq 0 then begin
        if n_elements(multi_pos) gt 0 then begin
           min_len = min([multi_xlen, multi_ylen])
           charsize = 5d * min_len
        endif else charsize = 2
     endif else charsize = charsize_in
     font = 1
     
     if n_elements(multi_pos) eq 0 then begin
        window, window_num, xsize = xsize, ysize = ysize
        pson, file = plotfile, /eps 
     endif
  
  endif else begin
     thick = 1
     if n_elements(charsize_in) eq 0 then charsize=1 else charsize = charsize_in
     if n_elements(multi_pos) eq 0 then begin
        if windowavailable(window_num) then begin 
           wset, window_num
           if !d.x_size ne xsize or !d.y_size ne ysize then make_win = 1 else make_win = 0
        endif else make_win = 1
        
        if make_win eq 1 then window, window_num, xsize = xsize, ysize = ysize
     endif
  endelse

 
  if n_elements(title) ne 0 then plot_title = title $
  else begin
     case type of
        'abs': plot_title = 'Magnitude in ' + plane_name + ' plane'
        'phase': plot_title = 'Phase in ' + plane_name + ' plane'
        'real': plot_title = 'Real part of ' + plane_name + ' plane'
        'imaginary': plot_title = 'Imaginary part of ' + plane_name + ' plane'
     endcase
  endelse

  if keyword_set(image_space) then begin
     case slice_axis of 
        0:begin
           plot_xtitle = 'theta y (radians)'
           plot_ytitle = plot_yname + ' (MHz)
        end
        1: begin
           plot_xtitle = 'theta x (radians)'
           plot_ytitle = plot_yname + ' (MHz)
        end
        2: begin
           plot_xtitle = 'theta x (radians)'
           plot_ytitle = 'theta y (radians)'
        end
     endcase
  endif else begin
     plot_xtitle = plot_xname + ' (Mpc!U-1!N)'
     if slice_axis lt 2 then plot_ytitle = plot_yname + ' (MHz)' $
     else plot_ytitle = plot_yname + ' (Mpc!U-1!N)'
  endelse

  plot_xarr = xarr_edges
  plot_yarr = yarr_edges

  if keyword_set(baseline_axis) then initial_title = '' else initial_title = plot_title

  cgplot, plot_xarr, plot_yarr, /nodata, xlog=xlog, ylog=ylog, xstyle=5, ystyle=5, title = initial_title, position = plot_pos, $
        xrange = minmax(plot_xarr), yrange = minmax(plot_yarr), thick = thick, charthick = charthick, xthick = xthick, $
          ythick = ythick, charsize = charsize, font = font, noerase = no_erase, background = background_color, $
          axiscolor = annotate_color
  cgimage, slice_plot_norm_rgb, /nointerp,/overplot,/noerase

  if keyword_set(mark_0) then begin
     cgplot, /overplot, plot_yarr * 0d, plot_yarr, color=annotate_color, psym=-0, thick = thick
     if slice_axis eq 2 then  cgplot, /overplot, plot_xarr, plot_xarr* 0d, color=annotate_color, psym=-3, thick = thick
  endif


  cgaxis, xaxis=0, xtick_get = xticks, xtitle = plot_xtitle, xrange = minmax(plot_xarr), $
          charthick = charthick, xthick = xthick, ythick = ythick, charsize = charsize, font = font, $
          xstyle = 1, color = annotate_color
  
  cgaxis, yaxis=0, ytick_get = yticks, ytitle = plot_ytitle, yrange = minmax(plot_yarr), $
          charthick = charthick, xthick = xthick, ythick = ythick, charsize = charsize, font = font, ystyle = 1, color = annotate_color
 
  if keyword_set(baseline_axis) then begin
     cgaxis, xaxis=1, xrange = minmax(plot_xarr * kperp_lambda_conv), xthick = xthick, $
             charthick = charthick, ythick = ythick, charsize = charsize, font = font, xstyle = 1, color = annotate_color

    if not keyword_set(no_title) then cgtext, xloc_title, yloc_title, plot_title, /normal, alignment=0.5, charsize=1.2 * charsize, $
                                               color = annotate_color, font = font
     cgtext, xloc_lambda, yloc_lambda, textoidl('(\lambda)', font = font), /normal, alignment=0.5, charsize=charsize, $
             color = annotate_color, font = font

  endif else $
     cgaxis, xaxis=1, xrange = minmax(plot_xarr), xtickv = xticks, xtickname = replicate(' ', n_elements(xticks)), $
             charthick = charthick, xthick = xthick, ythick = ythick, charsize = charsize, font = font, xstyle = 1, $
             color = annotate_color

  if keyword_set(baseline_axis) and slice_axis eq 2 then begin
     cgaxis, yaxis=1, yrange = minmax(plot_yarr * kperp_lambda_conv), xthick = xthick, $
             charthick = charthick, ythick = ythick, charsize = charsize, font = font, ystyle = 1, color = annotate_color

      cgtext, xloc2_lambda, yloc2_lambda, textoidl('(\lambda)', font = font), /normal, alignment=0.5, charsize=charsize, $
              color = annotate_color, font = font
 endif else $
     cgaxis, yaxis=1, yrange = minmax(plot_yarr), ytickv = yticks, ytickname = replicate(' ', n_elements(yticks)), $
             charthick = charthick, xthick = xthick, ythick = ythick, charsize = charsize, font = font, ystyle = 1, $
             color = annotate_color

  if type eq 'phase' then cb_range = data_range * 180/!dpi else cb_range = data_range


  cgcolorbar, color = annotate_color, /vertical, position = cb_pos, charsize = charsize, font = font, minrange = cb_range[0], $
              maxrange = cb_range[1], title = cb_title, minor=5

  if keyword_set(pub) and n_elements(multi_pos) eq 0 then begin
     psoff
     wdelete, window_num
  endif
  
  tvlct, r, g, b

end
