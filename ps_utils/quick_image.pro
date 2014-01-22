pro quick_image, image, xvals, yvals, data_range = data_range, xrange = xrange, yrange = yrange, log=log, xtitle = xtitle, ytitle = ytitle, title = title, $
    grey_scale = grey_scale, xlog = xlog, ylog = ylog, missing_value = missing_value, noerase = noerase, savefile = savefile, png = png, eps = eps
    
  if n_elements(savefile) gt 0 or keyword_set(png) or keyword_set(eps) then pub = 1 else pub = 0
  if pub eq 1 then begin
    if not (keyword_set(png) or keyword_set(eps)) then begin
      basename = cgRootName(savefile, directory=directory, extension=extension)
      
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
    if n_elements(savefile) eq 0 then begin
      if keyword_set(eps) then savefile = 'idl_quick_image.eps' else savefile = 'idl_quick_image'
      cd, current = current_dir
      print, 'no filename specified for quick_image output. Using ' + current_dir + path_sep() + savefile
    endif
    
    if keyword_set(png) and keyword_set(eps) then begin
      print, 'both eps and png cannot be set, using png'
      eps = 0
    endif
    
    if keyword_set(eps) then delete_ps = 0 else delete_ps = 1
  endif
  
  ;; precaution in case a slice is passed in but it still appears > 2d (ie shallow dimension)
  image = reform(image)
  if n_elements(size(image, /dim)) ne 2 then begin
    print, 'image must be 2 dimensional'
    return
  endif
  
  if max(abs(imaginary(image))) gt 0 then begin
    print, 'image is complex, showing real part'
    image = real_part(image)
  endif
  
  if keyword_set(missing_value) then begin
    good_locs = where(image ne missing_value)
    erase = 1
  endif else good_locs = indgen(n_elements(image))
  
  if keyword_set(log) then begin
    wh_low = where(image le 0, count_low, complement = wh_good, ncomplement = count_good)
    if count_low eq 0 then plot_image = alog10(image) else begin
      if count_good eq 0 then begin
        print, 'Entire image is 0 or negative -- log scaling will not work.'
        return
      endif else print, 'Warning: part of the image is 0 or negative, setting those values to a small value for log scaling'
      
      min_good = min(image[wh_good])
      
      plot_image = image*0.
      plot_image[wh_good] = alog10(image[wh_good])
      plot_image[wh_low] = alog10(min_good/10.)
    endelse
    
    if n_elements(data_range) eq 0 then data_range = 10^minmax(plot_image[good_locs])
    img_range = alog10(data_range)
    
    cb_log = 1
    ticks = loglevels(data_range)
    while n_elements(ticks) lt 2 do begin
      if n_elements(log_binsize) eq 0 then exp = floor(alog10(data_range[0])) else exp = exp-1
      
      nbins = ceil((data_range[1] - data_range[0]) / 10.^exp)+2
      if nbins gt 1 then begin
        temp = (findgen(nbins)+floor(data_range[0]/10.^exp))*10.^exp
        wh_inrange = where(temp ge data_range[0] and temp le data_range[1], count_inrange)
        if count_inrange ne 0 then ticks = temp[wh_inrange]
      endif
    endwhile
    divisions = n_elements(ticks)-1
    ticknames = number_formatter(ticks, format='(e13.2)', /print_exp)
    minor=0
  endif else begin
    plot_image = image
    if n_elements(data_range) eq 0 then data_range = minmax(plot_image[good_locs])
    img_range = data_range
    
    tickinterval = float(number_formatter((data_range[1]-data_range[0])/6., format = '(e13.0)'))
  endelse
  
  if n_elements(xvals) gt 1 then xrange = minmax(xvals)
  if n_elements(yvals) gt 1 then yrange = minmax(yvals)
  
  tvlct, r, g, b, /get
  if keyword_set(grey_scale) then cgloadct, 0, /reverse else cgloadct, 25, /brewer, /reverse
  
  if n_elements(xlog) ne 0 then axkeywords = create_struct('xlog', 1, 'xtickformat', 'exponent')
  if n_elements(ylog) ne 0 then $
    if n_elements(axkeywords) ne 0 then axkeywords = create_struct(axkeywords, 'ylog', 1, 'ytickformat', 'exponent') $
  else axkeywords = create_struct('ylog', 1, 'ytickformat', 'exponent')
  
  if n_elements(missing_value) gt 0 and not keyword_set(noerase) then cgerase
  
  if keyword_set(pub) then begin
    if n_elements(missing_value) gt 0 then begin
      alphabackgroundimage = cgsnapshot()
    endif
    cgps_open, savefile, /font, encapsulated=eps
  endif
  
  cgimage, plot_image, maxvalue = img_range[1], minvalue = img_range[0], position = [.15,.1,.8,.95], /axes, xrange = xrange, $
    yrange = yrange, xtitle = xtitle, ytitle = ytitle, title = title, axkeywords = axkeywords, missing_value = missing_value, noerase = noerase, $
    alphabackgroundimage = alphabackgroundimage
    
  cgcolorbar, range=data_range, position = [.92, .1,.95,.95], /vertical, ylog = cb_log, minor = minor, ticknames = ticknames, $
    divisions=divisions, ytickv = ticks, tickinterval=tickinterval, format='exponent'
    
    
  if keyword_set(pub) then cgps_close, png = png, delete_ps = delete_ps
  
  tvlct, r, g, b
  
  
end
