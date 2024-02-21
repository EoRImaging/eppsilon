pro binning_hist_plots, power, sim_noise, weights, mask_1dbin, power_1d, weights_1d, sim_noise_1d, $
    plot_log = plot_log, window_start = window_start, plotfilebase = plotfilebase, png = png, eps = eps, pdf = pdf

  if n_elements(plot_log) eq 0 then plot_log = 1

  if n_elements(plotfilebase) gt 0 or keyword_set(png) or keyword_set(eps) or keyword_set(pdf) then pub = 1 else pub = 0
  if pub eq 1 then begin
    if not (keyword_set(png) or keyword_set(eps) or keyword_set(pdf)) then begin
      basename = cgRootName(plotfilebase, directory=directory, extension=extension)

      case extension of
        'eps': eps=1
        'png': png=1
        'pdf': pdf=1
        '': png = 1
        else: begin
          print, 'Unrecognized extension, using png'
          png = 1
        end
      endcase

    endif

    if n_elements(plotfilebase) eq 0 and n_elements(multi_pos) eq 0 then begin
      plotfilebase = 'idl_binning_hist_plots'
      cd, current = current_dir
      print, 'no filename specified for binning_hist_plots output. Using ' + current_dir + path_sep() + plotfilebase
    endif

    if keyword_set(eps) then begin
      plot_exten = '.eps'
      delete_ps = 0
    endif else begin
      plot_exten = '.ps'
      delete_ps = 1
    endelse
  endif

  dims = size(power, /dimension)
  if n_elements(dims) eq 3 then begin
    binsize=0.25
    if keyword_set(plot_log) then yrange = [.1,1e4]
  endif else begin
    binsize=0.7
    if keyword_set(plot_log) then yrange = [.1,1e2]
  endelse

  noise_ratio = sim_noise*sqrt(weights)
  power_ratio = power*sqrt(weights)
  sigma_use = sqrt(1./weights)
  sigma_use[where(weights eq 0)] = 0

  nplots = max(mask_1dbin)
  npages = ceil(nplots/12.)
  nrow=[2,3]
  ncol = ceil(nplots/float(nrow*npages))
  layout_use = where(nrow*ncol - nplots eq min(nrow*ncol - nplots))
  nrow = nrow[layout_use[0]]
  ncol = ncol[layout_use[0]]

  if pub eq 1 then begin
    ps_aspect = nrow / float(ncol)
    if ps_aspect lt 1 then landscape = 1 else landscape = 0
    if (Keyword_Set(eps) and Keyword_Set(pdf) and landscape) then begin
      print, "Warning! eps forces plots to be portrait; run without eps to get landscape pdf plots."
    endif
    IF Keyword_Set(eps) THEN landscape = 0

    sizes = cgpswindow(LANDSCAPE=landscape, aspectRatio = ps_aspect, /sane_offsets)

    ; charthick = 3
    thick = 3
    xthick = 3
    ythick = 3
    symsize = 1

    temp = [sizes.xsize/ncol, sizes.ysize/nrow]
    charsize = 0.5d * mean(temp)

    legend_charsize = charsize/2.

    font = 1

    plotfile_power = plotfilebase + '_power'
    plotfile_ratio = plotfilebase + '_sigratio'
    if npages gt 1 then begin
      ;; multiple pages of plots, make separate files for each
      plotfile_power = plotfile_power + '_pg' + number_formatter(indgen(npages + 1)) + plot_exten
      plotfile_ratio = plotfile_ratio + '_pg' + number_formatter(indgen(npages + 1)) + plot_exten
    endif else begin
      plotfile_power += plot_exten
      plotfile_ratio += plot_exten
    endelse

  endif else begin
    charthick = 1
    thick = 1
    xthick = 1
    ythick = 1
    font = -1
    charsize = 2
    legend_charsize = charsize/2.
    symsize = 1
  endelse

  !p.multi=[0,ncol,nrow]
  if n_elements(window_start) then window_use = window_start-1 else window_use = 0
  page = 0

  for k=1, nplots do begin
    if (k-1) mod (nrow*ncol) eq 0 then begin
      window_use += 1
      page += 1

      if pub eq 0 then window,window_use, xsize=ncol*300, ysize=nrow*300 else begin
        if page gt 1 then begin
          if png then begin
            if pdf then delete_ps_use = 0 else delete_ps_use = delete_ps
            cgps_close, /png, delete_ps = delete_ps_use, density = 600
          endif
          if pdf then begin
            if not png then cgps_close
            cgps2pdf, plotfile_ratio[page-2], delete_ps=delete_ps
          endif
        endif

        cgps_open, plotfile_ratio[page-1], /font, encapsulated=eps, /nomatch, inches=sizes.inches, xsize=sizes.xsize, ysize=sizes.ysize, $
          xoffset=sizes.xoffset, yoffset=sizes.yoffset, landscape = landscape

      endelse
    endif

    wh_bin = where(mask_1dbin eq k, count_bin)
    if count_bin gt 0 then begin
      ref_dist = randomn(seed, count_bin)

      range = [-1,1]*max(abs([noise_ratio[wh_bin], power_ratio[wh_bin]]))
      quick_histplot, noise_ratio[wh_bin], binsize=binsize, range=range, $
        title = 'Bin #' + number_formatter(k) + ', npixels=' + number_formatter(count_bin), loghist = plot_log, yrange=yrange, $
        charthick = charthick, thick = thick, xthick = xthick, ythick = ythick, charsize = charsize, font = font, xrange = range, xstyle=1
      quick_histplot, power_ratio[wh_bin], binsize=binsize, range=range, /over, color='red', loghist = plot_log, thick = thick
      quick_histplot, ref_dist, binsize=binsize, range=range, /over, color='dark grey',loghist = plot_log, thick = thick/2.

      if keyword_set(plot_log) then value_yloc = 10^(!y.crange[1]*0.6) else value_yloc = !y.crange[1]*0.6

      if sim_noise_1d[k-1]*sqrt(weights_1d[k-1]) gt !x.crange[0] and sim_noise_1d[k-1]*sqrt(weights_1d[k-1]) lt !x.crange[1] then $
        cgplot, sim_noise_1d[k-1]*sqrt(weights_1d[k-1]), value_yloc, err_xhigh = 1, err_xlow = 1, psym=6, color='black', /over, $
        thick = thick/2., err_thick=thick/2., symsize = symsize/2. $
      else cgarrow, !x.crange[1]-1, value_yloc, !x.crange[1], value_yloc, /data, color='black', /solid, hsize = !d.x_size/200, thick = thick/2.
      if power_1d[k-1]*sqrt(weights_1d[k-1]) gt !x.crange[0] and power_1d[k-1]*sqrt(weights_1d[k-1]) lt !x.crange[1] then $
        cgplot, power_1d[k-1]*sqrt(weights_1d[k-1]), value_yloc, err_xhigh = 1, err_xlow = 1, psym=5, color='blue', /over, $
        thick = thick, err_thick=thick, symsize = symsize $
      else cgarrow, !x.crange[1]-1, value_yloc, !x.crange[1], value_yloc, /data, color='blue', /solid, hsize = !d.x_size/200, thick = thick
      al_legend, ['sim_noise/sigma', 'power/sigma', 'gaussian dist', '1d power/sigma'], textcolor=['black', 'red', 'dark grey', 'blue'], box=0, /right, $
        charthick = charthick, charsize = legend_charsize
    endif
  endfor
  if pub eq 1 then begin
    if png then begin
      if pdf then delete_ps_use = 0 else delete_ps_use = delete_ps
      cgps_close, /png, delete_ps = delete_ps_use, density = 600
    endif
    if pdf then begin
      if not png then cgps_close
      cgps2pdf, plotfile_ratio[page-1], delete_ps=delete_ps
    endif
  endif


  !p.multi=[0,ncol,nrow]
  page = 0

  for k=1, nplots do begin
    if (k-1) mod (nrow*ncol) eq 0 then begin
      window_use += 1
      page += 1

      if pub eq 0 then window,window_use, xsize=ncol*300, ysize=nrow*300 else begin
        if page gt 1 then begin
          if png then begin
            if pdf then delete_ps_use = 0 else delete_ps_use = delete_ps
            cgps_close, /png, delete_ps = delete_ps_use, density = 600
          endif
          if pdf then begin
            if not png then cgps_close
            cgps2pdf, plotfile_power[page-2], delete_ps=delete_ps
          endif
        endif

        cgps_open, plotfile_power[page-1], /font, encapsulated=eps, /nomatch, inches=sizes.inches, xsize=sizes.xsize, ysize=sizes.ysize, $
          xoffset=sizes.xoffset, yoffset=sizes.yoffset, landscape = landscape

      endelse
    endif

    wh_bin = where(mask_1dbin eq k, count_bin)
    if count_bin gt 0 then begin
      ref_dist = randomn(seed, count_bin) * sigma_use[wh_bin]
      binsize_use = binsize * mean(sigma_use[wh_bin])

      range = [-1,1]*max(abs([sim_noise[wh_bin], power[wh_bin]]))
      quick_histplot, sim_noise[wh_bin], binsize=binsize_use, range=range, $
        title = 'Bin #' + number_formatter(k) + ', npixels=' + number_formatter(count_bin), loghist = plot_log, yrange=yrange, $
        charthick = charthick, thick = thick, xthick = xthick, ythick = ythick, charsize = charsize, font = font, xrange = range, xstyle=1
      quick_histplot, power[wh_bin], binsize=binsize_use, range=range, /over, color='red', loghist = plot_log, thick = thick
      quick_histplot, ref_dist, binsize=binsize_use, range=range, /over, color='dark grey', loghist = plot_log, thick = thick/2.

      if keyword_set(plot_log) then value_yloc = 10^(!y.crange[1]*0.75) else value_yloc = !y.crange[1]*0.75

      if power_1d[k-1]*sqrt(weights_1d[k-1]) gt !x.crange[0] and power_1d[k-1]*sqrt(weights_1d[k-1]) lt !x.crange[1] then $
        cgplot, power_1d[k-1]*sqrt(weights_1d[k-1]), value_yloc, err_xhigh = sqrt(weights_1d[k-1]), $
        err_xlow = sqrt(weights_1d[k-1]), psym=5, color='blue', /over, thick = thick, err_thick = thick, symsize = symsize $
      else cgarrow, !x.crange[1]-sqrt(weights_1d[k-1]), value_yloc, !x.crange[1], value_yloc, /data, color='blue', /solid, hsize = !d.x_size/200, thick = thick
      al_legend, ['sim_noise', 'power', 'gaussian dist', '1d power'], textcolor=['black', 'red', 'dark grey', 'blue'], box=0, /right, $
        charthick = charthick, charsize = legend_charsize
    endif
  endfor
  if pub eq 1 then begin
    if png then begin
      if pdf then delete_ps_use = 0 else delete_ps_use = delete_ps
      cgps_close, /png, delete_ps = delete_ps_use, density = 600
    endif
    if pdf then begin
      if not png then cgps_close
      cgps2pdf, plotfile_ratio[page-1], delete_ps=delete_ps
    endif
  endif

  !p.multi=0

end
