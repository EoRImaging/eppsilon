pro fhd_data_plots, datafile, plot_path = plot_path, healpix=healpix, refresh_dft = refresh_dft, refresh_ps = refresh_ps, $
                    refresh_binning = refresh_binning, $
                    no_weighting = no_weighting, std_power = std_power, no_kzero = no_kzero, $
                    no_weighted_averaging = no_weighted_averaging, $
                    data_range = data_range, plot_weights = plot_weights, slice_nobin = slice_nobin, linear_kpar = linear_kpar, $
                    plot_uvf = plot_uvf, uvf_data_range = uvf_data_range, uvf_type = uvf_type, baseline_axis = baseline_axis, $
                    plot_wedge_line = plot_wedge_line, grey_scale = grey_scale, pub = pub

  ;; default to absolute value for uvf plots
  if keyword_set(plot_uvf) and n_elements(uvf_type) eq 0 then uvf_type = 'abs'

  ;; default to including baseline axis
  if n_elements(baseline_axis) eq 0 then baseline_axis = 1
  ;; default to linear kparallel bins
  if n_elements(linear_kpar) eq 0 then linear_kpar = 1

  ;; default to plot wedge line
  if n_elements(plot_wedge_line) eq 0 then plot_wedge_line=1

  if keyword_set(healpix) and keyword_set(refresh_dft) then refresh_ps = 1

  if n_elements(no_weighted_averaging) gt 0 and not keyword_set(no_weighted_averaging) and keyword_set(no_weighting) then $
     message, 'Cannot set no_weighted_averaging=0 and no_weighting=1'

  ftest = file_test(datafile)
  if datafile_test eq 0 then message, 'datafile not found'

  temp = strpos(datafile, '/', /reverse_search)
  froot = strmid(datafile, 0, temp+1)
  infilebase = strmid(datafile, temp+1)
  temp2 = strpos(infilebase, '.', /reverse_search)
  datafilebase = strmid(infilebase, 0, temp2)

  ;; froot = base_path('data') + 'fhd_ps_data/'
  ;; if keyword_set(healpix) then datafilebase = 'multi_freq_residuals_cube_healpix' else datafilebase ='multi_freq_residuals_cube_20k'
  ;; datafile = froot + datafilebase + '.sav'
  
  file_obj = obj_new('idl_savefile', datafile)
  varnames = file_obj->names()
  wh_obs = where(strlowcase(varnames) eq 'obs', count_obs)
  if count_obs ne 0 then file_obj->restore, 'obs' $
  else begin
     wh_obs = where(strlowcase(varnames) eq 'obs_arr', count_obs)
     if count_obs ne 0 then begin
        file_obj->restore, 'obs_arr'
        n_obs = n_elements(obs_arr)
        max_baseline_vals = dblarr(n_obs)
        for i=0, n_obs-1 do begin
           if abs((*obs_arr[i]).degpix - (*obs_arr[0]).degpix) gt 0 then message, 'inconsistent degpix values in obs_arr'
           if total(abs((*obs_arr[i]).fbin_i - (*obs_arr[0]).fbin_i)) gt 0 then message, 'inconsistent fbin_i values in obs_arr'
           if total(abs((*obs_arr[i]).freq - (*obs_arr[0]).freq)) gt 0 then message, 'inconsistent freq values in obs_arr'
           max_baseline_vals[i] = (*obs_arr[i]).max_baseline
        endfor

        wh_max = where(max_baseline_vals eq max(max_baseline_vals))
        obs = *obs_arr[wh_max[0]]
     endif else message, 'no obs or obs_arr in datafile'
  endelse
  if keyword_set(healpix) then file_obj->restore, 'nside'
  obj_destroy, file_obj
  
  n_freqbins = max(obs.fbin_i) + 1
  frequencies = dblarr(n_freqbins)
  for i=0, n_freqbins-1 do begin
     wh = where(obs.fbin_i eq i, count)
     if count eq 0 then stop
     
     frequencies[i] = mean(obs.freq[wh]) / 1e6 ;; in MHz
  endfor
  
  obs_tags = tag_names(obs)
  wh_baseline = where(strlowcase(obs_tags) eq 'max_baseline', count_whbaseline)
  ;; the max baseline in the obs structure is given in wavelengths, need to convert using the maximum frequency
  if count_whbaseline ne 0 then max_baseline = 3e8/max(obs.freq)*obs.max_baseline else max_baseline = 342.497
  
  fadd = ''
  if keyword_set(std_power) then fadd = fadd + '_sp'
  if keyword_set(no_weighting) then fadd = fadd + '_nowt'
  
  fadd_2dbin = ''
  wt_fadd_2dbin = ''
  ;;if keyword_set(fill_holes) then fadd_2dbin = fadd_2dbin + '_nohole'
  if keyword_set(no_weighted_averaging) then fadd_2dbin = fadd_2dbin + '_nowtave'
  if keyword_set(no_kzero) then begin
     fadd_2dbin = fadd_2dbin + '_nok0'
     wt_fadd_2dbin = wt_fadd_2dbin + '_nok0'
  endif
  if keyword_set(linear_kpar) then begin
     fadd_2dbin = fadd_2dbin + '_linkpar'
     wt_fadd_2dbin = wt_fadd_2dbin + '_linkpar'
  endif

  data_varnames = ['DIRTY_XX_CUBE', 'RES_XX_CUBE', 'MODEL_XX_CUBE', 'DIRTY_YY_CUBE', 'RES_YY_CUBE', 'MODEL_YY_CUBE']
  weight_varnames = ['WEIGHTS_XX_CUBE', 'WEIGHTS_YY_CUBE']
  weight_labels = ['XX', 'YY']
  weight_ind = [intarr(3), intarr(3)+1]
  wt_file_labels = '_weights_' + strlowcase(weight_labels[weight_ind])
  file_labels = ['_dirty_xx', '_res_xx', '_model_xx', '_dirty_yy', '_res_yy', '_model_yy']
  titles = ['dirty XX', 'residual XX', 'model XX', 'dirty YY', 'residual YY', 'model YY']

  ;; data_varnames = ['DIRTY_XX_CUBE', 'RES_XX_CUBE', 'MODEL_XX_CUBE']
  ;; weight_varnames = ['WEIGHTS_XX_CUBE']
  ;; weight_labels = ['XX']
  ;; weight_ind = [intarr(3)]
  ;; wt_file_labels = '_weights_' + strlowcase(weight_labels[weight_ind])
  ;; file_labels = ['_dirty_xx', '_res_xx', '_model_xx']
  ;; titles = ['dirty XX', 'residual XX', 'model XX']


  savefilebase = froot + datafilebase + file_labels + fadd
  savefiles_2d = savefilebase + fadd_2dbin + '_2dkpower.idlsave'
  test_save = file_test(savefiles_2d) *  (1 - file_test(savefiles_2d, /zero_length))

  nfiles = n_elements(savefiles_2d)
  for i=0, nfiles-1 do begin
     if test_save[i] eq 0 or keyword_set(refresh_ps) or keyword_set(refresh_binning) then begin
        if keyword_set(healpix) then begin
           pixelfile = datafile
           pixel_varname = 'hpx_inds'
           hpx_dftsetup_savefile = froot + 'multi_freq_residuals_cube_healpix_dftsetup.idlsave'
           weight_savefilebase = froot + datafilebase + wt_file_labels + fadd
           
           weight_refresh = intarr(nfiles)
           if keyword_set(refresh_dft) then begin
              temp = sort(weight_ind)
              weight_refresh[temp[where(weight_ind[temp]-shift(weight_ind[temp],1) ne 0)]] = 1
           endif

           fhd_3dps, datafile, data_varnames[i], datafile, weight_varnames[weight_ind[i]], frequencies, max_baseline, /healpix, $
                     nside = nside, pixelfile = pixelfile, pixelvar = pixel_varname, hpx_dftsetup_savefile = hpx_dftsetup_savefile, $
                     savefilebase = savefilebase[i], weight_savefilebase = weight_savefilebase[i], refresh = refresh_ps, $
                     dft_refresh_data=refresh_dft, dft_refresh_weight=weight_refresh[i], $
                     no_weighting = no_weighting, std_power = std_power, no_kzero = no_kzero, $
                     no_weighted_averaging = no_weighted_averaging, /quiet
        endif else $
           fhd_3dps, datafile, data_varnames[i], datafile, weight_varnames[weight_ind[i]], frequencies, max_baseline, $
                     degpix=obs.degpix, savefilebase = savefilebase[i], refresh = refresh_ps, no_weighting = no_weighting, $
                     std_power = std_power, no_kzero = no_kzero, no_weighted_averaging = no_weighted_averaging, /quiet
     endif
  endfor


  restore, savefiles_2d[0]
  wh_good_kperp = where(total(power, 2) gt 0, count)
  if count eq 0 then stop
  kperp_plot_range = [min(kperp_edges[wh_good_kperp]), max(kperp_edges[wh_good_kperp+1])]
  kperp_plot_range = [6e-3, min([max(kperp_edges[wh_good_kperp+1]),1.5e-1])]

  plot_folder = 'fhd_data/'
  plotfile_path = base_path('plots') + 'power_spectrum/' + plot_folder
  plotfile_base = plotfile_path + datafilebase + fadd

  plot_fadd = ''
  if keyword_set(grey_scale) then plot_fadd = plot_fadd + '_grey'

  plotfile_2d = plotfile_base + fadd_2dbin + '_kspace_power' + plot_fadd + '.eps'

  if not keyword_set(slice_nobin) then slice_fadd = '_binned' else slice_fadd = ''
  yslice_plotfile = plotfile_base + '_xz_plane' + plot_fadd + slice_fadd + '.eps'
  xslice_plotfile = plotfile_base + '_yz_plane' + plot_fadd + slice_fadd + '.eps'
  zslice_plotfile = plotfile_base + '_xy_plane' + plot_fadd + slice_fadd + '.eps'

  if keyword_set(plot_wedge_line) then begin
     z0_freq = 1420.40 ;; MHz
     redshifts = z0_freq/frequencies - 1
     mean_redshift = mean(redshifts)

     cosmology_measures, mean_redshift, wedge_factor = wedge_factor
     source_dist = 20d * !dpi / 180d
     wedge_amp = wedge_factor * source_dist
  endif else wedge_amp = 0d

                                ;if n_elements(data_range) eq 0 then data_range = [1e13, 1e20]

  if nfiles eq 1 then begin
     kpower_2d_plots, savefiles_2d, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
                      data_range = data_range, pub = pub, plotfile = plotfile_2d, title = titles, grey_scale = grey_scale, $
                      plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp, baseline_axis = baseline_axis
  endif else begin

     nrow = ceil(sqrt(nfiles))
     ncol = ceil(nfiles / double(nrow))
     positions = fltarr(4, ncol*nrow)
     
     row_val = reverse(reform(rebin(indgen(nrow), nrow, ncol), ncol*nrow))
     col_val = reform(rebin(reform(indgen(ncol), 1, ncol), nrow, ncol), ncol*nrow)
     
     
     if keyword_set(pub) then begin
        xmargin = 0.025
        ymargin = 0.025
     endif else begin
        xmargin = 0.0125
        ymargin = 0.0125
     endelse
     
     positions[0,*] = col_val/double(ncol)+xmargin
     positions[1,*] = row_val/double(nrow)+ymargin
     positions[2,*] = (col_val+1)/double(ncol)-xmargin
     positions[3,*] = (row_val+1)/double(nrow)-ymargin
     
     max_ysize = 1200
     max_xsize = 1800
     if nfiles gt 9 then begin
        multi_aspect =0.9 
        xsize = round((max_ysize/nrow) * ncol/multi_aspect)
        ysize = max_ysize
     endif else begin
        if keyword_set(baseline_axis) then multi_aspect = 0.75 else multi_aspect =0.5
        xsize = round((max_ysize/nrow) * ncol/multi_aspect)
        ysize = max_ysize

        if xsize gt max_xsize then begin
           factor = double(max_xsize) / xsize
           xsize = max_xsize
           ysize = round(ysize * factor)
        endif

     endelse
     window_num = 1
     if windowavailable(window_num) then begin 
        wset, window_num
        if !d.x_size ne xsize or !d.y_size ne ysize then make_win = 1 else make_win = 0
     endif else make_win = 1
     if make_win eq 1 then window, window_num, xsize = xsize, ysize = ysize
     erase
     
     if keyword_set(pub) then begin
        pson, file = plotfile_2d, /eps
     endif
     
     for i=0, nfiles-1 do begin        
        kpower_2d_plots, savefiles_2d[i], multi_pos = positions[*,i], multi_aspect = multi_aspect, $
                         kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, data_range = data_range, pub = pub, $
                         title = titles[i], grey_scale = grey_scale, plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp, $
                         baseline_axis = baseline_axis
     endfor
     
     if keyword_set(pub) then begin
        psoff
        wdelete, window_num
     endif

  endelse

  weight_plotfile_2d = plotfile_path + 'weights/' + datafilebase + '_weights' + fadd + wt_fadd_2dbin + '_kspace' + $
                       plot_fadd + '.eps'
  nweights = n_elements(weight_labels)

  if nweights eq 1 then begin
     kpower_2d_plots, savefiles_2d[0], /plot_weights, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
                      pub = pub, plotfile = weight_plotfile_2d, title = 'Weights ' + weight_labels, window_num = 2, $
                      grey_scale = grey_scale, plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp, baseline_axis = baseline_axis
  endif else begin
     
     nrow = ceil(sqrt(nweights))
     ncol = ceil(nweights / double(nrow))
     positions = fltarr(4, ncol*nrow)

     col_val = reverse(reform(rebin(indgen(ncol), ncol, nrow), ncol*nrow))
     row_val = reform(rebin(reform(indgen(nrow), 1, nrow), ncol, nrow), ncol*nrow)
     
     
     if keyword_set(pub) then begin
        xmargin = 0.025
        ymargin = 0.025
     endif else begin
        xmargin = 0.0125
        ymargin = 0.0125
     endelse
     
     positions[0,*] = col_val/double(ncol)+xmargin
     positions[1,*] = row_val/double(nrow)+ymargin
     positions[2,*] = (col_val+1)/double(ncol)-xmargin
     positions[3,*] = (row_val+1)/double(nrow)-ymargin
     
     max_ysize = 1200
     max_xsize = 1800
     if nfiles gt 9 then begin
        multi_aspect =0.9 
        xsize = round((max_ysize/nrow) * ncol/multi_aspect)
        ysize = max_ysize
     endif else begin
        if keyword_set(baseline_axis) then multi_aspect = 0.75 else multi_aspect =0.5
        xsize = round((max_ysize/nrow) * ncol/multi_aspect)
        ysize = max_ysize

        if xsize gt max_xsize then begin
           factor = double(max_xsize) / xsize
           xsize = max_xsize
           ysize = round(ysize * factor)
        endif

     endelse
     window_num = 2
     if windowavailable(window_num) then begin 
        wset, window_num
        if !d.x_size ne xsize or !d.y_size ne ysize then make_win = 1 else make_win = 0
     endif else make_win = 1
     if make_win eq 1 then window, window_num, xsize = xsize, ysize = ysize
     erase
     
     if keyword_set(pub) then begin
        pson, file = weight_plotfile_2d, /eps
     endif
     
     for i=0, nweights-1 do begin
        wh = where(weight_ind eq i, count)
        if count eq 0 then stop
        kpower_2d_plots, savefiles_2d[wh[i]], /plot_weights, multi_pos = positions[*,i], multi_aspect = multi_aspect, $
                         kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, pub = pub, $
                         title = 'Weights ' + weight_labels[i], grey_scale = grey_scale, plot_wedge_line = plot_wedge_line, $
                         wedge_amp = wedge_amp, baseline_axis = baseline_axis
     endfor
     
     if keyword_set(pub) then begin
        psoff
        wdelete, window_num
     endif
  endelse



end
