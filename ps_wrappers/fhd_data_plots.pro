pro fhd_data_plots, datafile, save_path = save_path, plot_path = plot_path, healpix=healpix, pol_inc = pol_inc, type_inc = type_inc, $
                    refresh_dft = refresh_dft, dft_fchunk = dft_fchunk, refresh_ps = refresh_ps, refresh_binning = refresh_binning, $
                    no_weighting = no_weighting, std_power = std_power, no_kzero = no_kzero, $
                    no_weighted_averaging = no_weighted_averaging, data_range = data_range, plot_weights = plot_weights, $
                    slice_nobin = slice_nobin, log_kpar = log_kpar, log_kperp = log_kperp, kpar_bin = kpar_bin, $
                    kperp_bin = kperp_bin, log_k1d = log_k1d, k1d_bin = k1d_bin, plot_uvf = plot_uvf, $
                    uvf_data_range = uvf_data_range, uvf_type = uvf_type, baseline_axis = baseline_axis, $
                    plot_wedge_line = plot_wedge_line, grey_scale = grey_scale, pub = pub

  ;; default to absolute value for uvf plots
  if keyword_set(plot_uvf) and n_elements(uvf_type) eq 0 then uvf_type = 'abs'

  ;; default to including baseline axis
  if n_elements(baseline_axis) eq 0 then baseline_axis = 1

  ;; default to plot wedge line
  if n_elements(plot_wedge_line) eq 0 then plot_wedge_line=1

  if keyword_set(healpix) and keyword_set(refresh_dft) then refresh_ps = 1
  if keyword_set(refresh_ps) then refresh_binning = 1

  if n_elements(no_weighted_averaging) gt 0 and not keyword_set(no_weighted_averaging) and keyword_set(no_weighting) then $
     message, 'Cannot set no_weighted_averaging=0 and no_weighting=1'

  if n_elements(pol_inc) eq 0 then pol_inc = ['xx', 'yy']
  pol_enum = ['xx', 'yy']
  npol = n_elements(pol_inc)
  pol_num = intarr(npol)
  for i=0, npol-1 do begin
     wh = where(pol_enum eq pol_inc[i], count)
     if count eq 0 then message, 'pol ' + pol_inc[i] + ' not recognized.'
     pol_num[i] = wh[0]
  endfor
  pol_inc = pol_enum[pol_num[uniq(pol_num, sort(pol_num))]]

  if n_elements(type_inc) eq 0 then type_inc = ['dirty', 'model', 'res']
  type_enum = ['dirty', 'model', 'res']
  ntype = n_elements(type_inc)
  type_num = intarr(ntype)
  for i=0, ntype-1 do begin
     wh = where(type_enum eq type_inc[i], count)
     if count eq 0 then message, 'type ' + type_inc[i] + ' not recognized.'
     type_num[i] = wh[0]
  endfor
  type_inc = type_enum[type_num[uniq(type_num, sort(type_num))]]

  datafile_test = file_test(datafile)
  if datafile_test eq 0 then message, 'datafile not found'

  temp = strpos(datafile, '/', /reverse_search)
  infilebase = strmid(datafile, temp+1)
  temp2 = strpos(infilebase, '.', /reverse_search)
  datafilebase = strmid(infilebase, 0, temp2)

  if n_elements(save_path) ne 0 then froot = save_path else froot = strmid(datafile, 0, temp+1)
 
  file_obj = obj_new('idl_savefile', datafile)
  varnames = file_obj->names()
  wh_obs = where(strlowcase(varnames) eq 'obs', count_obs)
  if count_obs ne 0 then file_obj->restore, 'obs' $
  else begin
     wh_obs = where(strlowcase(varnames) eq 'obs_arr', count_obs)
     if count_obs ne 0 then begin
        file_obj->restore, 'obs_arr'
        if size(obs_arr,/type) eq 10 then begin
           n_obs = n_elements(obs_arr)

           max_baseline_vals = dblarr(n_obs)
           obs_radec_vals = dblarr(n_obs, 2)
           zen_radec_vals = dblarr(n_obs, 2)
           for i=0, n_obs-1 do begin
              if abs((*obs_arr[i]).degpix - (*obs_arr[0]).degpix) gt 0 then message, 'inconsistent degpix values in obs_arr'
              if total(abs((*obs_arr[i]).fbin_i - (*obs_arr[0]).fbin_i)) gt 0 then message, 'inconsistent fbin_i values in obs_arr'
              if total(abs((*obs_arr[i]).freq - (*obs_arr[0]).freq)) gt 0 then message, 'inconsistent freq values in obs_arr'
              max_baseline_vals[i] = (*obs_arr[i]).max_baseline
              obs_radec_vals[i, *] = [(*obs_arr[i]).obsra, (*obs_arr[i]).obsdec]
              zen_radec_vals[i, *] = [(*obs_arr[i]).zenra, (*obs_arr[i]).zendec]              
           endfor

           max_baseline_lambda = max(max_baseline_vals)
           
           degpix = (*obs_arr[0]).degpix
           freq = (*obs_arr[0]).freq
           fbin_i = (*obs_arr[0]).fbin_i

       endif else begin
           n_obs = n_elements(obs_arr)
           max_baseline_lambda = max(obs_arr.max_baseline)

           obs_radec_vals = [[obs_arr.obsra],[obs_arr.obsdec]]
           zen_radec_vals = [[obs_arr.zenra],[obs_arr.zendec]]
           theta_vals = sqrt((obs_radec_vals[*,0] - zen_radec_vals[*,0])^2d + (obs_radec_vals[*,1] - zen_radec_vals[*,1])^2d)
           max_theta = max(theta_vals)

           n_freq_vals = obs_arr.n_freq
           if total(abs(obs_arr.n_freq - obs_arr[0].n_freq)) ne 0 then message, 'inconsistent number of frequencies in obs_arr'
           n_freq = obs_arr[0].n_freq

           degpix_vals = obs_arr.degpix
           if total(abs(obs_arr.degpix - obs_arr[0].degpix)) ne 0 then message, 'inconsistent degpix values in obs_arr'
           degpix = obs_arr[0].degpix

           obs_tags = tag_names(obs_arr)
           wh_freq = where(strlowcase(obs_tags) eq 'freq', count_freq)
           if count_freq ne 0 then begin
              freq_vals = obs_arr.freq
              fbin_i_vals = obs_arr.fbin_i
           endif else begin
              freq_vals = dblarr(n_freq, n_obs)
              fbin_i_vals = fltarr(n_freq, n_obs)
              for i=0, n_obs-1 do begin
                 freq_vals[*,i] = (*obs_arr[i].bin).freq
                 fbin_i_vals[*,i] = (*obs_arr[i].bin).fbin_i
              endfor
           endelse
           if total(abs(freq_vals - rebin(freq_vals[*,0], n_freq, n_obs))) ne 0 then message, 'inconsistent freq values in obs_arr'
           freq = freq_vals[*,0]

           if total(abs(fbin_i_vals - rebin(fbin_i_vals[*,0], n_freq, n_obs))) ne 0 then $
              message, 'inconsistent fbin_i values in obs_arr'
           fbin_i = fbin_i_vals[*,0]
        endelse

       theta_vals = angle_difference(obs_radec_vals[*,1], obs_radec_vals[*,0], zen_radec_vals[*,1], zen_radec_vals[*,0], $
                                     /degree, /nearest)
       max_theta = max(theta_vals)
       
     endif else message, 'no obs or obs_arr in datafile'
  endelse
  if keyword_set(healpix) then file_obj->restore, 'nside'
  obj_destroy, file_obj
  
  n_freqbins = max(fbin_i) + 1
  frequencies = dblarr(n_freqbins)
  for i=0, n_freqbins-1 do begin
     wh = where(fbin_i eq i, count)
     if count eq 0 then stop
     
     frequencies[i] = mean(freq[wh]) / 1e6 ;; in MHz
  endfor
  
  ;; the max baseline in the obs structure is given in wavelengths, need to convert using the maximum frequency
  max_baseline = 3e8/max(freq)*max_baseline_lambda ;;else max_baseline = 342.497
  
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
  if keyword_set(log_kpar) then begin
     fadd_2dbin = fadd_2dbin + '_logkpar'
     wt_fadd_2dbin = wt_fadd_2dbin + '_logkpar'
  endif
  if keyword_set(log_kperp) then begin
     fadd_2dbin = fadd_2dbin + '_logkperp'
     wt_fadd_2dbin = wt_fadd_2dbin + '_logkperp'
  endif

  fadd_1dbin = ''
  if keyword_set(log_k) then fadd_1dbin = fadd_1dbin + '_logk'

  n_cubes = npol*ntype
  type_pol_str = strarr(n_cubes)
  for i=0, npol-1 do type_pol_str[ntype*i:i*ntype+ntype-1] = type_inc + '_' + pol_inc[i]
  data_varnames = strupcase(type_pol_str + '_cube')
  weight_varnames = strupcase('weights_' + pol_inc + '_cube')
  weight_labels = strupcase(pol_inc)
  weight_ind = intarr(n_cubes)
  for i=0, npol-1 do weight_ind[ntype*i:i*ntype+ntype-1] = i
  wt_file_labels = '_weights_' + strlowcase(weight_labels[weight_ind])
  file_labels = '_' + strlowcase(type_pol_str)
  titles = strarr(n_cubes)
  for i=0, npol-1 do titles[ntype*i:i*ntype+ntype-1] = type_inc + ' ' + pol_inc[i]

  savefilebase = froot + datafilebase + file_labels + fadd
  savefiles_2d = savefilebase + fadd_2dbin + '_2dkpower.idlsave'
  test_save_2d = file_test(savefiles_2d) *  (1 - file_test(savefiles_2d, /zero_length))

  savefiles_1d = savefilebase + fadd_1dbin + '_1dkpower.idlsave'
  test_save_1d = file_test(savefiles_1d) *  (1 - file_test(savefiles_1d, /zero_length))

  if keyword_set(refresh_binning) then begin
     test_save_2d = test_save_2d*0
     test_save_1d = test_save_1d*0
  endif

  for i=0, n_cubes-1 do begin

     ;; if binsizes are specified, check that binsize is right
     if (n_elements(kperp_bin) ne 0 or n_elements(kpar_bin) ne 0) and test_save_2d[i] gt 0 then begin
        file_obj = obj_new('idl_savefile', savefiles_2d[i])
        if n_elements(kpar_bin) ne 0 then begin
           kpar_bin_req = kpar_bin
           file_obj->restore, 'kpar_bin'
           if abs(kpar_bin - kpar_bin_req) gt 0. then test_save_2d[i]=0
           kpar_bin = kpar_bin_req
        endif
        if test_save_2d[i] gt 0 and n_elements(kpar_bin) ne 0 then begin
           kperp_bin_req = kperp_bin
           file_obj->restore, 'kperp_bin'
           if abs(kperp_binsize - kperp_bin_req) gt 0. then test_save_2d[i]=0
           kperp_bin = kperp_bin_req
        endif
        obj_destroy, file_obj
     endif

     if n_elements(k1d_bin) ne 0 and test_save_1d[i] gt 0 then begin
        file_obj = obj_new('idl_savefile', savefiles_1d[i])
        k_bin_req = k_bin
        file_obj->restore, 'k_bin'
        if abs(k_bin - k_bin_req) gt 0. then test_save_1d[i]=0
        k_bin = k_bin_req
     endif

     test = test_save_2d[i] * test_save_1d[i]

     if test eq 0 then begin
        if keyword_set(healpix) then begin
           pixelfile = datafile
           pixel_varname = 'hpx_inds'
           hpx_dftsetup_savefile = froot + datafilebase + '_dftsetup.idlsave'
           weight_savefilebase = froot + datafilebase + wt_file_labels + fadd
           
           weight_refresh = intarr(n_cubes)
           if keyword_set(refresh_dft) then begin
              temp = weight_ind[uniq(weight_ind, sort(weight_ind))]
              for j=0, n_elements(temp)-1 do weight_refresh[(where(weight_ind eq temp[j]))[0]] = 1
           endif

           fhd_3dps, datafile, data_varnames[i], datafile, weight_varnames[weight_ind[i]], frequencies, max_baseline, /healpix, $
                     nside = nside, pixelfile = pixelfile, pixelvar = pixel_varname, hpx_dftsetup_savefile = hpx_dftsetup_savefile, $
                     savefilebase = savefilebase[i], weight_savefilebase = weight_savefilebase[i], refresh = refresh_ps, $
                     dft_refresh_data=refresh_dft, dft_refresh_weight=weight_refresh[i], dft_fchunk = dft_fchunk, $
                     no_weighting = no_weighting, std_power = std_power, no_kzero = no_kzero, log_kpar = log_kpar, $
                     log_kperp = log_kperp, kpar_bin = kpar_bin, kperp_bin = kperp_bin, $
                     no_weighted_averaging = no_weighted_averaging, /quiet
        endif else $
           fhd_3dps, datafile, data_varnames[i], datafile, weight_varnames[weight_ind[i]], frequencies, max_baseline, $
                     degpix=degpix, savefilebase = savefilebase[i], refresh = refresh_ps, no_weighting = no_weighting, $
                     std_power = std_power, no_kzero = no_kzero, log_kpar = log_kpar, log_kperp = log_kperp, kpar_bin = kpar_bin, $
                     kperp_bin = kperp_bin, no_weighted_averaging = no_weighted_averaging, /quiet
     endif
  endfor


  restore, savefiles_2d[0]
  wh_good_kperp = where(total(power, 2) gt 0, count)
  if count eq 0 then stop
  ;;kperp_plot_range = [min(kperp_edges[wh_good_kperp]), max(kperp_edges[wh_good_kperp+1])]
  kperp_plot_range = [6e-3, min([max(kperp_edges[wh_good_kperp+1]),1.1e-1])]
  

  if n_elements(plot_path) ne 0 then plotfile_path = plot_path else plotfile_path = froot
  ;; plot_folder = 'fhd_data/'
  ;; plotfile_path = base_path('plots') + 'power_spectrum/' + plot_folder
  plotfile_base = plotfile_path + datafilebase + [file_labels, wt_file_labels[uniq(weight_ind, sort(weight_ind))]] + fadd

  plot_fadd = ''
  if keyword_set(grey_scale) then plot_fadd = plot_fadd + '_grey'

  plotfiles_2d = plotfile_base + fadd_2dbin + '_2dkpower' + plot_fadd + '.eps'
  plotfile_1d = plotfile_path + datafilebase + fadd + fadd_1dbin + '_1dkpower' + '.eps'

  if not keyword_set(slice_nobin) then slice_fadd = '_binned' else slice_fadd = ''
  yslice_plotfile = plotfile_base + '_xz_plane' + plot_fadd + slice_fadd + '.eps'
  xslice_plotfile = plotfile_base + '_yz_plane' + plot_fadd + slice_fadd + '.eps'
  zslice_plotfile = plotfile_base + '_xy_plane' + plot_fadd + slice_fadd + '.eps'

  if keyword_set(plot_wedge_line) then begin
     z0_freq = 1420.40 ;; MHz
     redshifts = z0_freq/frequencies - 1
     mean_redshift = mean(redshifts)

     cosmology_measures, mean_redshift, wedge_factor = wedge_factor
     ;; assume 20 degrees from pointing center to first null
     source_dist = 20d * !dpi / 180d
     fov_amp = wedge_factor * source_dist

     ;; calculate angular distance to horizon
     horizon_amp = wedge_factor * ((max_theta+90d) * !dpi / 180d)

     wedge_amp = [fov_amp, horizon_amp]
  endif else wedge_amp = 0d

  nplots = n_cubes + npol
  
  savefiles_2d_use = strarr(nplots)
  plotfiles_2d_use = strarr(nplots)
  plot_weights = intarr(nplots)
  plot_titles = strarr(nplots)
  for i=0, npol-1 do begin
     savefiles_2d_use[i*(ntype+1):i*(ntype+1)+ntype-1] = savefiles_2d[i*ntype:i*ntype+ntype-1]
     savefiles_2d_use[i*(ntype+1)+ntype] = savefiles_2d[i*ntype]
     plotfiles_2d_use[i*(ntype+1):i*(ntype+1)+ntype-1] = plotfiles_2d[i*ntype:i*ntype+ntype-1]
     plotfiles_2d_use[i*(ntype+1)+ntype] = plotfiles_2d_wt[i]
     plot_weights[i*(ntype+1)+ntype]=1
     plot_titles[i*(ntype+1):i*(ntype+1)+ntype-1] = titles[i*ntype:i*ntype+ntype-1]
     plot_titles[i*(ntype+1)+ntype] = 'weights ' + weight_labels[i]
  endfor

  if keyword_set(pub) then begin
     for i=0, nplots-1 do begin
        
        if plot_weights[i] eq 0 then $
           kpower_2d_plots, savefiles_2d_plot[i], /pub, plotfile = plotfiles_2d[i], $
                            kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, data_range = data_range, $
                            title = plot_titles[i], grey_scale = grey_scale, plot_wedge_line = plot_wedge_line, $
                            wedge_amp = wedge_amp, baseline_axis = baseline_axis $
        else kpower_2d_plots, savefiles_2d_plot[i], /plot_weights, plotfile = plotfiles_2d[i],$
                              kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, /pub, $
                              title = plot_titles[i], grey_scale = grey_scale, plot_wedge_line = plot_wedge_line, $
                              wedge_amp = wedge_amp, baseline_axis = baseline_axis
     endfor

  endif else begin

     nrow = ntype + 1 ;; ntype + 1 for weights
     ncol = npol
     positions = fltarr(4, ncol*nrow)
     
     row_val = reverse(reform(rebin(indgen(nrow), nrow, ncol), ncol*nrow))
     col_val = reform(rebin(reform(indgen(ncol), 1, ncol), nrow, ncol), ncol*nrow)
     
     positions[0,*] = col_val/double(ncol)
     positions[1,*] = row_val/double(nrow)
     positions[2,*] = (col_val+1)/double(ncol)
     positions[3,*] = (row_val+1)/double(nrow)
     
     max_ysize = 1200
     max_xsize = 1800
     if n_cubes gt 9 then begin
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
     
     multi_size = [xsize/ncol, ysize/nrow]
     for i=0, nplots-1 do begin     
        if plot_weights[i] eq 0 then $
           kpower_2d_plots, savefiles_2d_plot[i], multi_pos = positions[*,i], multi_size = multi_size, $
                            kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, data_range = data_range, $
                            title = plot_titles[i], grey_scale = grey_scale, plot_wedge_line = plot_wedge_line, $
                            wedge_amp = wedge_amp, baseline_axis = baseline_axis $
        else kpower_2d_plots, savefiles_2d_plot[i], multi_pos = positions[*,i], multi_size = multi_size, /plot_weights, $
                              kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
                              title = plot_titles[i], grey_scale = grey_scale, plot_wedge_line = plot_wedge_line, $
                              wedge_amp = wedge_amp, baseline_axis = baseline_axis
     endfor


     ;; for i=0, npol-1 do begin
     ;;    wh_pol = where(weight_ind eq i, count_pol)
     ;;    if count_pol gt 0 then begin
     ;;       file_arr = savefiles_1d[wh_pol]
     ;;       kpower_1d_plots, file_arr, window_num = i+2, pub = pub, colors = colors, names = titles[wh_pol], delta = delta
     ;;    endif
     ;; endfor

  endelse

  file_arr = savefiles_1d
  kpower_1d_plots, file_arr, window_num = i+2, colors = colors, names = titles, delta = delta, pub = pub, plotfile = plotfile_1d

end
