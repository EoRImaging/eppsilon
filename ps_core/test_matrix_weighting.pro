pro test_matrix_weighting, save_filebase = save_filebase, covar_use = covar_use, bp_edge_val = bp_edge_val, $
    plot_covar = plot_covar, use_test_signal = use_test_signal, bp_edges_high = bp_edges_high, $
    png = png, eps = eps, pdf = pdf
    
  if n_elements(save_filebase) gt 0 or keyword_set(png) or keyword_set(eps) or keyword_set(pdf) then pub = 1 else pub = 0
  
  if pub eq 1 then begin
    save_path = base_path('plots') + 'single_use/'
    if n_elements(save_filebase) gt 0 then savefile = save_path + save_filebase $
    else savefile = save_path + 'covar_weight_test'
    
    if keyword_set(png) and keyword_set(eps) and keyword_set(pdf) then begin
      print, 'only one of eps, pdf and png can be set, using png'
      eps = 0
      png = 1
    endif
    
    if not keyword_set(png) and not keyword_set(pdf) and not keyword_set(eps) then png = 1
    
    if keyword_set(png) or keyword_set(pdf) then delete_ps = 1 else if keyword_set(eps) then delete_ps = 0
    
  endif
  
  if n_elements(covar_use) eq 0 then covar_use = 'both'
  covar_types = ['instrument', 'foreground', 'both']
  wh_covar = where(covar_types eq covar_use, count_covar)
  if count_covar eq 0 then message, 'covar_use not recognized, options are: ' + covar_types
  
  if n_elements(bp_edge_val) eq 0 then bp_edge_val = 1e4;!values.d_infinity
  
  delta_f = [0.08, 0.16, 0.16]
  n_freq = [384, 384, 192]
  n_trials = n_elements(delta_f)
  
  delta_eta = 1./(n_freq*delta_f)
  
  norm_err_factor_range = fltarr(n_trials, 2)
  peak_norm = fltarr(n_trials)
  peak_norm_meas = fltarr(n_trials)
  peak_ratio = fltarr(n_trials, 5)
  power_ratio = fltarr(n_trials, 5)
  peak_power = fltarr(n_trials, 2)
  unflagged_frac = fltarr(n_trials)
  for i=0, n_trials-1 do begin
  
    frequencies = findgen(n_freq[i])*delta_f[i] + 130 ;; in MHz
    
    test_signal = sin(36*!dpi*frequencies/(n_freq[i]*delta_f[i]))
    
    test_signal_ft = shift(fft(test_signal), n_freq[i]/2)*n_freq[i]*delta_f[i]
    
    identity = diag_matrix([fltarr(n_freq[i])+1d])
    
    ft_matrix = fft(identity, dimension=1)*n_freq[i]*delta_f[i]
    inv_ft_matrix = fft(identity, dimension=1, /inverse)*delta_eta[i]
    temp = matrix_multiply(ft_matrix, inv_ft_matrix)
    
    
    test_signal_ft2 = shift(matrix_multiply(ft_matrix, test_signal), n_freq[i]/2)
    test_signal_ft3 = shift(reform(matrix_multiply(transpose(test_signal), ft_matrix, /btranspose)), n_freq[i]/2)
    ;print, minmax(test_signal_ft-test_signal_ft2)
    ;print, minmax(test_signal_ft-test_signal_ft3)
    
    test_signal2 = matrix_multiply(inv_ft_matrix, shift(test_signal_ft2, n_freq[i]/2))
    test_signal3 = fft(shift(test_signal_ft, n_freq[i]/2), /inverse)*delta_eta[i]
    ;print, minmax(test_signal - test_signal2)
    ;print, minmax(test_signal - test_signal3)
    
    eta_vals = indgen(n_freq[i])*delta_eta[i] - n_freq[i]*delta_eta[i]
    
    eta_var = [1e4, fltarr(n_freq[i]-1)+1.]
    if keyword_set(use_test_signal) then eta_var = shift(test_signal_ft * 10. / max(test_signal_ft), n_freq[i]/2)*1e4
    ;eta_var = sin(4.*!dpi*frequencies/(n_freq[i]*delta_f[i]))
    covar_eta_fg = diag_matrix(shift(eta_var, n_freq[i]/2))
    
    ;; needs to be an inverse FFT
    covar_f_fg = matrix_multiply(inv_ft_matrix, matrix_multiply(shift(covar_eta_fg, n_freq[i]/2, n_freq[i]/2), conj(inv_ft_matrix), /btranspose))
    covar_eta_fg_check = shift(matrix_multiply(ft_matrix, matrix_multiply(covar_f_fg, conj(ft_matrix), /btranspose)), n_freq[i]/2, n_freq[i]/2)
    ;print, minmax(covar_eta_fg - covar_eta_fg_check)
    
    large_val = bp_edge_val
    freq_var = fltarr(n_freq[i])+1
    n_coarse = 24
    wh_edge = where(indgen(n_freq[i]) mod (n_freq[i]/n_coarse) eq 0 or indgen(n_freq[i]) mod (n_freq[i]/n_coarse) eq n_freq[i]/n_coarse-1, count_edges)
    if count_edges gt 0 then freq_var[wh_edge] = large_val
    
    covar_f_inst = diag_matrix(freq_var)
    covar_eta_inst = shift(matrix_multiply(ft_matrix, matrix_multiply(covar_f_inst, conj(ft_matrix), /btranspose)), n_freq[i]/2, n_freq[i]/2)
    
    ;signal_power = shift([5e6, fltarr(n_freq[i]-1)+1.], n_freq[i]/2)
    ;signal_power = abs(test_signal_ft)
    
    ;signal_eta = randomn(seed, n_freq[i])*sqrt(signal_power) + complex(0,1)*randomn(seed, n_freq[i])*sqrt(signal_power)
    ;signal = fft(shift(signal_eta, n_freq[i]/2), /inverse)*delta_eta[i]
    
    signal = fft([5e6, fltarr(n_freq[i]-1)+1.],/inverse)*delta_eta[i]
    if keyword_set(use_test_signal) then signal = test_signal
    
    signal = signal + randomn(seed, n_freq[i])*max(signal)/100.
    signal_nobp = signal
    
    if keyword_set(bp_edges_high) then signal_weights = freq_var $
    else signal_weights = 1./freq_var
    
    signal_weights_ft = shift(fft(signal_weights), n_freq[i]/2)*n_freq[i]*delta_f[i]
    signal = signal*signal_weights
    
    signal_ft = shift(fft(signal), n_freq[i]/2)*n_freq[i]*delta_f[i]
    signal_nobp_ft = shift(fft(signal_nobp), n_freq[i]/2)*n_freq[i]*delta_f[i]
    
    case covar_use of
      'instrument': covar_f = covar_f_inst
      'foreground': covar_f = covar_f_fg
      'both': covar_f = covar_f_fg+covar_f_inst
    endcase
    
    covar_eta = shift(matrix_multiply(ft_matrix, matrix_multiply(covar_f, conj(ft_matrix), /btranspose)), n_freq[i]/2, n_freq[i]/2)
    
    inv_covar_f = la_invert(covar_f)
    inv_covar_eta = shift(matrix_multiply(ft_matrix, matrix_multiply(inv_covar_f, conj(ft_matrix), /btranspose)), n_freq[i]/2, n_freq[i]/2)
    ;inv_covar_eta_inst = shift(matrix_multiply(ft_matrix, matrix_multiply(1./covar_f, conj(ft_matrix), /btranspose)), n_freq[i]/2, n_freq[i]/2)
    
    wt_signal = matrix_multiply(inv_covar_f, signal_nobp) * delta_f[i]
    ;wt_signal = matrix_multiply(inv_covar_f, signal) * delta_f[i]
    
    wt_signal_ft = shift(fft(wt_signal), n_freq[i]/2)*n_freq[i]*delta_f[i]
    
    norm = total(abs(inv_covar_eta)^2.,2)*delta_eta[i]
    temp = shift(abs(fft(total(abs(inv_covar_f)^2.,2))), n_freq[i]/2)*n_freq[i]*delta_f[i]
    
    wt_power = abs(wt_signal_ft)^2./norm
    wt_power_norm = wt_power * n_freq[i] / total((signal_weights)^2.)
    power = abs(signal_ft)^2./(n_freq[i]*delta_f[i])
    power_norm = abs(signal_ft)^2./(total(signal_weights^2.)*delta_f[i])
    power_nobp = abs(signal_nobp_ft)^2./(n_freq[i]*delta_f[i])
    
    ;print, minmax(wt_power/power)
    ;print, minmax(power/wt_power)
    ;wh_peak = (where(signal_power eq max(signal_power)))[0]
    wh_peak = (where(power_nobp eq max(power_nobp)))[0]
    peak_power[i,*] = [power_nobp[wh_peak], wt_power[wh_peak]]
    
    
    ;norm_err_factor_range[i,*] = minmax(wt_power/power)
    ratio_names = ['flag/no flag', 'flag norm/no flag',  'weight/flag', 'weight/no flag', 'weight norm/no flag']
    
    peak_ratio[i, 0] = power[wh_peak]/power_nobp[wh_peak]
    peak_ratio[i, 1] = power_norm[wh_peak]/power_nobp[wh_peak]
    peak_ratio[i, 2] = wt_power[wh_peak]/power[wh_peak]
    peak_ratio[i, 3] = wt_power[wh_peak]/power_nobp[wh_peak]
    peak_ratio[i, 4] = wt_power_norm[wh_peak]/power_nobp[wh_peak]
    
    power_ratio[i, 0] = total(power)/total(power_nobp)
    power_ratio[i, 1] = total(power_norm)/total(power_nobp)
    power_ratio[i, 2] = total(wt_power)/total(power)
    power_ratio[i, 3] = total(wt_power)/total(power_nobp)
    power_ratio[i, 4] = total(wt_power_norm)/total(power_nobp)
    unflagged_frac[i] = 1.-count_edges/float(n_freq[i])
  endfor
  
  if keyword_set(plot_covar) then begin
    if pub eq 1 then covar_savefile = savefile
    
    start_multi_params = {ncol:2, nrow:2, ordering:'row'}
    quick_image, covar_f, frequencies, frequencies, start_multi_params = start_multi_params, multi_pos = multi_pos, $
      xtitle = 'frequency (MHz)', ytitle = 'frequency (MHz)', title = 'freq covariance', $
      savefile = covar_savefile, png = png, eps = eps, pdf = pdf
    quick_image, covar_eta, eta_vals, eta_vals, multi_pos = multi_pos[*,1], /noerase, $
      xtitle = 'eta (1/MHz)', ytitle = 'eta (1/MHz)', title = 'eta covariance', $
      savefile = covar_savefile, png = png, eps = eps, pdf = pdf
    quick_image, inv_covar_f, frequencies, frequencies, multi_pos = multi_pos[*,2], /noerase, $
      xtitle = 'frequency (MHz)', ytitle = 'frequency (MHz)', title = 'freq inverse covariance', $
      savefile = covar_savefile, png = png, eps = eps, pdf = pdf
    quick_image, inv_covar_eta, eta_vals, eta_vals, multi_pos = multi_pos[*,3], /noerase, $
      xtitle = 'eta (1/MHz)', ytitle = 'eta (1/MHz)', title = 'eta inverse covariance', $
      savefile = covar_savefile, png = png, eps = eps, pdf = pdf
      
    if keyword_set(pub) then begin
      cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density = 600
    endif
  endif
  
  if pub eq 1 then power_savefile = savefile+'_power'
  
  if keyword_set(pub) then begin
    font = 1
    charsize = 1.25
    thick = 3
    xthick = 3
    ythick = 3
  endif else begin
    font = -1
    charsize = 1.2
  endelse
  
  nrow = 2
  ncol = 2
  
  col_val = reform(rebin(indgen(ncol), ncol, nrow), ncol*nrow)
  row_val = reverse(reform(rebin(reform(indgen(nrow), 1, nrow), ncol, nrow), ncol*nrow))
  
  multi_pos = fltarr(4,nrow*ncol)
  multi_pos[0,*] = col_val/double(ncol)
  multi_pos[1,*] = row_val/double(nrow)
  multi_pos[2,*] = (col_val+1)/double(ncol)
  multi_pos[3,*] = (row_val+1)/double(nrow)
  
  
  margin = [0.2, 0.2, 0.2, 0.2]
  plot_pos = [margin[0], margin[1], (1-margin[2]), (1-margin[3])]
  
  xlen = multi_pos[2,*] - multi_pos[0,*]
  ylen = multi_pos[3,*] - multi_pos[1,*]
  pos_use = [xlen * plot_pos[0] + multi_pos[0,*], ylen * plot_pos[1] + multi_pos[1,*], $
    xlen * plot_pos[2] + multi_pos[0,*], ylen * plot_pos[3] + multi_pos[1,*]]
    
  ; cgplot, frequencies, signal, position=pos_use[*,0], xtitle = 'frequency (MHz)', yrange = minmax([signal, wt_signal]), $
  ;   title = 'Signal & Weighted Signal', xstyle=1, thick = thick, xthick = xthick, ythick = ythick, font = font, charsize = charsize
  ; cgplot, frequencies, wt_signal, /over, color='red', font = font, charsize = charsize, linestyle=2
    
    
  wh_pos = where([power, wt_power, power_nobp] gt 0, count_pos)
  if count_pos gt 0 then yrange = minmax(([power, wt_power, power_nobp])[wh_pos])
  
  if pub ne 1 then begin
    if windowavailable(2) then begin
      wset, 2
      if !d.x_size ne 600 or !d.y_size ne 600 then make_win = 1 else make_win = 0
    endif else make_win = 1
    if make_win eq 1 then window, 2, xsize = 600, ysize = 600
  endif
  
  cgplot, eta_vals, power, position=pos_use[*,0], /ylog, xtitle = 'eta (1/MHz)', yrange = yrange, $
    title = 'Power & Weighted Power', xstyle=1, font = font, charsize = charsize
    
  cgplot, eta_vals, wt_power, /over, color='red', linestyle = 2, font = font, charsize = charsize
  ;  cgplot, eta_vals, wt_power, position=pos_use[*,1], /noerase, /ylog, xtitle = 'eta (1/MHz)', yrange = yrange, $
  ;    title = 'Weighted Power', xstyle=1, font = font, charsize = charsize
  
  ;cgplot, eta_vals, wt_power/power, position=pos_use[*,1], /noerase, xtitle = 'eta (1/MHz)', yrange = [0.5,1.5], $
  ;  title = 'Weighted Power/Power', xstyle=1, font = font, charsize = charsize
  
  xrange = minmax(eta_vals[n_freq[2]*.45:n_freq[2]*.55])
  cgplot, eta_vals, power, position=pos_use[*,1], /noerase, xtitle = 'eta (1/MHz)', yrange = [.1,1.1]*yrange[1], $
    title = 'Power & Weighted Power', xstyle=1, font = font, charsize = charsize, xrange = xrange
  cgplot, eta_vals, wt_power, /over, color='red', linestyle = 2, font = font, charsize = charsize
  
  
  cgplot, eta_vals, power_nobp, position=pos_use[*,2], /noerase, /ylog, yrange = yrange, title= 'Power w/o bandpass', $
    xtitle = 'eta (1/MHz)', font = font, charsize = charsize, xstyle=1
  cgplot, eta_vals, wt_power, /over, color='red', linestyle = 2, font = font, charsize = charsize
  
  
  ;  cgplot, eta_vals, wt_power/power_nobp, position=pos_use[*,3], /noerase, xtitle = 'eta (1/MHz)', yrange = [0.5,1.5], $
  ;    title = 'Weighted Power/Power no bandpass', xstyle=1, font = font, charsize = charsize
  cgplot, eta_vals, power_nobp, position=pos_use[*,3], /noerase, xtitle = 'eta (1/MHz)', yrange = [.1,1.1]*yrange[1], $
    title = 'Power w/o bandpass', xstyle=1, font = font, charsize = charsize, xrange = xrange
  cgplot, eta_vals, wt_power, /over, color='red', linestyle = 2, font = font, charsize = charsize
  
  
  
  ;cgplot, eta_vals, norm, position=pos_use[*,3], /noerase, color='blue', yrange = minmax(norm), title= 'Normalization', $
  ;  xtitle = 'eta (1/MHz)', font = font, charsize = charsize, xstyle=1
  
  
  
  
  if keyword_set(pub) then begin
    cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density = 600
  endif
  
  ;print, norm_err_factor_range
  print, 'peak ratios:'
  for i=0, 4 do print, ratio_names[i], peak_ratio[*, i], format = '(a20, 3f9.3)'
  print, ''
  print, 'integral ratios:'
  for i=0, 4 do print, ratio_names[i], power_ratio[*, i], format = '(a20, 3f9.3)'
  print, ''
  print, 'unflagged fraction', unflagged_frac, format = '(a20, 3f9.3)'
;print, ''
;print, peak_power
;print, peak_norm
;print, peak_norm_meas
;print, peak_norm_meas/peak_norm
  
;print, peak_power
  
end