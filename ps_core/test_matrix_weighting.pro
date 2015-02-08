pro test_matrix_weighting, save_filebase = save_filebase, png = png, eps = eps, pdf = pdf, plot_covar = plot_covar

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
  
  delta_f = [0.08, 0.16, 0.16]
  n_freq = [384, 384, 192]
  n_trials = n_elements(delta_f)
  
  delta_eta = 1./(n_freq*delta_f)
  
  norm_err_factor_range = fltarr(n_trials, 2)
  peak_ratio = fltarr(n_trials, 2)
  peak_power = fltarr(n_trials, 2)
  unflagged_frac = fltarr(n_trials)
  for i=0, n_trials-1 do begin
  
    frequencies = findgen(n_freq[i])*delta_f[i] + 130 ;; in MHz
    
    test_signal = sin(12.*!dpi*frequencies/(n_freq[i]*delta_f[i]))
    
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
    
    eta_var = [10, fltarr(n_freq[i]-1)+1.]
    ;eta_var = shift(test_signal_ft * 10. / max(test_signal_ft), n_freq[i]/2)*1e4
    ;eta_var = sin(4.*!dpi*frequencies/(n_freq[i]*delta_f[i]))
    covar_eta1 = diag_matrix(shift(eta_var, n_freq[i]/2))
    
    ;; needs to be an inverse FFT
    covar_f1 = matrix_multiply(inv_ft_matrix, matrix_multiply(shift(covar_eta1, n_freq[i]/2, n_freq[i]/2), conj(inv_ft_matrix), /btranspose))
    covar_eta1_check = shift(matrix_multiply(ft_matrix, matrix_multiply(covar_f1, conj(ft_matrix), /btranspose)), n_freq[i]/2, n_freq[i]/2)
    ;print, minmax(covar_eta1 - covar_eta1_check)
    
    freq_var = fltarr(n_freq[i])+1
    wh_edge = where(indgen(n_freq[i]) mod (n_freq[i]/24) eq 0 or indgen(n_freq[i]) mod (n_freq[i]/24) eq n_freq[i]/24-1, count_edges)
    if count_edges gt 0 then freq_var[wh_edge] = 2;!values.d_infinity
    covar_f2 = diag_matrix(freq_var)
    covar_eta2 = shift(matrix_multiply(ft_matrix, matrix_multiply(covar_f2, conj(ft_matrix), /btranspose)), n_freq[i]/2, n_freq[i]/2)
    
    signal = shift(fft([5e4, fltarr(n_freq[i]-1)+1.]), n_freq[i]/2)*n_freq[i]*delta_f[i]
    ;signal = test_signal
    signal = signal + randomn(seed, n_freq[i])*max(signal)/100.
    signal_nobp = signal
    if count_edges gt 0 then signal[wh_edge] = 0.5*signal[wh_edge]
    
    signal_ft = shift(fft(signal), n_freq[i]/2)*n_freq[i]*delta_f[i]
    signal_nobp_ft = shift(fft(signal_nobp), n_freq[i]/2)*n_freq[i]*delta_f[i]
    
    
    covar_f = covar_f1+covar_f2
    ;covar_f = covar_f2
    covar_eta = shift(matrix_multiply(ft_matrix, matrix_multiply(covar_f, conj(ft_matrix), /btranspose)), n_freq[i]/2, n_freq[i]/2)
    
    inv_covar_f = la_invert(covar_f)
    inv_covar_eta = shift(matrix_multiply(ft_matrix, matrix_multiply(inv_covar_f, conj(ft_matrix), /btranspose)), n_freq[i]/2, n_freq[i]/2)
    ;inv_covar_eta2 = shift(matrix_multiply(ft_matrix, matrix_multiply(1./covar_f, conj(ft_matrix), /btranspose)), n_freq[i]/2, n_freq[i]/2)
    
    wt_signal = matrix_multiply(inv_covar_f, signal) * delta_f[i]
    ;wt_signal2 = matrix_multiply(1./covar_f, signal)
    wt_signal_ft = shift(fft(wt_signal), n_freq[i]/2)*n_freq[i]*delta_f[i]
    ;wt_signal2_ft = shift(fft(wt_signal2), n_freq[i]/2)*n_freq[i]*delta_f[i]
    
    norm = total(abs(inv_covar_eta)^2.,2)*delta_eta[i]
    ;norm2 = total(abs(inv_covar_eta2)^2.,2)*delta_eta[i]
    
    wt_power = abs(wt_signal_ft)^2./norm
    ;wt_power2 = abs(wt_signal2_ft)^2./norm2
    power = abs(signal_ft)^2./(n_freq[i]*delta_f[i])
    power_nobp = abs(signal_nobp_ft)^2./(n_freq[i]*delta_f[i])
    
    ;print, minmax(wt_power/power)
    ;print, minmax(power/wt_power)
    norm_err_factor_range[i,*] = minmax(wt_power/power)
    peak_power[i,*] = [max(power) ,max(wt_power)]
    peak_ratio[i, 0] = max(wt_power)/max(power)
    peak_ratio[i, 1] = max(wt_power)/max(power_nobp)
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
  print, peak_ratio
  print, unflagged_frac
;print, peak_power
  
end