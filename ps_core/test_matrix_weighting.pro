pro test_matrix_weighting, save_filebase = save_filebase, covar_use = covar_use, bp_edge_val = bp_edge_val, $
    plot_covar = plot_covar, use_test_signal = use_test_signal, mode_use = mode_use, no_vec = no_vec, reg_var = reg_var, $
    png = png, eps = eps, pdf = pdf, random_bp_locs = random_bp_locs
    
  if n_elements(save_filebase) gt 0 or keyword_set(png) or keyword_set(eps) or keyword_set(pdf) then pub = 1 else pub = 0
  
  if pub eq 1 then begin
    save_path = base_path('plots') + 'single_use/'
    if n_elements(save_filebase) gt 0 then savefile = save_path + save_filebase $
    else savefile = save_path + 'covar_weight_test'
    
    if not keyword_set(png) and not keyword_set(pdf) and not keyword_set(eps) then png = 1
    
    if keyword_set(eps) then delete_ps = 1 else delete_ps = 0
    
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
  
  peak_ratio = fltarr(n_trials, 5)
  power_ratio = fltarr(n_trials, 5)
  peak_power = fltarr(n_trials, 3)
  vec_wt_sum = fltarr(n_trials, 2)
  unflagged_frac = fltarr(n_trials)
  for i=0, n_trials-1 do begin
  
    frequencies = findgen(n_freq[i])*delta_f[i] + 130 ;; in MHz
    
    test_signal = sin(36*!dpi*frequencies/(n_freq[i]*delta_f[i]))
    
    test_signal_ft = shift(fft(test_signal), n_freq[i]/2)*n_freq[i]*delta_f[i]
    
    identity = diag_matrix([fltarr(n_freq[i])+1d])
    
    ft_matrix = fft(identity, dimension=1)*n_freq[i]*delta_f[i]
    inv_ft_matrix = fft(identity, dimension=1, /inverse)*delta_eta[i]
    ft_matrix_inv = la_invert(ft_matrix)
    temp = matrix_multiply(ft_matrix, inv_ft_matrix)
    
    test_signal_ft2 = shift(matrix_multiply(ft_matrix, test_signal), n_freq[i]/2)
    test_signal_ft3 = shift(reform(matrix_multiply(transpose(test_signal), ft_matrix, /btranspose)), n_freq[i]/2)
    ;print, minmax(test_signal_ft-test_signal_ft2)
    ;print, minmax(test_signal_ft-test_signal_ft3)
    
    test_signal2 = matrix_multiply(inv_ft_matrix, shift(test_signal_ft2, n_freq[i]/2))
    test_signal3 = fft(shift(test_signal_ft, n_freq[i]/2), /inverse)*delta_eta[i]
    ;print, minmax(test_signal - test_signal2)
    ;print, minmax(test_signal - test_signal3)
    
    eta_vals = indgen(n_freq[i])*delta_eta[i] - n_freq[i]*delta_eta[i]/2.
    if n_elements(mode_use) eq 0 then mode_use = 0
    
    max_eta_val = 1e4
    eta_var = fltarr(n_freq[i])
    if covar_use eq 'foreground' then eta_var = eta_var + 1 ;; ensure invertability
    eta_var[mode_use] = max_eta_val
    if keyword_set(use_test_signal) then eta_var = shift(test_signal_ft * 10. / max(test_signal_ft), n_freq[i]/2)*max_eta_val
    ;eta_var = sin(4.*!dpi*frequencies/(n_freq[i]*delta_f[i]))
    covar_eta_fg = diag_matrix(shift(eta_var, n_freq[i]/2))
    
    ;; needs to be an inverse FFT
    covar_f_fg = matrix_multiply(inv_ft_matrix, matrix_multiply(shift(covar_eta_fg, n_freq[i]/2, n_freq[i]/2), conj(inv_ft_matrix), /btranspose))
    covar_eta_fg_check = shift(matrix_multiply(ft_matrix, matrix_multiply(covar_f_fg, conj(ft_matrix), /btranspose)), n_freq[i]/2, n_freq[i]/2)
    ;print, minmax(covar_eta_fg - covar_eta_fg_check)
    
    large_val = bp_edge_val
    if n_elements(reg_var) eq 0 then reg_var = 1
    freq_var = fltarr(n_freq[i])+reg_var
    n_coarse = 24
    wh_edge = where(indgen(n_freq[i]) mod (n_freq[i]/n_coarse) eq 0 or indgen(n_freq[i]) mod (n_freq[i]/n_coarse) eq n_freq[i]/n_coarse-1, count_edges)
    if keyword_set(random_bp_locs) then begin
      bp_locs = round(randomu(seed, count_edges)*(n_freq[i]-1))
      bp_locs = bp_locs[sort(bp_locs)]
      bp_locs = bp_locs[uniq(bp_locs)]
      while n_elements(bp_locs) lt count_edges do begin
        bp_locs = [bp_locs, round(randomu(seed, count_edges-n_elements(bp_locs))*(n_freq[i]-1))]
        bp_locs = bp_locs[sort(bp_locs)]
        bp_locs = bp_locs[uniq(bp_locs)]
      endwhile
      if max(bp_locs) ge n_freq[i] or min(bp_locs) lt 0 then stop
      wh_edge = bp_locs
    endif
    
    if count_edges gt 0 then freq_var[wh_edge] = large_val
    
    ;signal_power = shift([5e6, fltarr(n_freq[i]-1)+1.], n_freq[i]/2)
    ;signal_power = abs(test_signal_ft)
    
    ;signal_eta = randomn(seed, n_freq[i])*sqrt(signal_power) + complex(0,1)*randomn(seed, n_freq[i])*sqrt(signal_power)
    ;signal = fft(shift(signal_eta, n_freq[i]/2), /inverse)*delta_eta[i]
    
    signal_in_ft = fltarr(n_freq[i]) ;; in eta space
    signal_in_ft[mode_use] = 5e6 ;/ n_elements(mode_use)
    signal = fft(signal_in_ft,/inverse)*delta_eta[i] ;; in frequency space
    if keyword_set(use_test_signal) then signal = test_signal
    if max(abs(imaginary(signal))) eq 0 then signal = real_part(signal)
    
    noise_init = randomn(seed, n_freq[i])*max(signal)/100.
    
    signal_nobp = signal + noise_init
    
    if covar_use ne 'foreground' then noise = noise_init * sqrt(freq_var)
    signal = signal + noise
    
    covar_f_inst = diag_matrix(freq_var)
    ;    covar_f_inst = diag_matrix(noise)
    covar_eta_inst = shift(matrix_multiply(ft_matrix, matrix_multiply(covar_f_inst, conj(ft_matrix), /btranspose)), n_freq[i]/2, n_freq[i]/2)
    
    ;    filename = base_path('data') + 'single_use/foreground_covar_eta.txt'
    ;    header = ['delta_f=', 'n_freq=', 'delta_eta=']+number_formatter([delta_f[i], n_freq[i], delta_eta[i]])
    ;    textfast, covar_eta_fg, header, file_path=filename, /write
    ;
    ;    filename = base_path('data') + 'single_use/instrument_covar_f.txt'
    ;    header = ['delta_f=', 'n_freq=', 'delta_eta=']+number_formatter([delta_f[i], n_freq[i], delta_eta[i]])
    ;    textfast, covar_f_inst, header, file_path=filename, /write
    
    if not keyword_set(no_vec) then begin
      vec_wt = 1./freq_var
      wh_inf = where(finite(freq_var) ne 1, count_inf)
      ;      vec_wt = 1./noise^2.
      ;      wh_inf = where(finite(noise) ne 1, count_inf)
      if count_inf gt 0 then vec_wt[wh_inf] = 0
      wt_signal_vec = signal * vec_wt
      if count_inf gt 0 then wt_signal_vec[wh_inf] = 0
    endif else begin
      wt_signal_vec = signal
      vec_wt = fltarr(n_freq[i]) + 1.
    endelse
    
    signal_ft = shift(fft(signal), n_freq[i]/2)*n_freq[i]*delta_f[i] ;; signal (with noise) in eta
    wt_signal_vec_ft = shift(fft(wt_signal_vec), n_freq[i]/2)*n_freq[i]*delta_f[i] ;; variance weighted signal in eta
    signal_nobp_ft = shift(fft(signal_nobp), n_freq[i]/2)*n_freq[i]*delta_f[i] ;; signal w/o bandpass edges in eta
    
    case covar_use of
      'instrument': covar_f = covar_f_inst
      'foreground': covar_f = covar_f_fg
      'both': covar_f = covar_f_fg+covar_f_inst
    endcase
    
    covar_eta = shift(matrix_multiply(ft_matrix, matrix_multiply(covar_f, conj(ft_matrix), /btranspose)), n_freq[i]/2, n_freq[i]/2)
    
    inv_covar_f = la_invert(covar_f)
    inv_covar_eta = shift(matrix_multiply(ft_matrix, matrix_multiply(inv_covar_f, conj(ft_matrix), /btranspose)), n_freq[i]/2, n_freq[i]/2)
    
    inv_covar_eta2 = la_invert(covar_eta)
    ;    print, max(abs(inv_covar_eta - inv_covar_eta2))
    
    ;    print, minmax(abs(ft_matrix_inv - transpose(conj(ft_matrix))))
    
    inv_covar_eta_trace = trace(inv_covar_eta)
    inv_covar_f_trace = trace(inv_covar_f)
    
    ;    print, abs(inv_covar_f_trace), abs(inv_covar_eta_trace)
    ;    print, (2*trace(abs(inv_covar_f)) - trace(abs(inv_covar_f)^2.))/n_freq[i]
    inv_covar_f_inst = la_invert(covar_f_inst)
    mean_inv_var_inst = mean(diag_matrix(abs(inv_covar_f_inst)))
    mean_inv_var = mean(diag_matrix(abs(inv_covar_f)))
    
    ;wt_signal = matrix_multiply(inv_covar_f, signal_nobp) * delta_f[i]
    temp_signal = signal
    wh_sig_inf = where(finite(signal) ne 1, count_sig_inf)
    if count_sig_inf gt 0 then if max(abs(inv_covar_f[wh_sig_inf, *])) eq 0 then temp_signal[wh_sig_inf] = 0 else stop
    wt_signal = matrix_multiply(inv_covar_f, temp_signal) * delta_f[i] ;; covariance weighted signal in frequency
    
    wt_signal_ft = shift(fft(wt_signal), n_freq[i]/2)*n_freq[i]*delta_f[i] ;; covariance weighted signal in eta
    
    norm = total(abs(inv_covar_eta)^2.,2)*delta_eta[i]  ;; standard covariance normalization
    temp = shift(abs(fft(total(abs(inv_covar_f)^2.,2))), n_freq[i]/2)*n_freq[i]*delta_f[i]
    
    wt_power = abs(wt_signal_ft)^2./norm ;; standard covariance weighted power
    power = abs(wt_signal_vec_ft)^2./(n_freq[i]*delta_f[i]) ;; wrong normalization variance weighted power
    power_norm = abs(wt_signal_vec_ft)^2./(total(vec_wt^2.)*delta_f[i]) ;; correct normalization variance weighted power
    power_nobp = abs(signal_nobp_ft)^2./(n_freq[i]*delta_f[i]) ;; power of signal w/o bandpass edges
    
    wh_peak = where(power_nobp eq max(power_nobp), count_npeak)
    wh_peak = wh_peak[0]
    
    ;; peak power for signal w/o bandpass, var weighted correct norm, covar weighted
    peak_power[i,*] = [power_nobp[wh_peak], power_norm[wh_peak], wt_power[wh_peak]]
    vec_wt_sum[i,*] = [total(vec_wt), total(vec_wt^2.)]
    
    ;; empirical additional normalization that returns correct covariance power
    norm2 = total(vec_wt)^2./(n_freq[i]*total(vec_wt^2.))
    ;norm2 = power_norm[wh_peak]/power_nobp[wh_peak]
    wt_power2 = wt_power/norm2 ;; correctly normalized covar power
    
    ratio_names = ['flag/no flag', 'flag norm/no flag',  'weight/flag', 'weight/no flag', 'weight norm2/no flag']
    
    peak_ratio[i, 0] = power[wh_peak]/power_nobp[wh_peak]
    peak_ratio[i, 1] = power_norm[wh_peak]/power_nobp[wh_peak]
    peak_ratio[i, 2] = wt_power[wh_peak]/power[wh_peak]
    peak_ratio[i, 3] = wt_power[wh_peak]/power_nobp[wh_peak]
    peak_ratio[i, 4] = wt_power2[wh_peak]/power_nobp[wh_peak]
    
    power_ratio[i, 0] = total(power)/total(power_nobp)
    power_ratio[i, 1] = total(power_norm)/total(power_nobp)
    power_ratio[i, 2] = total(wt_power)/total(power)
    power_ratio[i, 3] = total(wt_power)/total(power_nobp)
    power_ratio[i, 4] = total(wt_power2)/total(power_nobp)
    unflagged_frac[i] = 1.-count_edges/float(n_freq[i])
  endfor
  
  if keyword_set(plot_covar) then begin
    if pub eq 1 then covar_savefile = savefile
    
    start_multi_params = {ncol:2, nrow:2, ordering:'row'}
    
    covar_f_range = minmax(covar_f(where(finite(covar_f))))
    covar_f_range = 10^float([floor(alog10(covar_f_range[0])), ceil(alog10(covar_f_range[1]))])
    quick_image, abs(covar_f), frequencies, frequencies, start_multi_params = start_multi_params, multi_pos = multi_pos, $
      xtitle = 'frequency (MHz)', ytitle = 'frequency (MHz)', title = '|freq covariance|', $
      savefile = covar_savefile, png = png, eps = eps, pdf = pdf, data_range = covar_f_range, /log
    quick_image, abs(covar_eta), eta_vals, eta_vals, multi_pos = multi_pos[*,1], /noerase, $
      xtitle = 'eta (1/MHz)', ytitle = 'eta (1/MHz)', title = '|eta covariance|', $
      savefile = covar_savefile, png = png, eps = eps, pdf = pdf, /log
      
    inv_covar_f_range = minmax(abs(inv_covar_f[where(abs(inv_covar_f) gt 0)]))
    inv_covar_f_range = 10^float([floor(alog10(inv_covar_f_range[0])), ceil(alog10(inv_covar_f_range[1]))])
    quick_image, abs(inv_covar_f), frequencies, frequencies, multi_pos = multi_pos[*,2], /noerase, $
      xtitle = 'frequency (MHz)', ytitle = 'frequency (MHz)', title = '|freq inverse covariance|', $
      savefile = covar_savefile, png = png, eps = eps, pdf = pdf, /log, $
      data_range = inv_covar_f_range
      ;data_range = [-1,1]*inv_covar_f_range[1], data_min_abs = inv_covar_f_range[0], color_profile = 'sym_log'
      
    inv_covar_eta_range = minmax(abs(inv_covar_eta[where(abs(inv_covar_eta) gt 0)]))
    inv_covar_eta_range = 10^float([floor(alog10(inv_covar_eta_range[0])), ceil(alog10(inv_covar_eta_range[1]))])
    quick_image, abs(inv_covar_eta), eta_vals, eta_vals, multi_pos = multi_pos[*,3], /noerase, $
      xtitle = 'eta (1/MHz)', ytitle = 'eta (1/MHz)', title = '|eta inverse covariance|', $
      savefile = covar_savefile, png = png, eps = eps, pdf = pdf, /log, $
      data_range = inv_covar_eta_range
      ;data_range = [-1,1]*inv_covar_eta_range[1], data_min_abs = inv_covar_eta_range[0], color_profile = 'sym_log'
   
    if keyword_set(pub) then begin
      if png then begin
        if pdf then delete_ps_use = 0 else delete_ps_use = delete_ps
        cgps_close, /png, delete_ps = delete_ps_use, density = 600
      endif
      if pdf then begin
        if not png then cgps_close
        cgps2pdf, covar_savefile, delete_ps=delete_ps
      endif
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
  
  wh_pos = where([wt_power, power_nobp] gt 0, count_pos)
  if count_pos gt 0 then yrange2 = minmax(([wt_power, power_nobp])[wh_pos])
  
  
  if pub ne 1 then begin
    if windowavailable(2) then begin
      wset, 2
      if !d.x_size ne 600 or !d.y_size ne 600 then make_win = 1 else make_win = 0
    endif else make_win = 1
    if make_win eq 1 then window, 2, xsize = 600, ysize = 600
  endif
  
  cgplot, eta_vals, power_norm, position=pos_use[*,0], /ylog, xtitle = 'eta (1/MHz)', yrange = yrange, $
    title = 'Power & Weighted Power', xstyle=1, font = font, charsize = charsize
    
  cgplot, eta_vals, wt_power, /over, color='red', linestyle = 2, font = font, charsize = charsize
  ;  cgplot, eta_vals, wt_power, position=pos_use[*,1], /noerase, /ylog, xtitle = 'eta (1/MHz)', yrange = yrange, $
  ;    title = 'Weighted Power', xstyle=1, font = font, charsize = charsize
  
  ;cgplot, eta_vals, wt_power/power, position=pos_use[*,1], /noerase, xtitle = 'eta (1/MHz)', yrange = [0.5,1.5], $
  ;  title = 'Weighted Power/Power', xstyle=1, font = font, charsize = charsize
  
  ind_range = [n_freq[2]*(wh_peak/float(n_freq[2])-.05)>0, n_freq[2]*(wh_peak/float(n_freq[2])+.05)<(n_freq[2]-1)]
  xrange = minmax(eta_vals[ind_range[0]:ind_range[1]])
  cgplot, eta_vals, power_norm, position=pos_use[*,1], /noerase, xtitle = 'eta (1/MHz)', yrange = [.1,1.1]*yrange[1], $
    title = 'Power & Weighted Power', xstyle=1, font = font, charsize = charsize, xrange = xrange
  cgplot, eta_vals, wt_power, /over, color='red', linestyle = 2, font = font, charsize = charsize
  
  
  cgplot, eta_vals, power_nobp, position=pos_use[*,2], /noerase, /ylog, yrange = yrange, title= 'Power w/o bandpass', $
    xtitle = 'eta (1/MHz)', font = font, charsize = charsize, xstyle=1
  cgplot, eta_vals, wt_power, /over, color='red', linestyle = 2, font = font, charsize = charsize
  
  
  ;  cgplot, eta_vals, wt_power/power_nobp, position=pos_use[*,3], /noerase, xtitle = 'eta (1/MHz)', yrange = [0.5,1.5], $
  ;    title = 'Weighted Power/Power no bandpass', xstyle=1, font = font, charsize = charsize
  cgplot, eta_vals, power_nobp, position=pos_use[*,3], /noerase, xtitle = 'eta (1/MHz)', yrange = [.1,1.1]*yrange2[1], $
    title = 'Power w/o bandpass', xstyle=1, font = font, charsize = charsize, xrange = xrange
  cgplot, eta_vals, wt_power, /over, color='red', linestyle = 2, font = font, charsize = charsize
  
  
  
  ;cgplot, eta_vals, norm, position=pos_use[*,3], /noerase, color='blue', yrange = minmax(norm), title= 'Normalization', $
  ;  xtitle = 'eta (1/MHz)', font = font, charsize = charsize, xstyle=1
  
  
  
  
  if keyword_set(pub) then begin
    if png then begin
      if pdf then delete_ps_use = 0 else delete_ps_use = delete_ps
      cgps_close, /png, delete_ps = delete_ps_use, density = 600
    endif
    if pdf then begin
      if not png then cgps_close
      cgps2pdf, power_savefile, delete_ps=delete_ps
    endif
  endif
  
  print, 'peak ratios:'
  for i=0, 4 do print, ratio_names[i], round(peak_ratio[*, i]*1000)/1000., format = '(a20, 3f9.3)'
  print, ''
  print, 'integral ratios:'
  for i=0, 4 do print, ratio_names[i], round(power_ratio[*, i]*1000)/1000., format = '(a20, 3f9.3)'
  print, ''
  print, 'unflagged fraction', round(unflagged_frac*1000)/1000., format = '(a20, 3f9.3)'
  ;print, 'peak vec: ', peak_power[*,1], format = '(a20, 3e9.1)'
  ;print, 'peak input: ', peak_power[*,0], format = '(a20, 3e9.1)'
  print, 'sum vec wt/bp: ', vec_wt_sum[*,0]/n_freq, format = '(a20, 3f9.3)'
  print, 'sum (vec wt)^2/bp: ', vec_wt_sum[*,1]/n_freq, format = '(a20, 3f9.3)'
  
end
