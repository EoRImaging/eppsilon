function z_mpc, frequencies, hubble_param = hubble_param, f_delta = f_delta, even_freq = even_freq, $
    redshifts = redshifts, comov_dist_los = comov_dist_los, z_mpc_delta = z_mpc_delta
    
  n_freq = n_elements(frequencies)
  
  ;; check whether or not the frequencies are evenly spaced.
  freq_diff = frequencies - shift(frequencies, 1)
  freq_diff = freq_diff[1:*]
  
  z0_freq = 1420.40d ;; MHz
  redshifts = z0_freq/frequencies - 1d
  cosmology_measures, redshifts, comoving_dist_los = comov_dist_los, hubble_param = hubble_param
  
  comov_los_diff = comov_dist_los - shift(comov_dist_los, -1)
  comov_los_diff = comov_los_diff[0:n_elements(comov_dist_los)-2]
  
  if max(freq_diff-freq_diff[0]) gt 1e-12 then begin
    ;; frequencies are not evenly spaced, need to be careful about z_mpc_delta/mean
    even_freq = 0
    
    freq_diff_hist = histogram(freq_diff, binsize = min(freq_diff)*.1, locations=locs, reverse_indices = ri)
    if max(freq_diff_hist)/float(n_freq) lt .5 then begin
      message, 'frequency channel width distribution is strange, nominal channel width is unclear.'
    endif else begin
      peak_bin = (where(freq_diff_hist eq max(freq_diff_hist), count_peak))[0]
      if count_peak eq 1 then peak_diffs = freq_diff[ri[ri[peak_bin] : ri[peak_bin+1]-1]] $
      else message, 'frequency channel width distribution is strange, nominal channel width is unclear.'
      
      f_delta = mean(peak_diffs)
    endelse
    
    nominal_freqs = findgen(floor(((max(frequencies)-min(frequencies))/f_delta))+1)*f_delta + min(frequencies)
    nominal_z = z0_freq/nominal_freqs - 1
    cosmology_measures, nominal_z, comoving_dist_los = nominal_comov_dist_los
    nominal_comov_diffs = nominal_comov_dist_los - shift(nominal_comov_dist_los, -1)
    nominal_comov_diffs = nominal_comov_diffs[0:n_elements(nominal_comov_diffs)-2]
    
    z_mpc_delta = mean(nominal_comov_diffs)
    z_mpc_mean = mean(nominal_comov_dist_los)
    
  endif else begin
    even_freq = 1
    
    f_delta = freq_diff[0]
    z_mpc_delta = float(mean(comov_los_diff))
    z_mpc_mean = float(mean(comov_dist_los))
    n_kz = n_freq
    
  endelse
  
  return, z_mpc_mean
end