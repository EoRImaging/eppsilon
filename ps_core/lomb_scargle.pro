pro lomb_scargle, n_kx=n_kx, n_ky=n_ky, n_kz=n_kz, n_freq=n_freq, even_freq=even_freq, z_mpc_delta=z_mpc_delta, $
  comov_dist_los=comov_dist_los, kz_mpc=kz_mpc, sum_sigma2=sum_sigma2, contrib_n_freq=contrib_n_freq, n_val=n_val, $
  wh_pos=wh_pos, wh_neg=wh_neg, sigma2_1=sigma2_1, sigma2_2=sigma2_2, cos_theta=cos_theta, sin_theta=sin_theta

  ;; for new power calc, need cos2, sin2, cos*sin transforms
  n_unfolded_kz = n_elements(kz_mpc)
  covar_cos = fltarr(n_kx, n_ky, n_freq)
  covar_sin = fltarr(n_kx, n_ky, n_freq)
  covar_cross = fltarr(n_kx, n_ky, n_freq)

  ;; comov_dist_los goes from large to small z
  if keyword_set(even_freq) then begin
    z_relative = dindgen(n_freq)*z_mpc_delta
  endif else z_relative = reverse(comov_dist_los-min(comov_dist_los))
  freq_kz_arr = rebin(reform(kz_mpc, 1, n_unfolded_kz), n_freq, n_unfolded_kz) * $
    rebin(z_relative, n_freq, n_unfolded_kz)

  cos_arr = cos(freq_kz_arr)
  sin_arr = sin(freq_kz_arr)

  sum_sigma2 = reform(sum_sigma2, n_kx*n_ky, n_freq)
  ;; doing 2 FTs so need 2 factors of z_mpc_delta.
  ;; No multiplication by N b/c don't need to fix IDL FFT
  covar_cos = matrix_multiply(sum_sigma2, cos_arr^2d) * (z_mpc_delta)^2.
  covar_sin = matrix_multiply(sum_sigma2, sin_arr^2d) * (z_mpc_delta)^2.
  covar_cross = matrix_multiply(sum_sigma2, cos_arr*sin_arr) * (z_mpc_delta)^2.

  wh_0f = where(contrib_n_freq eq 0, count_0f)
  if count_0f gt 0 then begin
    covar_cos[wh_0f, *] = 0
    covar_sin[wh_0f, *] = 0
    covar_cross[wh_0f, *] = 0
  endif

  ;; reform to get back to n_kx, n_ky, n_kz dimensions
  covar_cos = reform(covar_cos, n_kx, n_ky, n_unfolded_kz)
  covar_sin = reform(covar_sin, n_kx, n_ky, n_unfolded_kz)
  covar_cross = reform(covar_cross, n_kx, n_ky, n_unfolded_kz)

  ;; drop pixels with less than 1/3 of the frequencies
  wh_fewfreq = where(contrib_n_freq lt ceil(n_freq/3d), count_fewfreq)
  if count_fewfreq gt 0 then begin
    mask_fewfreq = contrib_n_freq * 0 + 1
    mask_fewfreq[wh_fewfreq] = 0
    mask_fewfreq = rebin(temporary(mask_fewfreq), n_kx, n_ky, n_unfolded_kz)
    covar_cos = temporary(covar_cos) * mask_fewfreq
    covar_sin = temporary(covar_sin) * mask_fewfreq
    covar_cross = temporary(covar_cross) * mask_fewfreq
    undefine, mask_fewfreq
  endif
  
  undefine, sum_sigma2, freq_kz_arr, cos_arr, sin_arr

  ;; get rotation angle to diagonalize covariance block
  theta = atan(2.*covar_cross, covar_cos - covar_sin)/2.

  cos_theta = cos(theta)
  sin_theta = sin(theta)
  undefine, theta

  ;; rotate errors to get orthogonal distribution
  sigma2_cos = covar_cos*cos_theta^2. + 2.*covar_cross*cos_theta*sin_theta + $
    covar_sin*sin_theta^2.
  sigma2_sin = covar_cos*sin_theta^2. - 2.*covar_cross*cos_theta*sin_theta + $
    covar_sin*cos_theta^2.
    
  ;; fold errors in kz 
  a1_0 = sigma2_cos[*,*,where(n_val eq 0)]
  a1_n = (sigma2_cos[*, *, wh_pos] + sigma2_cos[*, *, reverse(wh_neg)])/2.
  sigma2_1 = dblarr(n_kx, n_ky, n_kz)
  sigma2_1[*, *, 0] = a1_0
  sigma2_1[*, *, 1:n_kz-1] = a1_n

  b1_n = (sigma2_sin[*, *, wh_pos] + sigma2_sin[*, *, reverse(wh_neg)])/2.
  sigma2_2 = dblarr(n_kx, n_ky, n_kz)
  sigma2_2[*, *, 1:n_kz-1] = b1_n
    
  undefine, covar_cos, covar_sin, covar_cross, a1_0, a1_n, b1_0, b1_n

  return

end
