pro save_1D_text, textfile, k_edges, power, weights, noise_expval, hubble_param, $
  noise, theory_delta2 = theory_delta2, sim_noise_power = sim_noise_power, $
  sim_noise_diff = sim_noise_diff, nfiles = nfiles, hinv = hinv

  compile_opt idl2, logical_predicate

  nrows = n_elements(k_edges)
  ncol = 4
  if n_elements(theory_delta2) gt 0 then ncol += 1
  if nfiles eq 2 then begin
    ncol += 1
    if n_elements(sim_noise_power) gt 0 then ncol += 2
  endif else begin
    if n_elements(sim_noise_power) gt 0 then ncol += 1
  endelse
  data = fltarr(ncol, nrows)

  sigma_vals = sqrt(1./weights)
  wh_wt0 = where(weights eq 0, count_wh_wt0)
  if count_wh_wt0 gt 0 then sigma_vals[wh_wt0]= !values.f_infinity

  if ~keyword_set(hinv) then begin
    data[0, *] = k_edges
    data[1, *] = [0, power]
    data[2, *] = [0, sigma_vals]
    data[3, *] = [0, noise_expval]
    header = ['k bin max (Mpc^-1)', 'power (mK^2 Mpc^3)', 'sigma (mK^2 Mpc^3)', 'expected noise (mK^2 Mpc^3)']
    next_col_ind = 4
    if n_elements(theory_delta2) gt 0 then begin
      data[next_col_ind, *] = [0, theory_delta2]
      header = [header, 'delta^2 (mK^2)']
      next_col_ind += 1
    endif
    
    if nfiles eq 2 then begin
      data[next_col_ind, *] = [0, noise]
      header = [header, 'observed noise (mK^2 Mpc^3)']
      next_col_ind += 1
      if n_elements(sim_noise_power) gt 0 then begin
        data[next_col_ind, *] = [0, sim_noise_power]
        header = [header, 'sim noise power (mK^2 Mpc^3)']
        next_col_ind += 1

        data[next_col_ind, *] = [0, sim_noise_diff]
        header = [header, 'sim noise diff (mK^2 Mpc^3)']
        next_col_ind += 1
      endif
    endif else begin
      if n_elements(sim_noise_power) gt 0 then begin
        data[next_col_ind, *] = [0, sim_noise_power]
        header = [header, 'sim noise power (mK^2 Mpc^3)']
        next_col_ind += 1
      endif
    endelse
  endif else begin
    data[0, *] = k_edges / hubble_param
    data[1, *] = [0, power] * (hubble_param)^3d
    data[2, *] = [0, sigma_vals] * (hubble_param)^3d
    data[3, *] = [0, noise_expval] * (hubble_param)^3d
    header = ['k bin max (h Mpc^-1)', 'power (mK^2 h^-3 Mpc^3)', 'sigma (mK^2 h^-3 Mpc^3)', 'expected noise (mK^2 h^-3 Mpc^3)']
    next_col_ind = 4
    if n_elements(theory_delta2) gt 0 then begin
      data[next_col_ind, *] = [0, theory_delta2]
      header = [header, 'delta^2 (mK^2)']
      next_col_ind += 1
    endif
    
    if nfiles eq 2 then begin
      data[next_col_ind, *] = [0, noise] * (hubble_param)^3d
      header = [header, 'observed noise (mK^2 h^-3 Mpc^3)']
      next_col_ind += 1
      if n_elements(sim_noise_power) gt 0 then begin
        data[next_col_ind, *] = [0, sim_noise_power] * (hubble_param)^3d
        header = [header, 'sim noise power (mK^2 h^-3 Mpc^3)']
        next_col_ind += 1

        data[next_col_ind, *] = [0, sim_noise_diff] * (hubble_param)^3d
        header = [header, 'sim noise diff (mK^2 h^-3 Mpc^3)']
        next_col_ind += 1
      endif
    endif else begin
      if n_elements(sim_noise_power) gt 0 then begin
        data[next_col_ind, *] = [0, sim_noise_power] * (hubble_param)^3d
        header = [header, 'sim noise power (mK^2 h^-3 Mpc^3)']
        next_col_ind += 1
      endif
    endelse
  endelse

  ;; add hash to indicate header line for Danny
  header[0] = '# ' + header[0]
  data_use = Strarr(ncol,nrows+1)
  data_use[*,0]=header
  data_use[*,1:*]=(data)

  delimiter=String(9B)
  format_code=String(format='("(",A,"(A,",A,A,A,"))")',Strn(ncol),'"',delimiter,'"')

  openw,unit,textfile,/Get_LUN
  printf,unit,format=format_code,data_use
  free_lun,unit


end
