function create_binning_1d_options, binning_1d_options = binning_1d_options, $
    wedge_angles = wedge_angles, wedge_amps = wedge_amps, wedge_names = wedge_names, $
    coarse_harm_width = coarse_harm_width, coarse_harm0 = coarse_harm0, $
    log_k = log_k, k_bin = k_bin, kpar_range_1dave = kpar_range_1dave, $
    kperp_range_1dave = kperp_range_1dave, kperp_range_lambda_1dave = kperp_range_lambda_1dave, $
    kx_range_1dave = kx_range_1dave, kx_range_lambda_1dave = kx_range_lambda_1dave, $
    ky_range_1dave = ky_range_1dave, ky_range_lambda_1dave = ky_range_lambda_1dave, $
    kperp_range_lambda_kparpower = kperp_range_lambda_kparpower, $
    kpar_range_kperppower = kpar_range_kperppower, return_new = return_new

  update_tags = list()
  update_values = list()
  if n_elements(binning_1d_options) eq 0 then begin
    ;; default to linear bins
    if n_elements(log_k) eq 0 then log_k = 0

    binning_1d_options = {log_k: log_k}
  endif else begin
    if n_elements(log_k) gt 0 then begin
      update_tags.add, 'log_k'
      update_values.add, log_k
    endif
  endelse

  if n_elements(kperp_range_1dave) gt 1 and n_elements(kperp_range_lambda_1dave) gt 1 then begin
    message, 'both kperp_range_1dave and kperp_range_lambda_1dave cannot be set'
  endif
  if n_elements(kx_range_1dave) gt 1 and n_elements(kx_range_lambda_1dave) gt 1 then begin
    message, 'both kx_range_1dave and kx_range_lambda_1dave cannot be set'
  endif
  if n_elements(ky_range_1dave) gt 1 and n_elements(ky_range_lambda_1dave) gt 1 then begin
    message, 'both ky_range_1dave and ky_range_lambda_1dave cannot be set'
  endif

  remove_tags = list()
  if n_elements(wedge_angles) gt 0 then begin
    update_tags.add, 'wedge_angles'
    update_values.add, wedge_angles
  endif
  if n_elements(wedge_amps) gt 0 then begin
    update_tags.add, 'wedge_amps'
    update_values.add, wedge_amps
  endif
  if n_elements(wedge_names) gt 0 then begin
    update_tags.add, 'wedge_names'
    update_values.add, wedge_names
  endif
  if n_elements(coarse_harm_width) gt 0 then begin
    update_tags.add, 'coarse_harm_width'
    update_values.add, coarse_harm_width
  endif
  if n_elements(coarse_harm0) gt 0 then begin
    update_tags.add, 'coarse_harm0'
    update_values.add, coarse_harm0
  endif

  if n_elements(k_bin) gt 0 then begin
    update_tags.add, 'k_bin'
    update_values.add, k_bin
  endif

  if n_elements(kpar_range_1dave) gt 0 then begin
    if n_elements(kpar_range_1dave) ne 2 then begin
      message, 'kpar_range_1dave must be a 2 element vector'
    endif
    update_tags.add, 'kpar_range_1dave'
    update_values.add, kpar_range_1dave
  endif

  if n_elements(kperp_range_1dave) gt 0 then begin
    if n_elements(kperp_range_1dave) ne 2 then begin
      message, 'kperp_range_1dave must be a 2 element vector'
    endif
    update_tags.add, 'kperp_range_1dave'
    update_values.add, kperp_range_1dave
  endif
  if n_elements(kperp_range_lambda_1dave) gt 0 then begin
    if n_elements(kperp_range_lambda_1dave) ne 2 then begin
      message, 'kperp_range_lambda_1dave must be a 2 element vector'
    endif
    update_tags.add, 'kperp_range_lambda_1dave'
    update_values.add, kperp_range_lambda_1dave
    if tag_exist(binning_1d_options, 'kperp_range_1dave') then begin
      remove_tags.add, 'kperp_range_1dave'
    endif
  endif

  if n_elements(kx_range_1dave) gt 0 then begin
    if n_elements(kx_range_1dave) ne 2 then begin
      message, 'kx_range_1dave must be a 2 element vector'
    endif
    if min(kx_range_1dave) lt 0 then begin
      message, 'kx_range_1dave values must be positive (interpreted as absolute values)'
    endif
    update_tags.add, 'kx_range_1dave'
    update_values.add, kx_range_1dave
  endif
  if n_elements(kx_range_lambda_1dave) gt 0 then begin
    if n_elements(kx_range_lambda_1dave) ne 2 then begin
      message, 'kx_range_lambda_1dave must be a 2 element vector'
    endif
    if min(kx_range_lambda_1dave) lt 0 then begin
      message, 'kx_range_lambda_1dave values must be positive (interpreted as absolute values)'
    endif
    update_tags.add, 'kx_range_lambda_1dave'
    update_values.add, kx_range_lambda_1dave
    if tag_exist(binning_1d_options, 'kx_range_1dave') then begin
      remove_tags.add, 'kx_range_1dave'
    endif
  endif

  if n_elements(ky_range_1dave) gt 0 then begin
    if n_elements(ky_range_1dave) ne 2 then begin
      message, 'ky_range_1dave must be a 2 element vector'
    endif
    if min(ky_range_1dave) lt 0 then begin
      message, 'ky_range_1dave values must be positive (interpreted as absolute values)'
    endif
    update_tags.add, 'ky_range_1dave'
    update_values.add, ky_range_1dave
  endif
  if n_elements(ky_range_lambda_1dave) gt 0 then begin
    if min(ky_range_lambda_1dave) lt 0 then begin
      message, 'ky_range_lambda_1dave values must be positive (interpreted as absolute values)'
    endif
    update_tags.add, 'ky_range_lambda_1dave'
    update_values.add, ky_range_lambda_1dave
    if tag_exist(binning_1d_options, 'ky_range_1dave') then begin
      remove_tags.add, 'ky_range_1dave'
    endif
  endif

  if n_elements(kperp_range_lambda_kparpower) gt 0 then begin
    if n_elements(kperp_range_lambda_kparpower) ne 2 then begin
      message, 'kperp_range_lambda_kparpower must be a 2 element vector'
    endif
    update_tags.add, 'kperp_range_lambda_kparpower'
    update_values.add, kperp_range_lambda_kparpower
  endif

  if n_elements(kpar_range_kperppower) gt 0 then begin
    if n_elements(kpar_range_kperppower) ne 2 then begin
      message, 'kpar_range_kperppower must be a 2 element vector'
    endif
    update_tags.add, 'kpar_range_kperppower'
    update_values.add, kpar_range_kperppower
  endif

  if isa(remove_tags) gt 0 then begin
    if (n_elements(remove_tags) eq 0) then begin
      ;; undefine doesn't work on lists for some reason.
      remove_tags = 0
      undefine, remove_tags
    endif else begin
      remove_tags = remove_tags.toarray()
    endelse
  endif

  if isa(update_tags) gt 0 then begin
    if (n_elements(update_tags) eq 0) then begin
      ;; undefine doesn't work on lists for some reason.
      update_tags = 0
      undefine, update_tags
    endif else begin
      update_tags = update_tags.toarray()
    endelse
  endif

  if isa(update_values) gt 0 then begin
    if (n_elements(update_values) eq 0) then begin
      ;; undefine doesn't work on lists for some reason.
      update_values = 0
      undefine, update_values
    endif
  endif

  if n_elements(remove_tags) gt 0 or n_elements(update_tags) gt 0 or $
    keyword_set(return_new) then begin

    new_binning_1d_options = update_option_struct(binning_1d_options, update_tags, $
      update_values, remove_tags = remove_tags, return_new = return_new)

    return, new_binning_1d_options
  endif else begin
    return, binning_1d_options
  endelse
end
