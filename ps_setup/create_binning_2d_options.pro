function create_binning_2d_options, binning_2d_options = binning_2d_options, $
    no_kzero = no_kzero, log_kpar = log_kpar, log_kperp = log_kperp, $
    kpar_bin = kpar_bin, kperp_bin = kperp_bin, return_new = return_new

  update_tags = list()
  update_values = list()
  if n_elements(binning_2d_options) eq 0 then begin
    ;; default to including kpar=0 bin
    if n_elements(no_kzero) eq 0 then no_kzero = 0

    ;; default to linear bins
    if n_elements(log_kpar) eq 0 then log_kpar = 0
    if n_elements(log_kperp) eq 0 then log_kperp = 0

    binning_2d_options = {no_kzero: no_kzero, log_kpar: log_kpar, $
      log_kperp: log_kperp}
  endif else begin
    if n_elements(no_kzero) gt 0 then begin
      update_tags.add, 'no_kzero'
      update_values.add, no_kzero
    endif
    if n_elements(log_kpar) gt 0 then begin
      update_tags.add, 'log_kpar'
      update_values.add, log_kpar
    endif
    if n_elements(log_kperp) gt 0 then begin
      update_tags.add, 'log_kperp'
      update_values.add, log_kperp
    endif
  endelse

  if n_elements(kpar_bin) gt 0 then begin
    update_tags.add, 'kpar_bin'
    update_values.add, kpar_bin
  endif
  if n_elements(kperp_bin) gt 0 then begin
    update_tags.add, 'kperp_bin'
    update_values.add, kperp_bin
  endif

  if n_elements(update_tags) gt 0 or keyword_set(return_new) then begin
    update_tags = update_tags.toarray()

    new_binning_2d_options = update_option_struct(binning_2d_options, update_tags, $
      update_values, return_new = return_new)

    return, new_binning_2d_options
  endif else begin
    return, binning_2d_options
  endelse
end
