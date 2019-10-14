function create_ps_options, ps_options = ps_options, $
    ave_removal = ave_removal, wt_cutoffs = wt_cutoffs, $
    wt_measures = wt_measures, spec_window_type = spec_window_type, $
    no_spec_window = no_spec_window, allow_beam_approx = allow_beam_approx, $
    dft_z_use = dft_z_use, std_power = std_power, no_wtd_avg = no_wtd_avg, $
    inverse_covar_weight = inverse_covar_weight, return_new = return_new

  ;; if inverse covariance weighted don't use spectral window
  if keyword_set(inverse_covar_weight) then no_spec_window = 1

  update_tags = list()
  update_values = list()
  if n_elements(ps_options) eq 0 then begin
    ;; default to ave removal
    if n_elements(ave_removal) eq 0 then ave_removal = 1

    ;; default to erroring if beam_integral isn't present
    if n_elements(allow_beam_approx) eq 0 then allow_beam_approx = 0

    ;; default to not turning off spectral windowing
    if n_elements(no_spec_window) eq 0 then no_spec_window = 0

    ;; default to using the true z's if frequencies aren't evenly spaced
    ;; NOTE not sure this is the best choice but it is what is already in the code
    ;; NOTE this isn't currently used to update filenames. It probably should
    ;; but then we'd need to move the code that checks for even freq. sampling
    ;; upstream
    if n_elements(dft_z_use) eq 0 then dft_z_use = 'true'

    ;; default to Lomb-Scargle power calc
    if n_elements(std_power) eq 0 then std_power = 0

    ;; default to weighted averaging
    if n_elements(no_wtd_avg) eq 0 then no_wtd_avg = 0

    ;; default to inverse variance not covariance weighting
    if n_elements(inverse_covar_weight) eq 0 then inverse_covar_weight = 0

    ;; default to blackman-harris spectral window
    if not keyword_set(no_spec_window) then begin
      if n_elements(spec_window_type) eq 0 then spec_window_type = 'Blackman-Harris'
    endif else undefine, spec_window_type

    ;; density correction defaults
    if n_elements(wt_cutoffs) eq 0 then begin
      ;; default to wt_cutoffs = 1, wt_measures = 'min'
      wt_cutoffs = 1
      wt_measures = 'min'
    endif else if n_elements(wt_measures) eq 0 then begin
      print, 'wt_cutoffs is specified but wt_measures is not. Defaulting wt_measures to "min".'
      wt_measures = strarr(n_elements(wt_cutoffs)) + 'min'
    endif

    ps_options = {ave_removal: ave_removal, std_power: std_power, no_wtd_avg: no_wtd_avg, $
      inverse_covar_weight:inverse_covar_weight}

  endif else begin
    if n_elements(ave_removal) gt 0 then begin
      update_tags.add, 'ave_removal'
      update_values.add, ave_removal
    endif
    if n_elements(std_power) gt 0 then begin
      update_tags.add, 'std_power'
      update_values.add, std_power
    endif
    if n_elements(inverse_covar_weight) gt 0 then begin
      update_tags.add, 'inverse_covar_weight'
      update_values.add, inverse_covar_weight
    endif
  endelse

  if n_elements(wt_cutoffs) gt 0 then begin
    update_tags.add, 'wt_cutoffs'
    update_values.add, wt_cutoffs
  endif
  if n_elements(wt_measures) gt 0 then begin
    update_tags.add, 'wt_measures'
    update_values.add, wt_measures
  endif

  if n_elements(spec_window_type) gt 0 then begin
    update_tags.add, 'spec_window_type'
    update_values.add, spec_window_type
  endif

  if n_elements(dft_z_use) gt 0 then begin
    wh_val = where(['true', 'regular'] eq dft_z_use, count_whval)
    if count_whval eq 0 then begin
      message, 'dft_z_use must be one of "true" or "regular".'
    endif
    update_tags.add, 'dft_z_use'
    update_values.add, dft_z_use
  endif

  if n_elements(update_tags) gt 0 or keyword_set(return_new) then begin
    update_tags = update_tags.toarray()

    new_ps_options = update_option_struct(ps_options, update_tags, update_values, $
      return_new = return_new)

    return, new_ps_options
  endif else begin
    return, ps_options
  endelse
end
