function create_refresh_options, refresh_options = refresh_options, $
    refresh_dft = refresh_dft, refresh_weight_dft = refresh_weight_dft, $
    refresh_beam = refresh_beam, refresh_kcube = refresh_kcube, refresh_ps = refresh_ps, $
    refresh_binning = refresh_binning, refresh_info = refresh_info, $
    return_new = return_new

  ;; cascade refreshes
  if keyword_set(refresh_dft) then refresh_weight_dft = 1
  if keyword_set(refresh_dft) then refresh_beam = 1
  if keyword_set(refresh_dft) or keyword_set(refresh_weight_dft) or $
    keyword_set(refresh_beam) then refresh_kcube = 1
  if keyword_set(refresh_kcube) then refresh_ps = 1
  if keyword_set(refresh_ps) then refresh_binning = 1

  if n_elements(refresh_options) eq 0 then begin
    ;; default to no refreshes
    if n_elements(refresh_dft) eq 0 then refresh_dft = 0
    if n_elements(refresh_weight_dft) eq 0 then refresh_weight_dft = 0
    if n_elements(refresh_beam) eq 0 then refresh_beam = 0
    if n_elements(refresh_kcube) eq 0 then refresh_kcube = 0
    if n_elements(refresh_ps) eq 0 then refresh_ps = 0
    if n_elements(refresh_binning) eq 0 then refresh_binning = 0
    if n_elements(refresh_info) eq 0 then refresh_info = 0

    refresh_options = {refresh_dft: refresh_dft, refresh_weight_dft: refresh_weight_dft, $
      refresh_beam: refresh_beam, refresh_kcube: refresh_kcube, refresh_ps: refresh_ps, $
      refresh_binning:refresh_binning, refresh_info:refresh_info}

    return, refresh_options
  endif else begin
    update_tags = list()
    update_values = list()

    if n_elements(refresh_dft) gt 0 then begin
      update_tags.add, 'refresh_dft'
      update_values.add, refresh_dft
    endif
    if n_elements(refresh_weight_dft) gt 0 then begin
      update_tags.add, 'refresh_weight_dft'
      update_values.add, refresh_weight_dft
    endif
    if n_elements(refresh_beam) gt 0 then begin
      update_tags.add, 'refresh_beam'
      update_values.add, refresh_beam
    endif
    if n_elements(refresh_kcube) gt 0 then begin
      update_tags.add, 'refresh_kcube'
      update_values.add, refresh_kcube
    endif
    if n_elements(refresh_ps) gt 0 then begin
      update_tags.add, 'refresh_ps'
      update_values.add, refresh_ps
    endif
    if n_elements(refresh_binning) gt 0 then begin
      update_tags.add, 'refresh_binning'
      update_values.add, refresh_binning
    endif
    if n_elements(refresh_info) gt 0 then begin
      update_tags.add, 'refresh_info'
      update_values.add, refresh_info
    endif

    if n_elements(update_tags) eq 0 and not keyword_set(return_new) then begin
        message, 'no update tags and return_new not set.'
    endif

    update_tags = update_tags.toarray()

    if n_elements(update_tags) gt 0 or keyword_set(return_new) then begin
      new_refresh_options = update_option_struct(refresh_options, update_tags, $
        update_values, return_new = return_new)
    endif

    return, new_refresh_options
  endelse

end
