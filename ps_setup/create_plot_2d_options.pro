function create_plot_2d_options, plot_2d_options = plot_2d_options,$
    plot_wedge_line = plot_wedge_line, $
    kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, $
    kperp_plot_range = kperp_plot_range, kperp_lambda_plot_range = kperp_lambda_plot_range, $
    kpar_plot_range = kpar_plot_range, baseline_axis = baseline_axis, $
    delay_axis = delay_axis, cable_length_axis = cable_length_axis, $
    data_range = data_range, sigma_range = sigma_range, nev_range = nev_range, $
    snr_range = snr_range, noise_range = noise_range, nnr_range = nnr_range, $
    slice_range = slice_range, color_type = color_type, return_new = return_new

  if keyword_set(delay_axis) and keyword_set(cable_length_axis) then begin
    message, 'only one of delay_axis and cable_length_axis can be set'
  endif

  update_tags = list()
  update_values = list()
  if n_elements(plot_2d_options) eq 0 then begin
    ;; default to plot wedge line
    if n_elements(plot_wedge_line) eq 0 then plot_wedge_line=1

    ;; default to log axes
    if n_elements(kperp_linear_axis) eq 0 then kperp_linear_axis=0
    if n_elements(kpar_linear_axis) eq 0 then kpar_linear_axis=0

    ;; default to including baseline axis & delay axis
    if n_elements(baseline_axis) eq 0 then baseline_axis = 1
    if n_elements(cable_length_axis) eq 0 then cable_length_axis = 0
    if n_elements(delay_axis) eq 0 then begin
      if keyword_set(cable_length_axis) then begin
        delay_axis = 0
      endif else begin
        delay_axis = 1
      endelse
    endif

    plot_2d_options = {plot_wedge_line: plot_wedge_line, $
      kperp_linear_axis: kperp_linear_axis, kpar_linear_axis: kpar_linear_axis, $
      baseline_axis: baseline_axis, delay_axis: delay_axis, $
      cable_length_axis: cable_length_axis}
  endif else begin
    if n_elements(delay_axis) gt 0 then begin
      if keyword_set(delay_axis) then cable_length_axis = 0
    endif
    if n_elements(cable_length_axis) gt 0 then begin
      if keyword_set(cable_length_axis) then delay_axis = 0
    endif

    if n_elements(plot_wedge_line) gt 0 then begin
      update_tags.add, 'plot_wedge_line'
      update_values.add, plot_wedge_line
    endif
    if n_elements(kperp_linear_axis) gt 0 then begin
      update_tags.add, 'kperp_linear_axis'
      update_values.add, kperp_linear_axis
    endif
    if n_elements(kpar_linear_axis) gt 0 then begin
      update_tags.add, 'kpar_linear_axis'
      update_values.add, kpar_linear_axis
    endif
    if n_elements(baseline_axis) gt 0 then begin
      update_tags.add, 'baseline_axis'
      update_values.add, baseline_axis
    endif
    if n_elements(delay_axis) gt 0 then begin
      update_tags.add, 'delay_axis'
      update_values.add, delay_axis
    endif
    if n_elements(cable_length_axis) gt 0 then begin
      update_tags.add, 'cable_length_axis'
      update_values.add, cable_length_axis
    endif
  endelse

  if n_elements(kperp_plot_range) gt 0 then begin
    if n_elements(kperp_plot_range) ne 2 then begin
      message, 'kperp_plot_range must be a 2 element vector'
    endif
    update_tags.add, 'kperp_plot_range'
    update_values.add, kperp_plot_range
  endif
  if n_elements(kperp_lambda_plot_range) gt 0 then begin
    if n_elements(kperp_lambda_plot_range) ne 2 then begin
      message, 'kperp_lambda_plot_range must be a 2 element vector'
    endif
    update_tags.add, 'kperp_lambda_plot_range'
    update_values.add, kperp_lambda_plot_range
  endif
  if n_elements(kpar_plot_range) gt 0 then begin
    if n_elements(kpar_plot_range) ne 2 then begin
      message, 'kpar_plot_range must be a 2 element vector'
    endif
    update_tags.add, 'kpar_plot_range'
    update_values.add, kpar_plot_range
  endif
  if n_elements(data_range) gt 0 then begin
    if n_elements(data_range) ne 2 then begin
      message, 'data_range must be a 2 element vector'
    endif
    update_tags.add, 'data_range'
    update_values.add, data_range
  endif
  if n_elements(sigma_range) gt 0 then begin
    if n_elements(sigma_range) ne 2 then begin
      message, 'sigma_range must be a 2 element vector'
    endif
    update_tags.add, 'sigma_range'
    update_values.add, sigma_range
  endif
  if n_elements(nev_range) gt 0 then begin
    if n_elements(nev_range) ne 2 then begin
      message, 'nev_range must be a 2 element vector'
    endif
    update_tags.add, 'nev_range'
    update_values.add, nev_range
  endif
  if n_elements(snr_range) gt 0 then begin
    if n_elements(snr_range) ne 2 then begin
      message, 'snr_range must be a 2 element vector'
    endif
    update_tags.add, 'snr_range'
    update_values.add, snr_range
  endif
  if n_elements(noise_range) gt 0 then begin
    if n_elements(noise_range) ne 2 then begin
      message, 'noise_range must be a 2 element vector'
    endif
    update_tags.add, 'noise_range'
    update_values.add, noise_range
  endif
  if n_elements(nnr_range) gt 0 then begin
    if n_elements(nnr_range) ne 2 then begin
      message, 'nnr_range must be a 2 element vector'
    endif
    update_tags.add, 'nnr_range'
    update_values.add, nnr_range
  endif
  if n_elements(slice_range) gt 0 then begin
    if n_elements(slice_range) ne 2 then begin
      message, 'slice_range must be a 2 element vector'
    endif
    update_tags.add, 'slice_range'
    update_values.add, slice_range
  endif

  if n_elements(color_type) gt 0 then begin
    color_type_enum = ['integer', 'log', 'linear']
    wh_color_type = where(color_type_enum eq color_type, count_type)
    if count_type eq 0 then begin
      message, 'color_type must be one of: ' + color_type_enum
    endif
    update_tags.add, 'color_type'
    update_values.add, color_type
  endif

  if n_elements(update_tags) gt 0 or keyword_set(return_new) then begin
    update_tags = update_tags.toarray()

    new_plot_2d_options = update_option_struct(plot_2d_options, update_tags, update_values, $
      return_new = return_new)

    return, new_plot_2d_options
  endif else begin
    return, plot_2d_options
  endelse
end
