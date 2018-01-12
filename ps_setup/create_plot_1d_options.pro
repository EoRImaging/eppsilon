function create_plot_1d_options, plot_1d_options = plot_1d_options, $
    range_1d = range_1d, plot_1d_delta = plot_1d_delta, $
    plot_1d_error_bars = plot_1d_error_bars, plot_1d_nsigma = plot_1d_nsigma, $
    plot_eor_1d = plot_eor_1d, plot_flat_1d = plot_flat_1d, $
    no_text_1d = no_text_1d, return_new = return_new

  update_tags = list()
  update_values = list()
  if n_elements(plot_1d_options) eq 0 then begin
    ;; default to power not delta
    if n_elements(plot_1d_delta) eq 0 then plot_1d_delta = 0

    ;; default to plotting thermal nose line, not error bars
    if n_elements(plot_1d_error_bars) eq 0 then plot_1d_error_bars=0

    ;; default to not plotting eor or flat power line
    if n_elements(plot_eor_1d) eq 0 then plot_eor_1d=0
    if n_elements(plot_flat_1d) eq 0 then plot_flat_1d=0

    ;; default to printing text on 1d plots
    if n_elements(no_text_1d) eq 0 then no_text_1d=0

    plot_1d_options = {plot_1d_delta: plot_1d_delta, plot_1d_error_bars: plot_1d_error_bars, $
      plot_eor_1d: plot_eor_1d, plot_flat_1d: plot_flat_1d, no_text_1d:no_text_1d}
  endif else begin
    if n_elements(plot_1d_delta) gt 0 then begin
      update_tags.add, 'plot_1d_delta'
      update_values.add, plot_1d_delta
    endif
    if n_elements(plot_1d_error_bars) gt 0 then begin
      update_tags.add, 'plot_1d_error_bars'
      update_values.add, plot_1d_error_bars
    endif
    if n_elements(plot_eor_1d) gt 0 then begin
      update_tags.add, 'plot_eor_1d'
      update_values.add, plot_eor_1d
    endif
    if n_elements(plot_flat_1d) gt 0 then begin
      update_tags.add, 'plot_flat_1d'
      update_values.add, plot_flat_1d
    endif
    if n_elements(no_text_1d) gt 0 then begin
      update_tags.add, 'no_text_1d'
      update_values.add, no_text_1d
    endif
  endelse

  if n_elements(range_1d) gt 0 then begin
    if n_elements(range_1d) ne 2 then begin
      message, 'range_1d must be a 2 element vector'
    endif
    update_tags.add, 'range_1d'
    update_values.add, range_1d
  endif
  if n_elements(plot_1d_nsigma) gt 0 then begin
    update_tags.add, 'plot_1d_nsigma'
    update_values.add, plot_1d_nsigma
  endif

  if n_elements(update_tags) gt 0 or keyword_set(return_new) then begin

    update_tags = update_tags.toarray()

    if n_elements(update_tags) gt 0 or keyword_set(return_new) then begin
      new_plot_1d_options = update_option_struct(plot_1d_options, update_tags, update_values, $
        return_new = return_new)
    endif

    return, new_plot_1d_options
  endif else begin
    return, plot_1d_options
  endelse
end
