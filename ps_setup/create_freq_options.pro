function create_freq_options, $
    freq_options = freq_options, $
    freq_ch_range = freq_ch_range, $
    freq_flags = freq_flags, $
    freq_mask = freq_mask, $
    freq_flag_name = freq_flag_name, $
    freq_flag_repeat = freq_flag_repeat, $
    freq_avg_factor = freq_avg_factor, $
    force_even_freqs = force_even_freqs, $
    freq_avg_bins = freq_avg_bins, $
    freq_bin_name = freq_bin_name, $
    frequencies = frequencies, $
    n_vis_freq = n_vis_freq, $
    beam_int = beam_int, $
    new_freq_mask = new_freq_mask

  if keyword_set(freq_flag_repeat) and n_elements(freq_flags) gt 0 then begin
    if size(freq_flag_repeat, /type) ne 2 and n_elements(freq_flag_repeat) ne 1 then begin
      message, "freq_flag_repeat must be a scalar integer specifying the number of " $
      + "times to repeat the freq_flags for the full flag array."
    endif

    ;; remake freq flag array
    n_inital = n_elements(freq_flags)
    n_total = n_inital * freq_flag_repeat
    freq_flags = reform(rebin(reform(freq_flags, n_inital, 1), n_inital, freq_flag_repeat, /sample), n_total)
  endif

  update_tags = list()
  update_values = list()

  if n_elements(freq_options) eq 0 then begin
    ;; default to no averaging
    if n_elements(freq_avg_factor) eq 0 then begin
      freq_avg_factor = 1
    endif

    ;; default to not forcing even frequency spacing after averaging
    if n_elements(force_even_freqs) eq 0 then begin
      force_even_freqs = 0
    endif

    freq_options = {freq_avg_factor: freq_avg_factor, force_even_freqs: force_even_freqs}
  endif else begin
    if n_elements(freq_avg_factor) gt 0 then begin
      if size(freq_avg_factor, /type) ne 2 then begin
        message, "freq_avg_factor must be an integer."
      endif
      update_tags.add, 'freq_avg_factor'
      update_values.add, freq_avg_factor
    endif
    if n_elements(force_even_freqs) gt 0 then begin
      update_tags.add, 'force_even_freqs'
      update_values.add, force_even_freqs
    endif

  endelse

  if n_elements(freq_ch_range) gt 0 then begin
    update_tags.add, 'freq_ch_range'
    update_values.add, freq_ch_range
  endif

  if n_elements(freq_flags) gt 0 then begin
    update_tags.add, 'freq_flags'
    update_values.add, freq_flags
  endif

  if n_elements(freq_mask) gt 0 then begin
    update_tags.add, 'freq_mask'
    update_values.add, freq_mask
  endif

  if n_elements(freq_flag_name) gt 0 then begin
    update_tags.add, 'freq_flag_name'
    update_values.add, freq_flag_name
  endif

  if n_elements(freq_avg_bins) gt 0 then begin
    if freq_options.freq_avg_factor gt 1 then begin
      message, "freq_avg_factor and freq_avg_bins cannot both be set"
    endif else if freq_options.force_even_freqs gt 0 then begin
      message, "cannot use force_even_freqs with freq_avg_bins"
    endif else begin
      update_tags.add, 'freq_avg_bins'
      update_values.add, freq_avg_bins
    endelse
  endif

  if n_elements(freq_bin_name) gt 0 then begin
    update_tags.add, 'freq_bin_name'
    update_values.add, freq_bin_name
  endif

  if n_elements(frequencies) gt 0 then begin
    update_tags.add, 'frequencies'
    update_values.add, frequencies
  endif

  if n_elements(n_vis_freq) gt 0 then begin
    update_tags.add, 'n_vis_freq'
    update_values.add, n_vis_freq
  endif

  if n_elements(beam_int) gt 0 then begin
    update_tags.add, 'beam_int'
    update_values.add, beam_int
  endif

  if n_elements(new_freq_mask) gt 0 then begin
    update_tags.add, 'new_freq_mask'
    update_values.add, new_freq_mask
  endif

  if n_elements(update_tags) gt 0 or keyword_set(return_new) then begin
    update_tags = update_tags.toarray()

    new_freq_options = update_option_struct($
        freq_options, update_tags, update_values, return_new = return_new $
    )

    return, new_freq_options
  endif else begin
    return, freq_options
  endelse
end
