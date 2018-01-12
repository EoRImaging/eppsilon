function create_uvf_options, uvf_options = uvf_options, $
    image_window_name = image_window_name, image_window_frac_size = image_window_frac_size, $
    delta_uv_lambda = delta_uv_lambda, max_uv_lambda = max_uv_lambda, $
    uv_avg = uv_avg, uv_img_clip = uv_img_clip, $
    dft_fchunk = dft_fchunk, no_dft_progress = no_dft_progress, return_new = return_new

  update_tags = list()
  update_values = list()
  if n_elements(uvf_options) eq 0 then begin
    ;; default to giving dft progress reports
    if n_elements(no_dft_progress) eq 0 then no_dft_progress = 0

    uvf_options = {no_dft_progress: no_dft_progress}
  endif else begin
    if n_elements(no_dft_progress) gt 0 then begin
      update_tags.add, 'no_dft_progress'
      update_values.add, no_dft_progress
    endif
  endelse

  if n_elements(image_window_name) gt 0 then begin
    update_tags.add, 'image_window_name'
    update_values.add, image_window_name
  endif
  if n_elements(image_window_frac_size) gt 0 then begin
    if image_window_frac_size gt 1 or image_window_frac_size lt 0 then begin
      print, 'image_window_frac_size must be a value between 0 and 1, using default values.'
      undefine, image_window_frac_size
    endif else begin
      update_tags.add, 'image_window_frac_size'
      update_values.add, image_window_frac_size
    endelse
  endif

  if n_elements(delta_uv_lambda) gt 0 then begin
    update_tags.add, 'delta_uv_lambda'
    update_values.add, delta_uv_lambda
  endif
  if n_elements(max_uv_lambda) gt 0 then begin
    update_tags.add, 'max_uv_lambda'
    update_values.add, max_uv_lambda
  endif

  if n_elements(uv_avg) gt 0 then begin
    update_tags.add, 'uv_avg'
    update_values.add, uv_avg
  endif
  if n_elements(uv_img_clip) gt 0 then begin
    update_tags.add, 'uv_img_clip'
    update_values.add, uv_img_clip
  endif

  if n_elements(dft_fchunk) gt 0 then begin
    update_tags.add, 'dft_fchunk'
    update_values.add, dft_fchunk
  endif

  update_tags = update_tags.toarray()

  if n_elements(update_tags) gt 0 or keyword_set(return_new) then begin
    uvf_options = update_option_struct(uvf_options, update_tags, update_values, $
      return_new = return_new)
  endif

  return, uvf_options
end
