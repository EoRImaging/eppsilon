function create_file_tags, freq_ch_range = freq_ch_range, freq_flags = freq_flags, $
    freq_flag_name = freq_flag_name, uvf_options = uvf_options, ps_options = ps_options

  window_type_list = ['Hann', 'Hamming', 'Blackman', 'Nutall', 'Blackman-Nutall', $
                      'Blackman-Harris', 'Blackman-Harris^2', 'Tukey', 'None']
  win_tag_list = ['hann', 'ham', 'blm', 'ntl', 'bn', 'bh', 'bh2', 'tk', '']
  if tag_exist(ps_options, 'spec_window_type') then begin
    wh_type = where(strlowcase(window_type_list) eq strlowcase(ps_options.spec_window_type), count_type)
    if count_type eq 0 then begin
      wh_type = where(strlowcase(win_tag_list) eq strlowcase(ps_options.spec_window_type), $
        count_type)
    endif
    if count_type eq 0 then begin
      message, 'Spectral window type not recognized.'
    endif else begin
      ps_options.spec_window_type = window_type_list[wh_type[0]]
      sw_tag = '_' + win_tag_list[wh_type[0]]
    endelse
  endif else sw_tag = ''

  if tag_exist(uvf_options, 'image_window_name') ne 0 then begin
    wh_type = where(strlowcase(window_type_list) eq strlowcase(uvf_options.image_window_name), $
      count_type)
    if count_type eq 0 then begin
      wh_type = where(strlowcase(win_tag_list) eq strlowcase(uvf_options.image_window_name), $
        count_type)
    endif
    if count_type eq 0 then message, 'Image window type not recognized.' else begin
      image_window_name = window_type_list[wh_type[0]]
      if window_type_list[wh_type] ne 'None' and tag_exist(uvf_options, 'image_window_frac_size') ne 0 then begin
        iw_size_tag = number_formatter(uvf_options.image_window_frac_size)
      endif else iw_size_tag = ''
      if image_window_name eq 'None' then begin
        iw_tag = ''
      endif else begin
        iw_tag = '_' + win_tag_list[wh_type[0]] + iw_size_tag
      endelse
    endelse
  endif else iw_tag = ''

  if n_elements(freq_ch_range) ne 0 then begin
    if min(freq_ch_range) lt 0 or max(freq_ch_range) - min(freq_ch_range) lt 3 then begin
      message, 'invalid freq_ch_range'
    endif
    fch_tag = '_ch' + number_formatter(min(freq_ch_range)) + '-' +$
      number_formatter(max(freq_ch_range))
  endif else fch_tag = ''

  if n_elements(freq_flags) ne 0 then begin
    if n_elements(freq_flag_name) eq 0 then begin
      freq_flag_name = ''
    endif else begin
      if size(freq_flag_name, /type) ne 7 then freq_flag_name = number_formatter(freq_flag_name)
    endelse
    flag_tag = '_flag' + freq_flag_name
  endif else flag_tag = ''
  fch_tag = fch_tag + flag_tag

  if tag_exist(uvf_options, 'delta_uv_lambda') then begin
    uv_tag = '_deluv' + number_formatter(uvf_options.delta_uv_lambda)
  endif else begin
    uv_tag = ''
  endelse
  ; if tag_exist(uvf_options, 'max_uv_lambda') then begin
  ;   uv_tag = uv_tag + '_maxuv' + number_formatter(uvf_options.max_uv_lambda)
  ; endif
  if tag_exist(uvf_options, 'uv_avg') then begin
    uv_tag = uv_tag + '_uvavg' + number_formatter(uvf_options.uv_avg)
  endif
  if tag_exist(uvf_options, 'uv_img_clip') then begin
    uv_tag = uv_tag + '_uvimgclip' + number_formatter(uvf_options.uv_img_clip)
  endif
  uvf_tag = uv_tag + iw_tag + fch_tag

  if ps_options.std_power then kcube_tag = '_stdp' else kcube_tag = ''
  if ps_options.inverse_covar_weight then kcube_tag = kcube_tag + '_invcovar'
  if ps_options.ave_removal then kcube_tag = kcube_tag + '_averemove'
  kcube_tag = kcube_tag + sw_tag

  power_tag = kcube_tag
  if ps_options.no_wtd_avg then power_tag = power_tag + '_nowtavg'

  file_tags = {uvf_tag: uvf_tag, kcube_tag: kcube_tag, power_tag: power_tag}

  return, file_tags
end
