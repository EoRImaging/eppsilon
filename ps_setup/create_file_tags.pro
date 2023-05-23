function create_file_tags, uvf_options = uvf_options, freq_options = freq_options, ps_options = ps_options

  window_type_list = ['Hann', 'Hamming', 'Blackman', 'Nutall', 'Blackman-Nutall', $
                      'Blackman-Harris', 'Blackman-Harris^2', 'Tukey', 'None']
  win_tag_list = ['han', 'ham', 'blm', 'ntl', 'bn', 'bh', 'bh2', 'tk', '']
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
      if win_tag_list[wh_type[0]] ne '' then begin
        sw_tag = '_sw' + win_tag_list[wh_type[0]]
      endif else begin
        sw_tag = ''
      endelse
    endelse
  endif else begin
    sw_tag = ''
  endelse

  if tag_exist(freq_options, 'freq_ch_range') then begin
    if min(freq_options.freq_ch_range) lt 0 or max(freq_options.freq_ch_range) - min(freq_options.freq_ch_range) lt 3 then begin
      message, 'invalid freq_ch_range'
    endif
    freq_tag = '_ch' + number_formatter(min(freq_options.freq_ch_range)) + '-' +$
      number_formatter(max(freq_options.freq_ch_range))
  endif else begin
    freq_tag = ''
  endelse

  if tag_exist(freq_options, 'freq_flag_name') then begin
    if size(freq_options.freq_flag_name, /type) ne 7 then begin
      freq_flag_name_use = number_formatter(freq_options.freq_flag_name)
    endif else begin
      freq_flag_name_use = freq_options.freq_flag_name
    endelse
    freq_tag = freq_tag + '_flag' + freq_flag_name_use
  endif else begin
    if tag_exist(freq_options, 'freq_flags') then begin
      freq_tag = freq_tag + '_flag'
    endif
  endelse

  if freq_options.freq_avg_factor gt 1 then begin
    freq_tag = freq_tag + '_freqave' + number_formatter(freq_options.freq_avg_factor)
    if freq_options.force_even_freqs ne 0 then begin
      freq_tag = freq_tag + '_evenfreqs'
    endif
  endif

  if tag_exist(uvf_options, 'delta_uv_lambda') then begin
    uv_tag = '_deluv' + number_formatter(uvf_options.delta_uv_lambda)
  endif else begin
    uv_tag = ''
  endelse
  if uvf_options.full_image then begin
    uv_tag = uv_tag + '_fullimg'
  endif
  if not uvf_options.image_clip then begin
    uv_tag = uv_tag + '_noimgclip'
  endif
  if tag_exist(uvf_options, 'max_uv_lambda') then begin
    uv_tag = uv_tag + '_maxuv' + number_formatter(uvf_options.max_uv_lambda)
  endif
  if tag_exist(uvf_options, 'uv_avg') then begin
    uv_tag = uv_tag + '_uvavg' + number_formatter(uvf_options.uv_avg)
  endif
  if tag_exist(uvf_options, 'uv_img_clip') then begin
    uv_tag = uv_tag + '_uvimgclip' + number_formatter(uvf_options.uv_img_clip)
  endif
  uvf_tag = uv_tag

  if ps_options.std_power then begin
    kcube_tag = '_stdp'
  endif else begin
    kcube_tag = ''
  endelse

  if ps_options.freq_dft then begin
    kcube_tag = kcube_tag + '_dft'
    ;; add a tag for using the regularly spaced approximate z rather than true z
    if ps_options.dft_z_use eq 'regular' then begin
      kcube_tag = kcube_tag + 'reg'
    endif
  endif
  if ps_options.inverse_covar_weight then begin
    kcube_tag = kcube_tag + '_invcovar'
  endif
  if ps_options.ave_removal then begin
    kcube_tag = kcube_tag + '_averemove'
  endif
  kcube_tag = kcube_tag + sw_tag

  power_tag = kcube_tag
  if ps_options.no_wtd_avg then begin
    power_tag = power_tag + '_nowtavg'
  endif

  file_tags = {uvf_tag: uvf_tag, freq_tag:freq_tag, kcube_tag: kcube_tag, power_tag: power_tag}

  return, file_tags
end
