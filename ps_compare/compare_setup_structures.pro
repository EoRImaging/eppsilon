pro compare_setup_structures, folder_names_in, obs_names_in, $
    ps_foldernames = ps_foldernames, version_test = version_test, $
    spec_window_types = spec_window_types, delta_uv_lambda = delta_uv_lambda, $
    max_uv_lambda = max_uv_lambda, full_image = full_image, $
    image_clip = image_clip, ave_removal = ave_removal, $
    freq_dft = freq_dft, dft_z_use = dft_z_use, std_power = std_power, $
    all_type_pol = all_type_pol, freq_ch_range = freq_ch_range, $
    freq_flags = freq_flags, freq_flag_name = freq_flag_name, $
    freq_flag_repeat = freq_flag_repeat, $
    diff_plot_path = diff_plot_path, diff_save_path = diff_save_path, $
    freq_avg_factor = freq_avg_factor, force_even_freqs = force_even_freqs, $
    folder_names = folder_names, obs_info = obs_info, uvf_options = uvf_options, $
    freq_options = freq_options, ps_options = ps_options


  if n_elements(folder_names_in) gt 2 then message, 'only 1 or 2 folder_names allowed'
  if n_elements(folder_names_in) eq 0 then message, 'at least 1 folder name must be specified'
  if n_elements(obs_names_in) gt 2 then message, 'only 1 or 2 obs_names_in allowed'
  
  folder_names = get_folder(folder_names_in, loc_name = loc_name,  rts = rts, $
    dirty_folder = dirty_folder)

  if keyword_set(version_test) and n_elements(ps_foldername) eq 0 $
      and n_elements(folder_names_in) eq 1 and n_elements(obs_names_in) lt 2 then begin
    git_info = git_info(ps_repository_dir())
    ps_foldernames = ['ps_master', 'ps_' + git_info.branch]
  endif

  obs_info = ps_filenames(folder_names, obs_names_in, dirty_folder = dirty_folder, $
    exact_obsnames = exact_obsnames, rts = rts, sim = sim,  uvf_input = uvf_input, $
    casa = casa, data_subdirs = data_subdirs, ps_foldernames = ps_foldernames, $
    save_paths = save_paths, plot_paths = plot_paths, refresh_info = refresh_info, $
    no_wtvar_rts = no_wtvar_rts)

  if n_elements(diff_plot_path) eq 0 then begin
    if n_elements(diff_save_path) gt 0 then begin
      diff_plot_path = diff_save_path + path_sep() + 'plots' + path_sep()
    endif
  endif

  wh_noinfo = where(obs_info.info_files eq '', count_noinfo)
  if count_noinfo gt 0 then message, 'Info files are not all present'

  if n_elements(delta_uv_lambda) gt 1 then message, 'only 1 delta_uv_lambda allowed'

  if n_elements(max_uv_lambda) lt 2 and n_elements(full_image) lt 2 $
     and n_elements(image_clip) lt 2 then begin

    uvf_options = create_uvf_options(delta_uv_lambda = delta_uv_lambda, $
      max_uv_lambda = max_uv_lambda, full_image = full_image, image_clip = image_clip)

  endif else begin
    case n_elements(max_uv_lambda) of
      0:
      1: begin
        mul0 = max_uv_lambda
        mul1 = max_uv_lambda
      end
      2: begin
        mul0 = max_uv_lambda[0]
        mul1 = max_uv_lambda[1]
      end
      else: message, 'only 1 or 2 max_uv_lambda values allowed'
    endcase

    case n_elements(full_image) of
      0:
      1: begin
        fi0 = full_image
        fi1 = full_image
      end
      2: begin
        fi0 = full_image[0]
        fi1 = full_image[1]
      end
    endcase

    case n_elements(image_clip) of
        0:
        1: begin
          ic0 = image_clip
          ic1 = image_clip
        end
        2: begin
          ic0 = image_clip[0]
          ic1 = image_clip[1]
        end
      else: message, 'only 1 or 2 image_clip values allowed'
    endcase

    uvf_options0 = create_uvf_options(delta_uv_lambda = delta_uv_lambda, $
      max_uv_lambda = mul0, full_image = fi0, image_clip = ic0)
    uvf_options1 = create_uvf_options(delta_uv_lambda = delta_uv_lambda, $
      max_uv_lambda = mul1, full_image = fi1, image_clip = ic1)
    
    uvf_options = list(uvf_options0, uvf_options1)
  endelse

  if n_elements(freq_flags) gt 0 then begin
    if size(freq_flags, /type) eq 11 then begin
      ;; this is a list of freq flags, should have 1 or 2 elements
      n_freq_flags = n_elements(freq_flags)
      if n_freq_flags gt 2 then begin
        message, "If freq_flags are specified as a list, the number of elements gives " $
        + "the number of sets of freq flags, which can only be 1 or 2."
      endif

      ;; if a single zero is passed in for a freq_flag set, turn it into a null so later logic works
      for flag_ind=0, n_freq_flags-1 do begin
        if n_elements(freq_flags[flag_ind]) eq 1 then begin
          if freq_flags[flag_ind] eq 0 then begin
            freq_flags[freq_ind] = !Null
          endif
        endif
      endfor
    endif else begin
      ;; this is an array, check the dimensionality to determine number of sets
      n_flag_dims = size(freq_flags,/n_dim)
      if n_flag_dims gt 2 then begin
        message, "If freq_flags are specified as an array, it can only be a 1 or 2 dimensional array."
      endif
      if n_flag_dims eq 2 then begin
        n_freq_flags = (size(freq_flags,/dim))[1]
      endif else begin
        n_freq_flags = 1
      endelse

      if n_freq_flags gt 2 then begin
        message, "If freq_flags are specified as an array, the second dimension gives " $
        + "the number of sets of freq flags, which can only be 1 or 2."
      endif
      ;; In this case, turn it into a list for unified indexing later
      if n_freq_flags eq 1 then begin
        freq_flags = list(freq_flags)
      endif else begin
        freq_flags = list(freq_flags[*, 0], freq_flags[*, 1])
      endelse
    endelse
  endif

  if size(freq_ch_range,/n_dim) lt 2 and n_elements(freq_flags) lt 2 $
    and n_elements(freq_flag_name) lt 2 and n_elements(freq_flag_repeat) lt 2 $
    and n_elements(freq_avg_factor) lt 2 and n_elements(force_even_freqs) lt 2 then begin

    freq_options = create_freq_options( $
      freq_ch_range = freq_ch_range, $
      freq_flags = freq_flags, $
      freq_flag_name = freq_flag_name, $
      freq_flag_repeat = freq_flag_repeat, $
      freq_avg_factor = freq_avg_factor, $
      force_even_freqs = force_even_freqs)
  endif else begin

    if n_elements(freq_ch_range) gt 0 then begin
      case size(freq_ch_range,/n_dim) of
        1: begin
          fc0 = freq_ch_range
          fc1 = freq_ch_range
        end
        2: begin
          if (size(freq_ch_range,/dim))[1] ne 2 then begin
            message, 'only 1 or 2 sets of freq_ch_range values allowed'
          endif
          fc0 = freq_ch_range[*, 0]
          fc1 = freq_ch_range[*, 1]
        end
        else: message, 'only 1 or 2 sets of freq_ch_range values allowed'
      endcase
    endif
    if n_elements(freq_flags) gt 0 then begin
      case n_elements(freq_flags) of
        1: begin
          ff0 = freq_flags
          ff1 = freq_flags
        end
        2: begin
          ff0 = freq_flags[0]
          ff1 = freq_flags[1]
        end
        else: message, 'only 1 or 2 sets of freq_flags values allowed'
      endcase
    endif
    case n_elements(freq_flag_name) of
      0:
      1: begin
        fn0 = freq_flag_name
        fn1 = freq_flag_name
      end
      2: begin
        fn0 = freq_flag_name[0]
        fn1 = freq_flag_name[1]
      end
      else: message, 'only 1 or 2 freq_flag_name values allowed'
    endcase

    case n_elements(freq_flag_repeat) of
      0:
      1: begin
        fr0 = freq_flag_repeat
        fr1 = freq_flag_repeat
      end
      2: begin
        fr0 = freq_flag_repeat[0]
        fr1 = freq_flag_repeat[1]
      end
      else: message, 'only 1 or 2 freq_flag_repeat values allowed'
    endcase

    case n_elements(freq_avg_factor) of
      0:
      1: begin
        fa0 = freq_avg_factor
        fa1 = freq_avg_factor
      end
      2: begin
        fa0 = freq_avg_factor[0]
        fa1 = freq_avg_factor[1]
      end
      else: message, 'only 1 or 2 freq_avg_factor values allowed'
    endcase

    case n_elements(force_even_freqs) of
      0:
      1: begin
        ef0 = force_even_freqs
        ef1 = force_even_freqs
      end
      2: begin
        ef0 = force_even_freqs[0]
        ef1 = force_even_freqs[1]
      end
      else: message, 'only 1 or 2 force_even_freqs values allowed'
    endcase

    freq_options1 = create_freq_options( $
      freq_ch_range = fc0, $
      freq_flags = ff0, $
      freq_flag_name = fn0, $
      freq_flag_repeat = fr0, $
      freq_avg_factor = fa0, $
      force_even_freqs = ef0)

    freq_options2 = create_freq_options( $
      freq_ch_range = fc1, $
      freq_flags = ff1, $
      freq_flag_name = fn1, $
      freq_flag_repeat = fr1, $
      freq_avg_factor = fa1, $
      force_even_freqs = ef1)

    freq_options = list(freq_options1, freq_options2)
  endelse

  if n_elements(ave_removal) lt 2 and n_elements(wt_cutoffs) lt 2 and $
    n_elements(wt_measures) lt 2 and n_elements(spec_window_types) lt 2 and $
    n_elements(freq_dft) lt 2 and n_elements(dft_z_use) lt 2 and $
    n_elements(freq_avg_factor) lt 2 and n_elements(force_even_freqs) lt 2 and $
    n_elements(std_power) lt 2 then begin

    ps_options = create_ps_options(ave_removal = ave_removal, $
    wt_cutoffs = wt_cutoffs, wt_measures = wt_measures, $
    spec_window_type = spec_window_types, $
    freq_dft = freq_dft, dft_z_use = dft_z_use, $
    std_power = std_power)

  endif else begin
    case n_elements(ave_removal) of
      0:
      1: begin
        ar0 = ave_removal
        ar1 = ave_removal
      end
      2: begin
        ar0 = ave_removal[0]
        ar1 = ave_removal[1]
      end
      else: message, 'only 1 or 2 ave_removal values allowed'
    endcase

    case n_elements(wt_cutoffs) of
      0:
      1: begin
        wtc0 = wt_cutoffs
        wtc1 = wt_cutoffs
      end
      2: begin
        wtc0 = wt_cutoffs[0]
        wtc1 = wt_cutoffs[1]
      end
      else: message, 'only 1 or 2 wt_cutoffs allowed'
    endcase

    case n_elements(wt_measures) of
      0:
      1: begin
        wtm0 = wt_measures
        wtm1 = wt_measures
      end
      2: begin
        wtm0 = wt_measures[0]
        wtm1 = wt_measures[1]
      end
      else: message, 'only 1 or 2 wt_measures allowed'
    endcase

    case n_elements(spec_window_types) of
      0:
      1: begin
        spw0 = spec_window_types
        spw1 = spec_window_types
      end
      2: begin
        spw0 = spec_window_types[0]
        spw1 = spec_window_types[1]
      end
      else: message, 'only 1 or 2 spec_window_types allowed'
    endcase

    case n_elements(freq_dft) of
      0:
      1: begin
        dft0 = freq_dft
        dft1 = freq_dft
      end
      2: begin
        dft0 = freq_dft[0]
        dft1 = freq_dft[1]
      end
      else: message, 'only 1 or 2 freq_dft values allowed'
    endcase

    case n_elements(dft_z_use) of
      0:
      1: begin
        dftz0 = dft_z_use
        dftz1 = dft_z_use
      end
      2: begin
        dftz0 = dft_z_use[0]
        dftz1 = dft_z_use[1]
      end
      else: message, 'only 1 or 2 dft_z_use values allowed'
    endcase

    case n_elements(std_power) of
      0:
      1: begin
        sp0 = std_power
        sp1 = std_power
      end
      2: begin
        sp0 = std_power[0]
        sp1 = std_power[1]
      end
      else: message, 'only 1 or 2 std_power values allowed'
    endcase

    ps_options0 = create_ps_options(ave_removal = ar0, wt_cutoffs = wtc0, $
      wt_measures = wtm0, spec_window_type = spw0, freq_dft = dft0, $
      dft_z_use = dftz0, std_power = sp0)

    ps_options1 = create_ps_options(ave_removal = ar1, wt_cutoffs = wtc1, $
      wt_measures = wtm1, spec_window_type = spw1, freq_dft = dft1, $
      dft_z_use = dftz1, std_power = sp1)

    ps_options = list(ps_options0, ps_options1)
  endelse

end