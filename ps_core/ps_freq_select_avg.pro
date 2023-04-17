pro ps_freq_select_avg, file_struct, n_vis_freq, refresh_options = refresh_options, $
  freq_options = freq_options, ps_options = ps_options, $
  metadata_only = metadata_only, file_check_only = file_check_only

  if n_elements(metadata_only) eq 0 then metadata_only = 0
  if n_elements(file_check_only) eq 0 then file_check_only = 0

  nfiles = n_elements(file_struct.datafile) 

  if tag_exist(file_struct, 'beam_int') then begin
      beam_int = file_struct.beam_int
  endif

  new_n_vis_freq = n_vis_freq
  original_freqs = file_struct.frequencies
  if tag_exist(freq_options, 'freq_flags') then begin
    original_mask = freq_options.freq_mask
  endif
  if tag_exist(freq_options, 'freq_ch_range') then begin
    original_freqs = original_freqs[min(freq_options.freq_ch_range):max(freq_options.freq_ch_range)]
    if nfiles eq 2 then begin
      new_n_vis_freq = new_n_vis_freq[*, min(freq_options.freq_ch_range):max(freq_options.freq_ch_range)]
    endif else begin
      new_n_vis_freq = new_n_vis_freq[min(freq_options.freq_ch_range):max(freq_options.freq_ch_range)]
    endelse
    if tag_exist(freq_options, 'freq_flags') then begin
      original_mask = original_mask[min(freq_options.freq_ch_range):max(freq_options.freq_ch_range)]
    endif
    if tag_exist(file_struct, 'beam_int') then begin
      if nfiles eq 2 then begin
        beam_int = beam_int[*, min(freq_options.freq_ch_range):max(freq_options.freq_ch_range)]
      endif else begin
        beam_int = beam_int[min(freq_options.freq_ch_range):max(freq_options.freq_ch_range)]
      endelse
    endif
  endif
  if tag_exist(freq_options, 'freq_flags') then begin
    if nfiles eq 2 then begin
        mask_use = rebin(reform(original_mask, 1, n_elements(original_freqs)), size(beam_int,/dimension), /sample)
    endif else begin
        mask_use = original_mask
    endelse
    new_n_vis_freq = new_n_vis_freq * mask_use
  endif
  if tag_exist(file_struct, 'beam_int') then begin
    if nfiles eq 2 then begin
      beam_int = total(beam_int * new_n_vis_freq, 2) / total(new_n_vis_freq, 2)
    endif else begin
      beam_int = total(beam_int * new_n_vis_freq) / total(new_n_vis_freq)
    endelse
  endif

  new_n_vis_freq_select = new_n_vis_freq

  ;; average in frequency
  if freq_options.freq_avg_factor gt 1 then begin
    ;; warn about freq_ch_range also being set
    if tag_exist(freq_options, 'freq_ch_range') then begin
      print, "both freq_avg_factor and freq_ch_range are set. The values in " $
        + "freq_ch_range are being interpreted as the original channel numbers."
    endif
    ;; check that factor divides evenly into number of frequencies
    if n_elements(original_freqs) mod freq_options.freq_avg_factor ne 0 then begin
      message, "freq_avg_factor must divide evenly into number of frequencies to be " $
        + "averaged (accounting for freq_ch_range if set)."
    endif

    avg_n_freqs = n_elements(original_freqs) / freq_options.freq_avg_factor
    if tag_exist(freq_options, 'freq_flags') and (freq_options.force_even_freqs eq 0) then begin
      print, "both freq_avg_factor and freq_flags are set. The values in " $
        + "freq_flags are being interpreted as the original channel numbers."
      ;; compute new frequencies accounting for flagging.
      frequencies = fltarr(avg_n_freqs)
      for fi=0, avg_n_freqs-1 do begin
        this_inds = findgen(freq_options.freq_avg_factor) + (fi * freq_options.freq_avg_factor)
        this_freqs = original_freqs[this_inds]
        this_mask = original_mask[this_inds]
        wh_unflagged = where(this_mask eq 1, count_unflagged)
        if count_unflagged gt 0 then begin
          ;; take the mean over unflagged channels
          frequencies[fi] = mean(this_freqs[wh_unflagged])
        endif else begin
          ;; all channels are flagged, take the mean over all channels
          frequencies[fi] = mean(this_freqs)
        endelse
      endfor
    endif else begin
      ;; do a simple average of the input frequencies.
      frequencies = reform(original_freqs, freq_options.freq_avg_factor, avg_n_freqs)
      frequencies = mean(frequencies, dimension=1)
    endelse
    ;; After averaging, if there are any fully flagged frequencies they should remain flagged
    if tag_exist(freq_options, 'freq_flags') then begin
      new_mask = reform(original_mask, freq_options.freq_avg_factor, avg_n_freqs)
      new_mask = max(new_mask, dimension=1)
    endif
    if nfiles eq 2 then begin
      new_n_vis_freq = reform(new_n_vis_freq, nfiles, freq_options.freq_avg_factor, avg_n_freqs)
      new_n_vis_freq = total(new_n_vis_freq, 2)
    endif else begin
      new_n_vis_freq = reform(new_n_vis_freq, freq_options.freq_avg_factor, avg_n_freqs)
      new_n_vis_freq = total(new_n_vis_freq, 1)
    endelse

    ;; require a freq_dft if the resulting frequencies are not regularly spaced
    ;; check frequency spacing to determine if a frequency DFT is required.
    if ( $
      not ps_options.freq_dft $
      and tag_exist(freq_options, 'freq_flags') $
      and freq_options.force_even_freqs eq 0 $
    ) then begin
      freq_diff = frequencies - shift(frequencies, 1)
      freq_diff = freq_diff[1:*]
      if max(abs(freq_diff-freq_diff[0])) gt 1e-12 then begin
        ps_options = create_ps_options(ps_options = ps_options, freq_dft = 1)
      endif
    endif
    freq_options = create_freq_options(freq_options = freq_options, $
      frequencies = frequencies, n_vis_freq = new_n_vis_freq, beam_int = beam_int)
    if tag_exist(freq_options, 'freq_flags') then begin
      freq_options = create_freq_options(freq_options = freq_options, new_freq_mask = new_mask)
    endif
  endif else begin ;; endif averaging
    if tag_exist(freq_options, 'freq_ch_range') then begin
      freq_options = create_freq_options(freq_options = freq_options, $
        frequencies = original_freqs, n_vis_freq = new_n_vis_freq, beam_int = beam_int)
      if tag_exist(freq_options, 'freq_flags') then begin
        freq_options = create_freq_options(freq_options = freq_options, new_freq_mask = original_mask)
      endif
    endif else begin
      ;; get here if there's only flagging but no channel selection or averaging
      freq_options = create_freq_options(freq_options = freq_options, $
        frequencies = original_freqs, n_vis_freq = new_n_vis_freq, beam_int = beam_int)
    endelse
  endelse

  if metadata_only then begin
    ;; don't do anything with the actual data cubes
    return
  endif

  for i=0, nfiles-1 do begin
    test_uvf = file_valid(file_struct.uvf_savefile[i])
    if not test_uvf then test_uvf = check_old_path(file_struct, 'uvf_savefile', index=i)
    if not test_uvf and file_check_only then begin
      message, "the uvf file (" + file_struct.uvf_savefile[i] + ") does not exist."
    endif

    test_wt_uvf = file_valid(file_struct.uvf_weight_savefile[i])
    if not test_wt_uvf then test_wt_uvf = check_old_path(file_struct, 'uvf_weight_savefile', index=i)
    if not test_wt_uvf and file_check_only then begin
      message, "the weight_uvf file (" + file_struct.uvf_weight_savefile[i] + ") does not exist."
    endif

    if tag_exist(file_struct, 'beam_savefile') then begin
      test_beam = file_valid(file_struct.beam_savefile[i])
      if not test_beam then test_beam = check_old_path(file_struct, 'beam_savefile', index=i)
      if not test_beam and file_check_only then begin
        message, "the beam file (" + file_struct.beam_savefile[i] + ") does not exist."
      endif
    endif else test_beam = 1

    if tag_exist(freq_options, 'freq_flags') then begin
      if test_uvf eq 1 then begin
        old_freq_mask = getvar_savefile(file_struct.uvf_savefile[i], 'freq_mask')
        if total(abs(old_freq_mask - freq_options.freq_mask)) ne 0 then begin
          if file_check_only then begin
            message, "the uvf file (" + file_struct.uvf_savefile[i] + ") has the wrong freq flagging."
          endif
          test_uvf = 0
        endif
      endif
      if test_wt_uvf eq 1 then begin
        old_freq_mask = getvar_savefile(file_struct.uvf_weight_savefile[i], 'freq_mask')
        if total(abs(old_freq_mask - freq_options.freq_mask)) ne 0 then begin
          if file_check_only then begin
            message, "the weight_uvf file (" + file_struct.uvf_weight_savefile[i] + ") has the wrong freq flagging."
          endif
          test_wt_uvf = 0
        endif
      endif
      if test_beam eq 1 then begin
        old_freq_mask = getvar_savefile(file_struct.beam_savefile[i], 'freq_mask')
        if total(abs(old_freq_mask - freq_options.freq_mask)) ne 0 then begin
          if file_check_only then begin
            message, "the beam file (" + file_struct.beam_savefile[i] + ") has the wrong freq flagging."
          endif
          test_beam = 0
        endif
      endif
      ;; This is the original freq_mask, which is what we want to save in files
      freq_mask = freq_options.freq_mask
    endif

    n_freq_orig = n_elements(file_struct.frequencies)
    if test_uvf eq 0 or refresh_options.refresh_freq_select_avg then begin

      ;; check for the full cube to avoid redoing the DFT
      full_uvf_file = file_struct.uvf_full_savefile[i]
      test_full_uvf = file_valid(full_uvf_file)
      if not test_full_uvf then test_full_uvf = check_old_path(file_struct, full_uvf_file, index=i)

      if test_full_uvf eq 0 then begin
        message, 'full uvf cube does not exist. expected full file: ' + full_uvf_file
      endif
      restore, full_uvf_file

      if tag_exist(freq_options, 'freq_flags') then begin
        data_cube = data_cube * rebin(reform(freq_options.freq_mask, 1, 1, n_freq_orig), $
          size(data_cube, /dimension), /sample)
      endif
      if tag_exist(freq_options, 'freq_ch_range') then begin
        data_cube = data_cube[*, *, min(freq_options.freq_ch_range):max(freq_options.freq_ch_range)]
      endif

      ;; frequency averaging (really summing because this is just the numerator)
      if freq_options.freq_avg_factor gt 1 then begin
        data_shape = size(data_cube, /dimensions)
        data_cube = reform(data_cube, data_shape[0], data_shape[1], freq_options.freq_avg_factor, avg_n_freqs)
        data_cube = total(data_cube, 3)
      endif

      if tag_exist(freq_options, 'freq_flags') then begin
        save, file = file_struct.uvf_savefile[i], kx_rad_vals, ky_rad_vals, $
          data_cube, freq_mask, uvf_git_hash
      endif else begin
        save, file = file_struct.uvf_savefile[i], kx_rad_vals, ky_rad_vals, $
          data_cube, uvf_git_hash
      endelse
      undefine, data_cube, uvf_git_hash


    endif ;; endif test_uvf eq 0 or refresh

    if test_wt_uvf eq 0 or refresh_options.refresh_freq_select_avg then begin

      ;; if this is a limited freq. range cube, check for the full cube to avoid redoing the DFT
      full_uvf_wt_file = file_struct.uvf_weight_full_savefile[i]
      test_full_wt_uvf = file_valid(full_uvf_wt_file)
      if not test_full_wt_uvf then test_full_wt_uvf = check_old_path(file_struct, 'uvf_weight_savefile', index=i)

      if test_full_uvf eq 0 then begin
        message, 'full uvf weight cube does not exist. expected full file: ' + full_uvf_wt_file
      endif
      restore, full_uvf_wt_file

      if tag_exist(freq_options, 'freq_flags') then begin
        flag_arr = rebin(reform(freq_options.freq_mask, 1, 1, n_freq_orig), size(weights_cube, /dimension), /sample)
        weights_cube = weights_cube * flag_arr
        variance_cube = variance_cube * flag_arr
      endif
      if tag_exist(freq_options, 'freq_ch_range') then begin
        weights_cube = weights_cube[*, *, min(freq_options.freq_ch_range):max(freq_options.freq_ch_range)]
        variance_cube = variance_cube[*, *, min(freq_options.freq_ch_range):max(freq_options.freq_ch_range)]
      endif

      ;; frequency averaging (really summing because these are just the numerator and denominators)
      if freq_options.freq_avg_factor gt 1 then begin
        weights_cube = reform(weights_cube, data_shape[0], data_shape[1], freq_options.freq_avg_factor, avg_n_freqs)
        weights_cube = total(weights_cube, 3)
        variance_cube = reform(variance_cube, data_shape[0], data_shape[1], freq_options.freq_avg_factor, avg_n_freqs)
        variance_cube = total(variance_cube, 3)
      endif

      if tag_exist(freq_options, 'freq_flags') then begin
        save, file = file_struct.uvf_weight_savefile[i], kx_rad_vals, $
          ky_rad_vals, weights_cube, variance_cube, freq_mask, uvf_wt_git_hash
      endif else begin
        save, file = file_struct.uvf_weight_savefile[i], kx_rad_vals, $
          ky_rad_vals, weights_cube, variance_cube, uvf_wt_git_hash
      endelse
      undefine, weights_cube, variance_cube, uvf_wt_git_hash

    endif ;; endif test_wt_uvf eq 0 or refresh

    if test_beam eq 0 or refresh_options.refresh_freq_select_avg then begin

      ;; check for the full cube to avoid redoing the pixel selection
      full_beam_file = file_struct.beam_full_savefile[i]
      test_full_beam = file_valid(full_beam_file)
      if not test_full_beam then test_full_beam = check_old_path(file_struct, full_beam_file, index=i)

      if test_full_beam eq 0 then begin
        message, 'full beam file does not exist. expected full file: ' + full_beam_file
      endif

      pixels_kept = getvar_savefile(file_struct.beam_full_savefile[i], "pixels")
      beam_git_hash = getvar_savefile(file_struct.beam_full_savefile[i], "beam_git_hash")
      git, repo_path = ps_repository_dir(), result=this_run_git_hash

      ;; also have to get full cube to do frequency selections before averaging.
      beam = getvar_savefile(file_struct.beamfile[i], file_struct.beamvar)
      if tag_exist(file_struct, 'nside') ne 0 then begin
        pixel_nums = getvar_savefile(file_struct.pixelfile[0], file_struct.pixelvar[0])
      endif

      match, pixel_nums, pixels_kept, pixel_nums_inds, subb, count=n_pixels_keep
      if n_pixels_keep ne n_elements(pixel_nums) then begin
        pixels = pixel_nums[pixel_nums_inds]
        beam = beam[pixel_nums_inds, *]
      endif else begin
        pixels = pixel_nums
      endelse

      if tag_exist(freq_options, 'freq_ch_range') then begin
        beam = beam[*, *, min(freq_options.freq_ch_range):max(freq_options.freq_ch_range)]
      endif

      if max(beam) le 1.1 then begin
        ;; beam is peak normalized to 1
        print, 'WARNING: It appears that this beam is from a very old run ' + $
          'of FHD  (max of beam^2 is ~1). If that is not the case, there ' + $
          'might be something unexpected happening'
        temp = beam * rebin(reform(new_n_vis_freq_select[i, *], 1, n_freq_orig), n_pixels_keep, n_freq_orig, /sample)
      endif else if max(beam) le file_struct.n_obs[i]*1.1 then begin
        ;; beam is peak normalized to 1 for each obs, then summed over obs so peak is ~ n_obs
        print, 'WARNING: It appears that this beam is from an old run ' + $
          'of FHD  (max of beam^2 is ~n_obs). If that is not the case, there ' + $
          'might be something unexpected happening'
        temp = (beam/file_struct.n_obs[i]) * rebin(reform(new_n_vis_freq_select[i, *], 1, n_freq_orig), $
          n_pixels_keep, n_freq_orig, /sample)
      endif else begin
        ;; beam is peak normalized to 1 then multiplied by nvis_freq for each obs & summed
        temp = beam
        if tag_exist(freq_options, 'freq_flags') then begin
          ;; original_mask has freq_ch_range applied but not freq averaging
          temp *= rebin(reform(original_mask, 1, n_elements(original_mask)), size(beam, /dimension), /sample)
        endif
      endelse

      avg_beam = total(temp, 2) / total(new_n_vis_freq_select[i, *])

      nside = file_struct.nside

      beam_git_hash = this_run_git_hash

      save, file=file_struct.beam_savefile[i], avg_beam, pixels, freq_mask, nside, beam_git_hash

    endif ;; endif test_beam eq 0 or refresh

  endfor ;; end for loop over files

end