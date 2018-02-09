pro make_std_plots, file_struct_arr, savefiles_2d, titles, $
    plotfile_struct = plotfile_struct, plot_sim_noise = plot_sim_noise, $
    window_num = window_num, vs_note = vs_note, wedge_amp = wedge_amp, $
    plot_options = plot_options, plot_2d_options = plot_2d_options

  nfiles = file_struct_arr[0].nfiles
  npol = max(file_struct_arr.pol_index) + 1
  ntype = max(file_struct_arr.type_index) + 1
  n_cubes = n_elements(file_struct_arr)

  if plot_options.pub and plot_options.individual_plots then begin
    for i=0, n_cubes-1 do begin
      kpower_2d_plots, savefiles_2d[i], $
        plotfile = plotfile_struct.power[i], $
        plot_options = plot_options, plot_2d_options = plot_2d_options, $
        title_prefix = titles[i], wedge_amp = wedge_amp
    endfor

    for i=0, npol -1 do begin
      cube_ind = ntype * i

      kpower_2d_plots, savefiles_2d[cube_ind], /plot_sigma, $
        plotfile = plotfile_struct.error[i], $
        plot_options = plot_options, plot_2d_options = plot_2d_options, $
        title_prefix = file_struct_arr[cube_ind].pol, note = plot_options.note + ' ' + vs_note, $
        wedge_amp = wedge_amp
    endfor

    for i=0, npol -1 do begin
      cube_ind = ntype * i

      kpower_2d_plots, savefiles_2d[cube_ind], /plot_exp_noise, $
        plotfile = plotfile_struct.noise_expval[i], $
        plot_options = plot_options, plot_2d_options = plot_2d_options, $
        title_prefix = file_struct_arr[cube_ind].pol, note = plot_options.note + ' ' + vs_note, $
        wedge_amp = wedge_amp
    endfor

    for i=0, n_cubes-1 do begin
      kpower_2d_plots, savefiles_2d[i], /snr, $
        plotfile = plotfile_struct.snr[i], $
        plot_options = plot_options, plot_2d_options = plot_2d_options, $
        title_prefix = titles[i], wedge_amp = wedge_amp
    endfor

    if nfiles eq 2 then begin
      for i=0, n_cubes-1 do begin
        kpower_2d_plots, savefiles_2d[i], /plot_noise, $
          plotfile = plotfile_struct.noise[i], $
          plot_options = plot_options, plot_2d_options = plot_2d_options, $
          title_prefix = titles[i], note = plot_options.note + ' ' + vs_note, $
          wedge_amp = wedge_amp
      endfor

      if keyword_set(plot_sim_noise) then begin
        for i=0, n_cubes-1 do begin
          kpower_2d_plots, savefiles_2d[i], /plot_sim_noise, $
            plotfile = plotfile_struct.sim_noise[i], $
            plot_options = plot_options, plot_2d_options = plot_2d_options, $
            title_prefix = titles[i], note = plot_options.note + ' ' + vs_note, $
            wedge_amp = wedge_amp
        endfor
      endif

      if keyword_set(plot_sim_noise) then begin
        for i=0, n_cubes-1 do begin
          kpower_2d_plots, savefiles_2d[i], /plot_simnoise_diff, $
            plotfile = plotfile_struct.sim_noise_diff[i], $
            plot_options = plot_options, plot_2d_options = plot_2d_options, $
            title_prefix = titles[i], note = plot_options.note + ' ' + vs_note, $
            wedge_amp = wedge_amp
        endfor
      endif

      for i=0, n_cubes-1 do begin
        kpower_2d_plots, savefiles_2d[i], /nnr, plotfile = plotfile_struct.nnr[i], $
          plot_options = plot_options, plot_2d_options = plot_2d_options, $
          title_prefix = titles[i], note = plot_options.note + ' ' + vs_note, $
          wedge_amp = wedge_amp
      endfor

      if keyword_set(plot_sim_noise) then begin
        for i=0, n_cubes-1 do begin
          kpower_2d_plots, savefiles_2d[i], /sim_nnr, plotfile = plotfile_struct.sim_nnr[i], $
            plot_options = plot_options, plot_2d_options = plot_2d_options, $
            title_prefix = titles[i], note = plot_options.note + ' ' + vs_note, $
            wedge_amp = wedge_amp
        endfor
      endif
    endif
  endif else begin
    if ntype gt 1 then begin
      cube_inds = indgen(ntype,npol)
      plot_cube_order = intarr(ntype,npol)
      for pol_i=0, npol-1 do begin
        wh_pol = where(file_struct_arr.pol_index eq pol_i, count_wh_pol)
        if count_wh_pol eq 0 then message, 'no cubes for pol_index = ' + pol_i
        plot_type_order = sort(file_struct_arr[wh_pol].type)
        plot_cube_order[*,pol_i] = cube_inds[plot_type_order,pol_i]
      endfor
    endif else plot_cube_order = indgen(npol)

    if plot_2d_options.kperp_linear_axis then begin
      ;; aspect ratio doesn't work out for kperp_linear with multiple rows
      ncol = ntype*npol
      nrow = 1
    endif else begin
      ncol = ntype
      nrow = npol
    endelse
    start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}

    window_num = window_num+1
    for i=0, n_cubes-1 do begin
      cube_i = plot_cube_order[i]
      if i gt 0 then  pos_use = positions[*,i]
      if i eq n_cubes-1 and tag_exist(plot_options, 'note') then begin
        note_use = plot_options.note
      endif else begin
        undefine, note_use
      endelse
      if plot_options.pub then begin
        plotfile_use = plotfile_struct.power
      endif else begin
        undefine, plotfile_use
      endelse

      kpower_2d_plots, savefiles_2d[cube_i], multi_pos = pos_use, $
        start_multi_params = start_multi_params, plotfile = plotfile_use, $
        plot_options = plot_options, plot_2d_options = plot_2d_options, $
        title_prefix = titles[cube_i], note = note_use, wedge_amp = wedge_amp, $
        window_num = window_num

      if i eq 0 then begin
        positions = pos_use
        undefine, start_multi_params
      endif

    endfor
    if plot_options.pub then begin
      cgps_close, png = plot_options.png, pdf = plot_options.pdf, delete_ps = plot_options.delete_ps, density = 600
    endif
    undefine, positions, pos_use

    if plot_2d_options.kperp_linear_axis then begin
      ;; aspect ratio doesn't work out for kperp_linear with multiple rows
      ncol = 2*npol
      nrow = 1
    endif else begin
      ncol = 2
      nrow = npol
    endelse
    start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}

    window_num = window_num+1
    for i=0, npol*2-1 do begin
      if i gt 0 then pos_use = positions[*,i]
      if i eq npol*2-1 and tag_exist(plot_options, 'note') then begin
        note_use = plot_options.note + ' ' + vs_note
      endif else begin
        undefine, note_use
      endelse
      if plot_options.pub then begin
        plotfile_use = plotfile_struct.error
      endif else begin
        undefine, plotfile_use
      endelse

      pol_ind = i / 2
      cube_ind = ntype * pol_ind

      if i mod 2 eq 0 then begin
        kpower_2d_plots, savefiles_2d[cube_ind], multi_pos = pos_use, $
          start_multi_params = start_multi_params, plotfile = plotfile_use, /plot_sigma, $
          plot_options = plot_options, plot_2d_options = plot_2d_options, $
          title_prefix = file_struct_arr[cube_ind].pol, wedge_amp = wedge_amp, $
          window_num = window_num
      endif else begin
        kpower_2d_plots, savefiles_2d[cube_ind], multi_pos = pos_use, $
          start_multi_params = start_multi_params, /plot_exp_noise, $
          plot_options = plot_options, plot_2d_options = plot_2d_options, $
          title_prefix = file_struct_arr[cube_ind].pol, note = note_use, $
          wedge_amp = wedge_amp
      endelse
      if i eq 0 then begin
        positions = pos_use
        undefine, start_multi_params
      endif
    endfor
    undefine, positions, pos_use
    if plot_options.pub then begin
      cgps_close, png = plot_options.png, pdf = plot_options.pdf, delete_ps = plot_options.delete_ps, density = 600
    endif


    ;; now plot SNR -- no separate sigma plots
    if plot_2d_options.kperp_linear_axis then begin
      ;; aspect ratio doesn't work out for kperp_linear with multiple rows
      ncol = ntype*npol
      nrow = 1
    endif else begin
      ncol = ntype
      nrow = npol
    endelse
    start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}

    window_num = window_num+1
    for i=0, n_cubes-1 do begin
      cube_i = plot_cube_order[i]
      if i gt 0 then  pos_use = positions[*,i]
      if i eq n_cubes-1 and tag_exist(plot_options, 'note') then begin
        note_use = plot_options.note
      endif else begin
        undefine, note_use
      endelse
      if plot_options.pub then plotfile_use = plotfile_struct.snr else undefine, plotfile_use

      kpower_2d_plots, savefiles_2d[cube_i], /snr, multi_pos = pos_use, $
        start_multi_params = start_multi_params, plotfile = plotfile_use, $
        plot_options = plot_options, plot_2d_options = plot_2d_options, $
        title_prefix = titles[cube_i], note = note_use, $
        wedge_amp = wedge_amp, window_num = window_num

      if i eq 0 then begin
        positions = pos_use
        undefine, start_multi_params
      endif
    endfor
    undefine, positions, pos_use
    if plot_options.pub then begin
      cgps_close, png = plot_options.png, pdf = plot_options.pdf, delete_ps = plot_options.delete_ps, density = 600
    endif

    if keyword_set(plot_sim_noise) then begin
      window_num = window_num+1
      start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}

      for i=0, n_cubes-1 do begin
        cube_i = plot_cube_order[i]
        if i gt 0 then  pos_use = positions[*,i]
        if i eq n_cubes-1 and tag_exist(plot_options, 'note') then begin
          note_use = plot_options.note + ' ' + vs_note
        endif else begin
          undefine, note_use
        endelse
        if plot_options.pub then begin
          plotfile_use = plotfile_struct.sim_noise
        endif else begin
          undefine, plotfile_use
        endelse

        kpower_2d_plots, savefiles_2d[cube_i], /plot_sim_noise, multi_pos = pos_use, $
          start_multi_params = start_multi_params, plotfile = plotfile_use, $
          plot_options = plot_options, plot_2d_options = plot_2d_options, $
          title_prefix = titles[cube_i], note = note_use, $
          wedge_amp = wedge_amp, window_num = window_num

        if i eq 0 then begin
          positions = pos_use
          undefine, start_multi_params
        endif
      endfor
      undefine, positions, pos_use
      if plot_options.pub then begin
        cgps_close, png = plot_options.png, pdf = plot_options.pdf, delete_ps = plot_options.delete_ps, density = 600
      endif

      window_num = window_num+1
      start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}

      for i=0, n_cubes-1 do begin
        cube_i = plot_cube_order[i]
        if i gt 0 then  pos_use = positions[*,i]
        if i eq n_cubes-1 and tag_exist(plot_options, 'note') then begin
          note_use = plot_options.note + ' ' + vs_note
        endif else begin
          undefine, note_use
        endelse
        if plot_options.pub then begin
          plotfile_use = plotfile_struct.sim_snr
        endif else begin
          undefine, plotfile_use
        endelse

        kpower_2d_plots, savefiles_2d[cube_i], /sim_snr, multi_pos = pos_use, $
          start_multi_params = start_multi_params, plotfile = plotfile_use, $
          plot_options = plot_options, plot_2d_options = plot_2d_options, $
          title_prefix = titles[cube_i], note = note_use, $
          wedge_amp = wedge_amp, window_num = window_num

        if i eq 0 then begin
          positions = pos_use
          undefine, start_multi_params
        endif
      endfor
      undefine, positions, pos_use
      if plot_options.pub then begin
        cgps_close, png = plot_options.png, pdf = plot_options.pdf, delete_ps = plot_options.delete_ps, density = 600
      endif
    endif

    if nfiles eq 2 then begin

      window_num = window_num+1
      start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}
      undefine, pos_use

      for i=0, n_cubes-1 do begin
        cube_i = plot_cube_order[i]
        if i gt 0 then  pos_use = positions[*,i]
        if i eq n_cubes-1 and tag_exist(plot_options, 'note') then begin
          note_use = plot_options.note + ' ' + vs_note
        endif else begin
          undefine, note_use
        endelse
        if plot_options.pub then begin
          plotfile_use = plotfile_struct.noise
        endif else begin
          undefine, plotfile_use
        endelse

        kpower_2d_plots, savefiles_2d[cube_i], /plot_noise, multi_pos = pos_use, $
          start_multi_params = start_multi_params, plotfile = plotfile_use, $
          plot_options = plot_options, plot_2d_options = plot_2d_options, $
          title_prefix = titles[cube_i], note = note_use, $
          wedge_amp = wedge_amp, window_num = window_num

        if i eq 0 then begin
          positions = pos_use
          undefine, start_multi_params
        endif
      endfor
      undefine, positions, pos_use
      if plot_options.pub then begin
        cgps_close, png = plot_options.png, pdf = plot_options.pdf, delete_ps = plot_options.delete_ps, density = 600
      endif

      if keyword_set(plot_sim_noise) then begin
        window_num = window_num+1
        start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}

        for i=0, n_cubes-1 do begin
          cube_i = plot_cube_order[i]
          if i gt 0 then  pos_use = positions[*,i]
          if i eq n_cubes-1 and tag_exist(plot_options, 'note') then begin
            note_use = plot_options.note + ' ' + vs_note
          endif else begin
            undefine, note_use
          endelse
          if plot_options.pub then begin
            plotfile_use = plotfile_struct.sim_noise_diff
          endif else begin
            undefine, plotfile_use
          endelse

          kpower_2d_plots, savefiles_2d[cube_i], /plot_simnoise_diff, $
            multi_pos = pos_use, start_multi_params = start_multi_params, $
            plotfile = plotfile_use, $
            plot_options = plot_options, plot_2d_options = plot_2d_options, $
            title_prefix = titles[cube_i], note = note_use, $
            wedge_amp = wedge_amp, window_num = window_num

          if i eq 0 then begin
            positions = pos_use
            undefine, start_multi_params
          endif
        endfor
        undefine, positions, pos_use
        if plot_options.pub then begin
          cgps_close, png = plot_options.png, pdf = plot_options.pdf, delete_ps = plot_options.delete_ps, density = 600
        endif
      endif

      window_num = window_num+1
      start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}

      for i=0, n_cubes-1 do begin
        cube_i = plot_cube_order[i]
        if i gt 0 then  pos_use = positions[*,i]
        if i eq n_cubes-1 and tag_exist(plot_options, 'note') then begin
          note_use = plot_options.note + ' ' + vs_note
        endif else begin
          undefine, note_use
        endelse
        if plot_options.pub then begin
          plotfile_use = plotfile_struct.nnr
        endif else begin
          undefine, plotfile_use
        endelse

        kpower_2d_plots, savefiles_2d[cube_i], /nnr, multi_pos = pos_use, $
          start_multi_params = start_multi_params, plotfile = plotfile_use, $
          plot_options = plot_options, plot_2d_options = plot_2d_options, $
          title_prefix = titles[cube_i], note = note_use, $
          wedge_amp = wedge_amp, window_num = window_num

        if i eq 0 then begin
          positions = pos_use
          undefine, start_multi_params
        endif
      endfor
      undefine, positions, pos_use
      if plot_options.pub then begin
        cgps_close, png = plot_options.png, pdf = plot_options.pdf, delete_ps = plot_options.delete_ps, density = 600
      endif

      window_num = window_num+1
      start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}

      if keyword_set(plot_sim_noise) then begin
        for i=0, n_cubes-1 do begin
          cube_i = plot_cube_order[i]
          if i gt 0 then  pos_use = positions[*,i]
          if i eq n_cubes-1 and tag_exist(plot_options, 'note') then begin
            note_use = plot_options.note + ' ' + vs_note
          endif else begin
            undefine, note_use
          endelse
          if plot_options.pub then begin
            plotfile_use = plotfile_struct.sim_nnr
          endif else begin
            undefine, plotfile_use
          endelse

          kpower_2d_plots, savefiles_2d[cube_i], /sim_nnr, multi_pos = pos_use, $
            start_multi_params = start_multi_params, plotfile = plotfile_use, $
            plot_options = plot_options, plot_2d_options = plot_2d_options, $
            title_prefix = titles[cube_i], note = note_use, $
            wedge_amp = wedge_amp, window_num = window_num

          if i eq 0 then begin
            positions = pos_use
            undefine, start_multi_params
          endif
        endfor
        undefine, positions, pos_use
        if plot_options.pub then begin
          cgps_close, png = plot_options.png, pdf = plot_options.pdf, delete_ps = plot_options.delete_ps, density = 600
        endif
      endif

    endif
  endelse
end
