pro make_std_plots, file_struct_arr, savefiles_2d, titles, $
    plotfile_struct = plotfile_struct, plot_sim_noise = plot_sim_noise, $
    window_num = window_num, pub = pub, individual_plots = individual_plots, $
    png = png, eps = eps, pdf = pdf, kperp_plot_range = kperp_plot_range, $
    kpar_plot_range = kpar_plot_range, data_range = data_range, $
    sigma_range = sigma_range, nev_range = nev_range, snr_range = snr_range, $
    noise_range = noise_range, nnr_range = nnr_range, note = note, vs_note = vs_note, $
    plot_wedge_line = plot_wedge_line, hinv = hinv, wedge_amp = wedge_amp, $
    baseline_axis = baseline_axis, delay_axis = delay_axis, $
    cable_length_axis = cable_length_axis, $
    kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis

  nfiles = file_struct_arr[0].nfiles
  npol = max(file_struct_arr.pol_index) + 1
  ntype = max(file_struct_arr.type_index) + 1
  n_cubes = n_elements(file_struct_arr)

  if keyword_set(pub) and keyword_set(individual_plots) then begin
    for i=0, n_cubes-1 do kpower_2d_plots, savefiles_2d[i], png = png, eps = eps, $
      pdf = pdf, plotfile = plotfile_struct.power[i], $
      kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
      data_range = data_range, title_prefix = titles[i], note = note, $
      plot_wedge_line = plot_wedge_line, hinv = hinv, wedge_amp = wedge_amp, $
      baseline_axis = baseline_axis, delay_axis = delay_axis, $
      cable_length_axis = cable_length_axis, $
      kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis

    for i=0, npol -1 do kpower_2d_plots, savefiles_2d[i], /plot_sigma, png = png, $
      eps = eps, pdf = pdf, plotfile = plotfile_struct.error[i], $
      kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
      data_range = sigma_range, title_prefix = file_struct_arr[i].pol, $
      note = note + ' ' + vs_note, plot_wedge_line = plot_wedge_line, hinv = hinv, $
      wedge_amp = wedge_amp, baseline_axis = baseline_axis, delay_axis = delay_axis, $
      cable_length_axis = cable_length_axis, $
      kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis

    for i=0, npol -1 do kpower_2d_plots, savefiles_2d[i], /plot_exp_noise, $
      png = png, eps = eps, pdf = pdf, plotfile = plotfile_struct.noise_expval[i],$
      kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
      data_range = nev_range, title_prefix = file_struct_arr[i].pol, $
      note = note + ' ' + vs_note, plot_wedge_line = plot_wedge_line, hinv = hinv, $
      wedge_amp = wedge_amp, baseline_axis = baseline_axis, delay_axis = delay_axis, $
      cable_length_axis = cable_length_axis, $
      kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis

    for i=0, n_cubes-1 do kpower_2d_plots, savefiles_2d[i], png = png, eps = eps, $
      pdf = pdf, /snr, plotfile = plotfile_struct.snr[i], $
      kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
      data_range = snr_range, title_prefix = titles[i], note = note, $
      plot_wedge_line = plot_wedge_line, hinv = hinv, wedge_amp = wedge_amp, $
      baseline_axis = baseline_axis, delay_axis = delay_axis, $
      cable_length_axis = cable_length_axis, $
      kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis

    if nfiles eq 2 then begin
      for i=0, n_cubes-1 do kpower_2d_plots, savefiles_2d[i], png = png, $
        eps = eps, pdf = pdf, /plot_noise, plotfile = plotfile_struct.noise[i], $
        kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
        data_range = noise_range, title_prefix = titles[i], note = note + ' ' + vs_note, $
        plot_wedge_line = plot_wedge_line, hinv = hinv, wedge_amp = wedge_amp, $
        baseline_axis = baseline_axis, delay_axis = delay_axis, $
        cable_length_axis = cable_length_axis, $
        kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis

      if keyword_set(plot_sim_noise) then $
        for i=0, n_cubes-1 do kpower_2d_plots, savefiles_2d[i], png = png, $
          eps = eps, pdf = pdf, /plot_sim_noise, plotfile = plotfile_struct.sim_noise[i], $
          kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
          data_range = noise_range, title_prefix = titles[i], note = note + ' ' + vs_note, $
          plot_wedge_line = plot_wedge_line, hinv = hinv, wedge_amp = wedge_amp, $
          baseline_axis = baseline_axis, delay_axis = delay_axis, $
          cable_length_axis = cable_length_axis, $
          kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis

      if keyword_set(plot_sim_noise) then $
        for i=0, n_cubes-1 do kpower_2d_plots, savefiles_2d[i], png = png, $
          eps = eps, pdf = pdf, /plot_simnoise_diff, plotfile = plotfile_struct.sim_noise_diff[i], $
          kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
          data_range = noise_range, title_prefix = titles[i], note = note + ' ' + vs_note, $
          plot_wedge_line = plot_wedge_line, hinv = hinv, wedge_amp = wedge_amp, $
          baseline_axis = baseline_axis, delay_axis = delay_axis, $
          cable_length_axis = cable_length_axis, $
          kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis

      for i=0, n_cubes-1 do kpower_2d_plots, savefiles_2d[i], png = png, eps = eps, $
        pdf = pdf, /nnr, plotfile = plotfile_struct.nnr[i], $
        kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
        data_range = nnr_range, title_prefix = titles[i], note = note + ' ' + vs_note, $
        plot_wedge_line = plot_wedge_line, hinv = hinv, wedge_amp = wedge_amp, $
        baseline_axis = baseline_axis, delay_axis = delay_axis, $
        cable_length_axis = cable_length_axis, $
        kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis

      if keyword_set(plot_sim_noise) then $
        for i=0, n_cubes-1 do kpower_2d_plots, savefiles_2d[i], png = png, $
          eps = eps, pdf = pdf, /sim_nnr, plotfile = plotfile_struct.sim_nnr[i], $
          kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
          data_range = nnr_range, title_prefix = titles[i], note = note + ' ' + vs_note, $
          plot_wedge_line = plot_wedge_line, hinv = hinv, $
          wedge_amp = wedge_amp, baseline_axis = baseline_axis, delay_axis = delay_axis, $
          cable_length_axis = cable_length_axis, $
          kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis
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

    if keyword_set(kperp_linear_axis) then begin
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
      if i eq n_cubes-1 and n_elements(note) gt 0 then begin
        note_use = note
      endif else begin
        undefine, note_use
      endelse
      if keyword_set(pub) then begin
        plotfile_use = plotfile_struct.power
      endif else begin
        undefine, plotfile_use
      endelse

      kpower_2d_plots, savefiles_2d[cube_i], multi_pos = pos_use, $
        start_multi_params = start_multi_params, png = png, eps = eps, pdf = pdf, $
        plotfile = plotfile_use, kperp_plot_range = kperp_plot_range, $
        kpar_plot_range = kpar_plot_range, data_range = data_range, $
        title_prefix = titles[cube_i], note = note_use, plot_wedge_line = plot_wedge_line, $
        hinv = hinv, wedge_amp = wedge_amp, baseline_axis = baseline_axis, $
        delay_axis = delay_axis, cable_length_axis = cable_length_axis, $
        kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, $
        window_num = window_num

      if i eq 0 then begin
        positions = pos_use
        undefine, start_multi_params
      endif

    endfor
    if keyword_set(pub) then begin
      cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density = 600
    endif
    undefine, positions, pos_use

    if keyword_set(kperp_linear_axis) then begin
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
      if i eq npol*2-1 and n_elements(note) gt 0 then begin
        note_use = note + ' ' + vs_note
      endif else begin
        undefine, note_use
      endelse
      if keyword_set(pub) then begin
        plotfile_use = plotfile_struct.error
      endif else begin
        undefine, plotfile_use
      endelse

      pol_ind = i / 2

      if i mod 2 eq 0 then begin
        kpower_2d_plots, savefiles_2d[pol_ind], multi_pos = pos_use, $
        start_multi_params = start_multi_params, png = png, eps = eps, pdf = pdf, $
        plotfile = plotfile_use, /plot_sigma, data_range = sigma_range, $
        kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
        title_prefix = file_struct_arr[pol_ind].pol, plot_wedge_line = plot_wedge_line, $
        wedge_amp = wedge_amp, hinv = hinv, baseline_axis = baseline_axis, $
        delay_axis = delay_axis, cable_length_axis = cable_length_axis, $
        kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, $
        window_num = window_num
      endif else begin
        kpower_2d_plots, savefiles_2d[pol_ind], multi_pos = pos_use, $
        start_multi_params = start_multi_params, png = png, eps = eps, pdf = pdf, $
        /plot_exp_noise, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
        data_range = nev_range, title_prefix = file_struct_arr[pol_ind].pol, $
        note = note_use, plot_wedge_line = plot_wedge_line, hinv = hinv, $
        wedge_amp = wedge_amp, baseline_axis = baseline_axis, delay_axis = delay_axis, $
        cable_length_axis = cable_length_axis, $
        kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis
      endelse
      if i eq 0 then begin
        positions = pos_use
        undefine, start_multi_params
      endif
    endfor
    undefine, positions, pos_use
    if keyword_set(pub) then begin
      cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density = 600
    endif


    ;; now plot SNR -- no separate sigma plots
    if keyword_set(kperp_linear_axis) then begin
      ;; aspect ratio doesn't work out for kperp_linear with multiple rows
      ncol = ntype*npol
      nrow = 1
    endif else begin
      ncol = ntype
      nrow = npol
    endelse
    start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}

    window_num = window_num+1
    ;;snr_range = [1e0, 1e6]
    for i=0, n_cubes-1 do begin
      cube_i = plot_cube_order[i]
      if i gt 0 then  pos_use = positions[*,i]
      if i eq n_cubes-1 and n_elements(note) gt 0 then begin
        note_use = note
      endif else begin
        undefine, note_use
      endelse
      if keyword_set(pub) then plotfile_use = plotfile_struct.snr else undefine, plotfile_use

      kpower_2d_plots, savefiles_2d[cube_i], /snr, multi_pos = pos_use, $
        start_multi_params = start_multi_params, png = png, eps = eps, pdf = pdf, $
        plotfile = plotfile_use, kperp_plot_range = kperp_plot_range, $
        kpar_plot_range = kpar_plot_range, data_range = snr_range, $
        title_prefix = titles[cube_i], note = note_use, $
        plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp, $
        hinv = hinv, baseline_axis = baseline_axis, delay_axis = delay_axis, $
        cable_length_axis = cable_length_axis, kperp_linear_axis = kperp_linear_axis, $
        kpar_linear_axis = kpar_linear_axis, window_num = window_num
      if i eq 0 then begin
        positions = pos_use
        undefine, start_multi_params
      endif
    endfor
    undefine, positions, pos_use
    if keyword_set(pub) then begin
      cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density = 600
    endif

    if keyword_set(plot_sim_noise) then begin
      window_num = window_num+1
      start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}

      for i=0, n_cubes-1 do begin
        cube_i = plot_cube_order[i]
        if i gt 0 then  pos_use = positions[*,i]
        if i eq n_cubes-1 and n_elements(note) gt 0 then begin
          note_use = note + ' ' + vs_note
        endif else begin
          undefine, note_use
        endelse
        if keyword_set(pub) then begin
          plotfile_use = plotfile_struct.sim_noise
        endif else begin
          undefine, plotfile_use
        endelse

        kpower_2d_plots, savefiles_2d[cube_i], /plot_sim_noise, multi_pos = pos_use, $
          start_multi_params = start_multi_params, png = png, eps = eps, pdf = pdf, $
          plotfile = plotfile_use, kperp_plot_range = kperp_plot_range, $
          kpar_plot_range = kpar_plot_range, data_range = sigma_range, $
          title_prefix = titles[cube_i], note = note_use, $
          plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp, hinv = hinv, $
          baseline_axis = baseline_axis, delay_axis = delay_axis, $
          cable_length_axis = cable_length_axis, kperp_linear_axis = kperp_linear_axis, $
          kpar_linear_axis = kpar_linear_axis, window_num = window_num
        if i eq 0 then begin
          positions = pos_use
          undefine, start_multi_params
        endif
      endfor
      undefine, positions, pos_use
      if keyword_set(pub) then begin
        cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density = 600
      endif

      window_num = window_num+1
      start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}

      for i=0, n_cubes-1 do begin
        cube_i = plot_cube_order[i]
        if i gt 0 then  pos_use = positions[*,i]
        if i eq n_cubes-1 and n_elements(note) gt 0 then begin
          note_use = note + ' ' + vs_note
        endif else begin
          undefine, note_use
        endelse
        if keyword_set(pub) then begin
          plotfile_use = plotfile_struct.sim_snr
        endif else begin
          undefine, plotfile_use
        endelse

        kpower_2d_plots, savefiles_2d[cube_i], /sim_snr, multi_pos = pos_use,
          start_multi_params = start_multi_params, png = png, eps = eps, pdf = pdf, $
          plotfile = plotfile_use, kperp_plot_range = kperp_plot_range, $
          kpar_plot_range = kpar_plot_range, data_range = nnr_range, $
          title_prefix = titles[cube_i], note = note_use, $
          plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp, hinv = hinv, $
          baseline_axis = baseline_axis, delay_axis = delay_axis, $
          cable_length_axis = cable_length_axis, kperp_linear_axis = kperp_linear_axis, $
          kpar_linear_axis = kpar_linear_axis, window_num = window_num
        if i eq 0 then begin
          positions = pos_use
          undefine, start_multi_params
        endif
      endfor
      undefine, positions, pos_use
      if keyword_set(pub) then begin
        cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density = 600
      endif
    endif

    if nfiles eq 2 then begin

      window_num = window_num+1
      start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}
      undefine, pos_use

      for i=0, n_cubes-1 do begin
        cube_i = plot_cube_order[i]
        if i gt 0 then  pos_use = positions[*,i]
        if i eq n_cubes-1 and n_elements(note) gt 0 then begin
          note_use = note + ' ' + vs_note
        endif else begin
          undefine, note_use
        endelse
        if keyword_set(pub) then begin
          plotfile_use = plotfile_struct.noise
        endif else begin
          undefine, plotfile_use
        endelse

        kpower_2d_plots, savefiles_2d[cube_i], /plot_noise, multi_pos = pos_use, $
          start_multi_params = start_multi_params, png = png, eps = eps, pdf = pdf, $
          plotfile = plotfile_use, kperp_plot_range = kperp_plot_range, $
          kpar_plot_range = kpar_plot_range, data_range = noise_range, $
          title_prefix = titles[cube_i], note = note_use, $
          plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp, hinv = hinv, $
          baseline_axis = baseline_axis, delay_axis = delay_axis, $
          cable_length_axis = cable_length_axis, kperp_linear_axis = kperp_linear_axis, $
          kpar_linear_axis = kpar_linear_axis, window_num = window_num
        if i eq 0 then begin
          positions = pos_use
          undefine, start_multi_params
        endif
      endfor
      undefine, positions, pos_use
      if keyword_set(pub) then begin
        cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density = 600
      endif

      if keyword_set(plot_sim_noise) then begin
        window_num = window_num+1
        start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}

        for i=0, n_cubes-1 do begin
          cube_i = plot_cube_order[i]
          if i gt 0 then  pos_use = positions[*,i]
          if i eq n_cubes-1 and n_elements(note) gt 0 then begin
            note_use = note + ' ' + vs_note
          endif else begin
            undefine, note_use
          endelse
          if keyword_set(pub) then begin
            plotfile_use = plotfile_struct.sim_noise_diff
          endif else begin
            undefine, plotfile_use
          endelse

          kpower_2d_plots, savefiles_2d[cube_i], /plot_simnoise_diff, $
            multi_pos = pos_use, start_multi_params = start_multi_params, $
            png = png, eps = eps, pdf = pdf, plotfile = plotfile_use, $
            kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
            data_range = noise_range, title_prefix = titles[cube_i], note = note_use, $
            plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp, hinv = hinv, $
            baseline_axis = baseline_axis, delay_axis = delay_axis, $
            cable_length_axis = cable_length_axis, kperp_linear_axis = kperp_linear_axis, $
            kpar_linear_axis = kpar_linear_axis, window_num = window_num
          if i eq 0 then begin
            positions = pos_use
            undefine, start_multi_params
          endif
        endfor
        undefine, positions, pos_use
        if keyword_set(pub) then begin
          cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density = 600
        endif
      endif

      window_num = window_num+1
      start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}

      for i=0, n_cubes-1 do begin
        cube_i = plot_cube_order[i]
        if i gt 0 then  pos_use = positions[*,i]
        if i eq n_cubes-1 and n_elements(note) gt 0 then begin
          note_use = note + ' ' + vs_note
        endif else begin
          undefine, note_use
        endelse
        if keyword_set(pub) then begin
          plotfile_use = plotfile_struct.nnr
        endif else begin
          undefine, plotfile_use
        endelse

        kpower_2d_plots, savefiles_2d[cube_i], /nnr, multi_pos = pos_use, $
          start_multi_params = start_multi_params, png = png, eps = eps, pdf = pdf, $
          plotfile = plotfile_use, kperp_plot_range = kperp_plot_range, $
          kpar_plot_range = kpar_plot_range, data_range = nnr_range, $
          title_prefix = titles[cube_i], note = note_use, $
          plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp, hinv = hinv, $
          baseline_axis = baseline_axis, delay_axis = delay_axis, $
          cable_length_axis = cable_length_axis, kperp_linear_axis = kperp_linear_axis, $
          kpar_linear_axis = kpar_linear_axis, window_num = window_num
        if i eq 0 then begin
          positions = pos_use
          undefine, start_multi_params
        endif
      endfor
      undefine, positions, pos_use
      if keyword_set(pub) then begin
        cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density = 600
      endif

      window_num = window_num+1
      start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}

      if keyword_set(plot_sim_noise) then begin
        for i=0, n_cubes-1 do begin
          cube_i = plot_cube_order[i]
          if i gt 0 then  pos_use = positions[*,i]
          if i eq n_cubes-1 and n_elements(note) gt 0 then begin
            note_use = note + ' ' + vs_note
          endif else begin
            undefine, note_use
          endelse
          if keyword_set(pub) then begin
            plotfile_use = plotfile_struct.sim_nnr
          endif else begin
            bundefine, plotfile_use
          endelse

          kpower_2d_plots, savefiles_2d[cube_i], /sim_nnr, multi_pos = pos_use, $
            start_multi_params = start_multi_params, png = png, eps = eps, pdf = pdf, $
            plotfile = plotfile_use, kperp_plot_range = kperp_plot_range, $
            kpar_plot_range = kpar_plot_range, data_range = nnr_range, $
            title_prefix = titles[cube_i], note = note_use, $
            plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp, hinv = hinv, $
            baseline_axis = baseline_axis, delay_axis = delay_axis, $
            cable_length_axis = cable_length_axis, kperp_linear_axis = kperp_linear_axis, $
            kpar_linear_axis = kpar_linear_axis, window_num = window_num
          if i eq 0 then begin
            positions = pos_use
            undefine, start_multi_params
          endif
        endfor
        undefine, positions, pos_use
        if keyword_set(pub) then begin
          cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density = 600
        endif
      endif

    endif
  endelse
end
