;; Program to make 2D k-space plots. Either supply an idlsave file with the following variables:
;;     power, noise, weights, noise_expval, kperp_edges, kpar_edges, kperp_bin, kpar_bin, kperp_lambda_conv, delay_params, hubble_param
;; or supply them as keywords with the same name.
;;
;; Note that all input values are assumed to be in real Mpc, to convert to h^-1 Mpc, supply hubble_param and set hinv keyword (defaulted to set)
;;
;; Keywords:
;;  Inputs (can alternatively be supplied in power_savefile)
;;   power: a 2d power array in kperp, kpar space (either linearly or exponentially binned)
;;   noise: a 2d measured noise array in same space as power
;;   weights: a 2d array of the inverse variances in same space as power
;;   noise_expval: a 2d expected noise array in same space as power
;;   kperp_edges: a vector of the bin edges for kperp (should have one more element than the 1st power dimension)
;;   kpar_edges: a vector of the bin edges for kpar (should have one more element than the 2nd power dimension)
;;   kperp_lambda_conv: conversion factor to convert between kperp and lambda (only needed if plotting baseline axis)
;;   delay_params: a 2 element vector [delay_delta, delay_max] used if plotting delay axis
;;   hubble_param: value of h used in calculations, only needed if hinv is set
;;
;;   hinv: if set, pull h out of arrays and k values so everything is in h^-1 Mpc rather than real Mpc
;;
;;  Optionally, a plot_2d_options structure (created with create_plot_2d_options.pro)
;;    and/or a plot_options structure (created with create_plot_options.pro)
;;    can be passed which sets several of the keywords. If the keywords are set separately,
;;    the keyword values supercede the structure values.
;;
;;  Plot type flags: set only one of the following to dictate what to plot, power will be plotted if none are set
;;    plot_weights: inverse variance plot
;;    plot_noise: measured noise plot
;;    plot_sigma: sqrt(variance) = 1/sqrt(weights) plot
;;    plot_exp_noise: expected noise plot
;;    snr: power to sigma ratio plot
;;    nnr: measured noise to expected noise ratio plot
;;    plot_noise: measured noise plot
;;    plot_noise: measured noise plot
;;    pwr_ratio: ratio of powers in 2 savefiles or treat array supplied to power keyword as a ratio
;;
;;  Saving plots keywords
;;   plotfile: name of file to save plot to. If a recognized extension is present, plot will be saved as that type of file. Otherwise a png type will be used
;;   png: flag to save plot as a png. Overridden by plotfile extension if present
;;   eps: flag to save plot as encapsulated postscript. Overridden by plotfile extension if present
;;   pdf: flag to save plot as a pdf. Overridden by plotfile extension if present
;;
;;  Plot window & multi plot keywords
;;   start_multi_params: structure with tags {ncol, nrow, ordering} to indicate the start of a multi plot.
;;     If set, plot positions for future plots on the same page are returned in multi_pos
;;   multi_pos: either an output array of dimension [4, nplots] (if start_multi_params is set) containing plot positions for all plots
;;     or a 4 element vector input containing the plot positions for the current plot
;;     Only one of start_multi_params and multi_pos can be set.
;;   window_num: which idl window number to plot to. Window will be created if it doesn't exist, with a size based on the plot aspect ratio
;;
;;  Plot text keywords
;;   full_title: full plot title. Nothing will be added
;;   title_prefix: start of plot file, descriptors of plot type will be added
;;   no_title: enforce no plot title
;;   note: text to be written at bottom right corner of plot
;;   no_units: don't print units on color bar
;;   charsize: parameter to control text size on plots.
;;
;;  Plot axes and lines keywords
;;   baseline_axis: if set, plot a baseline axis at the top of the plot. Requires kperp_lambda_conv to be present
;;   delay_axis: if set, plot a delay axis on right hand side of the plot. Requires delay_params to be present
;;   cable_length_axis: if set, plot a cable-length reflection axis on right hand side of the plot. Requires delay_params to be present
;;   kperp_linear_axis: if set, use a linear axis for kperp. Either log or linear axes can be used with log or linearly binned data
;;   kpar_linear_axis: if set, use a linear axis for kpar. Either log or linear axes can be used with log or linearly binned data
;;   plot_wedge_line: if set, plot line to indicate wedge & EoR window
;;   wedge_amp: vector wedge amplitudes to plot lines at (ie 1st null, horizon)
;;
;;  Plot appearance keywords
;;   kperp_plot_range: range of kperp values to include [kperp_min, kperp_max]
;;   kpar_plot_range: range of kpar values to include [kpar_min, kpar_max]
;;   data_range: color bar range. plotted array will be clipped to these values.
;;   data_min_abs: minimum non-zero magnitude to plot if sym_log color profile is selected.
;;   color_profile: type of color bar, passed to log_color_calc. allowed values are:['log_cut', 'sym_log', 'abs'] default is 'log_cut'
;;   log_cut_val: minimum value to limit color array to, defaults to data_range[0] or the minimum positive value in the array.
;;     Only applies if color_profile is 'log_cut'.
;;   norm_2D: if set, normalize plot by norm_factor
;;   norm_factor: factor to use for normalizing plots if norm_2d is set. If not set, max of plotted array will be used and passed back in this keyword
;;   cb_size: color bar size (width) in normal units.
;;   margin: 4 element vector of margin sizes (x0, y0, x1, y1) in normal coordinates
;;   cb_margin: 2 element vector of margins on left & right of color bar in normal coordinates

pro kpower_2d_plots, power_savefile, power = power, noise_meas = noise_meas, $
  weights = weights, noise_expval = noise_expval, $
  kperp_edges = kperp_edges, kpar_edges = kpar_edges, kperp_bin = kperp_bin, kpar_bin = kpar_bin, $
  kperp_lambda_conv = kperp_lambda_conv, delay_params = delay_params, hubble_param = hubble_param, $
  plot_options = plot_options, plot_2d_options = plot_2d_options, $
  multi_pos = multi_pos, start_multi_params = start_multi_params, window_num = window_num, $
  plot_weights = plot_weights, plot_noise = plot_noise, plot_sim_noise = plot_sim_noise, $
  plot_simnoise_diff = plot_simnoise_diff, plot_sigma = plot_sigma, $
  plot_exp_noise = plot_exp_noise, $
  snr = snr, nnr = nnr, sim_nnr = sim_nnr, sim_snr = sim_snr, pwr_ratio = pwr_ratio, $
  plot_bin = plot_bin, bin_savefile = bin_savefile, $
  bin_contour = bin_contour, contour_levels = contour_levels, $
  plot_1d_noisefrac = plot_1d_noisefrac, $
  kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
  force_kperp_axis_range = force_kperp_axis_range, force_kpar_axis_range = force_kpar_axis_range, $
  data_range = data_range, data_min_abs = data_min_abs, $
  color_profile = color_profile, log_cut_val = log_cut_val, $
  color_type = color_type, invert_colorbar = invert_colorbar, $
  plotfile = plotfile, png = png, eps = eps, pdf = pdf, $
  no_title = no_title, full_title = full_title, title_prefix = title_prefix, note = note, $
  norm_2d = norm_2d, norm_factor = norm_factor, $
  wedge_amp = wedge_amp, plot_wedge_line = plot_wedge_line, $
  baseline_axis = baseline_axis, delay_axis = delay_axis, $
  cable_length_axis = cable_length_axis, kperp_linear_axis = kperp_linear_axis, $
  kpar_linear_axis = kpar_linear_axis, no_units = no_units, hinv = hinv, $
  charsize = charsize_in, cb_size = cb_size_in, margin=margin_in, $
  cb_margin = cb_margin_in, label_lt_0_oncb = label_lt_0_oncb

  if total([keyword_set(plot_weights), keyword_set(plot_sigma), $
  keyword_set(plot_exp_noise), keyword_set(plot_noise), $
    keyword_set(snr), keyword_set(nnr), keyword_set(pwr_ratio), $
    keyword_set(plot_sim_noise), keyword_set(plot_simnoise_diff), $
    keyword_set(sim_nnr), keyword_set(sim_snr)]) gt 1 then $
    message, 'only one of [plot_noise, plot_sim_noise, plot_simnoise_diff, ' + $
      'plot_sigma, plot_exp_noise, plot_weights, snr, nnr, sim_nnr, sim_snr] ' + $
      'keywords can be set'

  if n_elements(plot_options) gt 0 then begin
    if n_elements(hinv) eq 0 then hinv = plot_options.hinv
    if n_elements(note) eq 0 and tag_exist(plot_options, 'note') then note = plot_options.note
    if n_elements(png) eq 0 then png = plot_options.png
    if n_elements(eps) eq 0 then eps = plot_options.eps
    if n_elements(pdf) eq 0 then pdf = plot_options.pdf
  endif

  if n_elements(plot_2d_options) gt 0 then begin
    if n_elements(plot_wedge_line) eq 0 then plot_wedge_line = plot_2d_options.plot_wedge_line
    if n_elements(kperp_linear_axis) eq 0 then kperp_linear_axis = plot_2d_options.kperp_linear_axis
    if n_elements(kpar_linear_axis) eq 0 then kpar_linear_axis = plot_2d_options.kpar_linear_axis
    if n_elements(baseline_axis) eq 0 then baseline_axis = plot_2d_options.baseline_axis
    if n_elements(delay_axis) eq 0 then delay_axis = plot_2d_options.delay_axis
    if n_elements(cable_length_axis) eq 0 then cable_length_axis = plot_2d_options.cable_length_axis

    if n_elements(kperp_plot_range) eq 0 and tag_exist(plot_2d_options, 'kperp_plot_range') then begin
      kperp_plot_range = plot_2d_options.kperp_plot_range
    endif
    if n_elements(kpar_plot_range) eq 0 and tag_exist(plot_2d_options, 'kpar_plot_range') then begin
      kpar_plot_range = plot_2d_options.kpar_plot_range
    endif

    if n_elements(color_type) eq 0 and tag_exist(plot_2d_options, 'color_type') then begin
      color_type = plot_2d_options.color_type
    endif

    if n_elements(data_range) eq 0 then begin
      case 1 of
        keyword_set(plot_sigma): range_tag_use = 'sigma_range'
        keyword_set(plot_exp_noise): range_tag_use = 'nev_range'
        keyword_set(snr): range_tag_use = 'snr_range'
        keyword_set(plot_noise): range_tag_use = 'noise_range'
        keyword_set(plot_sim_noise): range_tag_use = 'noise_range'
        keyword_set(plot_simnoise_diff): range_tag_use = 'noise_range'
        keyword_set(nnr): range_tag_use = 'nnr_range'
        keyword_set(sim_nnr): range_tag_use = 'nnr_range'
        else: range_tag_use = 'data_range'
      endcase
      if tag_exist(plot_2d_options, range_tag_use) then begin
        tags = strlowcase(tag_names(plot_2d_options))
        tindex = (where(strcmp(tags, range_tag_use) EQ 1, count_tag))[0]
        if count_tag eq 0 then message, 'something went wrong with tag matching'
        data_range = plot_2d_options.(tindex)
      endif
    endif
  endif

  if keyword_set(delay_axis) and keyword_set(cable_length_axis) then begin
    message, 'Only one of delay_axis and cable_length_axis can be set'
  endif

  if n_elements(plotfile) gt 0 or keyword_set(png) or keyword_set(eps) or keyword_set(pdf) then begin
    pub = 1
  endif else begin
    pub = 0
  endelse

  if pub eq 1 then begin
    if not (keyword_set(png) or keyword_set(eps) or keyword_set(pdf)) then begin
      basename = cgRootName(plotfile, directory=directory, extension=extension)

      case extension of
        'eps': eps=1
        'png': png=1
        'pdf': pdf=1
        '': png = 1
        else: begin
          print, 'Unrecognized extension, using png'
          png = 1
        end
      endcase

    endif
    if n_elements(plotfile) eq 0 and n_elements(multi_pos) eq 0 then begin
      plotfile = 'idl_kpower_2d_plots'
      cd, current = current_dir
      print, 'no filename specified for kpower_2d_plots output. Using ' + $
        current_dir + path_sep() + plotfile
    endif

    n_ext_set = 0
    if keyword_set(png) then n_ext_set +=1
    if keyword_set(pdf) then n_ext_set +=1
    if keyword_set(eps) then n_ext_set +=1
    if n_ext_set gt 1 then begin
      print, 'only one of eps, png, pdf keywords can be set, using png'
      eps = 0
      pdf = 0
    endif

    if keyword_set(png) then begin
      plot_exten = '.png'
      delete_ps = 1
    endif else if keyword_set(pdf) then begin
      plot_exten = '.pdf'
      delete_ps = 1
    endif else if keyword_set(eps) then begin
      plot_exten = '.eps'
      delete_ps = 0
    endif
  endif


  if n_elements(window_num) eq 0 then window_num = 1

  if n_elements(start_multi_params) gt 0 and n_elements(multi_pos) gt 0 then begin
    message, 'If start_multi_params are passed, multi_pos cannot be passed ' + $
      'because then it is used as an output to pass back the positions for future plots.'
  endif

  if n_elements(multi_pos) gt 0 then begin
    if n_elements(multi_pos) ne 4 then begin
      message, 'multi_pos must be a 4 element plot position vector'
    endif
    if max(multi_pos) gt 1 or min(multi_pos) lt 0 then begin
      message, 'multi_pos must be in normalized coordinates (between 0 & 1)'
    endif
    if multi_pos[2] le multi_pos[0] or multi_pos[3] le multi_pos[1] then begin
      message, 'In multi_pos, x1 must be greater than x0 and y1 must be greater than y0 '
    endif
  endif


  if keyword_set(pwr_ratio) then begin
    if n_elements(power_savefile) gt 2 or n_elements(power_savefile) eq 1 then begin
      message, 'Only 2 files can be specified in power_savefile if pwr_ratio keyword is set'
    endif
  endif else begin
    if n_elements(power_savefile) gt 1 then begin
      message, 'Only 1 file can be specified in power_savefile unless pwr_ratio keyword is set'
    endif
  endelse

  if keyword_set(pwr_ratio) then begin
    if n_elements(power) eq 0 then begin
      kperp_edges = getvar_savefile(power_savefile[0], 'kperp_edges')
      if total(abs(kperp_edges - getvar_savefile(power_savefile[1], 'kperp_edges'))) ne 0 then begin
        message, 'kperp_edges do not match in savefiles'
      endif
      kpar_edges = getvar_savefile(power_savefile[0], 'kpar_edges')
      if total(abs(kpar_edges - getvar_savefile(power_savefile[1], 'kpar_edges'))) ne 0 then begin
        message, 'kpar_edges do not match in savefiles'
      endif
      kperp_bin = getvar_savefile(power_savefile[0], 'kperp_bin')
      if total(abs(kperp_bin - getvar_savefile(power_savefile[1], 'kperp_bin'))) ne 0 then begin
        message, 'kperp_bin do not match in savefiles'
      endif
      kpar_bin = getvar_savefile(power_savefile[0], 'kpar_bin')
      if total(abs(kpar_bin - getvar_savefile(power_savefile[1], 'kpar_bin'))) ne 0 then begin
        message, 'kpar_bin do not match in savefiles'
      endif
      kperp_lambda_conv = getvar_savefile(power_savefile[0], 'kperp_lambda_conv')
      if total(abs(kperp_lambda_conv - getvar_savefile(power_savefile[1], 'kperp_lambda_conv'))) ne 0 then begin
        message, 'kperp_lambda_conv do not match in savefiles'
      endif
      delay_params = getvar_savefile(power_savefile[0], 'delay_params')
      if total(abs(delay_params - getvar_savefile(power_savefile[1], 'delay_params'))) ne 0 then begin
        message, 'delay_params do not match in savefiles'
      endif
      hubble_param = getvar_savefile(power_savefile[0], 'hubble_param')
      if total(abs(hubble_param - getvar_savefile(power_savefile[1], 'hubble_param'))) ne 0 then begin
        message, 'hubble_param do not match in savefiles'
      endif

      power1 = getvar_savefile(power_savefile[0], 'power')
      power2 = getvar_savefile(power_savefile[1], 'power')

      power_use = power1 / power2
      wh0 = where(power2 eq 0, count0)
      if count0 gt 0 then power_use[wh0] = 0
    endif else power_use = power

    plot_type = 'power_ratio'

  endif else if n_elements(power_savefile) gt 0 then begin
    restore, power_savefile
    if n_elements(noise) then begin
      ;; backwards compatibility: noise_meas used to just be called noise
      noise_meas = noise
      undefine, noise
    endif
  endif

  if keyword_set(snr) then begin
    if n_elements(weights) eq 0 then message, 'weights is undefined'
    power_use = power * sqrt(weights)

    plot_type = 'snr'
  endif

  if keyword_set(plot_weights) then begin
    if n_elements(weights) eq 0 then message, 'weights is undefined'
    power_use = weights
    plot_type = 'weight'
  endif

  if keyword_set(plot_sigma) then begin
    if n_elements(weights) eq 0 then message, 'weights is undefined'
    power_use = 1/sqrt(weights)
    wh_wt0 = where(weights eq 0, count_wt0)
    if count_wt0 gt 0 then power_use[wh_wt0 ] = 0
    plot_type = 'sigma'
  endif

  if keyword_set(plot_exp_noise) then begin
    if n_elements(noise_expval) eq 0 then message, 'noise_expval is undefined'
    power_use = noise_expval
    plot_type = 'exp_noise'
  endif

  if keyword_set(plot_noise) then begin
    if n_elements(noise_meas) eq 0 then message, 'noise is undefined'
    power_use = noise_meas
    plot_type = 'noise_meas'
  endif

  if keyword_set(plot_sim_noise) then begin
    if n_elements(sim_noise) eq 0 then message, 'sim_noise is undefined'
    power_use = sim_noise
    plot_type = 'sim_noise'
  endif

  if keyword_set(plot_simnoise_diff) then begin
    if n_elements(sim_noise_diff) eq 0 then message, 'sim_noise_diff is undefined'
    power_use = sim_noise_diff
    plot_type = 'sim_noise_diff'
  endif

  if keyword_set(nnr) then begin
    if n_elements(noise_meas) eq 0 then message, 'noise_meas is undefined'
    if n_elements(noise_expval) eq 0 then message, 'noise_expval is undefined'
    power_use = noise_meas / noise_expval
    wh_err0 = where(noise_expval eq 0, count_err0)
    if count_err0 gt 0 then power_use[wh_err0] = 0

    plot_type = 'nnr'
  endif

  if keyword_set(sim_snr) then begin
    if n_elements(sim_noise) eq 0 then message, 'sim_noise is undefined in this file'
    if n_elements(weights) eq 0 then message, 'weights is undefined'
    power_use = abs(sim_noise) * sqrt(weights)
    plot_type = 'sim_snr'
  endif

  if keyword_set(sim_nnr) then begin
    if n_elements(sim_noise_diff) eq 0 then message, 'sim_noise_diff is undefined in this file'
    if n_elements(noise_expval) eq 0 then message, 'noise_expval is undefined'
    power_use = sim_noise_diff / noise_expval
    wh_err0 = where(noise_expval eq 0, count_err0)
    if count_err0 gt 0 then power_use[wh_err0] = 0

    plot_type = 'sim_nnr'
  endif

  if keyword_set(plot_bin) then begin
    power_use = bin_1to2d_ave
    plot_type = 'bin'
  endif

  if keyword_set(plot_1d_noisefrac) then begin
    power_use = noise_frac_1to2d
    plot_type = 'noise_frac'
  endif

  if n_elements(plot_type) eq 0 then begin
    plot_type = 'power'
    power_use = power
  endif

  if max(abs(imaginary(power_use))) gt 0 then begin
    print, 'data is complex, showing real part'
    power_use = real_part(power_use)
  endif


  if keyword_set(bin_contour) then begin
    if plot_type ne 'bin' then begin
      if n_elements(bin_savefile) eq 0 then begin
        message, 'bin_savefile must be supplied if bin_contour is set'
      endif

      kperp_edges = getvar_savefile(power_savefile, 'kperp_edges')
      if total(abs(kperp_edges - getvar_savefile(bin_savefile, 'kperp_edges'))) ne 0 then begin
        message, 'bin_savefile kperp_edges do not match power_savefile'
      endif
      kpar_edges = getvar_savefile(power_savefile, 'kpar_edges')
      if total(abs(kpar_edges - getvar_savefile(bin_savefile, 'kpar_edges'))) ne 0 then begin
        message, 'bin_savefile kpar_edges do not match power_savefile'
      endif
      kperp_bin = getvar_savefile(power_savefile, 'kperp_bin')
      if total(abs(kperp_bin - getvar_savefile(bin_savefile, 'kperp_bin'))) ne 0 then begin
        message, 'bin_savefile kperp_bin does not match power_savefile'
      endif
      kpar_bin = getvar_savefile(power_savefile, 'kpar_bin')
      if total(abs(kpar_bin - getvar_savefile(bin_savefile, 'kpar_bin'))) ne 0 then begin
        message, 'bin_savefile kpar_bin does not match power_savefile'
      endif
      kperp_lambda_conv = getvar_savefile(power_savefile, 'kperp_lambda_conv')
      if total(abs(kperp_lambda_conv - getvar_savefile(bin_savefile, 'kperp_lambda_conv'))) ne 0 then begin
        message, 'bin_savefile kperp_lambda_conv does not match power_savefile'
      endif
      delay_params = getvar_savefile(power_savefile, 'delay_params')
      if total(abs(delay_params - getvar_savefile(bin_savefile, 'delay_params'))) ne 0 then begin
        message, 'bin_savefile delay_params does not match power_savefile'
      endif
      hubble_param = getvar_savefile(power_savefile, 'hubble_param')
      if total(abs(hubble_param - getvar_savefile(bin_savefile, 'hubble_param'))) ne 0 then begin
        message, 'bin_savefile hubble_param does not match power_savefile'
      endif
    endif else begin
      if n_elements(bin_savefile) eq 0 then bin_savefile = power_savefile
    endelse

    bin = getvar_savefile(bin_savefile, 'bin_1to2d')
  endif

  n_kperp = n_elements(kperp_edges) - 1
  n_kpar = n_elements(kpar_edges) -1

  kperp_edges_use = kperp_edges
  kpar_edges_use = kpar_edges

  ;; Check whether binning is log or not
  log_bins = [0, 0]
  kperp_diffs = (kperp_edges_use - shift(kperp_edges_use, 1))[1:*]
  if kperp_edges_use[0] eq 0 or (abs(mean(kperp_diffs)-kperp_diffs[0]) gt 1e-5) then begin
    kperp_diffs = kperp_diffs[1:*] ;; lowest bin may have lower edge set to zero
  endif
  if total(abs(kperp_diffs - kperp_diffs[0])) gt n_kperp*2e-5 then log_bins[0] = 1
  kpar_diffs = (kpar_edges_use - shift(kpar_edges_use, 1))[1:*]
  if kpar_edges_use[0] eq 0 or (abs(mean(kpar_diffs)-kpar_diffs[0]) gt 1e-5) then begin
    kpar_diffs = kpar_diffs[1:*] ;; lowest bin may have lower edge set to zero
  endif
  if total(abs(kpar_diffs - kpar_diffs[0])) gt n_kpar*2e-5 then log_bins[1] = 1

  if n_elements(kpar_bin) gt 0 then kpar_bin_use = kpar_bin else begin
    ;; calculate kpar_bin from kpar_edges, don't use lowest bin because it might have lower edge set to zero
    if log_bins[1] eq 0 then kpar_bin = kpar_edges_use[2] - kpar_edges_use[1] $
    else kpar_bin = alog10(kpar_edges_use[2]) - alog10(kpar_edges_use[1])
    kpar_bin_use = kpar_bin
  endelse

  if keyword_set(hinv) then begin
    kperp_edges_use = kperp_edges_use / hubble_param
    kpar_edges_use = kpar_edges_use / hubble_param
    kpar_bin_use = kpar_bin_use / hubble_param
  endif

  if n_elements(kperp_plot_range) eq 0 then kperp_plot_range = minmax(kperp_edges_use)
  if n_elements(kpar_plot_range) eq 0 then kpar_plot_range = minmax(kpar_edges_use)

  if n_elements(plot_2d_options) gt 0 then begin
    if not tag_exist(plot_2d_options, 'kperp_plot_range') then begin
      plot_2d_options = create_plot_2d_options(plot_2d_options=plot_2d_options, $
        kperp_plot_range = kperp_plot_range)
    endif
    if not tag_exist(plot_2d_options, 'kpar_plot_range') then begin
      plot_2d_options = create_plot_2d_options(plot_2d_options=plot_2d_options, $
        kpar_plot_range = kpar_plot_range)
    endif
  endif

  units_str = ''
  case plot_type of
    'power': begin
      if keyword_set(hinv) then begin
        power_use = power_use * (hubble_param)^3d
        units_str = textoidl(' (mK^2 !8h!X^{-3} Mpc^3)', font = font)
      endif else units_str = textoidl(' (mK^2 Mpc^3)', font = font)

      if n_elements(color_type) eq 0 then color_type = 'log'
      plot_title = textoidl('P_k', font = font)
      if pub then plotfile_add = '_2dkpower' + plot_exten
    end
    'weight': begin
      units_str = ''
      plot_title = 'Weights'
      if n_elements(color_type) eq 0 then color_type = 'log'
      if pub then plotfile_add = '_2dweight' + plot_exten
    end
    'sigma': begin
      if keyword_set(hinv) then begin
        power_use = power_use * (hubble_param)^3d
        units_str = textoidl(' (mK^2 !8h!X^{-3} Mpc^3)', font = font)
      endif else units_str = textoidl(' (mK^2 Mpc^3)', font = font)

      if n_elements(color_type) eq 0 then color_type = 'log'
      plot_title = 'Error (sigma)'
      if pub then plotfile_add = '_2dsigma' + plot_exten
    end
    'exp_noise': begin
      if keyword_set(hinv) then begin
        power_use = power_use * (hubble_param)^3d
        units_str = textoidl(' (mK^2 !8h!X^{-3} Mpc^3)', font = font)
      endif else units_str = textoidl(' (mK^2 Mpc^3)', font = font)

      if n_elements(color_type) eq 0 then color_type = 'log'
      plot_title = 'Expected Noise'
      if pub then plotfile_add = '_2dnoise_expval' + plot_exten
    end
    'snr': begin
      units_str = ''
      if n_elements(color_type) eq 0 then color_type = 'log'
      plot_title = 'SNR (' + textoidl('P_k/\sigma', font = font) + ')'
      if pub then plotfile_add = '_2dsnr' + plot_exten
    end
    'noise_meas': begin
      if keyword_set(hinv) then begin
        power_use = power_use * (hubble_param)^3d
        units_str = textoidl(' (mK^2 !8h!X^{-3} Mpc^3)', font = font)
      endif else units_str = textoidl(' (mK^2 Mpc^3)', font = font)

      if n_elements(color_type) eq 0 then color_type = 'log'
      plot_title = 'Observed Noise'
      if pub then plotfile_add = '_2dnoise' + plot_exten
    end
    'sim_noise': begin
      if keyword_set(hinv) then begin
        power_use = power_use * (hubble_param)^3d
        units_str = textoidl(' (mK^2 !8h!X^{-3} Mpc^3)', font = font)
      endif else units_str = textoidl(' (mK^2 Mpc^3)', font = font)

      if n_elements(color_type) eq 0 then color_type = 'log'
      color_profile = 'abs'
      plot_title = '|Simulated Noise (cross)|'
      if pub then plotfile_add = '_2dsimnoise' + plot_exten
    end
    'sim_noise_diff': begin
      if keyword_set(hinv) then begin
        power_use = power_use * (hubble_param)^3d
        units_str = textoidl(' (mK^2 !8h!X^{-3} Mpc^3)', font = font)
      endif else units_str = textoidl(' (mK^2 Mpc^3)', font = font)

      if n_elements(color_type) eq 0 then color_type = 'log'
      plot_title = 'Simulated Noise Diff'
      if pub then plotfile_add = '_2dsimnoisediff' + plot_exten
    end
    'nnr': begin
      units_str = ''
      if n_elements(color_type) eq 0 then color_type = 'log'
      plot_title = 'Noise Ratio (' + textoidl('N_O/N_E', font = font) + ')'
      if pub then plotfile_add = '_2dnnr' + plot_exten
    end
    'sim_snr': begin
      units_str = ''
      if n_elements(color_type) eq 0 then color_type = 'log'
      ;color_profile = 'abs'
      plot_title = '|Sim Noise SNR| (|' + textoidl('Sim Noise/\sigma', font = font) + '|)'
      if pub then plotfile_add = '_2dsimsnr' + plot_exten
    end
    'sim_nnr': begin
      units_str = ''
      if n_elements(color_type) eq 0 then color_type = 'log'
      plot_title = 'Sim Noise diff Ratio (' + textoidl('N_{sim}/N_E', font = font) + ')'
      if pub then plotfile_add = '_2dsimnnr' + plot_exten
    end
    'power_ratio': begin
      units_str = ''
      if n_elements(color_type) eq 0 then color_type = 'log'
      plot_title = 'Power Ratio'
      if pub then plotfile_add = '_2dkpower' + plot_exten
    end
    'bin': begin
      units_str = '1D bin number * fill fraction'
      if n_elements(color_type) eq 0 then color_type = 'integer'
      plot_title = ''
      if pub then plotfile_add = '_1to2dbin' + plot_exten
    end
    'noise_frac': begin
      units_str = 'Noise fraction in 1d bin'
      if n_elements(color_type) eq 0 then color_type = 'linear'
      plot_title = ''
      if pub then plotfile_add = '_1to2dnoisefrac' + plot_exten
    end
  endcase

  if pub then begin
    if n_elements(plotfile) eq 0 then begin
      plotfile = strsplit(power_savefile[0], '.idlsave', /regex, /extract) + plotfile_add
    endif else begin
      if strcmp(strmid(plotfile, strlen(plotfile)-4), plot_exten, /fold_case) eq 0 then begin
        plotfile = plotfile + plot_exten
      endif
    endelse
  endif

  if keyword_set(no_units) then units_str = ''

  wh_kperp_inrange = where(kperp_edges_use[0:n_kperp-1] ge kperp_plot_range[0] and $
    kperp_edges_use[1:*] le kperp_plot_range[1], n_kperp_plot)
  wh_kpar_inrange = where(kpar_edges_use[0:n_kpar-1] ge kpar_plot_range[0] and $
    kpar_edges_use[1:*] le kpar_plot_range[1], n_kpar_plot)

  if n_kperp_plot eq 0 or n_kpar_plot eq 0 then message, 'No data in plot k range'

  if n_kperp_plot ne n_kperp then begin
    power_use = power_use[wh_kperp_inrange, *]
    if keyword_set(bin_contour) then bin = bin[wh_kperp_inrange, *]
    temp = [wh_kperp_inrange, wh_kperp_inrange[n_kperp_plot-1]+1]
    kperp_edges_use =kperp_edges_use[temp]
    n_kperp = n_kperp_plot
  endif

  lin_delay_kpar_slope = (delay_params[1] - delay_params[0])/(max(kpar_edges_use) - kpar_bin_use)
  lin_delay_kpar_intercept = delay_params[0] / (lin_delay_kpar_slope * kpar_bin_use)
  linear_delay_edges = lin_delay_kpar_slope * kpar_edges_use + lin_delay_kpar_intercept

  log_delay_kpar_slope = (alog10(delay_params[1]) - alog10(delay_params[0]))/(alog10(max(kpar_edges_use)) - alog10(kpar_bin_use))
  log_delay_kpar_intercept = alog10(delay_params[0]) / (log_delay_kpar_slope * alog10(kpar_bin_use))
  log_delay_edges = 10^(log_delay_kpar_slope * alog10(kpar_edges_use) + log_delay_kpar_intercept)

  if n_kpar_plot ne n_kpar then begin
    power_use = power_use[*, wh_kpar_inrange]
    if keyword_set(bin_contour) then bin = bin[*, wh_kpar_inrange]

    temp = [wh_kpar_inrange, wh_kpar_inrange[n_kpar_plot-1]+1]

    kpar_edges_use = kpar_edges_use[temp]

    linear_delay_edges = linear_delay_edges[temp]
    log_delay_edges = log_delay_edges[temp]

    n_kpar = n_kpar_plot
  endif

  if max(abs(power_use)) eq 0 then begin
    print, 'power is entirely zero.'
    no_plot = 1
  endif

  if keyword_set(norm_2d) then begin
    if n_elements(norm_factor) eq 0 then norm_factor = 1/max(power_use)
    power_use = power_use * norm_factor
  endif

  ;; check whether we need to have varying bin sizes (log/linear options don't match)
  log_axes = [1,1]
  if keyword_set(kperp_linear_axis) then log_axes[0] = 0
  if keyword_set(kpar_linear_axis) then log_axes[1] = 0

  if total(abs(log_bins-log_axes)) ne 0 then begin
    ;; need to make a new image array with varying bin sizes
    if log_bins[0] ne log_axes[0] then begin
      ;; kperp direction
      if log_bins[0] eq 0 then begin
        ;; linear binning, log axes
        wh_kperp0 = where(kperp_edges_use le 0, count_kperp0, complement = wh_kperp_good)
        if count_kperp0 gt 1 then stop

        kperp_log_edges = alog10(kperp_edges_use)
        if count_kperp0 eq 1 then begin
          kperp_log_diffs = (kperp_log_edges[1:*] - shift(kperp_log_edges[1:*], 1))[1:*]
          kperp_log_diffs = [kperp_log_diffs[0], kperp_log_diffs]
          kperp_log_edges[wh_kperp0] = kperp_log_edges[wh_kperp0+1] - kperp_log_diffs[wh_kperp0]
        endif else begin
          kperp_log_diffs = (kperp_log_edges - shift(kperp_log_edges, 1))[1:*]
        endelse

        image_kperp_delta = min(kperp_log_diffs)
        kperp_bin_widths = round(kperp_log_diffs / image_kperp_delta)
      endif else begin
        ;; log binning, linear axes
        kperp_diffs = (kperp_edges_use[1:*] - shift(kperp_edges_use[1:*], 1))[1:*]
        image_kperp_delta = min(kperp_diffs)/2d
        kperp_bin_widths = round(kperp_diffs / image_kperp_delta)
      endelse
      rebin_x = 1
    endif else begin
      ;; axes and binning agree
      if log_axes[0] eq 1 then kperp_log_edges = alog10(kperp_edges_use)
      rebin_x = 0
    endelse

    if log_bins[1] ne log_axes[1] then begin
      ;; kpar direction
      if log_bins[1] eq 0 then begin
        ;; linear binning, log axes
        wh_kpar0 = where(kpar_edges_use le 0, count_kpar0, complement = wh_kpar_good)
        if count_kpar0 gt 1 then stop

        kpar_log_edges = alog10(kpar_edges_use)
        delay_log_edges = alog10(linear_delay_edges)
        if count_kpar0 eq 1 then begin
          kpar_log_diffs = (kpar_log_edges[1:*] - shift(kpar_log_edges[1:*], 1))[1:*]
          kpar_log_diffs = [kpar_log_diffs[0], kpar_log_diffs]
          kpar_log_edges[wh_kpar0] = kpar_log_edges[wh_kpar0+1] - kpar_log_diffs[wh_kpar0]

          delay_log_diffs = (delay_log_edges[1:*] - shift(delay_log_edges[1:*], 1))[1:*]
          delay_log_diffs = [delay_log_diffs[0], delay_log_diffs]
          delay_log_edges[wh_kpar0] = delay_log_edges[wh_kpar0+1] - delay_log_diffs[wh_kpar0]

        endif else kpar_log_diffs = (kpar_log_edges - shift(kpar_log_edges, 1))[1:*]

        image_kpar_delta = min(kpar_log_diffs)/2d
        kpar_bin_widths = round(kpar_log_diffs / image_kpar_delta)
      endif else begin
        ;; log binning, linear axes
        kpar_diffs = (kpar_edges_use[1:*] - shift(kpar_edges_use[1:*], 1))[1:*]
        image_kpar_delta = min(kpar_diffs)/2d
        kpar_bin_widths = round(kpar_diffs / image_kpar_delta)

        delay_diffs = (log_delay_edges[1:*] - shift(log_delay_edges[1:*], 1))[1:*]
        delay_edges = findgen(n_kpar)*delay_diffs + log_delay_edges[0]
      endelse
      rebin_y = 1
    endif else begin
      ;; axes agree
      if log_bins[1] eq 1 then begin
        kpar_log_edges = alog10(kpar_edges_use)
        delay_log_edges = log_delay_edges
      endif else delay_edges = linear_delay_edges
      rebin_y = 0
    endelse

    if rebin_x eq 1 then begin
      ;; now get width for each input bin in image array
      nkperp_image = total(kperp_bin_widths)

      h_kperp = histogram(total(kperp_bin_widths,/cumulative)-1,binsize=1, min=0, $
        reverse_indices=ri_kperp)
      undefine, h_kperp
      kperp_inds = ri_kperp[0:nkperp_image-1]-ri_kperp[0]
    endif else begin
      nkperp_image = n_kperp
      kperp_inds = indgen(nkperp_image)
    endelse

    if rebin_y eq 1 then begin
      nkpar_image = total(kpar_bin_widths)

      h_kpar = histogram(total(kpar_bin_widths,/cumulative)-1,binsize=1, min=0, $
        reverse_indices=ri_kpar)
      undefine, h_kpar
      kpar_inds = rebin(reform(ri_kpar[0:nkpar_image-1]-ri_kpar[0], 1, nkpar_image), $
        nkperp_image, nkpar_image)
    endif else begin
      nkpar_image = n_kpar
      kpar_inds = rebin(reform(indgen(nkpar_image), 1, nkpar_image), nkperp_image, nkpar_image)
    endelse

    kperp_inds = rebin(kperp_inds, nkperp_image, nkpar_image)
    power_plot = power_use[kperp_inds, kpar_inds]
    if keyword_set(bin_contour) then begin
      bin_plot = bin[kperp_inds, kpar_inds]
      if log_axes[0] eq 0 then begin
        bin_kperp = findgen(nkperp_image)*(max(kperp_edges_use) - min(kperp_edges_use))/nkperp_image + min(kperp_edges_use)
      endif else begin
        bin_kperp = 10^(findgen(nkperp_image)*(max(kperp_log_edges) - min(kperp_log_edges))/nkperp_image + min(kperp_log_edges))
      endelse
      if log_axes[1] eq 0 then begin
        bin_kpar = findgen(nkpar_image)*(max(kpar_edges_use) - min(kpar_edges_use))/nkpar_image + min(kpar_edges_use)
      endif else begin
        bin_kpar = 10^(findgen(nkpar_image)*(max(kpar_log_edges) - min(kpar_log_edges))/nkpar_image + min(kpar_log_edges))
      endelse
    endif

    ;; now expand array in any non-rebinned direction to prevent interpolation
    if rebin_x eq 0 or nkperp_image lt 15 or nkperp_image lt 0.1*nkpar_image then begin
      power_plot = congrid(power_plot, nkperp_image*20, nkpar_image)
    endif
    if rebin_y eq 0 or nkpar_image lt 15  or nkpar_image lt 0.1*nkperp_image then begin
      power_plot = congrid(power_plot, nkperp_image, nkpar_image*20)
    endif
    if keyword_set(bin_contour) then begin
      if rebin_x eq 0 or nkperp_image lt 15 then begin
        bin_plot = congrid(bin_plot, nkperp_image*10, nkpar_image)
        if log_axes[0] eq 0 then begin
          bin_kperp = findgen(nkperp_image*10)*(max(kperp_edges_use) - min(kperp_edges_use))/nkperp_image + min(kperp_edges_use)
        endif else begin
          bin_kperp = 10^(findgen(nkperp_image*10)*(max(kperp_log_edges) - min(kperp_log_edges))/(nkperp_image*10) + min(kperp_log_edges))
        endelse
      endif
      if rebin_y eq 0or nkpar_image lt 15 then begin
        bin_plot = congrid(bin_plot, nkperp_image, nkpar_image*10)
        if log_axes[1] eq 0 then begin
          bin_kpar = findgen(nkpar_image*10)*(max(kpar_edges_use) - min(kpar_edges_use))/nkpar_image + min(kpar_edges_use)
        endif else begin
          bin_kpar = 10^(findgen(nkpar_image*10)*(max(kpar_log_edges) - min(kpar_log_edges))/(nkpar_image*10) + min(kpar_log_edges))
        endelse
      endif
    endif
  endif else begin
    ;; axes & binning agree for both directions
    ;; expand image array to prevent interpolation in postscript
    power_plot = congrid(power_use, n_kperp*10, n_kpar*10)
    if log_axes[0] eq 1 then kperp_log_edges = alog10(kperp_edges_use)
    if log_axes[1] eq 1 then begin
      kpar_log_edges = alog10(kpar_edges_use)
      delay_log_edges = log_delay_edges
    endif else delay_edges = linear_delay_edges

    if keyword_set(bin_contour) then begin
      bin_plot = congrid(bin, n_kperp*10, n_kpar*10)
      if log_axes[0] eq 0 then begin
        bin_kperp = findgen(nkperp_image*10)*(max(kperp_edges_use) - min(kperp_edges_use))/nkperp_image + min(kperp_edges_use)
      endif else begin
        bin_kperp = 10^(findgen(nkperp_image*10)*(max(kperp_log_edges) - min(kperp_log_edges))/nkperp_image + min(kperp_log_edges))
      endelse
      if log_axes[1] eq 0 then begin
        bin_kpar = findgen(nkpar_image*10)*(max(kpar_edges_use) - min(kpar_edges_use))/nkpar_image + min(kpar_edges_use)
      endif else begin
        bin_kpar = 10^(findgen(nkpar_image*10)*(max(kpar_log_edges) - min(kpar_log_edges))/nkpar_image + min(kpar_log_edges))
      endelse
    endif
  endelse


  tvlct, r, g, b, /get

  background_color = 'white'
  annotate_color = 'black'

  if not keyword_set(no_plot) then begin
    case color_type of
      'integer': begin
        color_range = [0, ceil(max(power_plot))]
        n_colors = color_range[1] - color_range[0] + 1

        if n_colors gt 256 then begin
          print, 'Too many bins to color them all differently.'
          n_colors = 256
        endif

        power_log_norm = power_plot

        max_bin = ceil(max(power_plot))
        cb_tick_size = ceil(max_bin/8.)
        cb_tick_vals = indgen(8) * cb_tick_size

        wh_large = where(cb_tick_vals gt max_bin, count_large, complement = wh_good)
        if count_large gt 0 then cb_tick_vals = cb_tick_vals[wh_good]
        cb_ticknames = number_formatter(cb_tick_vals)
        cb_ticks = cb_tick_vals + 0.5

        if max_bin gt max(cb_tick_vals) then begin
          cb_tick_vals = [cb_tick_vals, max_bin]
          cb_ticknames = [cb_ticknames, ' ']
          cb_ticks = [cb_ticks, color_range[1]+1]
        endif

        if not keyword_set(invert_colorbar) then begin
          cgloadct, 25, /brewer, /reverse, BOTTOM = 0, NCOLORS = n_colors, clip = [0, 235]
        endif else begin
          cgloadct, 25, /brewer, BOTTOM = 0, NCOLORS = n_colors, clip = [20, 255]
        endelse

      end
      'log': begin
        log_color_calc, power_plot, power_log_norm, cb_ticks, cb_ticknames, $
          color_range, n_colors, data_range = data_range, color_profile = color_profile, $
          log_cut_val = log_cut_val, min_abs = data_min_abs, oob_low = oob_low, $
          label_lt_0 = label_lt_0_oncb, invert_colorbar = invert_colorbar
        end
      'linear': begin
        if n_elements(data_range) eq 0 then data_range = minmax(power_plot)
        color_range = [0, 255]
        n_colors = color_range[1] - color_range[0] + 1

        if n_elements(color_profile) gt 0 then begin
          if color_profile eq 'sym_log' then begin
            color_profile_use = 'sym_lin'
          endif
        endif else begin
          color_profile_use = 'lin'
        endelse

        if color_profile_use eq 'sym_lin' then begin
          data_range = [-1, 1] * max(abs(data_range))

          n_pos_neg_colors = floor((n_colors-1)/2)
          zero_color = n_pos_neg_colors
          neg_color_range = zero_color-1 + [-1*(n_pos_neg_colors-1), 0]
          pos_color_range = zero_color+1 + [0,n_pos_neg_colors-1]

          if (n_pos_neg_colors*2. + 1) lt n_colors then begin
            ndiff = n_colors - (n_pos_neg_colors*2. + 1)
            n_colors = n_colors - ndiff
            color_range[1] = color_range[1] - ndiff
          endif

          power_log_norm = power_plot * 0.
          wh_pos = where(power_plot gt 0d, count_pos)
          if count_pos gt 0 then begin
            power_log_norm[wh_pos] = power_plot[wh_pos]*(pos_color_range[1] - pos_color_range[0]) / $
              (data_range[1]-0) + pos_color_range[0]
          endif
          wh_neg = where(power_plot lt 0d, count_neg)
          if count_neg gt 0 then begin
            power_log_norm[wh_neg] = (power_plot[wh_neg]-data_range[0])*(neg_color_range[1] - neg_color_range[0]) / $
              (0-data_range[0]) + neg_color_range[0]
          endif
          wh_zero = where(power_plot eq 0, count_zero)
          if count_zero gt 0 then begin
            power_log_norm[wh_zero] = zero_color
          endif

        endif else begin
          power_log_norm = (power_plot-data_range[0])*(color_range[1] - color_range[0]) / $
            (data_range[1]-data_range[0]) + color_range[0]
        endelse

        wh_low = where(power_plot lt data_range[0], count_low)
        if count_low gt 0 then power_log_norm[wh_low] = color_range[0]
        wh_high = where(power_plot gt data_range[1], count_high)
        if count_high gt 0 then power_log_norm[wh_high] = color_range[1]

        cb_tick_size = (data_range[1] - data_range[0])/8.
        cb_tick_vals = findgen(12) * cb_tick_size + data_range[0]

        wh_large = where((cb_tick_vals - data_range[1]) gt 0.1, count_large, complement = wh_good)
        if count_large gt 0 then cb_tick_vals = cb_tick_vals[wh_good]

        if abs(alog10(data_range[1] - data_range[0])) > 3 then begin
          format = '(e8.1)'
        endif else begin
          format = '(f9.5)'
        endelse

        cb_ticknames = number_formatter(cb_tick_vals, format = format, /print_exp)
        cb_ticks = (cb_tick_vals - cb_tick_vals[0])*(color_range[1] - color_range[0]) / $
          (data_range[1]-data_range[0])

        if data_range[1] gt max(cb_tick_vals) then begin
          cb_tick_vals = [cb_tick_vals, data_range[1]]
          cb_ticknames = [cb_ticknames, ' ']
          cb_ticks = [cb_ticks, color_range[1]]
        endif

        if color_profile_use eq 'sym_lin' then begin
          if keyword_set(invert_colorbar) then begin
            cgLoadCT, 13, /brewer, /reverse, clip=[20, 220], bottom=0, ncolors=n_pos_neg_colors
            cgloadct, 0, clip = [255, 255], bottom = zero_color, ncolors = 1
            cgLoadCT, 16, /brewer, clip=[20, 220], bottom=zero_color+1, ncolors=n_pos_neg_colors
          endif else begin
            cgLoadCT, 16, /brewer, /reverse, clip=[20, 220], bottom=0, ncolors=n_pos_neg_colors
            cgloadct, 0, clip = [255, 255], bottom = zero_color, ncolors = 1
            cgLoadCT, 13, /brewer, clip=[20, 220], bottom=zero_color+1, ncolors=n_pos_neg_colors
          endelse
        endif else begin
          if not keyword_set(invert_colorbar) then begin
            cgloadct, 25, /brewer, /reverse, BOTTOM = 0, NCOLORS = n_colors, clip = [0, 235]
          endif else begin
            cgloadct, 25, /brewer, BOTTOM = 0, NCOLORS = n_colors, clip = [20, 255]
          endelse
        endelse
      end
    endcase
  endif

  if keyword_set(bin_contour) then begin
    if n_elements(contour_levels) then begin
      if size(contour_levels, /type) then begin
        if contour_levels eq 'all' then levels = indgen(ceil(max(bin_plot))) + 1 $
        else message, 'contour_levels is an unrecognized string'
      endif else levels = contour_levels
    endif else begin
      levels = (indgen(6)+1) * ceil(max(bin_plot)/6.)
      if min(levels) gt 1 then levels = [1, levels]
      wh_large = where(levels gt max(power_plot), count_large, complement = wh_good)
      if count_large gt 0 then levels = levels[wh_good]
    endelse
  endif

  if keyword_set(force_kperp_axis_range) or keyword_set(force_kpar_axis_range) then begin
    ;; we want to draw axes outside the range of the plotted data
    image_relative_pos = float([0, 0, 1, 1])

    if keyword_set(force_kperp_axis_range) then begin
      if min(force_kperp_axis_range) gt min(kperp_edges_use) $
          or max(force_kperp_axis_range) lt max(kperp_edges_use) then begin
        message, 'force_kperp_axis_range is used for specifying axes outside the ' + $
          'range of the plotted data. To limit the plotted data range, use kperp_plot_range'
      endif

      if log_axes[0] eq 0 then begin
        image_relative_pos[0] = (min(kperp_edges_use) - min(force_kperp_axis_range))/(max(force_kperp_axis_range) - min(force_kperp_axis_range))
        image_relative_pos[2] = (max(kperp_edges_use) - min(force_kperp_axis_range))/(max(force_kperp_axis_range) - min(force_kperp_axis_range))
      endif else begin
        log_kperp_axis = alog10(force_kperp_axis_range)
        image_relative_pos[0] = (min(kperp_log_edges) - min(log_kperp_axis))/(max(log_kperp_axis) - min(log_kperp_axis))
        image_relative_pos[2] = (max(kperp_log_edges) - min(log_kperp_axis))/(max(log_kperp_axis) - min(log_kperp_axis))
      endelse

    endif


    if keyword_set(force_kpar_axis_range) then begin
      if min(force_kpar_axis_range) gt min(kpar_edges_use) $
          or max(force_kpar_axis_range) lt max(kpar_edges_use) then begin
        message, 'force_kpar_axis_range is used for specifying axes outside the ' + $
          'range of the plotted data. To limit the plotted data range, use kpar_plot_range'
      endif

      if log_axes[0] eq 0 then begin
        image_relative_pos[1] = (min(kpar_edges_use) - min(force_kpar_axis_range))/(max(force_kpar_axis_range) - min(force_kpar_axis_range))
        image_relative_pos[3] = (max(kpar_edges_use) - min(force_kpar_axis_range))/(max(force_kpar_axis_range) - min(force_kpar_axis_range))
      endif else begin
        log_kpar_axis = alog10(force_kpar_axis_range)
        image_relative_pos[1] = (min(kpar_log_edges) - min(log_kpar_axis))/(max(log_kpar_axis) - min(log_kpar_axis))
        image_relative_pos[3] = (max(kpar_log_edges) - min(log_kpar_axis))/(max(log_kpar_axis) - min(log_kpar_axis))
      endelse

    endif
  endif

  screen_size = get_screen_size()
  max_xsize = screen_size[0]
  max_ysize = screen_size[1]
  base_size = 600
  if n_elements(multi_pos) eq 4 then begin
    ;; work out positions scaled to the area allowed in multi_pos with proper aspect ratio
    multi_xlen = (multi_pos[2]-multi_pos[0])
    multi_ylen = (multi_pos[3]-multi_pos[1])
    multi_center = [multi_pos[0] + multi_xlen/2d, multi_pos[1] + multi_ylen/2d]

    multi_size = [!d.x_vsize*multi_xlen, !d.y_vsize*multi_ylen]
  endif

  ;; Work out plot & colorbar positions
  ;; in units of plot area (incl. margins)
  if n_elements(cb_size_in) eq 0 then cb_size = 0.025 else cb_size = cb_size_in
  if n_elements(margin_in) lt 4 then begin
    margin = [0.2, 0.2, 0.02, 0.1]
    if keyword_set(baseline_axis) and not keyword_set(no_title) then margin[3] = 0.15
    if keyword_set(delay_axis) or keyword_set(cable_length_axis) then begin
      if keyword_set(kpar_linear_axis) then begin
        ;; need more space because not just powers of 10 labels
        margin[2] = 0.2
      endif else begin
        margin[2] = 0.1
      endelse
    endif
  endif else margin = margin_in

  if n_elements(cb_margin_in) lt 2 then begin
    cb_margin = [0.2, 0.02]
  endif else cb_margin = cb_margin_in

  plot_pos = [margin[0], margin[1], (1-cb_margin[1]-cb_size-cb_margin[0]-margin[2]), $
    (1-margin[3])]
  cb_pos = [(1-cb_margin[1]-cb_size), margin[1], (1-cb_margin[1]), (1-margin[3])]

  plot_len = [plot_pos[2]-plot_pos[0], plot_pos[3] - plot_pos[1]]
  if min(plot_len) le 0 then stop

  plot_aspect = (plot_pos[3] - plot_pos[1]) / (plot_pos[2] - plot_pos[0])

  if log_axes[0] eq 0 then begin
    if keyword_set(force_kperp_axis_range) then begin
      kperp_length = max(force_kperp_axis_range) - min(force_kperp_axis_range)
    endif else begin
      kperp_length = max(kperp_edges_use) - min(kperp_edges_use)
    endelse
  endif else begin
    if keyword_set(force_kperp_axis_range) then begin
      kperp_length = max(log_kperp_axis) - min(log_kperp_axis)
    endif else begin

      if min(kperp_edges_use) lt 0 then begin
        wh_zero = where(kperp_edges_use eq 0, n_zero)
        if n_zero ne 0 then stop

        wh_pos = where(kperp_edges_use ge 0, n_pos, complement = wh_neg, ncomplement = n_neg)
        if n_neg gt 0 then neg_leng = max(alog10((-1)*kperp_edges_use[wh_neg])) - min(alog10((-1)*kperp_edges_use[wh_neg]))
        if n_pos gt 0 then pos_leng = max(alog10(kperp_edges_use[wh_pos])) - min(alog10(kperp_edges_use[wh_pos]))

        kperp_length = neg_leng + pos_leng + pos_leng/(n_pos-1)

      endif else begin
        kperp_length = max(kperp_log_edges) - min(kperp_log_edges)
      endelse
    endelse
  endelse

  if log_axes[1] eq 0 then begin
    if keyword_set(force_kpar_axis_range) then begin
      kpar_length = max(force_kpar_axis_range) - min(force_kpar_axis_range)
    endif else begin
      kpar_length = alog10(max(kpar_edges_use)) - alog10(min(kpar_edges_use[where(kpar_edges_use gt 0)]))
    endelse
  endif else begin
    if keyword_set(force_kpar_axis_range) then begin
      kpar_length = max(log_kpar_axis) - min(log_kpar_axis)
    endif else begin
      kpar_length = max(kpar_log_edges) - min(kpar_log_edges)
    endelse
  endelse

  data_aspect = float(kpar_length / kperp_length)

  aspect_ratio =  data_aspect /plot_aspect

  if keyword_set(pub) then begin
    ;; prevent crazy aspect ratios -- impossible to make decent looking plots.
    max_factor = 2
    if aspect_ratio gt max_factor then begin
      aspect_ratio = max_factor
    endif
    if aspect_ratio lt 1/float(max_factor) then begin
      aspect_ratio = 1/float(max_factor)
    endif
  endif

  if aspect_ratio gt 1 then begin
    y_factor = aspect_ratio
    x_factor = 1.
  endif else begin
    y_factor = 1.
    x_factor = 1./aspect_ratio
  endelse

  if n_elements(multi_pos) eq 4 or n_elements(start_multi_params) gt 0 then begin
    if n_elements(start_multi_params) gt 0 then begin
      ;; calculate desired window size and positions for all plots
      ncol = start_multi_params.ncol
      nrow = start_multi_params.nrow

      multi_pos = fltarr(4, ncol*nrow)

      if tag_exist(start_multi_params, 'ordering') eq 0 then ordering = 'row' $
      else ordering = start_multi_params.ordering

      case ordering of
        'col': begin
          ;; col-major values
          col_val = reform(rebin(reform(indgen(ncol), 1, ncol), nrow, ncol), ncol*nrow)
          row_val = reverse(reform(rebin(indgen(nrow), nrow, ncol), ncol*nrow))
        end
        'row': begin
          ;; row-major values
          col_val = reform(rebin(indgen(ncol), ncol, nrow), ncol*nrow)
          row_val = reverse(reform(rebin(reform(indgen(nrow), 1, nrow), ncol, nrow), ncol*nrow))
        end
        else: message, 'unrecognized ordering value in start_multi_params, use "col" or "row" '
      endcase

      multi_pos[0,*] = col_val/double(ncol)
      multi_pos[1,*] = row_val/double(nrow)
      multi_pos[2,*] = (col_val+1)/double(ncol)
      multi_pos[3,*] = (row_val+1)/double(nrow)

      ;; define window size based on aspect ratio
      base_size_use = base_size
      xsize = round(base_size * x_factor * double(ncol))
      ysize = round(base_size * y_factor * double(nrow))
      if not keyword_set(pub) then begin
        while (ysize gt max_ysize) or (xsize gt max_xsize) do begin
          if base_size_use gt 100 then begin
            base_size_use = base_size_use - 100
          endif else begin
            base_size_use = base_size_use * .75
          endelse
          xsize = round(base_size_use * x_factor * double(ncol))
          ysize = round(base_size_use * y_factor * double(nrow))
        endwhile
      endif

      ;; if pub is set, start ps output
      if keyword_set(pub) then begin
        ps_aspect = (y_factor * float(nrow)) / (x_factor * float(ncol))

        if ps_aspect lt 1 then landscape = 1 else landscape = 0
        IF Keyword_Set(eps) THEN landscape = 0
        sizes = cgpswindow(LANDSCAPE=landscape, aspectRatio = ps_aspect, /sane_offsets)

        cgps_open, plotfile, /font, encapsulated=eps, /nomatch, inches=sizes.inches, $
          xsize=sizes.xsize, ysize=sizes.ysize, xoffset=sizes.xoffset, $
          yoffset=sizes.yoffset, landscape = landscape

      endif else begin
        ;; make or set window
        if windowavailable(window_num) then begin
          wset, window_num
          if !d.x_size ne xsize or !d.y_size ne ysize then begin
            make_win = 1
          endif else begin
            make_win = 0
          endelse
        endif else make_win = 1
        if make_win eq 1 then window, window_num, xsize = xsize, ysize = ysize
        cgerase, background_color
      endelse

      ;; calculate multi_size & multi x/ylen not calculated earlier
      multi_xlen = (multi_pos[2,0]-multi_pos[0,0])
      multi_ylen = (multi_pos[3,0]-multi_pos[1,0])
      multi_center = [multi_pos[0,0] + multi_xlen/2d, multi_pos[1,0] + multi_ylen/2d]

      multi_size = [!d.x_vsize*multi_xlen, !d.y_vsize*multi_ylen]

      multi_pos_use = multi_pos[*,0]
    endif else multi_pos_use = multi_pos

    base_size_use = mean(round([!d.x_size*multi_xlen/x_factor, !d.y_size*multi_ylen/y_factor]))

    multi_aspect = multi_size[1]/float(multi_size[0])

    new_aspect = aspect_ratio/multi_aspect
    if new_aspect gt 1 then begin
      y_factor = 1.
      x_factor = 1/new_aspect
    endif else begin
      y_factor = new_aspect
      x_factor = 1.
    endelse

    new_xlen = multi_xlen*x_factor
    new_ylen = multi_ylen*y_factor
    new_multi = [multi_center[0] - new_xlen/2d, multi_center[1] - new_ylen*y_factor/2d, $
      multi_center[0] + new_xlen/2d, multi_center[1] + new_ylen*y_factor/2d]

    new_pos = [new_xlen * plot_pos[0] + new_multi[0], new_ylen * plot_pos[1] + new_multi[1], $
      new_xlen * plot_pos[2] + new_multi[0], new_ylen * plot_pos[3] + new_multi[1]]

    new_cb_pos = [new_xlen * cb_pos[0] + new_multi[0], new_ylen * cb_pos[1] + new_multi[1], $
      new_xlen * cb_pos[2] + new_multi[0], new_ylen * cb_pos[3] + new_multi[1]]

    plot_pos = new_pos
    cb_pos = new_cb_pos

    no_erase = 1
  endif else begin
    base_size_use = base_size
    xsize = round(base_size_use * x_factor)
    ysize = round(base_size_use * y_factor)

    if keyword_set(pub) then begin
      ps_aspect = y_factor / x_factor

      if ps_aspect lt 1 then landscape = 1 else landscape = 0
      IF Keyword_Set(eps) THEN landscape = 0
      sizes = cgpswindow(LANDSCAPE=landscape, aspectRatio = ps_aspect, /sane_offsets)

      cgps_open, plotfile, /font, encapsulated=eps, /nomatch, inches=sizes.inches, $
        xsize=sizes.xsize, ysize=sizes.ysize, xoffset=sizes.xoffset, $
        yoffset=sizes.yoffset, landscape = landscape

    endif else begin
      while (ysize gt max_ysize) or (xsize gt max_xsize) do begin
        base_size_use = base_size_use - 100
        xsize = round(base_size_use * x_factor)
        ysize = round(base_size_use * y_factor)
      endwhile


      if windowavailable(window_num) then begin
        wset, window_num
        if !d.x_size ne xsize or !d.y_size ne ysize then begin
          make_win = 1
        endif else begin
          make_win = 0
        endelse
      endif else make_win = 1
      if make_win eq 1 then window, window_num, xsize = xsize, ysize = ysize
      cgerase, background_color
    endelse

    no_erase = 0
  endelse

  if not keyword_set(no_title) then begin
    xloc_title = (plot_pos[2] - plot_pos[0])/2. + plot_pos[0]
    if n_elements(multi_pos) gt 0 then begin
      yloc_title = plot_pos[3] + 0.6* (multi_pos_use[3]-plot_pos[3])
    endif else begin
      yloc_title = plot_pos[3] + 0.6* (1-plot_pos[3])
    endelse
  endif

  if n_elements(multi_pos) gt 0 then begin
    xloc_lambda = plot_pos[0] - 0.15* (plot_pos[0]-multi_pos_use[0])
    yloc_lambda = plot_pos[3] + 0.15* (multi_pos_use[3]-plot_pos[3])


    xloc_delay = plot_pos[2] + 0.20 * (multi_pos_use[2]-plot_pos[2])
    yloc_delay = plot_pos[1] - 0.20 * (plot_pos[1]-multi_pos_use[1])

    xloc_note = .99*multi_pos_use[2]
    yloc_note = multi_pos_use[1] + 0.1* (plot_pos[1]-multi_pos_use[1])
  endif else begin
    xloc_lambda = plot_pos[0] - 0.1* (plot_pos[0]-0)
    yloc_lambda = plot_pos[3] + 0.1* (1-plot_pos[3])

    xloc_delay = plot_pos[2] + 0.20*(1-plot_pos[2])
    yloc_delay = plot_pos[1] - 0.2 * (plot_pos[1]-0)

    xloc_note = .99
    yloc_note = 0 + 0.1* (plot_pos[1]-0)
  endelse

  if keyword_set(pub) then begin
    charthick = 3
    thick = 3
    xthick = 3
    ythick = 3
    if n_elements(charsize_in) eq 0 then begin
      if n_elements(multi_pos) gt 0 then begin
        charsize = 1.2d * (mean(multi_size)/10000.)
      endif else charsize = 2
    endif else charsize = charsize_in

    font = 1

    DEVICE, /ISOLATIN1
    perp_char = '!9' + String("136B) + '!X' ;"


  endif else begin
    charthick = 1
    thick = 1
    xthick = 1
    ythick = 1
    font = -1
    if n_elements(charsize_in) eq 0 then begin
      if n_elements(multi_pos) gt 0 then begin
        ;charsize = 0.8d * (multi_size[0]/float(base_size_use))
        charsize = base_size_use / 250.
      endif else charsize = 2
    endif else charsize = charsize_in

    perp_char = '!9' + string(120B) + '!X'

  endelse


  if keyword_set(force_kperp_axis_range) then begin
    if log_axes[0] eq 1 then begin
      plot_kperp = 10^log_kperp_axis
    endif else begin
      plot_kperp = force_kperp_axis_range
    endelse
  endif else begin
    if log_axes[0] eq 1 then begin
      plot_kperp = 10^kperp_log_edges
    endif else begin
      plot_kperp = kperp_edges_use
    endelse
  endelse

  if keyword_set(force_kpar_axis_range)then begin
    if log_axes[1] eq 1 then begin
      plot_kpar = 10^log_kpar_axis
      plot_delay = 10^interpol(delay_log_edges, kpar_log_edges, log_kpar_axis)
    endif else begin
      plot_kpar = force_kpar_axis_range
      plot_delay = interpol(delay_edges, kpar_edges_use, force_kpar_axis_range)
    endelse
  endif else begin
    if log_axes[1] eq 1 then begin
      plot_kpar = 10^kpar_log_edges
      plot_delay = 10^delay_log_edges
    endif else begin
      plot_kpar = kpar_edges_use
      plot_delay = delay_edges
    endelse
  endelse
  cable_index_ref = 0.81
  ;; delay is in ns, factor of 2 to account for reflection bounce
  plot_cable_length = plot_delay * cable_index_ref * 0.3/2.

  ;; if plot title includes sigma need to replace 'sigma' with appropriate character
  ;; (textoidl has to be called after cgps_open)
  if keyword_set(plot_sigma) then begin
    plot_title = repstr(plot_title, 'sigma', textoidl('\sigma', font=font))
  endif

  if keyword_set(title_prefix) then begin
    plot_title = title_prefix + ' ' + plot_title
  endif
  if n_elements(full_title) ne 0 then plot_title = full_title
  if keyword_set(no_title) then undefine, plot_title

  if keyword_set (hinv) then xtitle = textoidl('k_{perp} (!8h!X Mpc^{-1})', font = font) $
  else xtitle = textoidl('k_{perp} (Mpc^{-1})', font = font)
  xtitle = repstr(xtitle, 'perp', perp_char)
  if keyword_set (hinv) then ytitle = textoidl('k_{||} (!8h!X Mpc^{-1})', font = font) $
  else ytitle = textoidl('k_{||} (Mpc^{-1})', font = font)

  if keyword_set(no_title) or keyword_set(baseline_axis) then begin
    initial_title = ''
  endif else begin
    initial_title = plot_title
  endelse

  if log_axes[0] eq 1 then xtickformat = 'exponent' else begin
    nticks = 4

    log_size = round(min(alog10(plot_kperp)))
    xtick_width = round((max(plot_kperp) - min(plot_kperp))/(nticks-1.)/(10.^log_size))*(10.^log_size)
    xtick_start = round(min(plot_kperp)/xtick_width)*(10.^log_size)
    xticks_in = round((dindgen(nticks)*xtick_width + xtick_start)/(10.^log_size))*(10.^log_size)

    x_nticks = nticks-1
    n_minor = 4

    if keyword_set(baseline_axis) then begin
      if keyword_set(hinv) then begin
        baseline_range = minmax(plot_kperp * hubble_param * kperp_lambda_conv)
      endif else begin
        baseline_range = minmax(plot_kperp* kperp_lambda_conv)
      endelse

      log_size2 = round(min(alog10(baseline_range)))
      xtick_width2 = round((max(baseline_range) - min(baseline_range))/(nticks-1.)/(10.^log_size2))*(10.^log_size2)
      xtick_start2 = round(min(baseline_range)/xtick_width2)*(10.^log_size2)
      xticks_in2 = round((dindgen(nticks)*xtick_width2 + xtick_start2)/(10.^log_size2))*(10.^log_size2)
    endif

  endelse
  if log_axes[1] eq 1 then ytickformat = 'exponent'

  if keyword_set(no_plot) then return

  axkeywords = {xlog: log_axes[0], ylog: log_axes[1], xstyle: 5, ystyle: 5, thick: thick, charthick: charthick, xthick: xthick, $
    ythick: ythick, charsize: charsize, font: font}

  if keyword_set(force_kperp_axis_range) or keyword_set(force_kpar_axis_range) then begin
    plot_pos_use = plot_pos

    plot_pos_use[0] = plot_pos[0] + image_relative_pos[0] * (plot_pos[2] - plot_pos[0])
    plot_pos_use[2] = plot_pos[0] + image_relative_pos[2] * (plot_pos[2] - plot_pos[0])
    plot_pos_use[1] = plot_pos[1] + image_relative_pos[1] * (plot_pos[3] - plot_pos[1])
    plot_pos_use[3] = plot_pos[1] + image_relative_pos[3] * (plot_pos[3] - plot_pos[1])

    if log_axes[0] eq 1 then begin
      data_plot_kperp = 10^log_kperp_axis
    endif else begin
      data_plot_kperp = force_kpar_axis_range
    endelse

    cgimage, power_log_norm, /nointerp, xrange = minmax(data_plot_kperp), $
      yrange = minmax(plot_kpar), title=initial_title, position = plot_pos_use, $
      noerase = no_erase, color = annotate_color, background = background_color, $
      axkeywords = axkeywords, /axes

    cgplot, fltarr(10), /nodata, xrange = minmax(plot_kperp), yrange = minmax(plot_kpar), $
      /noerase, title=initial_title, position = plot_pos, color = annotate_color, $
      background = background_color, xlog = log_axes[0], ylog = [1], xstyle = 5, $
      ystyle = 5, thick = thick, charthick = charthick, xthick = xthick, $
      ythick = ythick, charsize = charsize, font = font

  endif else begin
    cgimage, power_log_norm, /nointerp, xrange = minmax(plot_kperp), $
      yrange = minmax(plot_kpar), title=initial_title, position = plot_pos, $
      noerase = no_erase, color = annotate_color, background = background_color, $
      axkeywords = axkeywords, /axes
  endelse

  if keyword_set(bin_contour) then begin
    if keyword_set(png) then begin
      c_thick = [3,1,1,1,1]
    endif else begin
      c_thick = [2,1,1,1,1]
    endelse
    cgcontour, bin_plot, bin_kperp, bin_kpar, levels = levels, $
      thick = thick, charsize = charsize, font = font, label=0, /onimage, c_thick = c_thick
  endif

  if keyword_set(plot_wedge_line) then begin
    n_lines = n_elements(wedge_amp)
    sorted_amp = reverse(wedge_amp[sort(wedge_amp)])
    if n_lines gt 1 then begin
      linestyles = [0, 2, 1]
    endif else begin
      linestyles=2
    endelse

    for i=0, n_lines-1 do begin
      cgplot, /overplot, plot_kperp, plot_kperp * sorted_amp[i], color = annotate_color, $
      thick = thick+1, psym=-0, linestyle = linestyles[i mod n_elements(linestyles)]
    endfor
  endif

  cgaxis, xaxis=0, xtick_get = xticks, xtickv = xticks_in, xticks = x_nticks, $
    xminor=n_minor, xrange = minmax(plot_kperp), xtitle = xtitle, charthick = charthick, $
    xthick = xthick, ythick = ythick, charsize = charsize, font = font, $
    xtickformat = xtickformat, xstyle = 1, color = annotate_color

  cgaxis, yaxis=0, ytick_get = yticks, ytitle = ytitle, yrange = minmax(plot_kpar), $
    charthick = charthick, xthick = xthick, ythick = ythick, charsize = charsize, font = font, $
    ytickformat = ytickformat, ystyle = 1, color = annotate_color
  if keyword_set(baseline_axis) then begin
    ;; baselines don't care about hinv -- take it back out.
    if keyword_set(hinv) then begin
      baseline_range = minmax(plot_kperp * hubble_param * kperp_lambda_conv)
    endif else begin
      baseline_range = minmax(plot_kperp* kperp_lambda_conv)
    endelse

    if keyword_set(no_title) then begin
      xtitle = textoidl('(\lambda)', font = font)
    endif else begin
      undefine, xtitle
    endelse

    cgaxis, xaxis=1, xtickv = xticks_in2, xticks = x_nticks, xminor=n_minor, $
      xrange = baseline_range, xtickformat = xtickformat, xthick = xthick, xtitle = xtitle, $
      charthick = charthick, ythick = ythick, charsize = charsize, font = font, $
      xstyle = 1, color = annotate_color

    if not keyword_set(no_title) then begin
      cgtext, xloc_title, yloc_title, plot_title, /normal, alignment=0.5, $
        charsize=1.2 * charsize, color = annotate_color, font = font
      cgtext, xloc_lambda, yloc_lambda, textoidl('(\lambda)', font = font), $
        /normal, alignment=0.5, charsize=charsize, color = annotate_color, font = font
    endif
  endif else $
    cgaxis, xaxis=1, xrange = minmax(plot_kperp), xtickv = xticks, $
      xtickname = replicate(' ', n_elements(xticks)), $
      charthick = charthick, xthick = xthick, ythick = ythick, charsize = charsize, $
      font = font, xstyle = 1, color = annotate_color
  if keyword_set(delay_axis) or keyword_set(cable_length_axis) then begin

    if keyword_set(delay_axis) then begin
      yrange_use = minmax(plot_delay)
      units_text = '(ns)'
    endif else begin
      yrange_use = minmax(plot_cable_length)
      units_text = '(cbl m)'
    endelse

    cgaxis, yaxis=1, yrange = yrange_use, ytickformat = ytickformat, charthick = charthick, $
      xthick = xthick, ythick = ythick, charsize = charsize, font = font, $
      ystyle = 1, color = annotate_color

    cgtext, xloc_delay, yloc_delay, units_text, /normal, alignment=0.5, charsize=charsize*0.9, $
      color = annotate_color, font = font

  endif else begin
    cgaxis, yaxis=1, yrange = minmax(plot_kpar), ytickv = yticks, $
      ytickname = replicate(' ', n_elements(yticks)), charthick = charthick, $
      xthick = xthick, ythick = ythick, charsize = charsize, font = font, ystyle = 1, $
      color = annotate_color
  endelse

  if n_elements(note) ne 0 then begin
    if keyword_set(pub) then begin
      char_factor = 0.75
    endif else begin
      char_factor = 1
    endelse
    cgtext, xloc_note, yloc_note, note, /normal, alignment=1, $
      charsize = char_factor*charsize, color = annotate_color, font = font
  endif

  cgcolorbar, color = annotate_color, /vertical, position = cb_pos, $
    bottom = color_range[0], ncolors = n_colors, minor = 0, ticknames = cb_ticknames, $
    ytickv = cb_ticks, yticks = n_elements(cb_ticks) -1, title = units_str, $
    charsize = charsize, font = font, oob_low = oob_low

  if n_elements(oob_low) gt 0 and not keyword_set(label_lt_0_oncb) then begin
    cgtext, cb_pos[2], cb_pos[1]-(cb_pos[3]-cb_pos[1])*.05, '<0', /normal, $
    alignment=1, charsize = charsize, color = annotate_color, font = font
  endif

  if keyword_set(pub) and n_elements(multi_pos) eq 0 then begin
    cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density=600
  endif

  tvlct, r, g, b

end
