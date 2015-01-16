pro hellebore_wrapper, folder_name, obs_range, rts = rts, casa = casa, version = version, refresh_dft = refresh_dft, refresh_ps = refresh_ps, $
    refresh_binning = refresh_binning, refresh_info = refresh_info, refresh_beam = refresh_beam, dft_ian = dft_ian, $
    pol_inc = pol_inc, sim = sim, freq_ch_range = freq_ch_range, freq_flag_name = freq_flag_name, $
    no_spec_window = no_spec_window, spec_window_type = spec_window_type, std_power = std_power, no_wtd_avg = no_wtd_avg, $
    cut_image = cut_image, individual_plots = individual_plots, plot_filebase = plot_filebase, png = png, eps = eps, pdf = pdf, $
    plot_slices = plot_slices, slice_type = slice_type, uvf_plot_type = uvf_plot_type, plot_stdset = plot_stdset, $
    plot_kpar_power = plot_kpar_power, plot_kperp_power = plot_kperp_power, plot_k0_power = plot_k0_power, $
    kperp_range_1dave = kperp_range_1dave, kpar_range_1dave = kpar_range_1dave, uv_avg = uv_avg, uv_img_clip = uv_img_clip,$
    kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, $
    t32 = t32, set_data_ranges = set_data_ranges, plot_ranges = plot_ranges, slice_range = slice_range, cube_power_info = cube_power_info
    
  if n_elements(folder_name) eq 0 then message, 'folder name is required'
  if n_elements(folder_name) gt 1 then message, 'Only one folder_name can be supplied'
  
  
  if n_elements(obs_range) gt 0 then begin
    if size(obs_range,/type) eq 7 then begin
      if n_elements(obs_range) gt 1 then $
        message, 'obs_range can be specified as a single string to use as the name or as a 2 element obsid range'
      obs_name = obs_range
    endif else begin
      if n_elements(obs_range) gt 2 then message, 'obs_range can be specified as a single string to use as the name or as a 2 element obsid range'
      if n_elements(obs_range) eq 2 then obs_name = number_formatter(obs_range[0]) + '-' + number_formatter(obs_range[1]) else begin
        obs_name = number_formatter(obs_range[0])
      endelse
    endelse
  endif
  
  obs_info = hellebore_filenames(folder_name, obs_name, sim = sim, rts = rts, casa = casa, refresh_info = refresh_info)
  if tag_exist(obs_info, 'uvf_input') then uvf_input=obs_info.uvf_input
  
  if n_elements(plot_ranges) gt 0 then begin
    if tag_exist(plot_ranges, 'sigma') then sigma_range = plot_ranges.sigma
    if tag_exist(plot_ranges, 'noise_exp') then nev_range = plot_ranges.noise_exp
    if tag_exist(plot_ranges, 'power') then data_range = plot_ranges.power
    if tag_exist(plot_ranges, 'noise_ratio') then nnr_range = plot_ranges.noise_ratio
    if tag_exist(plot_ranges, 'snr') then snr_range = plot_ranges.snr
    if tag_exist(plot_ranges, 'noise') then noise_range = plot_ranges.noise
    
    if tag_exist(plot_ranges, 'power_1d') then range_1d = plot_ranges.power_1d
    
    if tag_exist(plot_ranges, 'kperp') then kperp_plot_range = plot_ranges.kperp
    if tag_exist(plot_ranges, 'kpar') then kpar_plot_range = plot_ranges.kpar
  endif
  
  if keyword_set(rts) then begin
    if keyword_set(casa) then message, 'rts and casa keywords cannot both be set'
    
    if obs_info.info_files[0] ne '' then datafile = obs_info.info_files[0] else $
      if obs_info.cube_files[0] ne '' then datafile = obs_info.cube_files[0] else $
      datafile = rts_fits2idlcube(obs_info.datafiles.(0), obs_info.weightfiles.(0), obs_info.variancefiles.(0), $
      pol_inc, save_path = obs_info.folder_names[0]+path_sep(), refresh = refresh_dft)
      
    plot_filebase = obs_info.rts_types[0] + '_' + obs_info.obs_names[0]
    
    if n_elements(set_data_ranges) eq 0 then set_data_ranges = 1
    if keyword_set(set_data_ranges) then begin
      if n_elements(sigma_range) eq 0 then sigma_range = [1e22, 1e26]
      if n_elements(nev_range) eq 0 then nev_range = [5e22, 5e26]
      
      if n_elements(data_range) eq 0 then data_range = [1e14, 1e20]
      if n_elements(nnr_range) eq 0 then nnr_range = [5e-9,5e-7]
      if n_elements(snr_range) eq 0 then snr_range = [1e-10, 1e-3]
      
      if n_elements(noise_range) eq 0 then noise_range = [2e15, 1e18]
    endif
    
  endif else if keyword_set(casa) then begin
  
    if obs_info.info_files[0] ne '' then datafile = obs_info.info_files[0] else $
      datafile = casa_fits2idlcube(obs_info.cube_files.(0), obs_info.weightfiles.(0), obs_info.variancefiles.(0), $
      pol_inc, save_path = obs_info.folder_names[0]+path_sep(), refresh = refresh_dft)
      
    plot_filebase = obs_info.casa_types[0] ;+ '_' + obs_info.obs_names[0]
    
    if n_elements(set_data_ranges) eq 0 then set_data_ranges = 0
    if keyword_set(set_data_ranges) then begin
      if n_elements(sigma_range) eq 0 then sigma_range = [1e15, 1e17]
      if n_elements(nev_range) eq 0 then nev_range = [1e16, 1e18]
      
      if n_elements(data_range) eq 0 then data_range = [1e6, 2e13]
      if n_elements(nnr_range) eq 0 then nnr_range = [1e-8, 1e-6]
      if n_elements(snr_range) eq 0 then snr_range = [1e-10, 1e-2]
      
      if n_elements(noise_range) eq 0 then noise_range = [1e9, 1e12]
    endif
    
  endif else begin
  
    if obs_info.info_files[0] ne '' then datafile = obs_info.info_files[0] else datafile = obs_info.cube_files.(0)
    
    if tag_exist(obs_info, 'beam_files') then beamfiles = obs_info.beam_files
    
    plot_filebase = obs_info.folder_basenames[0] + '_' + obs_info.obs_names[0]
    note = obs_info.fhd_types[0]
    
    if keyword_set(sim) then begin
      plot_eor_1d=1
      if n_elements(range_1d) eq 0 then range_1d = [1e5, 1e7]
    endif
    
    if n_elements(set_data_ranges) eq 0 and not keyword_set(sim) then set_data_ranges = 1
    
    if keyword_set(set_data_ranges) then begin
      if keyword_set(sim) then begin
        if n_elements(sigma_range) eq 0 then sigma_range = [7e9, 6e14]
        if n_elements(nev_range) eq 0 then nev_range = [2e10, 4e14]
        
        if n_elements(data_range) eq 0 then data_range = [2e2, 2e6]
        if n_elements(nnr_range) eq 0 then nnr_range = [1e-14, 8e-10]
        if n_elements(snr_range) eq 0 then snr_range = [5e-10, 1e-5]
        
        if n_elements(noise_range) eq 0 then noise_range = nev_range
      endif else begin
        if keyword_set(obs_info.integrated[0]) then begin
          if n_elements(sigma_range) eq 0 then sigma_range = [2e5, 2e9]
          if n_elements(nev_range) eq 0 then nev_range = [2e6, 2e10]
        endif else begin
          if n_elements(sigma_range) eq 0 then sigma_range = [1e4, 2e6]
          if n_elements(nev_range) eq 0 then nev_range = [5e4, 2e7]
        endelse
        
        if n_elements(data_range) eq 0 then data_range = [1e3, 1e15]
        if n_elements(nnr_range) eq 0 then nnr_range = [1e-1, 1e1]
        if n_elements(snr_range) eq 0 then snr_range = [1e-5, 1e7]
        
        if n_elements(noise_range) eq 0 then noise_range = nev_range
      endelse
    endif
    
    if n_elements(freq_flag_name) ne 0 then begin
      case freq_flag_name of
        'low': freq_flags = [0,7,8,15,16,23,24,31,32,39,40,47,48,55,56,63,64,71,72,79,80,87,88,95,$
          96,103,104,111,112,119,120,127,128,135,136,143,144,151,152,159,160,167,168,175,176,183,184,191]
          
        'edge': freq_flags = [0,1,6,7,8,9,14,15,16,17,22,23,24,25,30,31,32,33,38,39,40,41,46,47,48,$
          49,54,55,56,57,62,63,64,65,70,71,72,73,78,79,80,81,86,87,88,89,94,95,96,97,102,103,104,105,$
          110,111,112,113,118,119,120,121,126,127,128,129,134,135,136,137,142,143,144,145,150,151,152,$
          153,158,159,160,161,166,167,168,169,174,175,176,177,182,183,184,185,190,191]
          
        else: message, 'freq_flag_name not recognized'
      endcase
    endif
    
    
  endelse
  
  ;; dft_fchunk applies only to Healpix datasets (it's ignored otherwise) and it specifies how many frequencies to process
  ;;   through the dft at once. This keyword allows for trade-offs between memory use and speed.
  ;; The optimum setting varies by computer and the speed is NOT a linear function of this parameter
  ;;   (it's not even monotonic) so some experimenting is required. The memory required is approximately linear --
  ;;   the higher this value the more memory required.
  ;; The maximum value of this parameter is the number of frequency slices in the cube
  ;;   (if it's set too large it will be reduced to the maximum)
  
  dft_fchunk = 1
  
  ;; save_path specifies a location to save the power spectrum files.
  ;; This is also where the code looks for intermediate save files to avoid re-running code.
  ;; If this is parameter is not set, the files will be saved in the same directory as the datafile.
  
  ;; plot_path specifies a location to save plot files.
  ;; If this parameter is not set, the plots will be saved in the same directory as the datafile.
  
  ;  if keyword_set(rts) then std_savepath = base_path('data') + 'rts_data/' $
  ;  else if keyword_set(sim) then std_savepath = base_path('data') + 'fhd_sim_data/' else std_savepath = base_path('data') + 'fhd_ps_data/'
  ;
  ;  if n_elements(save_path) gt 0 then begin
  ;    pos = strpos(save_path, std_savepath)
  ;    if pos ne -1 then save_path_ext = strmid(save_path, pos + strlen(std_savepath)) else save_path_ext = ''
  ;    if keyword_set(rts) then plot_path = base_path('plots') + 'power_spectrum/rts_data/' + save_path_ext else $
  ;      if keyword_set(sim) then plot_path = base_path('plots') + 'power_spectrum/fhd_sim/' + save_path_ext else $
  ;      plot_path = base_path('plots') + 'power_spectrum/fhd_data/' + save_path_ext
  ;  endif else begin
  save_path = obs_info.save_paths[0]
  plot_path = obs_info.plot_paths[0]
  ;  endelse
  
  
  ;; the following sets the save_path to a 'psv'+version directory inside the datafile[0] directory and
  ;; creates the directory if it doesn't exist
  if n_elements(version) gt 0 then begin
    save_path = save_path + 'psv' + number_formatter(version) + path_sep()
    plot_path = plot_path + 'psv' + number_formatter(version) + path_sep()
  endif
  
  if not file_test(save_path, /directory) then file_mkdir, save_path
  if not file_test(plot_path, /directory) then file_mkdir, plot_path
  
  
  ;; savefilebase specifies a base name to use for the save files
  
  ;; plot_filebase specifies a base name to use for the plot files
  
  ;; freq_ch_range specifies which frequency channels to include in the power spectrum.
  ;; Fewer number of channels makes the dfts faster
  
  ;; pol_inc specifies which polarizations to generate the power spectra for.
  
  ;; cut_image keyword only applies to Healpix datasets. It allows for limiting the field of view in the
  ;; image plane to match calculated k-modes (centered on image center).
  ;; Currently defaults to on. Set equal to 0 to turn it off, 1 to turn it on
  
  ;; There are 3 refresh flags to indicate that various stages should be re-calculated
  ;;   (rather than using previous save files if they exist).
  ;; If an early stage is recalculated, all subsequent stages will also be recalculated
  ;; The earliest stage is refresh_dft, which is only used for Healpix datasets (it's ignored otherwise)
  ;; The next stage is refresh_ps and the last stage is refresh_binning.
  ;; To set any of these flags, set them equal to 1 (true)
  
  ;; options for spectral windowing:
  ;; available window funtions are: ['Hann', 'Hamming', 'Blackman', 'Nutall', 'Blackman-Nutall', 'Blackman-Harris']
  ;; Default is to use Blackman-Harris, for no spectral windowing set no_spec_window = 1
  ;; To use another window type use the spec_window_type keyword, eg spec_window_type = 'hann'
  
  ;; options for binning:
  ;; log_kperp, log_kpar and log_k1d are flags: set to 1 (true) for logarithmic bins
  ;; kperp_bin, kpar_bin and k1d_bin take scalar values to control bin sizes.
  ;;   (The actual binsize for linear binning and the log binsize for log binning -- bins per decade = 1/binsize)
  
  
  ;; options for plotting:
  ;; kperp_linear_axis is a flag, set to 1 to use a linear kperp axis (default is log axis)
  ;; kpar_linear_axis is a flag, set to 1 to use a linear kpar axis (default is log axis)
  ;; data_range specifies the min & max value of the signal colorbar (values outside that range are clipped to those values)
  ;; sigma_range, nev_range, snr_range, noise_range, nnr_range control the other colorbar ranges
  ;; baseline_axis is a flag (defaulted to true) to mark baseline; length along top axis of 2d plots (set to 0 to turn off)
  ;; delay_axis is a flag (defaulted to true) to mark delay time along right axis of 2d plots (set to 0 to turn off)
  ;; hinv is a flag (defaulted to true) to use h^-1 Mpc rather than physical Mpc in plot units (set to 0 to turn off)
  ;; plot_wedge_line is a flag (defaulted to true) to plot a line marking the wedge (both horizon & FoV) (set to 0 to turn off)
  ;; png & eps are flags to make save plots as png or eps files rather than displaying to the screen
  
  
  fhd_data_plots, datafile, beamfiles = beamfiles, plot_path = plot_path, plot_filebase = plot_filebase, save_path = save_path, savefilebase = savefilebase, $
    pol_inc = pol_inc, rts = rts, casa = casa, dft_fchunk=dft_fchunk, $
    refresh_dft = refresh_dft, refresh_ps = refresh_ps, refresh_binning = refresh_binning, refresh_info = refresh_info, refresh_beam = refresh_beam, $
    freq_ch_range = freq_ch_range, freq_flags = freq_flags, freq_flag_name = freq_flag_name, $
    no_spec_window = no_spec_window, spec_window_type = spec_window_type, std_power = std_power, no_wtd_avg = no_wtd_avg, $
    sim = sim, cut_image = cut_image, dft_ian = dft_ian, uvf_input = uvf_input, uv_avg = uv_avg, uv_img_clip = uv_img_clip, $
    log_kpar = log_kpar, log_kperp = log_kperp, kpar_bin = kpar_bin, kperp_bin = kperp_bin, $
    log_k1d = log_k1d, k1d_bin = k1d_bin, kperp_range_1dave = kperp_range_1dave, kpar_range_1dave = kpar_range_1dave,$
    kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
    plot_slices = plot_slices, slice_type = slice_type, uvf_plot_type = uvf_plot_type, plot_stdset = plot_stdset, $
    plot_kpar_power = plot_kpar_power, plot_kperp_power = plot_kperp_power, plot_k0_power = plot_k0_power, $
    data_range = data_range, sigma_range = sigma_range, nev_range = nev_range, snr_range = snr_range, noise_range = noise_range, nnr_range = nnr_range, $
    range_1d = range_1d, slice_range = slice_range, $
    baseline_axis = baseline_axis, delay_axis = delay_axis, hinv = hinv, $
    plot_wedge_line = plot_wedge_line, plot_eor_1d = plot_eor_1d, $
    individual_plots = individual_plots, note = note, png = png, eps = eps, pdf = pdf, cube_power_info = cube_power_info
    
    
  if not keyword_set(set_data_ranges) and keyword_set(plot_stdset) then begin
    print, 'data_range used: ', number_formatter(data_range, format = '(e7.1)')
    print, 'sigma_range used: ', number_formatter(sigma_range, format = '(e7.1)')
    print, 'nev_range used: ', number_formatter(nev_range, format = '(e7.1)')
    if n_elements(nnr_range) gt 0 then print, 'nnr_range used: ', number_formatter(nnr_range, format = '(e7.1)')
    print, 'snr_range used: ', number_formatter(snr_range, format = '(e7.1)')
    if n_elements(noise_range) gt 0 then print, 'noise_range used: ', number_formatter(noise_range, format = '(e7.1)')
  endif
  
  if keyword_set(plot_stdset) then begin
    if n_elements(plot_ranges) eq 0 then begin
      plot_ranges = create_struct('sigma', sigma_range, 'noise_exp', nev_range, 'power', data_range, 'snr', snr_range, $
        'kperp', kperp_plot_range, 'kpar', kpar_plot_range)
      if n_elements(nnr_range) gt 0 then plot_ranges = create_struct(plot_ranges, 'noise_ratio', nnr_range)
      if n_elements(noise_range) gt 0 then plot_ranges = create_struct(plot_ranges, 'noise', noise_range)
      if n_elements(range_1d) gt 0 then plot_ranges = create_struct(plot_ranges, 'power_1d', range_1d)
    endif else begin
      if not tag_exist(plot_ranges, 'sigma') then plot_ranges = create_struct(plot_ranges, 'sigma', sigma_range) else plot_ranges.sigma = sigma_range
      if not tag_exist(plot_ranges, 'noise_exp') then plot_ranges = create_struct(plot_ranges, 'noise_exp', nev_range) else plot_ranges.noise_exp = nev_range
      if not tag_exist(plot_ranges, 'power') then plot_ranges = create_struct(plot_ranges, 'power', data_range) else plot_ranges.power = data_range
      if not tag_exist(plot_ranges, 'snr') then plot_ranges = create_struct(plot_ranges, 'snr', snr_range) else plot_ranges.snr = snr_range
      if not tag_exist(plot_ranges, 'power_1d') then plot_ranges = create_struct(plot_ranges, 'power_1d', range_1d) else plot_ranges.power_1d = range_1d
      if n_elements(nnr_range) gt 0 then if not tag_exist(plot_ranges, 'noise_ratio') then $
        plot_ranges = create_struct(plot_ranges, 'noise_ratio', nnr_range) else plot_ranges.noise_ratio = nnr_range
      if n_elements(noise_range) gt 0 then if not tag_exist(plot_ranges, 'noise') then $
        plot_ranges = create_struct(plot_ranges, 'noise', noise_range) else plot_ranges.noise = noise_range
        
      if not tag_exist(plot_ranges, 'kperp') then plot_ranges = create_struct(plot_ranges, 'kperp', kperp_plot_range) else plot_ranges.kperp = kperp_plot_range
      if not tag_exist(plot_ranges, 'kpar') then plot_ranges = create_struct(plot_ranges, 'kpar', kpar_plot_range) else plot_ranges.kpar = kpar_plot_range
    endelse
  endif
  
end
