pro hellebore_wrapper, folder_name, rts = rts, version = version, refresh_dft = refresh_dft, refresh_ps = refresh_ps, $
    refresh_binning = refresh_binning, refresh_info = refresh_info, dft_ian = dft_ian, $
    pol_inc = pol_inc, sim = sim, freq_ch_range = freq_ch_range, freq_flag_name = freq_flag_name, $
    no_spec_window = no_spec_window, spec_window_type = spec_window_type, std_power = std_power, noise_sim = noise_sim, $
    cut_image = cut_image, individual_plots = individual_plots, plot_filebase = plot_filebase, png = png, eps = eps, $
    plot_slices = plot_slices, slice_type = slice_type, $
    kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, t32 = t32, set_data_ranges = set_data_ranges
    
  if keyword_set(rts) then begin
    ;    froot = base_path('data') + 'rts_data/test2/'
    ;
    ;    data_dir = froot + 'BdaggerV/'
    ;    datafiles = file_search(data_dir + '*.fits')
    ;
    ;    weights_dir = froot + 'Bdagger1/'
    ;    weightfiles = file_search(weights_dir + '*.fits')
    ;
    ;    variance_dir = froot + 'BdaggerB/'
    ;    variancefiles = file_search(variance_dir + '*.fits')
    ;
    ;    datafile =  rts_fits2idlcube(datafiles, weightfiles, variancefiles, pol_inc, save_path = froot)
  
    datafile = file_search(base_path('data') + 'rts_data/wellington_data2/*idlcube.idlsave')
    
  endif else if keyword_set(sim) then begin
    datafile = base_path('data') + 'fhd_sim_data/fhd_v300/Healpix/Sim_obs_' + ['even','odd']+ '_cube.sav'
    
  endif else begin
  
    if n_elements(folder_name) eq 0 then folder_name = base_path('data') + 'fhd_ps_data/128T_cubes/aug23_3hr_first/'
    
    ;; check for folder existence, otherwise look for common folder names to figure out full path. If none found, try base_path('data') + 'fhd_ps_data/128T_cubes/'
    folder_test = file_test(folder_name, /directory)
    if folder_test eq 0 then begin
      pos_fhd_data = strpos(folder_name, 'fhd_ps_data')
      if pos_fhd_data gt -1 then begin
        test_name = base_path('data') + strmid(folder_name, pos_fhd_data)
        folder_test = file_test(test_name, /directory)
        if folder_test eq 1 then folder_name = test_name
      endif
    endif
    if folder_test eq 0 then begin
      pos_fhd_128 = strpos(folder_name, '128T_cubes')
      if pos_fhd_128 gt -1 then begin
        test_name = base_path('data') + 'fhd_ps_data/' + strmid(folder_name, pos_fhd_128)
        folder_test = file_test(test_name, /directory)
        if folder_test eq 1 then folder_name = test_name
      endif
    endif
    if folder_test eq 0 then begin
      test_name = base_path('data') + 'fhd_ps_data/128T_cubes/' + folder_name
      folder_test = file_test(test_name, /directory)
      if folder_test eq 1 then folder_name = test_name
    endif
    
    if folder_test eq 0 then message, 'folder not found'
    
    fhd_type = file_basename(folder_name)
    
    if n_elements(obs_range) gt 0 then begin
      if size(obs_range,/type) eq 7 then begin
        if n_elements(obs_range) gt 1 then $
          message, 'obs_range can be specified as a single string to use as the name or as a 2 element obsid range'
        obs_name = obs_range
        obs_range = long(strsplit(obs_name, '-', /extract))
        if n_elements(obs_range) eq 1 then obs_name_single = obs_name
      endif else begin
        if n_elements(obs_range) gt 2 then message, 'obs_range can be specified as a single string to use as the name or as a 2 element obsid range'
        if n_elements(obs_range) eq 2 then obs_name = number_formatter(obs_range[0]) + '-' + number_formatter(obs_range[1]) else begin
          obs_name = number_formatter(obs_range[0]) + '-' + number_formatter(obs_range[0])
          obs_name_single = number_formatter(obs_range[0])
        endelse
      endelse
    endif else begin
      obs_name = ''
      obs_name_single = ''
    endelse
    
    ;; first look for integrated info files with names like Combined_obs_...
    info_file = file_search(folder_name + '/Combined_obs_*info*', count = n_infofile)
    if n_infofile gt 0 then begin
      if obs_name eq '' then begin
        if n_infofile gt 1 then print, 'More than 1 info files found, using first one'
        datafile = info_file[0]
        obs_name = stregex(datafile, '[0-9]+-[0-9]+', /extract)
        obs_range = long(strsplit(obs_name, '-', /extract))
      endif else begin
        if n_infofile gt 1 then message, 'More than one info file found with given obs_range'
        datafile = info_file[0]
        obs_name = stregex(datafile, '[0-9]+-[0-9]+', /extract)
        obs_range = long(strsplit(obs_name, '-', /extract))
      endelse
      
    endif else if n_elements(obs_range) lt 2 then begin
      ;; then look for single obs info files
      info_file = file_search(folder_name + '/' + obs_name_single + '*info*', count = n_infofile)
      if n_infofile gt 0 then begin
        info_basename = file_basename(info_file)
        if obs_name eq '' then begin
          if n_infofile gt 1 then print, 'More than 1 info files found, using first one'
          datafile = info_file[0]
          obs_name = stregex(info_basename[0], '[0-9]+', /extract)
          obs_range = long(obs_name)
        endif else begin
          if n_infofile gt 1 then message, 'More than one info file found with given obs_range'
          datafile = info_file[0]
          obs_name = stregex(info_basename[0], '[0-9]+', /extract)
          obs_range = long(obs_name)
        endelse
        
      endif
    endif
    
    if n_infofile eq 0 then begin
      ;; first look for integrated cube files with names like Combined_obs_...
      cube_files = file_search(folder_name + '/Combined_obs_' + obs_name + '*_cube.sav', count = n_cubefiles)
      if n_cubefiles gt 0 then begin
        if obs_name eq '' then begin
          obs_name_arr = stregex(cube_files, '[0-9]+-[0-9]+', /extract)
          wh_first = where(obs_name_arr eq obs_name_arr[0], count_first)
          if count_first lt n_elements(cube_files) then $
            print, 'More than one obs_range found, using first range (' + obs_name_arr[0] + ', ' + number_formatter(count_first) + ' files)'
          if count_first gt 2 then message, 'More than two cubes found with first obs_range'
          datafile = cube_files[wh_first]
          obs_name = obs_name_arr[0]
          obs_range = long(strsplit(obs_name, '-', /extract))
        endif else begin
          if n_elements(cube_files) gt 2 then message, 'More than two cubes found with given obs_range'
          
          obs_name_arr = stregex(cube_files, '[0-9]+-[0-9]+', /extract)
          if obs_name_arr[1] ne obs_name_arr[0] then message, 'Cube files do not have the same obs ranges.'
          obs_name = obs_name_arr[0]
          obs_range = long(strsplit(obs_name, '-', /extract))
        endelse
        
      endif else if n_elements(obs_range) lt 2 then begin
        ;; then look for single obs cube files
        cube_files = file_search(folder_name + '/' + obs_name + '*_cube.sav', count = n_cubefiles)
        if n_cubefiles gt 0 then begin
          cube_basename = file_basename(cube_files)
          if obs_name eq '' then begin
            obs_name_arr = stregex(cube_basename, '[0-9]+', /extract)
            wh_first = where(obs_name_arr eq obs_name_arr[0], count_first)
            if count_first lt n_elements(cube_files) then $
              print, 'More than one obs_range found, using first range (' + obs_name_arr[0] + ', ' + number_formatter(count_first) + ' files)'
            if count_first gt 2 then message, 'More than two cubes found with first obs_range'
            datafile = cube_files[wh_first]
            obs_name = obs_name_arr[0]
            obs_range = long(obs_name)
          endif else begin
            if n_elements(cube_files) gt 2 then message, 'More than two cubes found with given obs_range'
            
            datafile = cube_files
            obs_name_arr = stregex(cube_basename, '[0-9]+', /extract)
            if obs_name_arr[1] ne obs_name_arr[0] then message, 'Cube files do not have the same obs ranges.'
            obs_name = obs_name_arr[0]
            obs_range = long(obs_name)
          endelse
          
        endif
      endif
    endif
    
    if n_elements(datafile) eq 0 then message, 'No cube or info files found in folder ' + folder_name
    
    
    if n_elements(obs_range) eq 1 then integrated = 0 else if obs_range[1] - obs_range[0] gt 0 then integrated = 1 else integrated = 0
    
    plot_filebase = fhd_type + '_' + obs_name
    
    if n_elements(set_data_ranges) eq 0 then set_data_ranges = 1
    if keyword_set(set_data_ranges) then begin
      if keyword_set(integrated) then sigma_range = [2e0, 2e2] else sigma_range = [1e2, 2e4]
      if keyword_set(integrated) then nev_range = [5e0, 2e3] else nev_range = [5e2, 2e5]
      
      data_range = [1e-2, 1e8]
      nnr_range = [1e-1, 1e1]
      snr_range = [1e-4, 1e6]
      
      noise_range = nev_range
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
  
  ;; the following sets the save_path to a 'psv'+version directory inside the datafile[0] directory and
  ;; creates the directory if it doesn't exist
  if n_elements(version) gt 0 then begin
    save_path = file_dirname(datafile[0], /mark_directory) + 'psv' + number_formatter(version) + path_sep()
    if not file_test(save_path, /directory) then file_mkdir, save_path
  endif
  
  if keyword_set(rts) then std_savepath = base_path('data') + 'rts_data/' $
  else if keyword_set(sim) then std_savepath = base_path('data') + 'fhd_sim_data/' else std_savepath = base_path('data') + 'fhd_ps_data/'
  
  if n_elements(save_path) gt 0 then begin
    pos = strpos(save_path, std_savepath)
    if pos ne -1 then save_path_ext = strmid(save_path, pos + strlen(std_savepath)) else save_path_ext = ''
  endif else begin
    pos = strpos(file_dirname(datafile[0], /mark_directory), std_savepath)
    if pos ne -1 then save_path_ext = strmid(file_dirname(datafile[0], /mark_directory), pos + strlen(std_savepath)) $
    else save_path_ext = ''
  endelse
  
  ;; savefilebase specifies a base name to use for the save files
  
  
  ;; plot_path specifies a location to save plot files.
  ;; If this parameter is not set, the plots will be saved in the same directory as the datafile.
  if keyword_set(rts) then plot_path = base_path('plots') + 'power_spectrum/rts_data/' + save_path_ext $
  else $
    if keyword_set(sim) then plot_path = base_path('plots') + 'power_spectrum/fhd_sim/' + save_path_ext else $
    plot_path = base_path('plots') + 'power_spectrum/fhd_data/' + save_path_ext
    
  if not file_test(plot_path, /directory) then file_mkdir, plot_path
  
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
  ;; grey_scale is a flag to use a black/white color scale rather than the default color scale
  ;; png & eps are flags to make save plots as png or eps files rather than displaying to the screen
  
  
  fhd_data_plots, datafile, dft_fchunk=dft_fchunk, plot_path = plot_path, plot_filebase = plot_filebase, save_path = save_path, savefilebase = savefilebase, $
    pol_inc = pol_inc, rts = rts, $
    refresh_dft = refresh_dft, refresh_ps = refresh_ps, refresh_binning = refresh_binning, refresh_info = refresh_info, $
    freq_ch_range = freq_ch_range, freq_flags = freq_flags, freq_flag_name = freq_flag_name, $
    no_spec_window = no_spec_window, spec_window_type = spec_window_type, std_power = std_power, $
    noise_sim = noise_sim, cut_image = cut_image, dft_ian = dft_ian, $
    log_kpar = log_kpar, log_kperp = log_kperp, kpar_bin = kpar_bin, kperp_bin = kperp_bin, $
    log_k1d = log_k1d, k1d_bin = k1d_bin, kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, $
    plot_slices = plot_slices, slice_type = slice_type, $
    data_range = data_range, sigma_range = sigma_range, nev_range = nev_range, snr_range = snr_range, noise_range = noise_range, nnr_range = nnr_range, $
    baseline_axis = baseline_axis, delay_axis = delay_axis, hinv = hinv, $
    plot_wedge_line = plot_wedge_line, grey_scale = grey_scale, individual_plots = individual_plots, png = png, eps = eps
    
    
  if not keyword_set(set_data_ranges) then begin
    print, 'data_range used: ', number_formatter(data_range, format = '(e7.1)')
    print, 'sigma_range used: ', number_formatter(sigma_range, format = '(e7.1)')
    print, 'nev_range used: ', number_formatter(nev_range, format = '(e7.1)')
    print, 'nnr_range used: ', number_formatter(nnr_range, format = '(e7.1)')
    print, 'snr_range used: ', number_formatter(snr_range, format = '(e7.1)')
    print, 'noise_range used: ', number_formatter(noise_range, format = '(e7.1)')
  endif
end
