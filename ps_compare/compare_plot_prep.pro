pro compare_plot_prep, folder_names, obs_info, ps_foldernames = ps_foldernames, $
    cube_types, pols, comp_type, compare_files, $
    uvf_options = uvf_options, ps_options = ps_options, $
    plot_options = plot_options, plot_2d_options = plot_2d_options, $
    plot_slices = plot_slices, slice_type = slice_type, fadd_2dbin = fadd_2dbin, $
    freq_ch_range = freq_ch_range, $
    plot_filebase = plot_filebase, save_path = save_path, savefilebase = savefilebase, $
    axis_type_1d = axis_type_1d, full_compare = full_compare

  comp_type_enum = ['diff', 'diff_ratio', 'ratio']
  wh_comp_type = where(comp_type_enum eq comp_type, count_comp_type)
  if count_comp_type eq 0 then begin
    message, 'comp_type not recognized, must be one of: ' + strjoin(comp_type_enum, ' ,')
  endif
  if comp_type eq 'ratio' and keyword_set(full_compare) then comp_type = 'diff_ratio'

  if n_elements(obs_info.info_files) gt 2 then message, 'Only 1 or 2 info_files can be used'

  if n_elements(ps_options) eq 2 and tag_exist(ps_options[0], 'spec_window_type') then begin
    if (ps_options[0].spec_window_type ne ps_options[1].spec_window_type) then begin

      type_list = ['Hann', 'Hamming', 'Blackman', 'Nutall', 'Blackman-Nutall', $
        'Blackman-Harris', 'Blackman-Harris^2', 'None']
      sw_tag_list = ['hann', 'ham', 'blm', 'ntl', 'bn', 'bh', 'bh2', '']
      sw_tags = strarr(2)
      for i=0, 1 do begin
        wh_type = where(strlowcase(type_list) eq strlowcase(ps_options[i].spec_window_type), $
          count_type)
        if count_type eq 0 then begin
          wh_type = where(strlowcase(sw_tag_list) eq strlowcase(ps_options[i].spec_window_type), $
            count_type)
        endif
        if count_type eq 0 then begin
          message, 'Spectral window type not recognized.'
        endif else begin
          ps_options[i].spec_window_type = type_list[wh_type[0]]
          if ps_options[i].spec_window_type eq 'None' then begin
            sw_tags[i] = '_nosw'
          endif else begin
            sw_tags[i] = '_' + sw_tag_list[wh_type[0]]
          endelse
        endelse
      endfor
    endif
  endif else sw_tags = ''

  if n_elements(uvf_options) eq 2 and tag_exist(uvf_options[0], 'image_window_name') then begin
    if (uvf_options[0].image_window_name ne uvf_options[1].image_window_name) $
      or (uvf_options[0].image_window_frac_size ne uvf_options[1].image_window_frac_size) then begin

      type_list = ['Tukey', 'None']
      iw_tag_list = ['tk', '']

      iw_tag = strarr(2)
      iw_size_tag = strarr(2)
      for i=0,1 do begin
        wh_type = where(strlowcase(type_list) eq strlowcase(uvf_options[i].image_window_name), $
          count_type)
        if count_type eq 0 then begin
          wh_type = where(strlowcase(iw_tag_list) eq strlowcase(uvf_options[i].image_window_name), $
            count_type)
        endif
        if count_type eq 0 then begin
          message, 'Image window type not recognized.'
        endif else begin
          uvf_options[i].image_window_name = type_list[wh_type[0]]
          if type_list[wh_type] ne 'None' and tag_exist(uvf_options[i], 'image_window_frac_size') then begin
           iw_size_tag[i] = number_formatter(uvf_options[i].image_window_frac_size)
          endif else begin
            iw_size_tag[i] = ''
          endelse
          if uvf_options[i].image_window_name eq 'None' then begin
            iw_tag[i] = ''
          endif else begin
            iw_tag[i] = '_' + iw_tag_list[wh_type[0]] + iw_size_tag[i]
          endelse
        endelse
      endfor
    endif
  endif
  if n_elements(iw_tag) eq 0 then begin
    iw_tag = ''
    iw_size_tag = ''
  endif

  ;; density correction defaults & file naming for 2D & 1D files
  if n_elements(ps_options) eq 2 then begin
    if ps_options[0].wt_cutoffs ne ps_options[1].wt_cutoffs and $
        abs(ps_options[1].wt_cutoffs - ps_options[1].wt_cutoffs) le 1e-3 then begin

      ps_options[1].wt_cutoffs = ps_options[1].wt_cutoffs
    endif

    if (ps_options[0].wt_cutoffs ne ps_options[1].wt_cutoffs $
        or ps_options[0].wt_measures ne ps_options[1].wt_measures) then begin
      n_wtcuts = 2

      if comp_type eq 'diff' then begin
        message, 'no more than one wt_cutoff and wt_measure can be set for difference plots'
      endif
    endif else begin
      n_wtcuts = 1
    endelse
  endif else begin
    n_wtcuts = 1
  endelse

  kperp_density_names = strarr(n_wtcuts)
  wh_cutoff0 = where(ps_options.wt_cutoffs eq 0, count_cutoff0, complement = wh_cutoff_n0, $
    ncomplement = count_cutoff_n0)
  wh_std = where(ps_options.wt_cutoffs eq 1 and ps_options.wt_measures eq 'min', count_std)

  if count_cutoff0 gt 0 then kperp_density_names[wh_cutoff0] = '_nodensitycorr'
  if count_cutoff_n0 gt 0 then begin
    kperp_density_names[wh_cutoff_n0] = '_kperp_density_' + ps_options.wt_measures[wh_cutoff_n0] $
      + '_gt' + number_formatter(ps_options.wt_cutoffs[wh_cutoff_n0])
  endif
  if count_std gt 0 then kperp_density_names[wh_std] = '_dencorr'

  if n_elements(freq_ch_range) ne 0 then begin
    if min(freq_ch_range) lt 0 or max(freq_ch_range) - min(freq_ch_range) lt 3 then begin
      message, 'invalid freq_ch_range'
    endif
    fch_tag = '_ch' + number_formatter(min(freq_ch_range)) + '-' + $
      number_formatter(max(freq_ch_range))
  endif else fch_tag = ''

  n_diffs = max([n_elements(obs_info.info_files), n_elements(cube_types), $
    n_elements(pols), n_elements(ps_options), n_elements(uvf_options)])
  if n_diffs gt 2 then message, 'only 1 or 2 each of [folder_names, ps_foldernames, ' + $
    'obs_names, cube_types, pols, spec_window_types, wt_cutoffs, ave_removal, ' + $
    'image_window_name, image_window_frac_size] allowed'

  if (n_elements(obs_info.info_files) eq 2 or n_elements(ps_options) eq 2 $
    or n_elements(uvf_options) eq 2) and n_elements(cube_types) eq 0 $
    and n_elements(pols) eq 0 and n_elements(full_compare) eq 0 then full_compare = 1

  if keyword_set(full_compare) and n_elements(info_files) eq 1 and n_elements(ps_options) eq 1 $
      and n_elements(uvf_options) eq 1 then begin

    message, 'full_compare can only be set if one of [folder names, obs names, ' + $
      'spec windows, wt_cutoffs/measures, ave_removal values, image window name/size] is length 2.'
  endif

  if comp_type eq 'ratio' and keyword_set(full_compare) then comp_type = 'diff_ratio'

  if not keyword_set(full_compare) then begin
    if n_elements(cube_types) eq 0 then if n_diffs eq 1 then cube_types = ['dirty', 'res'] else cube_types = 'res'
    n_diffs = max([n_elements(obs_info.info_files), n_elements(cube_types), n_elements(pols), $
      n_elements(ps_options), n_elements(uvf_options)])

    if n_elements(pols) eq 0 then if n_diffs eq 1 then pols=['xx', 'yy'] else pols = 'xx'
    n_diffs = max([n_elements(obs_info.info_files), n_elements(cube_types), n_elements(pols), $
      n_elements(ps_options), n_elements(uvf_options)])

    if n_diffs eq 1 then message, 'at least one of [folder_names, ps_foldernames, ' + $
      'obs_names, cube_types, pols, spec_window_types, wt_cutoffs, ave_removal, ' + $
      'image_window_name, image_window_frac_size] must be a 2 element vector'

    case comp_type of
      'diff': n_cubes = 1
      'ratio': n_sets=1
      'diff_ratio': n_sets=2
    endcase
  endif else begin
    case comp_type of
      'diff': undefine, cube_types, pols
      'diff_ratio': begin
        pols = ['xx', 'yy']
        if n_elements(cube_types) ne 2 then cube_types = ['res', 'dirty']

        n_sets=4
      end
    endcase
  endelse

  max_file = n_elements(obs_info.info_files)-1
  max_type = n_elements(cube_types)-1
  max_pol = n_elements(pols)-1
  max_ps = n_elements(ps_options) - 1
  max_uvf = n_elements(uvf_options) - 1

  if n_elements(axis_type_1d) eq 0 then axis_type_1d = 'sym_log'

  if n_elements(folder_names) eq 2 or n_elements(ps_foldernames) eq 2 then begin
    if n_elements(save_path) eq 0 then save_path = obs_info.diff_save_path
    plot_options = create_plot_options(plot_options = plot_options, $
      note = obs_info.diff_note)
    if not tag_exist(plot_options, 'plot_path') then if tag_exist(obs_info, 'diff_plot_path') then begin
      plot_options = create_plot_options(plot_options = plot_options, $
        plot_path = obs_info.diff_plot_path)
    endif else begin
      plot_options = create_plot_options(plot_options = plot_options, $
        plot_path = save_path + path_sep() + 'plots' + path_sep())
    endelse
  endif else begin
    if n_elements(save_path) eq 0 then save_path = obs_info.save_paths[0]
    plot_options = create_plot_options(plot_options = plot_options, $
      note = obs_info.fhd_types[0])
    if not tag_exist(plot_options, 'plot_path') then begin
      plot_options = create_plot_options(plot_options = plot_options, $
        plot_path = obs_info.plot_paths[0])
    endif
  endelse

  case comp_type of
  'diff': op_str = ' minus '
  'diff_ratio': op_str = ' minus '
  'ratio': op_str = ' over '
  endcase

  if n_elements(ps_options) eq 2 then begin
    if (ps_options[0].spec_window_type ne ps_options[1].spec_window_type) then begin
      plot_options.note = plot_options.note + ' ' + spec_window_types[0] + $
        op_str + spec_window_types[1]
    endif

    if (ps_options[0].ave_removal ne ps_options[1].ave_removal) then begin
      ave_removal_str = ['no_ave_removal', 'ave_removal']
      plot_options.note = plot_options.note + ' ' + ave_removal_str[ave_removal[0]] + $
        op_str + ave_removal_str[ave_removal[1]]
      ar_tag_list = ['', '_averemoval']
      ar_tags = [ar_tag_list[ave_removal[0]], ar_tag_list[ave_removal[1]]]
    endif
  endif
  if n_elements(ar_tags) eq 0 then ar_tags = ''

  if n_wtcuts eq 2 then begin
    density_tags = kperp_density_names
  endif else begin
    density_tags = ''
  endelse

  if n_elements(iw_tag) gt 1 then begin
    plot_options.note = plot_options.note + ' ' + iw_tag[0] + op_str + iw_tag[1]
  endif

  file_struct_arr1 = fhd_file_setup(obs_info.info_files[0], $
    uvf_options = uvf_options[0], ps_options = ps_options[0], $
    freq_ch_range = freq_ch_range)
  if n_elements(obs_info.info_files) eq 2 or max_ps eq 1 or max_uvf eq 1 then begin
    file_struct_arr2 = fhd_file_setup(obs_info.info_files[max_file], $
    uvf_options = uvf_options[max_uvf], ps_options = ps_options[max_ps],$
    freq_ch_range = freq_ch_range)
  endif else begin
    file_struct_arr2 = file_struct_arr1
  endelse
  type_pol_str1 = file_struct_arr1.type_pol_str
  type_pol_str2 = file_struct_arr2.type_pol_str

  if keyword_set(plot_slices) then begin
    if n_elements(slice_type) eq 0 then slice_type = 'kspace'
    slice_type_enum = ['raw', 'divided', 'kspace', 'sumdiff', 'weights']

    wh_slice_type = where(slice_type_enum eq slice_type, count_slice_type)
    if count_slice_type eq 0 then begin
      print, 'slice_type not recognized, using default'
      slice_type = 'kspace'
    endif

    if slice_type ne 'kspace' then message, 'only kspace slice difference plots are currently supported'

    slice_tags = ['xz', 'yz', 'xy']
    n_slices = n_elements(slice_tags)
  endif else n_slices = 1

  ;; get same parts of power_tag to add to plot file name
  if n_elements(file_struct_arr2) eq 0 then same_power_tag = file_struct_arr1[0].power_tag else begin
    if file_struct_arr1[0].power_tag eq file_struct_arr2[0].power_tag then same_power_tag = file_struct_arr1[0].power_tag else begin
      tag_arr1 = strsplit(file_struct_arr1[0].power_tag, '_',/extract)
      tag_arr2 = strsplit(file_struct_arr2[0].power_tag, '_',/extract)
      for i=0, n_elements(tag_arr1)-1 do begin
        wh_in2 = where(tag_arr2 eq tag_arr1[i], count_in2)
        if count_in2 gt 0 then begin
          if n_elements(same_arr) eq 0 then same_arr = tag_arr1[i] else same_arr = [same_arr, tag_arr1[i]]
        endif
      endfor
      if n_elements(same_arr) gt 0 then same_power_tag = '_' + strjoin(same_arr, '_')
    endelse
  endelse

  if n_wtcuts eq 1 then same_density_tag = kperp_density_names else same_density_tag = ''

  case comp_type of
    'diff': op_str = '_minus_'
    'diff_ratio': op_str = '_diffratio_'
    'ratio': op_str = '_over_'
  endcase

  if keyword_set(full_compare) then begin

    if n_elements(folder_names) eq 2 and $
        folder_names[0] ne folder_names[n_elements(folder_names)-1] then begin

      plot_filebase = obs_info.name_same_parts + same_power_tag + same_density_tag +'__' + $
        obs_info.name_diff_parts[0] + '_' + obs_info.obs_names[0] + $
        iw_tag[0] + sw_tags[0] + ar_tags[0] + density_tags[0] + op_str + $
        obs_info.name_diff_parts[1]  + '_' + obs_info.obs_names[1] + $
        iw_tag[n_elements(iw_tag)-1] + sw_tags[n_elements(sw_tags)-1] + $
        ar_tags[n_elements(ar_tags)-1] + density_tags[n_wtcuts-1]
    endif else begin
      plot_filebase = obs_info.folder_basenames[0] + same_power_tag + same_density_tag + '__' + $
        obs_info.obs_names[0] + iw_tag[0] + sw_tags[0] + ar_tags[0] + density_tags[0] + op_str + $
        obs_info.obs_names[max_file] + iw_tag[n_elements(iw_tag)-1] + sw_tags[n_elements(sw_tags)-1] + $
        ar_tags[n_elements(ar_tags)-1] + density_tags[n_wtcuts-1]
    endelse
  endif else begin

    if n_elements(folder_names) eq 1 then begin
      if n_elements(obs_info.obs_names) gt 1 then begin
        plot_filebase = obs_info.folder_basenames[0] + same_power_tag + same_density_tag + $
          '_' + obs_info.obs_names[0] + '_' + cube_types[0] + '_' + pols[0] + $
          iw_tag[0] + sw_tags[0] + ar_tags[0] + density_tags[0] +  op_str + $
          obs_info.obs_names[0] + '_' + cube_types[max_type] + '_' + pols[max_pol] + $
          iw_tag[n_elements(iw_tag)-1] + sw_tags[n_elements(sw_tags)-1] + $
          ar_tags[n_elements(ar_tags)-1] + density_tags[n_wtcuts-1]
      endif else begin
        if obs_info.integrated[0] eq 0 then plot_start = obs_info.folder_basenames[0] + '_' + obs_info.obs_names[0] else plot_start = obs_info.fhd_types[0]

        plot_filebase = plot_start + same_power_tag + same_density_tag + '_' + $
          cube_types[0] + '_' + pols[0] + $
          iw_tag[0] + sw_tags[0] + ar_tags[0] + density_tags[0] +  op_str + $
          cube_types[max_type] + '_' + pols[max_pol] + $
          iw_tag[n_elements(iw_tag)-1] + sw_tags[n_elements(sw_tags)-1] + $
          ar_tags[n_elements(ar_tags)-1] + density_tags[n_wtcuts-1]
      endelse
    endif else begin
      plot_filebase = obs_info.name_same_parts + same_power_tag + $
        same_density_tag + '__' + strjoin([obs_info.name_diff_parts[0], $
        cube_types[0], pols[0]], '_') + $
        iw_tag[0] + sw_tags[0] + ar_tags[0] + density_tags[0] +  op_str + $
        strjoin([obs_info.name_diff_parts[1], cube_types[max_type], pols[max_pol]], '_') + $
        iw_tag[n_elements(iw_tag)-1] + sw_tags[n_elements(sw_tags)-1] + $
        ar_tags[n_elements(ar_tags)-1] + density_tags[n_wtcuts-1]
    endelse
  endelse
  plot_filebase = plot_filebase + fch_tag

  if not file_test(plot_options.plot_path, /directory) then file_mkdir, plot_options.plot_path
  if not file_test(save_path, /directory) then file_mkdir, save_path

  if keyword_set(full_compare) then begin
    case comp_type of
      'diff': begin
        if n_elements(type_pol_str1) ne n_elements(type_pol_str2) then begin
          message, 'all_type_pol cannot be used with these folders, they contain ' + $
            'different number of types & pols'
        endif
        n_cubes = n_elements(type_pol_str1)

        for i=0, n_cubes-1 do begin
          temp = where(type_pol_str2 eq type_pol_str1[i], count_typepol)
          if count_typepol eq 0 then begin
            message, 'all_type_pol cannot be used with these folders, they contain '+ $
              'different sets of types & pols'
          endif
        endfor

        type_pol_str = [[type_pol_str1],[type_pol_str1]]

        titles = type_pol_str[*, 0]
      end
      'diff_ratio': begin
        type_pol_str = strarr(n_sets, 2)
        type_pol_str[0, *] = cube_types[0] + '_' + pols[0]
        type_pol_str[1, *] = cube_types[1] + '_' + pols[0]
        type_pol_str[2, *] = cube_types[0] + '_' + pols[max_pol]
        type_pol_str[3, *] = cube_types[1] + '_' + pols[max_pol]

        titles = strarr(2,3)
        for i=0, 1 do begin
          titles[i,*] = [type_pol_str[2*i,0] + '/' + type_pol_str[2*i+1,0], $
            type_pol_str[2*i,1] + '/' + type_pol_str[2*i+1,1], 'Ratio Difference']
        endfor
      end
    endcase
  endif else begin
    case comp_type of
      'diff': begin
        type_pol_str = strarr(1, 2)
        type_pol_str[0]= cube_types[0] + '_' + pols[0]
        type_pol_str[1] = cube_types[max_type] + '_' + pols[max_pol]

        titles = type_pol_str[0] + '-' + type_pol_str[1]
      end
      'diff_ratio': begin
        type_pol_str = strarr(n_sets, 2)

        type_pol_str[0, *] = cube_types + '_' + pols[0]
        type_pol_str[1, *] = cube_types + '_' + pols[max_pol]

        titles = [type_pol_str[0,0] + '/' + type_pol_str[0,1], $
          type_pol_str[1,0] + '/' + type_pol_str[1,1], 'Ratio Difference']
      end
      'ratio': begin
        type_pol_str = strarr(1, 2)
        type_pol_str[0]= cube_types[0] + '_' + pols[0]
        type_pol_str[1] = cube_types[max_type] + '_' + pols[max_pol]

        titles = type_pol_str[0] + '/' + type_pol_str[1]
      end
    endcase
  endelse

  if comp_type eq 'diff' then cube_set_dim = n_cubes else cube_set_dim = n_sets
  type_pol_locs = intarr(cube_set_dim, 2)
  for i=0, cube_set_dim-1 do begin
    wh_type_pol1 = where(file_struct_arr1.type_pol_str eq type_pol_str[i,0], count_type_pol)
    if count_type_pol eq 0 then begin
      message, 'requested type_pol not found: ' + type_pol_str[i,0] + ' not in ' + folder_names[0]
    endif else begin
      type_pol_locs[i,0] = wh_type_pol1
    endelse
    wh_type_pol2 = where(file_struct_arr2.type_pol_str eq type_pol_str[i,1], count_type_pol)
    if count_type_pol eq 0 then begin
      message, 'requested type_pol not found: ' + type_pol_str[i,1] + ' not in ' + folder_names[1]
    endif else begin
      type_pol_locs[i,1] = wh_type_pol2
    endelse
  endfor

  if tag_exist(file_struct_arr1, 'n_obs') then n_obs1 = file_struct_arr1[0].n_obs
  if tag_exist(file_struct_arr2, 'n_obs') then n_obs2 = file_struct_arr2[0].n_obs

  if n_elements(n_obs1) gt 0 and n_elements(n_obs2) gt 0 then begin
    if not tag_exist(plot_options, 'note') then begin
      plot_options = create_plot_options(plot_options = plot_options, $
        note = '(' + number_formatter(round(mean(file_struct_arr1[0].n_obs))) + ',' + $
        number_formatter(round(mean(file_struct_arr2[0].n_obs))) + ')')
    endif else begin
      plot_options.note = plot_options.note + ' (' + $
        number_formatter(round(mean(file_struct_arr1[0].n_obs))) + $
        ',' + number_formatter(round(mean(file_struct_arr2[0].n_obs))) + ')'
    endelse
  endif


  if keyword_set(plot_2d_options.plot_wedge_line) then begin
    z0_freq = 1420.40 ;; MHz

    freq_use = file_struct_arr1[0].frequencies
    if n_elements(freq_ch_range) ne 0 then begin
      if max(freq_ch_range) gt n_elements(freq_use)-1 then message, 'invalid freq_ch_range'
      freq_use = freq_use[freq_ch_range[0]:freq_ch_range[1]]
    endif

    redshifts = z0_freq/freq_use - 1 ;; frequencies will be identical if kx, ky, kz match
    mean_redshift = mean(redshifts)

    cosmology_measures, mean_redshift, wedge_factor = wedge_factor
    ;; assume 20 degrees from pointing center to first null
    source_dist = 20d * !dpi / 180d
    fov_amp = wedge_factor * source_dist

    ;; calculate angular distance to horizon
    max_theta = max([file_struct_arr1[0].max_theta, file_struct_arr2[0].max_theta])
    horizon_amp = wedge_factor * ((max_theta+90d) * !dpi / 180d)

    wedge_amp = [fov_amp, horizon_amp]

    wedge_1dbin_names = ['', '_no_fov_wedge', '_no_horizon_wedge']
  endif else begin
    wedge_amp = 0d
    wedge_1dbin_names = ''
  endelse


  if comp_type eq 'diff' then begin
    mid_savefile_2d = strarr(n_slices, cube_set_dim)
    if not keyword_set(plot_slices) then begin
      savefiles_1d = strarr(cube_set_dim, n_elements(wedge_1dbin_names))
      mid_savefile_3d = strarr(n_slices, cube_set_dim)
    endif
  endif
  input_savefile1 = strarr(n_slices, cube_set_dim)
  input_savefile2 = strarr(n_slices, cube_set_dim)

  if plot_options.pub then plotfiles_2d = strarr(n_slices)

  for slice_i=0, n_slices-1 do begin
    for cube_i=0, cube_set_dim-1 do begin

      type_pol1 = type_pol_str[cube_i, 0]
      type_pol2 = type_pol_str[cube_i, 1]

      if comp_type eq 'diff' then begin

        if n_elements(obs_info.info_files) eq 2 then begin
          fileparts_1 = strsplit(file_struct_arr1[0].general_filebase, '_', /extract)
          fileparts_2 = strsplit(file_struct_arr2[0].general_filebase, '_', /extract)
          match_test = strcmp(fileparts_1, fileparts_2)
          wh_diff = where(match_test eq 0, count_diff, complement = wh_same, $
            ncomplement = count_same)
        endif

        if n_elements(savefilebase) eq 0 then begin
          if n_elements(obs_info.info_files) eq 1 then begin
            savefilebase_use = file_struct_arr1[0].general_filebase + '_' + $
              type_pol1 + op_str + type_pol2
          endif else begin
            if count_diff eq 0 then begin
              savefilebase_use = file_struct_arr1[0].general_filebase + '_' + type_pol1
              if type_pol1 ne type_pol2 then begin
                savefilebase_use = savefilebase_use + '_' + type_pol1 + op_str + type_pol2
              endif
            endif else begin
              if count_same gt 0 then begin
                if type_pol1 ne type_pol2 then begin
                  savefilebase_use = strjoin(fileparts_1[wh_same], '_') + '__' + $
                    strjoin(fileparts_1[wh_diff]) + '_' + type_pol1 + op_str + $
                    strjoin(fileparts_2[wh_diff]) + '_' + type_pol2
                endif else begin
                  savefilebase_use = strjoin(fileparts_1[wh_same], '_') + '__' + $
                    strjoin(fileparts_1[wh_diff]) + op_str + $
                    strjoin(fileparts_2[wh_diff]) + '__' + type_pol1
                endelse
              endif else begin
                if type_pol1 ne type_pol2 then begin
                  savefilebase_use = file_struct_arr1[0].general_filebase + '_' + type_pol_str[0] + $
                    op_str + file_struct_arr2[0].general_filebase + '_' + type_pol_str[1]
                endif else begin
                  savefilebase_use = file_struct_arr1[0].general_filebase + op_str + $
                    file_struct_arr2[0].general_filebase + '__' + type_pol1
                endelse
              endelse
            endelse
          endelse
        endif else begin
          if type_pol1 eq type_pol2 then begin
            savefilebase_use = savefilebase + '_' + type_pol1
          endif else begin
            savefilebase_use = savefilebase + '_' + type_pol1 + op_str + type_pol2
          endelse
        endelse
      endif

      if keyword_set(plot_slices) then begin
        if comp_type eq 'diff' then begin
          mid_savefile_2d[slice_i, cube_i] = save_path + savefilebase_use + $
            kperp_density_names + '_power_' + slice_tags[slice_i] + '_plane.idlsave'
        endif

        slice_filebase1 = file_struct_arr1[type_pol_locs[cube_i, 0]].savefile_froot + $
          file_struct_arr1[type_pol_locs[cube_i, 0]].savefilebase + $
          file_struct_arr1[type_pol_locs[cube_i, 0]].power_tag
        input_savefile1[slice_i, cube_i] = slice_filebase1 + '_' + slice_tags[slice_i] + $
          '_plane.idlsave'

        slice_filebase2 = file_struct_arr2[type_pol_locs[cube_i, 1]].savefile_froot + $
          file_struct_arr2[type_pol_locs[cube_i, 1]].savefilebase + $
          file_struct_arr2[type_pol_locs[cube_i, 1]].power_tag
        input_savefile2[slice_i, cube_i] = slice_filebase2 + '_' + slice_tags[slice_i] + $
          '_plane.idlsave'

      endif else begin
        if comp_type eq 'diff' then begin
          mid_savefile_3d[slice_i, cube_i] = save_path + savefilebase_use + '_power.idlsave'
          mid_savefile_2d[slice_i, cube_i] = save_path + savefilebase_use + $
            kperp_density_names + '_2dkpower.idlsave'

          for wedge_i=0, n_elements(wedge_1dbin_names)-1 do begin
            savefiles_1d[cube_i, wedge_i] = save_path + savefilebase_use + $
              kperp_density_names + wedge_1dbin_names[wedge_i] + '_1dkpower.idlsave'
          endfor

          input_savefile1[slice_i, cube_i] = file_struct_arr1[type_pol_locs[cube_i, 0]].power_savefile
          if file_test(input_savefile1[slice_i, cube_i]) eq 0 then begin
            message, 'No power file for ' + type_pol1 + ' and info_file: ' + $
              obs_info.info_files[0]
          endif

          input_savefile2[slice_i, cube_i] = file_struct_arr2[type_pol_locs[cube_i, 1]].power_savefile
          if file_test(input_savefile2[slice_i, cube_i]) eq 0 then begin
            message, 'No power file for ' + type_pol2 + ' and info_file: ' + $
              obs_info.info_files[n_elements(obs_info.info_files)-1]
          endif
        endif else begin

          input_savefile1[slice_i, cube_i] = file_struct_arr1[type_pol_locs[cube_i, 0]].savefile_froot + $
            file_struct_arr1[type_pol_locs[cube_i, 0]].savefilebase + $
            file_struct_arr1[type_pol_locs[cube_i, 0]].power_tag + fadd_2dbin + $
            kperp_density_names[0] + '_2dkpower.idlsave'
          input_savefile2[slice_i, cube_i] = file_struct_arr2[type_pol_locs[cube_i, 1]].savefile_froot + $
            file_struct_arr2[type_pol_locs[cube_i, 1]].savefilebase + $
            file_struct_arr2[type_pol_locs[cube_i, 1]].power_tag + fadd_2dbin + $
            kperp_density_names[n_wtcuts-1] + '_2dkpower.idlsave'

        endelse
      endelse

      if plot_options.pub then begin
        if keyword_set(plot_slices) then begin
          plotfiles_2d[slice_i] = plot_options.plot_path + plot_filebase + '_power_' + $
            slice_tags[slice_i] + '_plane' + plot_exten
        endif else begin

          plotfiles_2d[slice_i] = plot_options.plot_path + plot_filebase + '_2dkpower' + plot_exten

          if comp_type eq 'diff' and cube_i eq 0 then begin
            plotfiles_1d = plot_options.plot_path + plot_filebase + ['', wedge_1dbin_names] + $
              '_1dkpower' + plot_exten
          endif
        endelse
      endif else plot_exten = ''

      kperp_lambda_conv = getvar_savefile(input_savefile1[slice_i, cube_i], 'kperp_lambda_conv')
      hubble_param = getvar_savefile(input_savefile1[slice_i, cube_i], 'hubble_param')
      if not tag_exist(plot_2d_options, 'kperp_plot_range') then begin
        kperp_plot_range = [5./kperp_lambda_conv, min([file_struct_arr1.kspan/2., $
          file_struct_arr1.max_baseline_lambda])/kperp_lambda_conv]

        if plot_options.hinv then begin
          kperp_plot_range = kperp_plot_range / hubble_param
        endif

        plot_2d_options = create_plot_2d_options(plot_2d_options = plot_2d_options, $
          kperp_plot_range = kperp_plot_range)
      endif
    endfor
  endfor

  compare_files = {input_savefile1:input_savefile1, input_savefile2:input_savefile2, $
    titles:titles, n_slices:n_slices, $
    wedge_amp:wedge_amp, kperp_density_names:kperp_density_names, $
    kperp_plot_range:kperp_plot_range}

  if n_elements(n_cubes) gt 0 then compare_files = create_struct(compare_files, 'n_cubes', n_cubes)
  if n_elements(n_sets) gt 0 then compare_files = create_struct(compare_files, 'n_sets', n_sets)

  if plot_options.pub then begin
    compare_files = create_struct(compare_files, 'plotfiles_2d', plotfiles_2d)
  endif
  if n_elements(mid_savefile_2d) gt 0 then begin
    compare_files = create_struct(compare_files, 'mid_savefile_2d', mid_savefile_2d)
  endif
  if n_elements(mid_savefile_3d) gt 0 then begin
    compare_files = create_struct(compare_files, 'mid_savefile_3d', mid_savefile_3d)
  endif
  if n_elements(savefiles_1d) gt 0 then begin
    compare_files = create_struct(compare_files, 'savefiles_1d', savefiles_1d)
  endif
  if n_elements(slice_tags) gt 0 then begin
    compare_files = create_struct(compare_files, 'slice_tags', slice_tags)
  endif
end
