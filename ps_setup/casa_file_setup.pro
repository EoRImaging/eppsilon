function casa_file_setup, filename, pol_inc, save_path = save_path, $
    refresh_info = refresh_info, weight_savefilebase = weight_savefilebase_in, $
    variance_savefilebase = variance_savefilebase_in, $
    freq_ch_range = freq_ch_range, freq_flags = freq_flags, freq_flag_name = freq_flag_name, $
    uvf_savefilebase = uvf_savefilebase_in, savefilebase = savefilebase_in, $
    uvf_options = uvf_options, ps_options = ps_options

  if n_elements(pol_inc) ne 0 then pol_inc_in = pol_inc

  ;; check to see if filename is an info file
  wh_info = strpos(filename, '_info.idlsave')
  if max(wh_info) gt -1 then begin
    if n_elements(filename) gt 1 then message, 'only 1 info file can be passed'
    info_file = filename

    restore, info_file

    if n_elements(save_path) ne 0 then froot = save_path $
    else info_filebase = cgRootName(info_file, directory=froot)

    if n_elements(pol_inc_in) ne 0 then begin
      match, pol_inc, pol_inc_in, suba, subb, count = count_pol
      if count_pol ne n_elements(pol_inc_in) then refresh_info = 1
    endif

    if keyword_set(refresh_info) then begin
      if n_elements(metadata_struct) gt 0 then datafile = metadata_struct.datafile $
      else message, 'refresh_info is set but info file does not contain metadata_struct'

      datafile_test = file_test(datafile)
      if min(datafile_test) eq 0 then begin
        ;; can't find datafile. try some other potential folders
        if n_elements(save_path) ne 0 then begin
          datafile_test = file_test(save_path + file_basename(datafile))
          if min(datafile_test) gt 0 then datafile = save_path + file_basename(datafile)
        endif

        if min(datafile_test) eq 0 then begin
          datafile_test = file_test(file_dirname(info_file, /mark_directory) + file_basename(datafile))
          if min(datafile_test) gt 0 then datafile = file_dirname(info_file, /mark_directory) + file_basename(datafile)
        endif

        if min(datafile_test) eq 0 then message, 'refresh_info is set but datafile cannot be found'
      endif

    endif
  endif

  if keyword_set(refresh_info) or max(wh_info) eq -1 then begin
    if n_elements(datafile) eq 0 then datafile = filename

    if n_elements(pol_inc) eq 0 then pol_inc = ['xx', 'yy']
    pol_enum = ['xx', 'yy']
    npol = n_elements(pol_inc)
    pol_num = intarr(npol)
    for i=0, npol-1 do begin
      wh = where(pol_enum eq pol_inc[i], count)
      if count eq 0 then message, 'pol ' + pol_inc[i] + ' not recognized.'
      pol_num[i] = wh[0]
    endfor
    pol_inc = pol_enum[pol_num[uniq(pol_num, sort(pol_num))]]


    nfiles = n_elements(datafile)

    if nfiles gt 2 then message, 'Only 1 or 2 datafiles supported'
    if nfiles eq 2 then if datafile[0] eq datafile[1] then begin
      print, 'datafiles are identical'
      datafile = datafile[0]
      nfiles = 1
    endif

    if n_elements(savefilebase_in) gt 1 then message, 'only one savefilebase allowed'

    if nfiles eq 1 then infile_label = '' else begin
      data_filebase = [cgRootName(datafile[0]), cgRootName(datafile[1])]

      fileparts_1 = strsplit(data_filebase[0], '_', /extract)
      fileparts_2 = strsplit(data_filebase[1], '_', /extract)
      match_test = strcmp(fileparts_1, fileparts_2)
      wh_diff = where(match_test eq 0, count_diff, complement = wh_same, ncomplement = count_same)

      if count_diff gt 0 then infile_label = [strjoin(fileparts_1[wh_diff]), strjoin(fileparts_2[wh_diff])] $
      else infile_label = strarr(2)
    endelse

    if n_elements(savefilebase_in) eq 0 then begin
      if nfiles eq 1 then begin
        data_filebase = cgRootName(datafile, directory=datafile_dir)
        if n_elements(save_path) ne 0 then froot = save_path $
        else froot = datafile_dir
        uvf_froot = froot

        general_filebase = data_filebase

      endif else begin
        if n_elements(save_path) ne 0 then froot = save_path $
        else data_filebase = cgRootName(datafile[0], directory=froot)
        data_filebase = [data_filebase, cgRootName(datafile[1])]

        if count_diff eq 0 then general_filebase = data_filebase[0] + '_joint' else begin
          if count_same gt 0 then general_filebase = strjoin(fileparts_1[wh_same], '_') + '__' + strjoin(fileparts_1[wh_diff]) $
            + '_' + strjoin(fileparts_2[wh_diff]) + '_joint' $
          else general_filebase = data_filebase[0] + data_filebase[1] + '_joint'
        endelse

        general_filebase = general_filebase

      endelse
    endif else begin
      if n_elements(save_path) gt 0 then froot = save_path else begin
        savefilebase_in_base = cgRootName(savefilebase_in[0], directory=froot)
        if froot eq '.' then data_filebase = cgRootName(datafile[0], directory=froot)
      endelse

      general_filebase = savefilebase_in_base
    endelse

    if not keyword_set(refresh_info) then begin
      ; check for pre-saved info file. If it exists, restore it and return structure. Otherwise construct structure.
      info_file = froot + general_filebase + '_info.idlsave'
      info_filetest = file_test(info_file)
      if info_filetest eq 1 then begin
        ;; check if info file has the right polarizations
        info_pol_inc = getvar_savefile(info_file, 'pol_inc')
        match2, pol_inc, info_pol_inc, suba, subb
        if min([suba, subb]) ge 0 then begin
          ;; looks good, restore & check for directory structure changes
          file_struct_arr = casa_file_setup(info_file, pol_inc, save_path = save_path, $
              refresh_info = refresh_info, weight_savefilebase = weight_savefilebase_in, $
              variance_savefilebase = variance_savefilebase_in, $
              freq_ch_range = freq_ch_range, freq_flags = freq_flags, freq_flag_name = freq_flag_name, $
              uvf_savefilebase = uvf_savefilebase_in, savefilebase = savefilebase_in, $
              uvf_options = uvf_options, ps_options = ps_options)

          return, file_struct_arr
        endif
      endif
    endif

    datafile_test = file_test(datafile)
    if min(datafile_test) eq 0 then message, 'datafile not found'

    void = getvar_savefile(datafile[0], names = varnames)

    type_inc = ['res']
    ntypes = n_elements(type_inc)
    ncubes = npol * ntypes
    type_pol_str = strarr(ncubes)
    for i=0, npol-1 do type_pol_str[ntypes*i:i*ntypes+(ntypes-1)] = type_inc + '_' + pol_inc[i]

    data_varname = pol_inc + '_data'
    weight_varname = pol_inc + '_weights'
    variance_varname = pol_inc + '_variances'

    if n_elements(weight_savefilebase_in) gt 0 and n_elements(weight_savefilebase_in) ne nfiles then $
      message, 'if weight_savefilebase is specified it must have the same number of elements as datafiles'
    if n_elements(uvf_savefilebase_in) gt 0 and n_elements(uvf_savefilebase_in) ne nfiles then $
      message, 'if uvf_savefilebase is specified it must have the same number of elements as datafiles'

    for j=0, nfiles-1 do begin
      void = getvar_savefile(datafile[j], names = varnames)

      for i=0, ncubes-1 do begin
        wh = where(strlowcase(varnames) eq strlowcase(data_varname[i]), count)
        if count eq 0 then message, data_varname[i] + ' is not present in datafile (datafile=' + datafile[j] + ')'

        data_size = getvar_savefile(datafile[j], data_varname[i], /return_size)
        if data_size[0] gt 0 then this_data_dims = data_size[1:data_size[0]]

        if i eq 0 and j eq 0 then data_dims = this_data_dims else if total(abs(this_data_dims - data_dims)) ne 0 then message, 'data dimensions in files do not match'
      endfor

      for i=0, npol-1 do begin
        wh = where(strlowcase(varnames) eq strlowcase(weight_varname[i]), count)

        if count eq 0 then message, weight_varname[i] + ' is not present in datafile (datafile=' + datafile[j] + ')'

        wt_size = getvar_savefile(datafile[j], weight_varname[i], /return_size)
        if wt_size[0] gt 0 then wt_dims = wt_size[1:wt_size[0]]
        if total(abs(wt_dims - data_dims)) ne 0 then message, 'weight and data dimensions in files do not match'

        wh = where(strlowcase(varnames) eq strlowcase(variance_varname[i]), count)

        if count eq 0 then begin
          print, variance_varname[i] + ' is not present in datafile (datafile=' + datafile[j] + '). Using weights^2 instead of variances'
          no_var = 1
        endif else begin
          no_var = 0
          var_size = getvar_savefile(datafile[j], variance_varname[i], /return_size)
          if var_size[0] gt 0 then var_dims = var_size[1:var_size[0]]
          if total(abs(var_dims - data_dims)) ne 0 then message, 'variance and data dimensions in files do not match'
        endelse
      endfor

      if j eq 0 then frequencies = getvar_savefile(datafile[j], 'frequencies') else $
        if total(abs(frequencies-getvar_savefile(datafile[j], 'frequencies'))) ne 0 then $
        message, 'frequencies do not match between datafiles'

      if j eq 0 then ra_vals = getvar_savefile(datafile[j], 'ra_vals') else $
        if total(abs(ra_vals-getvar_savefile(datafile[j], 'ra_vals'))) ne 0 then $
        message, 'ra_vals do not match between datafiles'

      if j eq 0 then dec_vals = getvar_savefile(datafile[j], 'dec_vals') else $
        if total(abs(dec_vals-getvar_savefile(datafile[j], 'dec_vals'))) ne 0 then $
        message, 'dec_vals do not match between datafiles'
    endfor

    ;; convert to MHz if in Hz
    if mean(frequencies) gt 1000. then frequencies = frequencies/1e6

    n_freq = n_elements(frequencies)

    ;; made up for now
    freq_resolution = (frequencies[1]-frequencies[0])*1e6;; native resolution of visibilities in Hz
    time_resolution = 2.;; native resolution of visibilities in s
    n_vis = 1e9 + fltarr(nfiles)

    ;; pointing offset from zenith (for calculating horizon distance for wedge line)
    max_theta = 0.

    ;; degpix
    degpix = mean(abs([ra_vals[1]-ra_vals[0], dec_vals[1]-dec_vals[0]]))

    ;; making up based on image size/resolution
    kspan = 1/(degpix * !pi/180.)
    kpix = kspan/mean([n_elements(ra_vals), n_elements(dec_vals)])
    max_baseline_lambda = kspan

    metadata_struct = {datafile:datafile, weightfile: datafile, variancefile:datafile, $
      cube_varname:data_varname, weight_varname:weight_varname, variance_varname:variance_varname, $
      frequencies:frequencies, freq_resolution:freq_resolution, time_resolution:time_resolution, $
      n_vis:n_vis, max_baseline_lambda:max_baseline_lambda, max_theta:max_theta, degpix:degpix, kpix:kpix, kspan:kspan, $
      general_filebase:general_filebase, type_pol_str:type_pol_str, infile_label:infile_label, ntypes:ntypes}

    if no_var then metadata_struct = create_struct(metadata_struct, 'no_var', 1)

    if n_elements(vis_noise) gt 0 then metadata_struct = create_struct(metadata_struct, 'vis_noise', vis_noise)

    save, filename = info_file, metadata_struct, pol_inc
  endif

  npol = n_elements(pol_inc)
  nfiles = n_elements(metadata_struct.datafile)
  ncubes = npol * metadata_struct.ntypes
  if tag_exist(metadata_struct, 'no_var') then no_var = 1 else no_var = 0

  if n_elements(freq_flags) ne 0 then begin
    freq_mask = intarr(n_elements(metadata_struct.frequencies)) + 1
    freq_mask[freq_flags] = 0

    wh_freq_use = where(freq_mask gt 0, count_freq_use)
    if count_freq_use lt 3 then message, 'Too many frequencies flagged'
  endif

  if n_elements(freq_flags) ne 0 then begin
    if min(freq_flags) lt 0 or max(freq_flags) gt n_elements(metadata_struct.frequencies) $
      then message, 'invalid freq_flags'
  endif

  file_tags = create_file_tags(freq_ch_range = freq_ch_range, freq_flags = freq_flags, $
    freq_flag_name = freq_flag_name, uvf_options = uvf_options, ps_options = ps_options)

  wt_file_label = '_weights_' + strlowcase(pol_inc)
  file_label = '_' + strlowcase(metadata_struct.type_pol_str)
  savefilebase = metadata_struct.general_filebase + file_tags.uvf_tag + file_label

  if n_elements(uvf_savefilebase_in) lt nfiles then begin
    if nfiles eq 1 then begin
      ;; if we're only dealing with one file and uvf_savefilebase isn't specified then use same base for uvf files
      uvf_savefilebase = general_filebase + file_tags.uvf_tag + file_label
    endif else begin
      ;; need 2 uvf files for each type/pol
      if n_elements(save_path) ne 0 then uvf_froot = replicate(save_path, nfiles, ncubes) else begin
        uvf_froot = strarr(nfiles, ncubes)
        ;; test datafile path to see if it exists. if not, use froot
        for i=0, nfiles-1 do begin
          datafile_path = file_dirname(metadata_struct.datafile[i], /mark_directory)
          if file_test(datafile_path, /directory) then uvf_froot[i, *] = datafile_path else uvf_froot[i, *] = froot
        endfor
      endelse
      uvf_savefilebase = strarr(nfiles, ncubes)
      uvf_label = strarr(nfiles, ncubes)
      for i=0, nfiles-1 do begin
        uvf_savefilebase[i, *] = cgRootName(metadata_struct.datafile[i]) + file_tags.uvf_tag + file_label
        uvf_label[i, *] = metadata_struct.infile_label[i] + file_label
      endfor
    endelse
  endif else begin
    if n_elements(save_path) gt 0 then uvf_froot = replicate(save_path, nfiles, ncubes) else begin
      temp = file_dirname(uvf_savefilebase_in, /mark_directory)
      uvf_froot = strarr(nfiles, ncubes)
      for i=0, nfiles-1 do begin
        if temp[i] ne '.' then uvf_froot[i,*] = temp[i] else begin
          datafile_path = file_dirname(metadata_struct.datafile[i], /mark_directory)
          if file_test(datafile_path, /directory) then uvf_froot[i, *] = datafile_path else uvf_froot[i, *] = froot
        endelse
      endfor
    endelse
    uvf_savefilebase = file_basename(uvf_savefilebase_in) + file_tags.uvf_tag + file_label
  endelse

  ;; add sw tag to general_filebase so that plotfiles havefile_tags.uvf_tag in them
  general_filebase = metadata_struct.general_filebase + file_tags.uvf_tag


  uvf_savefile = uvf_froot + uvf_savefilebase + '_uvf.idlsave'
  uf_savefile = uvf_froot + uvf_savefilebase + '_uf_plane.idlsave'
  vf_savefile = uvf_froot + uvf_savefilebase + '_vf_plane.idlsave'
  uv_savefile = uvf_froot + uvf_savefilebase + '_uv_plane.idlsave'

  uf_raw_savefile = uvf_froot + uvf_savefilebase + '_uf_plane_raw.idlsave'
  vf_raw_savefile = uvf_froot + uvf_savefilebase + '_vf_plane_raw.idlsave'
  uv_raw_savefile = uvf_froot + uvf_savefilebase + '_uv_plane_raw.idlsave'


  uf_sum_savefile = froot + savefilebase + '_sum_uf_plane.idlsave'
  vf_sum_savefile = froot + savefilebase + '_sum_vf_plane.idlsave'
  uv_sum_savefile = froot + savefilebase + '_sum_uv_plane.idlsave'
  uf_diff_savefile = froot + savefilebase + '_diff_uf_plane.idlsave'
  vf_diff_savefile = froot + savefilebase + '_diff_vf_plane.idlsave'
  uv_diff_savefile = froot + savefilebase + '_diff_uv_plane.idlsave'


  kcube_savefile = froot + savefilebase + file_tags.kcube_tag + '_kcube.idlsave'
  power_savefile = froot + savefilebase + file_tags.power_tag + '_power.idlsave'
  fits_power_savefile = froot + savefilebase + file_tags.power_tag + '_power.fits'

  if n_elements(weight_savefilebase_in) eq 0 then begin
    wt_base = cgrootname(metadata_struct.weightfile[0], directory = wt_froot)
    if file_test(wt_froot, /directory) eq 0 then wt_froot = froot

    if nfiles eq 1 then begin
      if n_elements(save_path) gt 0 then wt_froot = save_path

      weight_savefilebase = wt_base + file_tags.uvf_tag + wt_file_label
    endif else begin
      if n_elements(save_path) gt 0 then wt_froot = save_path else begin
        wt_froot = strarr(nfiles, npol)
        for i=0, nfiles-1 do begin
          wtfile_path = file_dirname(metadata_struct.weightfile[i], /mark_directory)
          if file_test(wtfile_path, /directory) then wt_froot[i, *] = wtfile_path else wt_froot[i, *] = froot
        endfor
      endelse

      weight_savefilebase = strarr(nfiles, npol)
      for i=0, nfiles-1 do weight_savefilebase[i, *] = cgrootname(metadata_struct.weightfile[i]) + file_tags.uvf_tag + wt_file_label
    endelse
  endif else begin
    if n_elements(save_path) gt 0 then wt_froot = save_path else begin
      temp = file_dirname(weight_savefilebase_in, /mark_directory)
      if nfiles eq 1 then begin
        if temp ne '.' then wt_froot = temp else begin
          wt_froot = file_dirname(metadata_struct.weightfile, /mark_directory)
          if file_test(wt_froot) eq 0 then wt_froot = froot
        endelse
      endif else begin
        wt_froot = strarr(nfiles, npol)
        for i=0, nfiles-1 do begin
          if temp[i] ne '.' then uvf_froot[i,*] = temp[i] else begin
            wtfile_path = file_dirname(metadata_struct.weightfile[i], /mark_directory)
            if file_test(wtfile_path, /directory) then wt_froot[i, *] = wtfile_path else wt_froot[i, *] = froot
          endelse
        endfor
      endelse
    endelse
    weight_savefilebase = file_basename(weight_savefilebase_in) + file_tags.uvf_tag + wt_file_label
  endelse

  uvf_weight_savefile = wt_froot + weight_savefilebase + '_uvf.idlsave'
  uf_weight_savefile = wt_froot + weight_savefilebase + '_uf_plane.idlsave'
  vf_weight_savefile = wt_froot + weight_savefilebase + '_vf_plane.idlsave'
  uv_weight_savefile = wt_froot + weight_savefilebase + '_uv_plane.idlsave'


  for i=0, ncubes-1 do begin
    pol_index = i / metadata_struct.ntypes
    type_index = i mod metadata_struct.ntypes

    res_uvf_inputfiles = strarr(nfiles, 2)
    res_uvf_varname = strarr(nfiles, 2)


    file_struct = {datafile:metadata_struct.datafile, weightfile:metadata_struct.weightfile, variancefile:metadata_struct.variancefile, $
      datavar:metadata_struct.cube_varname[i], weightvar:metadata_struct.weight_varname[pol_index], variancevar:metadata_struct.variance_varname[pol_index], $
      frequencies:metadata_struct.frequencies, freq_resolution:metadata_struct.freq_resolution, time_resolution:metadata_struct.time_resolution, $
      n_vis:metadata_struct.n_vis, max_baseline_lambda:metadata_struct.max_baseline_lambda, max_theta:metadata_struct.max_theta, $
      degpix:metadata_struct.degpix, kpix:metadata_struct.kpix, kspan:metadata_struct.kspan, $
      uf_savefile:uf_savefile[*,i], vf_savefile:vf_savefile[*,i], uv_savefile:uv_savefile[*,i], $
      uf_raw_savefile:uf_raw_savefile[*,i], vf_raw_savefile:vf_raw_savefile[*,i], $
      uv_raw_savefile:uv_raw_savefile[*,i], $
      uf_sum_savefile:uf_sum_savefile[i], vf_sum_savefile:vf_sum_savefile[i], $
      uv_sum_savefile:uv_sum_savefile[i], uf_diff_savefile:uf_diff_savefile[i], $
      vf_diff_savefile:vf_diff_savefile[i], uv_diff_savefile:uv_diff_savefile[i], $
      uf_weight_savefile:uf_weight_savefile[*, pol_index], vf_weight_savefile:vf_weight_savefile[*, pol_index], $
      uv_weight_savefile:uv_weight_savefile[*, pol_index], $
      kcube_savefile:kcube_savefile[i], power_savefile:power_savefile[i], fits_power_savefile:fits_power_savefile[i],$
      savefile_froot:froot, savefilebase:savefilebase[i], general_filebase:general_filebase, $
      weight_savefilebase:weight_savefilebase[*,pol_index], $
      res_uvf_inputfiles:res_uvf_inputfiles, res_uvf_varname:res_uvf_varname, $
      file_label:file_label[i], uvf_label:uvf_label[*,i], wt_file_label:wt_file_label[pol_index], $
      uvf_tag:file_tags.uvf_tag, kcube_tag:file_tags.kcube_tag, power_tag:file_tags.power_tag, type_pol_str:metadata_struct.type_pol_str[i]}

    file_struct = create_struct(file_struct, 'uvf_savefile', uvf_savefile[*,i], $
      'uvf_weight_savefile', uvf_weight_savefile[*, pol_index])

    if no_var then file_struct = create_struct(file_struct, 'no_var', 1)

    if n_elements(freq_mask) gt 0 then file_struct = create_struct(file_struct, 'freq_mask', freq_mask)

    if tag_exist(metadata_struct, 'vis_noise') gt 0 then file_struct = create_struct(file_struct, 'vis_noise', reform(metadata_struct.vis_noise[pol_index, *]))

    if i eq 0 then file_struct_arr = replicate(file_struct, ncubes) else file_struct_arr[i] = file_struct
  endfor

  return, file_struct_arr

end
