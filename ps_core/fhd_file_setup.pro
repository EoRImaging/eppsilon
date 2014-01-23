function fhd_file_setup, filename, pol_inc, weightfile = weightfile, variancefile = variancefile, pixelfile = pixelfile, $
    dirtyvar = dirtyvar, modelvar = modelvar, weightvar = weightvar, variancevar = variancevar, $
    pixelvar = pixelvar, save_path = save_path, image = image, dft_ian = dft_ian, $
    weight_savefilebase = weight_savefilebase_in, $
    uvf_savefilebase = uvf_savefilebase_in, savefilebase = savefilebase_in, $
    freq_ch_range = freq_ch_range, spec_window_type = spec_window_type, noise_sim = noise_sim
    
  ;; check to see if filename is an info file
  wh_info = strpos(filename, '_info.idlsave')
  if max(wh_info) eq -1 then datafile = filename else begin
    if n_elements(filename) gt 1 then message, 'only 1 info file can be passed'
    info_file = filename
    
    restore, info_file
    
    ;; check to see if datafile(s), uvf, kcube or power files exist
    files_test = file_test([reform(file_struct_arr.datafile, n_elements(file_struct_arr.datafile)), $
      reform(file_struct_arr.uvf_savefile, n_elements(file_struct_arr.uvf_savefile)), $
      reform(file_struct_arr.kcube_savefile, n_elements(file_struct_arr.kcube_savefile)), $
      reform(file_struct_arr.power_savefile, n_elements(file_struct_arr.power_savefile))])
    if max(files_test) eq 0 then begin
      ;; can't find any of the files
      ;; test to see if directories have changed (ie info_file created on different system)
      ;; if so set them to something reasonable
    
      froot_test = file_test(file_struct_arr.savefile_froot, /directory)
      wh_froot_fail = where(froot_test eq 0, count_froot_fail)
      uvf_froot_test = file_test(file_dirname(file_struct_arr.uvf_savefile), /directory)
      wh_uvf_froot_fail = where(uvf_froot_test eq 0, count_uvf_froot_fail)
      wt_froot_test = file_test(file_dirname(file_struct_arr.uvf_weight_savefile), /directory)
      wh_wt_froot_fail = where(wt_froot_test eq 0, count_wt_froot_fail)
      
      if n_elements(save_path) gt 0 then begin
        if count_froot_fail gt 0 then froot = save_path
        if count_uvf_froot_fail gt 0 then uvf_froot = save_path
        if count_wt_froot_fail gt 0 then wt_froot = save_path
      endif else begin
        void = cgrootname(info_file, directory = info_dir)
        
        if count_froot_fail gt 0 then begin
          if n_elements(savefilebase_in) then begin
            void = cgrootname(savefilebase_in, directory = savefilebase_dir)
            froot = savefilebase_dir
          endif else froot = info_dir
        endif
        
        if count_uvf_froot_fail gt 0 then begin
          if n_elements(uvf_savefilebase_in) then begin
            void = cgrootname(uvf_savefilebase_in, directory = uvf_savefilebase_dir)
            uvf_froot = uvf_savefilebase_dir
          endif else uvf_froot = info_dir
        endif
        
        if count_wt_froot_fail gt 0 then begin
          if n_elements(weight_savefilebase_in) then begin
            void = cgrootname(weight_savefilebase_in, directory = wt_savefilebase_dir)
            wt_froot = wt_savefilebase_dir
          endif else  wt_froot = info_dir
        endif
      endelse     
      
    endif
    
    if n_elements(froot) ne 0 then begin
      file_struct_arr.savefile_froot[*] = froot
      file_struct_arr.uf_sum_savefile = froot + file_basename(file_struct_arr.uf_sum_savefile)
      file_struct_arr.vf_sum_savefile = froot + file_basename(file_struct_arr.vf_sum_savefile)
      file_struct_arr.uv_sum_savefile = froot + file_basename(file_struct_arr.uv_sum_savefile)
      file_struct_arr.uf_diff_savefile = froot + file_basename(file_struct_arr.uf_diff_savefile)
      file_struct_arr.vf_diff_savefile = froot + file_basename(file_struct_arr.vf_diff_savefile)
      file_struct_arr.uv_diff_savefile = froot + file_basename(file_struct_arr.uv_diff_savefile)
      file_struct_arr.kcube_savefile = froot + file_basename(file_struct_arr.kcube_savefile)
      file_struct_arr.power_savefile = froot + file_basename(file_struct_arr.power_savefile)
      file_struct_arr.fits_power_savefile = froot + file_basename(file_struct_arr.fits_power_savefile)
    endif
    
    if n_elements(uvf_froot) ne 0 then begin
      file_struct_arr.uvf_savefile = uvf_froot + file_basename(file_struct_arr.uvf_savefile)
      file_struct_arr.uf_savefile = uvf_froot + file_basename(file_struct_arr.uf_savefile)
      file_struct_arr.vf_savefile = uvf_froot + file_basename(file_struct_arr.vf_savefile)
      file_struct_arr.uv_savefile = uvf_froot + file_basename(file_struct_arr.uv_savefile)
      file_struct_arr.uf_raw_savefile = uvf_froot + file_basename(file_struct_arr.uf_raw_savefile)
      file_struct_arr.vf_raw_savefile = uvf_froot + file_basename(file_struct_arr.vf_raw_savefile)
      file_struct_arr.uv_raw_savefile = uvf_froot + file_basename(file_struct_arr.uv_raw_savefile)
    endif
    
    if n_elements(wt_froot) ne 0 then file_struct_arr.uvf_weight_savefile = wt_froot + file_basename(file_struct_arr.uvf_weight_savefile)
    
    ;; check again (with new directories) to see if datafile(s), uvf, kcube or power files exist
    files_test = file_test([reform(file_struct_arr.datafile, n_elements(file_struct_arr.datafile)), $
      reform(file_struct_arr.uvf_savefile, n_elements(file_struct_arr.uvf_savefile)), $
      reform(file_struct_arr.kcube_savefile, n_elements(file_struct_arr.kcube_savefile)), $
      reform(file_struct_arr.power_savefile, n_elements(file_struct_arr.power_savefile))])
      
    if max(files_test) eq 0 then message, 'Cannot find any datafile(s), uvf, kcube or power files, please specify paths using keywords.'
    
    
    return, file_struct_arr
  endelse
  
  nfiles = n_elements(datafile)
  if nfiles gt 2 then message, 'only 1 or 2 datafiles is supported'
  if nfiles eq 2 then if datafile[0] eq datafile[1] then begin
    print, 'datafiles are identical'
    datafile = datafile[0]
    nfiles = 1
  endif
  
  if n_elements(savefilebase_in) gt 1 then message, 'only one savefilebase allowed'
  
  if n_elements(freq_ch_range) ne 0 then begin
    if min(freq_ch_range) lt 0 or max(freq_ch_range) - min(freq_ch_range) lt 3 then message, 'invalid freq_ch_range'
    fch_tag = '_ch' + number_formatter(min(freq_ch_range)) + '-' + number_formatter(max(freq_ch_range))
  endif else fch_tag = ''
  
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
      if n_elements(save_path) ne 0 then froot = save_path $
      else data_filebase = cgRootName(datafile, directory=froot)
      uvf_froot = froot
      
      general_filebase = data_filebase + fch_tag
      infile_label = ''
      
    endif else begin
      if n_elements(save_path) ne 0 then froot = save_path $
      else data_filebase = cgRootName(datafile[0], directory=froot)
      data_filebase = [data_filebase, cgRootName(datafile[1])]
      
      if count_diff eq 0 then general_filebase = data_filebase[0] + '_joint' else begin
        if count_same gt 0 then general_filebase = strjoin(fileparts_1[wh_same], '_') + '__' + strjoin(fileparts_1[wh_diff]) $
          + '_' + strjoin(fileparts_2[wh_diff]) + '_joint' $
        else general_filebase = data_filebase[0] + data_filebase[1] + '_joint'
      endelse
      
      general_filebase = general_filebase + fch_tag
      
    endelse
  endif else begin
    if n_elements(save_path) gt 0 then froot = save_path else begin
      savefilebase_in_base = cgRootName(savefilebase_in[0], directory=froot)
      if froot eq '.' then data_filebase = cgRootName(datafile[0], directory=froot)
    endelse
    
    general_filebase = savefilebase_in_base + fch_tag
  endelse
  
  ; check for pre-saved info file. If it exists, restore it and return structure. Otherwise construct structure.
  info_file = froot + general_filebase + '_info.idlsave'
  info_filetest = file_test(info_file)
  if info_filetest eq 1 then begin
    ;; check if info file has the right polarizations
    info_pol_inc = getvar_savefile(info_file, 'pol_inc')
    match2, pol_inc, info_pol_inc, suba, subb
    if min([suba, subb]) ge 0 then begin
      ;; looks good, restore & check for directory structure changes
      file_struct_arr = fhd_file_setup(info_file, save_path = save_path, weight_savefilebase = weight_savefilebase_in, $
        uvf_savefilebase = uvf_savefilebase_in, savefilebase = savefilebase_in)
        
      return, file_struct_arr
    endif
  endif
  
  datafile_test = file_test(datafile)
  if min(datafile_test) eq 0 then message, 'datafile not found'
  
  file_obj = obj_new('idl_savefile', datafile[0])
  varnames = file_obj->names()
  
  
  if keyword_set(noise_sim) then begin
    datafile = datafile[0]
    nfiles = 1
    
    type_inc = ['noisesim']
    ntypes = n_elements(type_inc)
    ncubes = npol * ntypes
    type_pol_str = strarr(ncubes)
    for i=0, npol-1 do type_pol_str[ntypes*i:i*ntypes+(ntypes-1)] = type_inc + '_' + pol_inc[i]
  endif else begin
    ;; check for the existence of dirty, model, residual cubes
    if max(strmatch(varnames, 'dirty*cube',/fold_case)) then begin
      if n_elements(dirtyvar) eq 0 then dirty_varname = strupcase('dirty_' + pol_inc + '_cube') else dirty_varname = dirtyvar
      if n_elements(dirty_varname) ne npol then $
        if n_elements(dirty_varname) eq 1 then dirty_varname = replicate(dirty_varname, npol) $
      else message, 'dirtyvar must be a scalar or have the same number of elements as pol_inc'
      
      type_inc = ['dirty']
      if npol gt 1 then cube_varname = transpose(dirty_varname) else cube_varname = dirty_varname
    endif
    
    if max(strmatch(varnames, 'model*cube',/fold_case)) then begin
      if n_elements(modelvar) eq 0 then model_varname = strupcase('model_' + pol_inc + '_cube') else model_varname = modelvar
      if n_elements(model_varname) ne npol then $
        if n_elements(model_varname) eq 1 then model_varname = replicate(model_varname, npol) $
      else message, 'modelvar must be a scalar or have the same number of elements as pol_inc'
      
      if n_elements(type_inc) eq 0 then begin
        type_inc = ['model']
        if npol gt 1 then cube_varname = transpose(model_varname) else cube_varname = model_varname
      endif else begin
        type_inc = [type_inc, 'model']
        if npol gt 1 then cube_varname = [cube_varname,transpose(model_varname)] else cube_varname = [cube_varname, model_varname]
      endelse
    endif
    
    if max(strmatch(varnames, 'res*cube',/fold_case)) then begin
      if n_elements(residualvar) eq 0 then residual_varname = strupcase('res_' + pol_inc + '_cube') else dirty_varname = residualvar
      if n_elements(residual_varname) ne npol then $
        if n_elements(residual_varname) eq 1 then residual_varname = replicate(residual_varname, npol) $
      else message, 'residualvar must be a scalar or have the same number of elements as pol_inc'
      
      if n_elements(type_inc) eq 0 then begin
        type_inc = ['res']
        if npol gt 1 then cube_varname = transpose(residual_varname) else cube_varname = residual_varname
      endif else begin
        type_inc = [type_inc, 'res']
        if npol gt 1 then cube_varname = [cube_varname,transpose(residual_varname)] else cube_varname = [cube_varname, residual_varname]
      endelse
    endif else if type_inc eq ['dirty', 'model'] then begin
      ;; residual can be constructed from dirty-model
      type_inc = [type_inc, 'res']
      if npol gt 1 then cube_varname = [cube_varname,transpose(strarr(npol))] else cube_varname = [cube_varname, '']
    endif
    ntypes = n_elements(type_inc)
    ncubes = npol * ntypes
    type_pol_str = strarr(ncubes)
    for i=0, npol-1 do type_pol_str[ntypes*i:i*ntypes+(ntypes-1)] = type_inc + '_' + pol_inc[i]
  endelse
  
  
  if n_elements(weightvar) eq 0 then weight_varname = strupcase('weights_' + pol_inc + '_cube') else weight_varname = weightvar
  if n_elements(weight_varname) ne npol then $
    if n_elements(weight_varname) eq 1 then weight_varname = replicate(weight_varname, npol) $
  else message, 'weightvar must be a scalar or have the same number of elements as pol_inc'
  if n_elements(variancevar) eq 0 then variance_varname = strupcase('variance_' + pol_inc + '_cube') else $
    variance_varname = variancevar
  if n_elements(variance_varname) ne npol then $
    if n_elements(variance_varname) eq 1 then variance_varname = replicate(variance_varname, npol) $
  else message, 'variancevar must be a scalar or have the same number of elements as pol_inc'
  
  wt_file_label = '_weights_' + strlowcase(pol_inc)
  file_label = '_' + strlowcase(type_pol_str)
  
  if n_elements(weightfile) eq 0 then weightfile = datafile $
  else if n_elements(weightfile) ne nfiles then message, 'weightfile must have the same number of elements as datafile'
  if n_elements(variancefile) eq 0 then variancefile = datafile $
  else if n_elements(variancefile) ne nfiles then message, 'variancefile must have the same number of elements as datafile'
  
  if n_elements(pixelfile) gt 0 and n_elements(pixelfile) ne nfiles then $
    message, 'pixelfile must have the same number of elements as datafile'
  if n_elements(pixelvar) gt 0 and n_elements(pixelvar) ne nfiles then $
    if n_elements(pixelvar) eq 1 then pixelvar = replicate(pixelvar, nfiles) $
  else message, 'pixelvar must be a scalar or have the same number of elements as datafile'
  
  
  if n_elements(weight_savefilebase_in) gt 0 and n_elements(weight_savefilebase_in) ne nfiles then $
    message, 'if weight_savefilebase is specified it must have the same number of elements as datafiles'
  if n_elements(uvf_savefilebase_in) gt 0 and n_elements(uvf_savefilebase_in) ne nfiles then $
    message, 'if uvf_savefilebase is specified it must have the same number of elements as datafiles'
    
  if n_elements(spec_window_type) ne 0 then begin
    type_list = ['Hann', 'Hamming', 'Blackman', 'Nutall', 'Blackman-Nutall', 'Blackman-Harris']
    sw_tag_list = ['hann', 'ham', 'blm', 'ntl', 'bn', 'bh']
    
    wh_type = where(strlowcase(type_list) eq strlowcase(spec_window_type), count_type)
    if count_type eq 0 then message, 'Spectral window type not recognized.' else begin
      spec_window_type = type_list[wh_type[0]]
      sw_tag = '_' + sw_tag_list[wh_type[0]]
    endelse
  endif else sw_tag = ''
  
  
  if keyword_set(dft_ian) then dft_label = '_ian' else dft_label = ''
  
  if n_elements(savefilebase_in) eq 0 or n_elements(uvf_savefilebase_in) lt nfiles then begin
    if nfiles eq 1 then begin
      if n_elements(savefilebase_in) eq 0 then savefilebase = general_filebase + file_label + sw_tag
      
      ;; if we're only dealing with one file and uvf_savefilebase isn't specified then use same base for uvf files
      if n_elements(uvf_savefilebase_in) eq 0 then uvf_savefilebase = general_filebase + fch_tag + file_label + dft_label
    endif else begin
      if n_elements(savefilebase_in) eq 0 then savefilebase = general_filebase + file_label + sw_tag
      
      if n_elements(uvf_savefilebase_in) eq 0 then begin
        ;; need 2 uvf files for each type/pol
        if n_elements(save_path) ne 0 then uvf_froot = replicate(save_path, nfiles, ncubes) else begin
          uvf_froot = strarr(nfiles, ncubes)
          for i=0, nfiles-1 do uvf_froot[i, *] = file_dirname(datafile[i], /mark_directory)
        endelse
        uvf_savefilebase = strarr(nfiles, ncubes)
        uvf_label = strarr(nfiles, ncubes)
        for i=0, nfiles-1 do begin
          uvf_savefilebase[i, *] = data_filebase[i] + fch_tag + file_label + dft_label
          uvf_label[i, *] = infile_label[i] + file_label
        endfor
      endif
    endelse
  endif
  
  if n_elements(savefilebase_in) eq 1 then savefilebase = file_basename(savefilebase_in) + fch_tag + file_label + sw_tag
  
  ;; add sw tag to general_filebase so that plotfiles have sw_tag in them
  general_filebase = general_filebase + sw_tag
  
  if n_elements(uvf_savefilebase_in) eq nfiles then begin
    if n_elements(save_path) gt 0 then uvf_froot = replicate(save_path, nfiles, ncubes) else begin
      temp = file_dirname(uvf_savefilebase_in, /mark_directory)
      uvf_froot = strarr(nfiles, ncubes)
      for i=0, nfiles-1 do if temp[i] ne '.' then uvf_froot[i,*] = temp[i] $
      else uvf_froot[i,*] = file_dirname(datafile[i], /mark_directory)
    endelse
    uvf_savefilebase = file_basename(uvf_savefilebase_in) + fch_tag + file_label + dft_label
  endif
  
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
  
  
  kcube_savefile = froot + savefilebase + '_kcube.idlsave'
  power_savefile = froot + savefilebase + '_power.idlsave'
  fits_power_savefile = froot + savefilebase + '_power.fits'
  
  if n_elements(weight_savefilebase_in) eq 0 then begin
    wt_base = cgrootname(weightfile[0], directory = wt_froot)
    
    if nfiles eq 1 then begin
      if n_elements(save_path) gt 0 then wt_froot = save_path
      
      weight_savefilebase = wt_base + fch_tag + wt_file_label
    endif else begin
      if n_elements(save_path) gt 0 then wt_froot = save_path else begin
        wt_froot = strarr(nfiles, npol)
        for i=0, nfiles-1 do wt_froot[i,*] = file_dirname(weightfile[i], /mark_directory)
      endelse
      
      weight_savefilebase = strarr(nfiles, npol)
      for i=0, nfiles-1 do weight_savefilebase[i, *] = cgrootname(weightfile[i]) + fch_tag + wt_file_label
    endelse
  endif else begin
    if n_elements(save_path) gt 0 then wt_froot = save_path else begin
      temp = file_dirname(weight_savefilebase_in, /mark_directory)
      if nfiles eq 1 then begin
        if temp ne '.' then wt_froot = temp else wt_froot = file_dirname(weightfile, /mark_directory)
      endif else begin
        wt_froot = strarr(nfiles, npol)
        for i=0, nfiles-1 do if temp[i] ne '.' then wt_froot[i,*] = temp[i] else $
          wt_froot[i,*] = file_dirname(weightfile[i], /mark_directory)
      endelse
    endelse
    weight_savefilebase = file_basename(weight_savefilebase_in) + fch_tag + wt_file_label
  endelse
  
  uvf_weight_savefile = wt_froot + weight_savefilebase + '_uvf.idlsave'
  
  for j=0, nfiles-1 do begin
    file_obj = obj_new('idl_savefile', datafile[j])
    varnames = file_obj->names()
    
    wh_nside = where(strlowcase(varnames) eq 'nside', count_nside)
    data_dims = file_obj->size(cube_varname[0], /dimensions)
    if count_nside gt 0 and n_elements(data_dims) eq 2 then this_healpix = 1 else this_healpix = 0
    
    if j gt 0 then if (this_healpix eq 1 and healpix eq 0) or (this_healpix eq 0 and healpix eq 1) then $
      message, 'One datafile is in healpix and the other is not.'
    if this_healpix eq 1 then begin
      if j eq 0 then file_obj->restore, 'nside' else begin
        nside1 = nside
        file_obj->restore, 'nside'
        if nside1 ne nside then message, 'nside parameter does not agree between datafiles'
        undefine, nside1
      endelse
      healpix = 1
      
      if n_elements(pixelfile) ne nfiles then begin
        if j eq 0 then pixelfile = datafile[j] else pixelfile = [pixelfile, datafile[j]]
      endif
      if j eq 0 then begin
        if n_elements(pixelvar) eq 0 then pixel_varname = replicate('hpx_inds', nfiles) else pixel_varname =  pixelvar
      endif
    endif else healpix = 0
    
    wh_obs = where(strlowcase(varnames) eq 'obs', count_obs)
    if count_obs ne 0 then file_obj->restore, 'obs' else begin
      wh_obs = where(strlowcase(varnames) eq 'obs_arr', count_obs)
      if count_obs ne 0 then begin
        file_obj->restore, 'obs_arr'
        if size(obs_arr,/type) eq 10 then begin
          n_obs = n_elements(obs_arr)
          
          max_baseline_vals = dblarr(n_obs)
          obs_radec_vals = dblarr(n_obs, 2)
          zen_radec_vals = dblarr(n_obs, 2)
          for i=0, n_obs-1 do begin
            if abs((*obs_arr[i]).degpix - (*obs_arr[0]).degpix) gt 0 then message, 'inconsistent degpix values in obs_arr'
            if abs((*obs_arr[i]).kpix - (*obs_arr[0]).kpix) gt 0 then message, 'inconsistent kpix values in obs_arr'
            if total(abs((*obs_arr[i]).freq - (*obs_arr[0]).freq)) gt 0 then message, 'inconsistent freq values in obs_arr'
            if abs((*obs_arr[i]).n_freq - (*obs_arr[0]).n_freq) gt 0 then message, 'inconsistent n_freq values in obs_arr'
            
            max_baseline_vals[i] = (*obs_arr[i]).max_baseline
            obs_radec_vals[i, *] = [(*obs_arr[i]).obsra, (*obs_arr[i]).obsdec]
            zen_radec_vals[i, *] = [(*obs_arr[i]).zenra, (*obs_arr[i]).zendec]
          endfor
          
          if j eq 0 then max_baseline_lambda = max(max_baseline_vals) $
          else max_baseline_lambda = max([max_baseline_lambda, max_baseline_vals])
          
          if j eq 0 then degpix = (*obs_arr[0]).degpix else if (*obs_arr[0]).degpix ne degpix then $
            message, 'degpix does not agree between datafiles'
          if j eq 0 then kpix = (*obs_arr[0]).kpix else if (*obs_arr[0]).kpix ne kpix then $
            message, 'kpix does not agree between datafiles'
            
          if j eq 0 then fhd_dim = (*obs_arr[0]).dimension else if (*obs_arr[0]).dimension ne fhd_dim then $
            message, 'dimension does not agree between datafiles'
          if j eq 0 then fhd_elem = (*obs_arr[0]).elements else if (*obs_arr[0]).elements ne fhd_elem then $
            message, 'elements does not agree between datafiles'
            
          if fhd_dim ne fhd_elem then message, 'fhd image is not square in x & y'
          kspan = kpix * fhd_dim
          
          
          if j eq 0 then freq = (*obs_arr[0]).freq else if total(abs(freq-(*obs_arr[0]).freq)) ne 0 then $
            message, 'frequencies do not agree between datafiles'
          if j eq 0 then n_freq = (*obs_arr[0]).n_freq else if (*obs_arr[0]).n_freq ne n_freq then $
            message, 'n_freq does not agree between datafiles'
          stop
        endif else begin
          n_obs = n_elements(obs_arr)
          
          if j eq 0 then max_baseline_lambda = max(obs_arr.max_baseline) $
          else max_baseline_lambda = max([max_baseline_lambda, obs_arr.max_baseline])
          
          obs_radec_vals = [[obs_arr.obsra],[obs_arr.obsdec]]
          zen_radec_vals = [[obs_arr.zenra],[obs_arr.zendec]]
          
          if total(abs(obs_arr.n_freq - obs_arr[0].n_freq)) ne 0 then message, 'inconsistent number of frequencies in obs_arr'
          if j eq 0 then n_freq = obs_arr[0].n_freq else if obs_arr[0].n_freq ne n_freq then $
            message, 'n_freq does not agree between datafiles'
            
          if total(abs(obs_arr.degpix - obs_arr[0].degpix)) ne 0 then message, 'inconsistent degpix values in obs_arr'
          if j eq 0 then degpix = obs_arr[0].degpix else if obs_arr[0].degpix ne degpix  then $
            message, 'degpix does not agree between datafiles'
            
          if total(abs(obs_arr.kpix - obs_arr[0].kpix)) ne 0 then message, 'inconsistent kpix values in obs_arr'
          if j eq 0 then kpix = obs_arr[0].kpix else if obs_arr[0].kpix ne kpix then $
            message, 'kpix does not agree between datafiles'
            
          if total(abs(obs_arr.dimension - obs_arr[0].dimension)) ne 0 then message, 'inconsistent dimension values in obs_arr'
          if j eq 0 then fhd_dim = obs_arr[0].dimension else if obs_arr[0].dimension ne fhd_dim then $
            message, 'dimension does not agree between datafiles'
            
          if total(abs(obs_arr.elements - obs_arr[0].elements)) ne 0 then message, 'inconsistent elements values in obs_arr'
          if j eq 0 then fhd_elem = obs_arr[0].elements else if obs_arr[0].elements ne fhd_elem then $
            message, 'elements does not agree between datafiles'
            
          if fhd_dim ne fhd_elem then message, 'fhd image is not square in x & y'
          kspan = kpix * fhd_dim
          
          obs_tags = tag_names(obs_arr)
          wh_freq = where(strlowcase(obs_tags) eq 'freq', count_freq)
          if count_freq ne 0 then freq_vals = obs_arr.freq else begin
            wh_bin = where(strlowcase(obs_tags) eq 'bin', count_bin)
            wh_base_info = where(strlowcase(obs_tags) eq 'baseline_info', count_base_info)
            if count_bin ne 0 or count_base_info then begin
              freq_vals = dblarr(n_freq, n_obs)
              if count_bin ne 0 then for i=0, n_obs-1 do freq_vals[*,i] = (*obs_arr[i].bin).freq $
              else for i=0, n_obs-1 do freq_vals[*,i] = (*obs_arr[i].baseline_info).freq
            endif else stop
          endelse
          if total(abs(freq_vals - rebin(freq_vals[*,0], n_freq, n_obs))) ne 0 then message, 'inconsistent freq values in obs_arr'
          if j eq 0 then freq = freq_vals[*,0] else if total(abs(freq - freq_vals[*,0])) ne 0 then $
            message, 'frequencies do not agree between datafiles'
          freq_resolution = freq[1]-freq[0]
          
          dt_vals = dblarr(n_obs)
          for i=0, n_obs-1 do begin
            times = (*obs_arr[i].baseline_info).jdate
            ;; only allow time resolutions of n*.5 sec
            dt_vals[i] = round(((times[1]-times[0])*24*3600)*2.)/2.
          endfor
          if total(abs(dt_vals - dt_vals[0])) ne 0 then message, 'inconsistent time averaging in obs_arr'
          if j eq 0 then time_resolution = dt_vals[0] else $
            if total(abs(time_resolution - dt_vals[0])) ne 0 then message, 'time averaging does not agree between datafiles'
            
        endelse
        
        theta_vals = angle_difference(obs_radec_vals[*,1], obs_radec_vals[*,0], zen_radec_vals[*,1], zen_radec_vals[*,0], $
          /degree, /nearest)
        if j eq 0 then max_theta = max(theta_vals) else max_theta = max([max_theta, theta_vals])
        
        if j eq 0 then n_vis = fltarr(nfiles)
        n_vis[j] = total(obs_arr.n_vis)
        
      endif else message, 'no obs or obs_arr in datafile'
    endelse
    
    wh_navg = where(strlowcase(varnames) eq 'n_avg', count_obs)
    if count_obs ne 0 then begin
      if j eq 0 then file_obj->restore, 'n_avg' else begin
        n_avg1 = n_avg
        file_obj->restore, 'n_avg'
        if n_avg1 ne n_avg then message, 'n_avg parameter does not agree between datafiles'
        undefine, n_avg1
      endelse
    endif else begin
      print, 'no n_avg present, assuming n_avg=32'
      if j eq 0 then n_avg = 32 else if n_avg ne 32 then message, 'n_avg parameter does not agree between datafiles'
    endelse
    
    obj_destroy, file_obj
  endfor
  
  ;; fix uvf savefiles for gridded uv
  if healpix eq 0 and not keyword_set(image) then begin
    uvf_savefile = datafile
    uvf_weight_savefile = weightfile
  endif
  
  n_freqbins = n_freq / n_avg
  frequencies = dblarr(n_freqbins)
  for i=0, n_freqbins-1 do begin
    frequencies[i] = mean(freq[i*n_avg:i*n_avg+(n_avg-1)]) / 1e6 ;; in MHz
  endfor
  
  for i=0, ncubes-1 do begin
    pol_index = i / ntypes
    type_index = i mod ntypes
    
    if keyword_set(noise_sim) then begin
      data_varname = ''
      res_uvf_inputfiles = strmid(uvf_savefile[i], 0, strpos(uvf_savefile[i], 'noisesim')) + 'dirty' + $
        strmid(uvf_savefile[i], strpos(uvf_savefile[i], 'noisesim')+strlen('noisesim'))
      res_uvf_varname = strarr(n_elements(res_uvf_inputfiles)) + 'data_cube'
    endif else begin
      data_varname = cube_varname[type_index, pol_index]
      if data_varname ne '' then begin
        res_uvf_inputfiles = strarr(nfiles,2)
        res_uvf_varname = strarr(nfiles,2)
      endif else begin
        if healpix or keyword_set(image) then begin
          res_uvf_inputfiles = uvf_savefile[*, ntypes*pol_index:ntypes*pol_index+1]
          res_uvf_varname = strarr(nfiles, 2) + 'data_cube'
        endif else begin
          res_uvf_inputfiles = strarr(nfiles,2)
          res_uvf_varname = strarr(nfiles,2)
          for j=0, nfiles-1 do begin
            res_uvf_inputfiles[j,*] = datafile[j]
            res_uvf_varname[j,*] = [cube_varname[0, pol_index], cube_varname[1, pol_index]]
          endfor
        endelse
      endelse
    endelse
    
    file_struct = {datafile: datafile, weightfile: weightfile, variancefile:variancefile, $
      datavar:data_varname, weightvar:weight_varname[pol_index], variancevar:variance_varname[pol_index], $
      frequencies:frequencies, freq_resolution:freq_resolution, time_resolution:time_resolution, $
      n_vis:n_vis, max_baseline_lambda:max_baseline_lambda, max_theta:max_theta, degpix:degpix, kpix:kpix, kspan:kspan, $
      uf_savefile:uf_savefile[*,i], vf_savefile:vf_savefile[*,i], uv_savefile:uv_savefile[*,i], $
      uf_raw_savefile:uf_raw_savefile[*,i], vf_raw_savefile:vf_raw_savefile[*,i], $
      uv_raw_savefile:uv_raw_savefile[*,i], $
      uf_sum_savefile:uf_sum_savefile[i], vf_sum_savefile:vf_sum_savefile[i], $
      uv_sum_savefile:uv_sum_savefile[i], uf_diff_savefile:uf_diff_savefile[i], $
      vf_diff_savefile:vf_diff_savefile[i], uv_diff_savefile:uv_diff_savefile[i], $
      kcube_savefile:kcube_savefile[i], power_savefile:power_savefile[i], fits_power_savefile:fits_power_savefile[i],$
      savefile_froot:froot, savefilebase:savefilebase[i], general_filebase:general_filebase, $
      weight_savefilebase:weight_savefilebase[*,pol_index], $
      res_uvf_inputfiles:res_uvf_inputfiles, res_uvf_varname:res_uvf_varname, $
      file_label:file_label[i], uvf_label:uvf_label[*,i], wt_file_label:wt_file_label[pol_index]}
      
    if healpix or keyword_set(image) then file_struct = create_struct(file_struct, 'uvf_savefile', uvf_savefile[*,i], $
      'uvf_weight_savefile', uvf_weight_savefile[*, pol_index])
      
    if healpix then file_struct = create_struct(file_struct, 'pixelfile', pixelfile, 'pixelvar', pixel_varname, 'nside', nside)
    
    if i eq 0 then file_struct_arr = replicate(file_struct, ncubes) else file_struct_arr[i] = file_struct
  endfor
  
  save, filename = info_file, file_struct_arr, pol_inc
  
  return, file_struct_arr
  
end
