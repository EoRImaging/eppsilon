function fhd_file_setup, filename, pol_inc, weightfile = weightfile, variancefile = variancefile, pixelfile = pixelfile, $
    dirtyvar = dirtyvar, modelvar = modelvar, weightvar = weightvar, variancevar = variancevar, $
    pixelvar = pixelvar, save_path = save_path, image = image, dft_ian = dft_ian, $
    weight_savefilebase = weight_savefilebase_in, $
    uvf_savefilebase = uvf_savefilebase_in, savefilebase = savefilebase_in, $
    freq_ch_range = freq_ch_range, freq_flags = freq_flags, freq_flag_name = freq_flag_name, $
    spec_window_type = spec_window_type, noise_sim = noise_sim, std_power = std_power, refresh_info = refresh_info
    
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
      else datafile = file_struct_arr[0].datafile
      
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
      
    endif else begin
    
      if n_elements(metadata_struct) eq 0 then begin
      
        npol = n_elements(pol_inc)
        
        if keyword_set(noise_sim) then begin
          type_inc = ['noisesim']
          ntypes = n_elements(type_inc)
          ncubes = npol * ntypes
          type_pol_str = strarr(ncubes)
          for i=0, npol-1 do type_pol_str[ntypes*i:i*ntypes+(ntypes-1)] = type_inc + '_' + pol_inc[i]
        endif else begin
          ;; check for the existence of dirty, model, residual cubes
          if max(strmatch(file_struct_arr.datavar, 'dirty*cube',/fold_case)) then begin
            if n_elements(dirtyvar) eq 0 then dirty_varname = strupcase('dirty_' + pol_inc + '_cube') else dirty_varname = dirtyvar
            if n_elements(dirty_varname) ne npol then $
              if n_elements(dirty_varname) eq 1 then dirty_varname = replicate(dirty_varname, npol) $
            else message, 'dirtyvar must be a scalar or have the same number of elements as pol_inc'
            
            type_inc = ['dirty']
            if npol gt 1 then cube_varname = transpose(dirty_varname) else cube_varname = dirty_varname
          endif
          
          if max(strmatch(file_struct_arr.datavar, 'model*cube',/fold_case)) then begin
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
          
          if max(strmatch(file_struct_arr.datavar, 'res*cube',/fold_case)) then begin
            if n_elements(residualvar) eq 0 then residual_varname = strupcase('res_' + pol_inc + '_cube') else residual_varname = residualvar
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
        undefine, type_inc, ncubes, dirty_varname, residual_varname, model_varname, cube_varname
        
        nfiles = n_elements(file_struct_arr[0].datafile)
        if nfiles eq 1 then infile_label = '' else begin
          data_filebase = [cgRootName(file_struct_arr[0].datafile[0]), cgRootName(file_struct_arr[0].datafile[1])]
          
          fileparts_1 = strsplit(data_filebase[0], '_', /extract)
          fileparts_2 = strsplit(data_filebase[1], '_', /extract)
          match_test = strcmp(fileparts_1, fileparts_2)
          wh_diff = where(match_test eq 0, count_diff, complement = wh_same, ncomplement = count_same)
          
          if count_diff gt 0 then infile_label = [strjoin(fileparts_1[wh_diff]), strjoin(fileparts_2[wh_diff])] $
          else infile_label = strarr(2)
          
          undefine, data_filebase, fileparts_1, fileparts_2, match_test, wh_diff, count_diff, wh_same, count_same
        endelse
        undefine, nfiles
        
        
        ;; test for sw_tag in general_filebase and remove it if present.
        sw_tag_list = ['hann', 'ham', 'blm', 'ntl', 'bn', 'bh']
        for i=0, n_elements(sw_tag_list)-1 do begin
          sw_pos = strpos(file_struct_arr[0].general_filebase, '_' + sw_tag_list[i], /reverse_search)
          if sw_pos ne -1 then begin
            general_filebase = strmid(file_struct_arr[0].general_filebase, 0, sw_pos)
            break
          endif
        endfor
        if n_elements(general_filebase) eq 0 then general_filebase = file_struct_arr[0].general_filebase
        undefine, sw_tag_list
        
        metadata_struct = {datafile:file_struct_arr[0].datafile, weightfile:file_struct_arr[0].weightfile, $
          variancefile:file_struct_arr[0].variancefile, cube_varname:reform(file_struct_arr.datavar, ntypes, npol), $
          weight_varname:file_struct_arr.weightvar, variance_varname:file_struct_arr.variancevar, $
          frequencies:file_struct_arr[0].frequencies, freq_resolution:file_struct_arr[0].freq_resolution, $
          time_resolution:file_struct_arr[0].time_resolution, n_vis:file_struct_arr[0].n_vis, $
          max_baseline_lambda:file_struct_arr[0].max_baseline_lambda, max_theta:file_struct_arr[0].max_theta, $
          degpix:file_struct_arr[0].degpix, kpix:file_struct_arr[0].kpix, kspan:file_struct_arr[0].kspan, $
          general_filebase:general_filebase, type_pol_str:type_pol_str, infile_label:infile_label, ntypes:ntypes}
          
        if tag_exist(file_struct_arr[0], 'nside') then begin
          healpix = 1
          metadata_struct = create_struct(metadata_struct, 'pixelfile', file_struct_arr[0].pixelfile, 'pixel_varname', $
            file_struct_arr[0].pixelvar, 'nside', file_struct_arr[0].nside)
        endif else healpix = 0
        
        if tag_exist(file_struct_arr[0], 'no_var') then begin
          metadata_struct = create_struct(metadata_struct, 'no_var', 1)
          no_var = 1
        endif else no_var = 0
        
        
        save, filename = info_file, metadata_struct, pol_inc
      endif
    endelse
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
    ;; look for polarization labels in the filename
    ;    pol_pos = stregex(datafile, '[x-y][x-y]')
    ;    if min(pol_pos) eq -1 then begin
    if nfiles gt 2 then message, 'Only 1 or 2 datafiles supported'
    if nfiles eq 2 then if datafile[0] eq datafile[1] then begin
      print, 'datafiles are identical'
      datafile = datafile[0]
      nfiles = 1
    endif
    ;      temp = strarr(nfiles, npol)
    ;      for i = 0, nfiles-1 do temp[i,*] = datafile[i]
    ;    endif else begin
    ;      if nfiles ne 2*npol and nfiles ne npol then message, 'Only 1 or 2 datafiles per polarization supported'
    ;      if nfiles eq npol then nfiles = 1 else nfiles = 2
    ;      temp = strarr(nfiles, npol)
    ;      for i=0, npol do begin
    ;        this_pol = where(strpos(datafile, pol_inc[i]) ne -1, count)
    ;        if count ne nfiles then message, 'Expected ' + numberformatter(nfiles) + ' for ' + pol_inc[i] + ' polarization, only found ' + count
    ;        temp[*, i] = datafile[this_pol]
    ;      endfor
    ;      for i=0, nfiles*npol do begin
    ;        wh_duplicate = where(temp eq temp[i], count_duplicate)
    ;        if count_duplicate gt 0 then message, 'multiple identical datafiles indentified.'
    ;      endfor
    ;      datafile = temp
    ;    endelse
    
    if nfiles eq 2 and keyword_set(noise_sim) then begin
      datafile=datafile[0, *]
      nfiles=1
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
          file_struct_arr = fhd_file_setup(info_file, save_path = save_path, weight_savefilebase = weight_savefilebase_in, $
            uvf_savefilebase = uvf_savefilebase_in, savefilebase = savefilebase_in)
            
          return, file_struct_arr
        endif
      endif
    endif
    
    datafile_test = file_test(datafile)
    if min(datafile_test) eq 0 then message, 'datafile not found'
    
    void = getvar_savefile(datafile[0], names = varnames)
    
    if keyword_set(noise_sim) then begin
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
      
      
      
    for j=0, nfiles-1 do begin
      void = getvar_savefile(datafile[j], names = varnames)
      
      wh_nside = where(strlowcase(varnames) eq 'nside', count_nside)
      if count_nside gt 0 then this_healpix = 1 else this_healpix = 0
      
      for i=0, ncubes-1 do begin
        wh = where(strlowcase(varnames) eq strlowcase(cube_varname[i]), count)
        if count eq 0 then message, cube_varname[i] + ' is not present in datafile (datafile=' + datafile[j] + ')'
        
        data_size = getvar_savefile(datafile[j], cube_varname[i], /return_size)
        if data_size[0] gt 0 then this_data_dims = data_size[1:data_size[0]]
        
        if i eq 0 and j eq 0 then data_dims = this_data_dims else if total(abs(this_data_dims - data_dims)) ne 0 then message, 'data dimensions in files do not match'
      endfor
      
      void = getvar_savefile(weightfile[j], names = wt_varnames)
      void = getvar_savefile(variancefile[j], names = var_varnames)
      for i=0, npol-1 do begin
        wh = where(strlowcase(wt_varnames) eq strlowcase(weight_varname[i]), count)
        
        if count eq 0 then message, weight_varname[i] + ' is not present in weightfile (weightfile=' + weightfile[j] + ')'
        
        wt_size = getvar_savefile(weightfile[j], weight_varname[i], /return_size)
        if wt_size[0] gt 0 then wt_dims = wt_size[1:wt_size[0]]
        if total(abs(wt_dims - data_dims)) ne 0 then message, 'weight and data dimensions in files do not match'
        
        wh = where(strlowcase(var_varnames) eq strlowcase(variance_varname[i]), count)
        
        if count eq 0 then begin
          print, variance_varname[i] + ' is not present in variancefile (variancefile=' + variancefile[j] + '). Using weights^2 instead of variances'
          no_var = 1
        endif else begin
          no_var = 0
          var_size = getvar_savefile(variancefile[j], variance_varname[i], /return_size)
          if var_size[0] gt 0 then var_dims = var_size[1:var_size[0]]
          if total(abs(var_dims - data_dims)) ne 0 then message, 'variance and data dimensions in files do not match'
        endelse
      endfor
      
      if j gt 0 then if (this_healpix eq 1 and healpix eq 0) or (this_healpix eq 0 and healpix eq 1) then $
        message, 'One datafile is in healpix and the other is not.'
      if this_healpix eq 1 then begin
        if j eq 0 then nside = getvar_savefile(datafile[j], 'nside') else begin
          nside1 = nside
          nside = getvar_savefile(datafile[j], 'nside')
          if nside1 ne nside then message, 'nside parameter does not agree between datafiles'
          undefine, nside1
        endelse
        healpix = 1
        
        if n_elements(pixelfile) ne nfiles then begin
          if j eq 0 then pixelfile = [datafile[j]] else pixelfile = [pixelfile, datafile[j]]
        endif
        if j eq 0 then begin
          if n_elements(pixelvar) eq 0 then pixel_varname = replicate('hpx_inds', nfiles) else pixel_varname =  pixelvar
        endif
        
        void = getvar_savefile(pixelfile[j], names = pix_varnames)
        wh = where(strlowcase(pix_varnames) eq strlowcase(pixel_varname[j]), count)
        
        if count eq 0 then message, pixel_varname[j] + ' is not present in pixelfile (pixelfile=' + pixelfile[j] + ')'
        
        pix_size = getvar_savefile(pixelfile[j], pixel_varname[j], /return_size)
        if pix_size[0] gt 0 then pix_dims = pix_size[1:pix_size[0]]
        if total(abs(pix_dims - data_dims[0])) ne 0 then message, 'Number of Healpix pixels does not match data dimension'
        
        
        
        
      endif else healpix = 0
      
      wh_obs = where(strlowcase(varnames) eq 'obs', count_obs)
      if count_obs ne 0 then obs_arr = getvar_savefile(datafile[j], 'obs') else begin
        wh_obs = where(strlowcase(varnames) eq 'obs_arr', count_obs)
        obs_arr = getvar_savefile(datafile[j], 'obs_arr')
      endelse
      
      if count_obs ne 0 then begin
      
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
          
        theta_vals = angle_difference(obs_radec_vals[*,1], obs_radec_vals[*,0], zen_radec_vals[*,1], zen_radec_vals[*,0], $
          /degree, /nearest)
        if j eq 0 then max_theta = max(theta_vals) else max_theta = max([max_theta, theta_vals])
        
        if j eq 0 then n_vis = fltarr(nfiles)
        n_vis[j] = total(obs_arr.n_vis)
        
        if tag_exist(obs_arr[0], 'vis_noise') then begin
          vis_noise_arr = fltarr([n_obs, size(*obs_arr[0].vis_noise, /dimension)])
          for i=0, n_obs-1 do vis_noise_arr[i, *, *] = *obs_arr[i].vis_noise
          vis_noise_arr = total(vis_noise_arr, 1)/n_obs
          noise_dims = size(vis_noise_arr, /dimensions)
          if noise_dims[0] ne npol or noise_dims[1] ne n_freq then message, 'vis_noise dimensions do not match npol, n_freq'
          if j eq 0 then vis_noise = vis_noise_arr else vis_noise = (vis_noise_arr + vis_noise*j)/(j+1.)
        endif
        
      endif else message, 'no obs or obs_arr in datafile'
      
      
      wh_navg = where(strlowcase(varnames) eq 'n_avg', count_obs)
      if count_obs ne 0 then begin
        if j eq 0 then n_avg = getvar_savefile(datafile[j], 'n_avg') else begin
          n_avg1 = n_avg
          n_avg = getvar_savefile(datafile[j], 'n_avg')
          if n_avg1 ne n_avg then message, 'n_avg parameter does not agree between datafiles'
          undefine, n_avg1
        endelse
      endif else begin
        print, 'no n_avg present, assuming n_avg=32'
        if j eq 0 then n_avg = 32 else if n_avg ne 32 then message, 'n_avg parameter does not agree between datafiles'
      endelse
      if healpix then data_nfreq = data_dims[1] else data_nfreq = data_dims[2]
      if n_freq/n_avg ne data_nfreq then message, 'number of frequencies does not match number of data slices'
      
      
      
    endfor
    
    ;; fix uvf savefiles for gridded uv
    if healpix eq 0 and not keyword_set(image) then begin
      uvf_savefile = datafile
      uvf_weight_savefile = weightfile
    endif
    
    n_freqbins = n_freq / n_avg
    frequencies = dblarr(n_freqbins)
    if n_elements(vis_noise) gt 0 then vis_noise_avg = fltarr(npol, n_freqbins)
    for i=0, n_freqbins-1 do begin
      frequencies[i] = mean(freq[i*n_avg:i*n_avg+(n_avg-1)]) / 1e6 ;; in MHz
      if n_elements(vis_noise) gt 0 then vis_noise_avg[*, i] = sqrt(total(vis_noise[*, i*n_avg:i*n_avg+(n_avg-1)]^2., 2))
    endfor
    if n_elements(vis_noise) gt 0 then vis_noise = vis_noise_avg
    
    metadata_struct = {datafile: datafile, weightfile: weightfile, variancefile:variancefile, $
      cube_varname:cube_varname, weight_varname:weight_varname, variance_varname:variance_varname, $
      frequencies:frequencies, freq_resolution:freq_resolution, time_resolution:time_resolution, $
      n_vis:n_vis, max_baseline_lambda:max_baseline_lambda, max_theta:max_theta, degpix:degpix, kpix:kpix, kspan:kspan, $
      general_filebase:general_filebase, type_pol_str:type_pol_str, infile_label:infile_label, ntypes:ntypes}
      
    if healpix then metadata_struct = create_struct(metadata_struct, 'pixelfile', pixelfile, 'pixel_varname', pixel_varname, 'nside', nside)
    
    if no_var then metadata_struct = create_struct(metadata_struct, 'no_var', 1)
    
    if n_elements(vis_noise) gt 0 then metadata_struct = create_struct(metadata_struct, 'vis_noise', vis_noise)
    
    save, filename = info_file, metadata_struct, pol_inc
  endif
  
  npol = n_elements(pol_inc)
  nfiles = n_elements(metadata_struct.datafile)
  ncubes = npol * metadata_struct.ntypes
  if tag_exist(metadata_struct, 'nside') then healpix = 1 else healpix = 0
  if tag_exist(metadata_struct, 'no_var') then no_var = 1 else no_var = 0
  
  type_list = ['Hann', 'Hamming', 'Blackman', 'Nutall', 'Blackman-Nutall', 'Blackman-Harris']
  sw_tag_list = ['hann', 'ham', 'blm', 'ntl', 'bn', 'bh']
  if n_elements(spec_window_type) ne 0 then begin
    wh_type = where(strlowcase(type_list) eq strlowcase(spec_window_type), count_type)
    if count_type eq 0 then message, 'Spectral window type not recognized.' else begin
      spec_window_type = type_list[wh_type[0]]
      sw_tag = '_' + sw_tag_list[wh_type[0]]
    endelse
  endif else sw_tag = ''
  
  if n_elements(freq_flags) ne 0 then begin
    freq_mask = intarr(n_elements(metadata_struct.frequencies)) + 1
    freq_mask[freq_flags] = 0
    
    wh_freq_use = where(freq_mask gt 0, count_freq_use)
    if count_freq_use lt 3 then message, 'Too many frequencies flagged'
    
  endif
  
  if n_elements(freq_ch_range) ne 0 then begin
    if min(freq_ch_range) lt 0 or max(freq_ch_range) - min(freq_ch_range) lt 3 then message, 'invalid freq_ch_range'
    fch_tag = '_ch' + number_formatter(min(freq_ch_range)) + '-' + number_formatter(max(freq_ch_range))
  endif else fch_tag = ''
  
  if n_elements(freq_flags) ne 0 then begin
    if min(freq_flags) lt 0 or max(freq_flags) gt n_elements(metadata_struct.frequencies) then message, 'invalid freq_ch_range'
    if n_elements(freq_flag_name) eq 0 then freq_flag_name = '' else $
      if size(freq_flag_name, /type) ne 7 then freq_flag_name = number_formatter(freq_flag_name)
    flag_tag = '_flag' + freq_flag_name
  endif else flag_tag = ''
  fch_tag = fch_tag + flag_tag
  
  if keyword_set(std_power) then power_tag = '_stdp' else power_tag = ''
  power_tag = power_tag + sw_tag
  
  if keyword_set(dft_ian) then dft_label = '_ian' else dft_label = ''
  
  wt_file_label = '_weights_' + strlowcase(pol_inc)
  file_label = '_' + strlowcase(metadata_struct.type_pol_str)
  savefilebase = metadata_struct.general_filebase + fch_tag + file_label
  
  if n_elements(uvf_savefilebase_in) lt nfiles then begin
    if nfiles eq 1 then begin
      ;; if we're only dealing with one file and uvf_savefilebase isn't specified then use same base for uvf files
      uvf_savefilebase = general_filebase + fch_tag + file_label + dft_label
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
        uvf_savefilebase[i, *] = cgRootName(metadata_struct.datafile[i]) + fch_tag + file_label + dft_label
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
    uvf_savefilebase = file_basename(uvf_savefilebase_in) + fch_tag + file_label + dft_label
  endelse
  
  ;; add sw tag to general_filebase so that plotfiles have fch_tag & sw_tag in them
  general_filebase = metadata_struct.general_filebase + fch_tag
  
  
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
  
  
  kcube_savefile = froot + savefilebase + power_tag + '_kcube.idlsave'
  power_savefile = froot + savefilebase + power_tag + '_power.idlsave'
  fits_power_savefile = froot + savefilebase + power_tag + '_power.fits'
  
  if n_elements(weight_savefilebase_in) eq 0 then begin
    wt_base = cgrootname(metadata_struct.weightfile[0], directory = wt_froot)
    if file_test(wt_froot, /directory) eq 0 then wt_froot = froot
    
    if nfiles eq 1 then begin
      if n_elements(save_path) gt 0 then wt_froot = save_path
      
      weight_savefilebase = wt_base + fch_tag + wt_file_label
    endif else begin
      if n_elements(save_path) gt 0 then wt_froot = save_path else begin
        wt_froot = strarr(nfiles, npol)
        for i=0, nfiles-1 do begin
          wtfile_path = file_dirname(metadata_struct.weightfile[i], /mark_directory)
          if file_test(wtfile_path, /directory) then wt_froot[i, *] = wtfile_path else wt_froot[i, *] = froot
        endfor
      endelse
      
      weight_savefilebase = strarr(nfiles, npol)
      for i=0, nfiles-1 do weight_savefilebase[i, *] = cgrootname(metadata_struct.weightfile[i]) + fch_tag + wt_file_label
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
    weight_savefilebase = file_basename(weight_savefilebase_in) + fch_tag + wt_file_label
  endelse
  
  uvf_weight_savefile = wt_froot + weight_savefilebase + '_uvf.idlsave'
  
  for i=0, ncubes-1 do begin
    pol_index = i / metadata_struct.ntypes
    type_index = i mod metadata_struct.ntypes
    
    if keyword_set(noise_sim) then begin
      data_varname = ''
      res_uvf_inputfiles = strmid(uvf_savefile[i], 0, strpos(uvf_savefile[i], 'noisesim')) + 'dirty' + $
        strmid(uvf_savefile[i], strpos(uvf_savefile[i], 'noisesim')+strlen('noisesim'))
      res_uvf_varname = strarr(n_elements(res_uvf_inputfiles)) + 'data_cube'
    endif else begin
      if metadata_struct.ntypes gt 1 then data_varname = metadata_struct.cube_varname[type_index, pol_index] else data_varname = metadata_struct.cube_varname[pol_index]
      if data_varname ne '' then begin
        res_uvf_inputfiles = strarr(nfiles,2)
        res_uvf_varname = strarr(nfiles,2)
      endif else begin
        if healpix or keyword_set(image) then begin
          res_uvf_inputfiles = uvf_savefile[*, metadata_struct.ntypes*pol_index:metadata_struct.ntypes*pol_index+1]
          res_uvf_varname = strarr(nfiles, 2) + 'data_cube'
        endif else begin
          res_uvf_inputfiles = strarr(nfiles,2)
          res_uvf_varname = strarr(nfiles,2)
          for j=0, nfiles-1 do begin
            res_uvf_inputfiles[j,*] = datafile[j]
            res_uvf_varname[j,*] = [metadata_struct.cube_varname[0, pol_index], metadata_struct.cube_varname[1, pol_index]]
          endfor
        endelse
      endelse
    endelse
    
    file_struct = {datafile:metadata_struct.datafile, weightfile:metadata_struct.weightfile, variancefile:metadata_struct.variancefile, $
      datavar:data_varname, weightvar:metadata_struct.weight_varname[pol_index], variancevar:metadata_struct.variance_varname[pol_index], $
      frequencies:metadata_struct.frequencies, freq_resolution:metadata_struct.freq_resolution, time_resolution:metadata_struct.time_resolution, $
      n_vis:metadata_struct.n_vis, max_baseline_lambda:metadata_struct.max_baseline_lambda, max_theta:metadata_struct.max_theta, $
      degpix:metadata_struct.degpix, kpix:metadata_struct.kpix, kspan:metadata_struct.kspan, $
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
      file_label:file_label[i], uvf_label:uvf_label[*,i], wt_file_label:wt_file_label[pol_index], $
      fch_tag:fch_tag, power_tag:power_tag, type_pol_str:metadata_struct.type_pol_str[i]}
      
    if healpix or keyword_set(image) then file_struct = create_struct(file_struct, 'uvf_savefile', uvf_savefile[*,i], $
      'uvf_weight_savefile', uvf_weight_savefile[*, pol_index])
      
    if healpix then file_struct = create_struct(file_struct, 'pixelfile', metadata_struct.pixelfile, 'pixelvar', $
      metadata_struct.pixel_varname, 'nside', metadata_struct.nside)
      
    if no_var then file_struct = create_struct(file_struct, 'no_var', 1)
    
    if n_elements(freq_mask) gt 0 then file_struct = create_struct(file_struct, 'freq_mask', freq_mask)
    
    if tag_exist(metadata_struct, 'vis_noise') gt 0 then file_struct = create_struct(file_struct, 'vis_noise', reform(metadata_struct.vis_noise[pol_index, *]))
    
    if i eq 0 then file_struct_arr = replicate(file_struct, ncubes) else file_struct_arr[i] = file_struct
  endfor
  
  return, file_struct_arr
  
end