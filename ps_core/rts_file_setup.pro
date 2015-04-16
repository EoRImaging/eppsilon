function rts_file_setup, filename, pol_inc, save_path = save_path, refresh_info = refresh_info, $
    weight_savefilebase = weight_savefilebase_in, variance_savefilebase = variance_savefilebase_in, $
    uvf_savefilebase = uvf_savefilebase_in, savefilebase = savefilebase_in, $
    spec_window_type = spec_window_type, delta_uv_lambda = delta_uv_lambda, max_uv_lambda = max_uv_lambda
    
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
    
    if not tag_exist(metadata_struct, 'nfiles') then begin
      ;; this is an old info file from before polarizations could be in different files.
      ;;Need to adjust the metadata_struct to the new format
    
      nfiles = n_elements(metadata_struct.datafile)
      npol = n_elements(pol_inc)
      n_freq = n_elements(metadata_struct.frequencies)
      
      datafile = strarr(npol, nfiles)
      weightfile = strarr(npol, nfiles)
      variancefile = strarr(npol, nfiles)
      beamfile = strarr(npol, nfiles)
      n_vis = fltarr(npol, nfiles)
      n_vis_freq = fltarr(npol, nfiles, n_freq)
      if tag_exist(metadata_struct, 'nside') then pixelfile = strarr(npol, nfiles)
      if tag_exist(metadata_struct, 'vis_noise') then vis_noise = fltarr(npol, nfiles, n_freq)
      for i=0, npol-1 do begin
        datafile[i,*] = metadata_struct.datafile
        weightfile[i,*] = metadata_struct.weightfile
        variancefile[i,*] = metadata_struct.variancefile
        if tag_exist(metadata_struct, 'beamfile') then beamfile[i,*] = metadata_struct.beamfile else beamfile[i,*] = metadata_struct.datafile
        n_vis[i,*] = metadata_struct.n_vis
        if tag_exist(metadata_struct, 'n_vis_freq') then begin
          nfvis_dims = size(metadata_struct.n_vis_freq, /dimension)
          if nfvis_dims[0] eq nfiles then n_vis_freq[i,*,*] = metadata_struct.n_vis_freq $
          else n_vis_freq[i,*,*] = rebin(reform(metadata_struct.n_vis_freq, 1, n_freq), nfiles, n_freq, /sample)
        endif else n_vis_freq[i,*,*] = rebin(metadata_struct.n_vis, nfiles, n_freq)/n_freq
        if tag_exist(metadata_struct, 'nside') then pixelfile[i,*] = metadata_struct.pixelfile
        if tag_exist(metadata_struct, 'vis_noise') then vis_noise = rebin(reform(metadata_struct.vis_noise[i,*], 1, n_freq), nfiles, n_freq, /sample)
      endfor
      
      
      type_inc = strarr(metadata_struct.ntypes)
      for i=0, metadata_struct.ntypes-1 do type_inc[i] = (strsplit(metadata_struct.type_pol_str[i], '_',/extract))[0]
      
      cube_varname = metadata_struct.cube_varname
      weight_varname = metadata_struct.weight_varname
      variance_varname = metadata_struct.variance_varname
      ;; if tag_exist(metadata_struct, 'beam_varname') then beam_varname = metadata_struct.beam_varname else beam_varname = strupcase('beam_' + pol_inc + '_cube')
      frequencies = metadata_struct.frequencies
      freq_resolution = metadata_struct.freq_resolution
      time_resolution = metadata_struct.time_resolution
      max_baseline_lambda = metadata_struct.max_baseline_lambda
      max_theta = metadata_struct.max_theta
      degpix = metadata_struct.degpix
      kpix = metadata_struct.kpix
      kspan = metadata_struct.kspan
      general_filebase = metadata_struct.general_filebase
      infile_label = metadata_struct.infile_label
      type_pol_str = transpose(reform(metadata_struct.type_pol_str, metadata_struct.ntypes, npol))
      if tag_exist(metadata_struct, 'n_obs') then n_obs = metadata_struct.n_obs else n_obs = lonarr(npol, nfiles) + 1
      if tag_exist(metadata_struct, 'nside') then begin
        pixel_varname = metadata_struct.pixel_varname
        nside = metadata_struct.nside
      endif
      if tag_exist(metadata_struct, 'no_var') then no_var = metadata_struct.no_var
      undefine, metadata_struct
      
      metadata_struct = {datafile: datafile, weightfile: weightfile, variancefile:variancefile, beamfile:beamfile, $
        cube_varname:cube_varname, weight_varname:weight_varname, variance_varname:variance_varname, $;; beam_varname:beam_varname, $
        frequencies:frequencies, freq_resolution:freq_resolution, time_resolution:time_resolution, $
        n_vis:n_vis, n_vis_freq:n_vis_freq, max_baseline_lambda:max_baseline_lambda, max_theta:max_theta, degpix:degpix, kpix:kpix, kspan:kspan, $
        general_filebase:general_filebase, infile_label:infile_label, type_pol_str:type_pol_str, type_inc:type_inc, n_obs:n_obs, pol_inc:pol_inc, nfiles:nfiles}
        
      if n_elements(nside) gt 0 then metadata_struct = create_struct(metadata_struct, 'pixelfile', pixelfile, 'pixel_varname', pixel_varname, 'nside', nside)
      
      if n_elements(no_var) gt 0 then metadata_struct = create_struct(metadata_struct, 'no_var', 1)
      
      if n_elements(vis_noise) gt 0 then metadata_struct = create_struct(metadata_struct, 'vis_noise', vis_noise)
      
      save, filename = info_file, metadata_struct
      
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
    
    ;; check whether datafiles have polarization identifiers in filename
    pol_exist = stregex(datafile, '[xy][xy]', /boolean, /fold_case)
    if max(pol_exist) gt 0 then begin
      wh_pol_exist = where(pol_exist gt 0, count_pol_exist, ncomplement = count_no_pol)
      if count_no_pol gt 0 then message, 'some datafiles have pol identifiers and some do not'
      
      pols = strlowcase(stregex(datafile, '[xy][xy]', /extract, /fold_case))
      if n_elements(weightfile) gt 0 then begin
        wt_pols = strlowcase(stregex(weightfile, '[xy][xy]', /extract, /fold_case))
        match, pols, wt_pols, suba, subb, count = count_pol_match
        if count_pol_match ne n_elements(pols) then message, 'weightfile polarizations must match datafile polarizations'
      endif
      if n_elements(variancefile) gt 0 then begin
        var_pols = strlowcase(stregex(variancefile, '[xy][xy]', /extract, /fold_case))
        match, pols, var_pols, suba, subb, count = count_pol_match
        if count_pol_match ne n_elements(pols) then message, 'variancefile polarizations must match datafile polarizations'
      endif
      if n_elements(beamfile) gt 0 then begin
        bm_pols = strlowcase(stregex(beamfile, '[xy][xy]', /extract, /fold_case))
        match, pols, bm_pols, suba, subb, count = count_pol_match
        if count_pol_match ne n_elements(pols) then message, 'beamfile polarizations must match datafile polarizations'
      endif
      
      pol_inc = pols[0]
      pol_num = intarr(n_elements(pols))
      for i=0, n_elements(pols)-1 do begin
        wh_pol = where(pol_inc eq pols[i], count_pol)
        if count_pol eq 1 then pol_num[i] = wh_pol[0] else begin
          pol_inc = [pol_inc, pols[i]]
          pol_num[i] = n_elements(pol_inc)-1
        endelse
      endfor
      
      npol = n_elements(pol_inc)
      nfiles = ceil(n_elements(datafile)/float(npol))
      if nfiles gt 2 then message, 'Only 1 or 2 datafiles supported per pol'
      if n_elements(datafile) ne npol*nfiles then message, 'The same number of files must be included per pol'
      temp = strarr(npol, nfiles)
      temp_wt = strarr(npol, nfiles)
      temp_var = strarr(npol, nfiles)
      ;; temp_bm = strarr(npol, nfiles)
      for i=0, npol-1 do begin
        wh_pol =  where(pol_num eq i, count_pol)
        if count_pol ne nfiles then message, 'The same number of files must be included per pol'
        temp[i,*] = datafile[wh_pol]
        if n_elements(weightfile) gt 0 then temp_wt[i,*] = weightfile[wh_pol] else temp_wt = temp
        if n_elements(variancefile) gt 0 then temp_var[i,*] = variancefile[wh_pol] else temp_var = temp
        ;; if n_elements(beamfile) gt 0 then temp_bm[i,*] = beamfile[wh_pol] else temp_bm = temp
        if n_elements(pixelfile) gt 0 then temp_pix[i,*] = pixelfile[wh_pol] else temp_pix = temp
      endfor
      datafile = temporary(temp)
      weightfile = temporary(temp_wt)
      variancefile = temporary(temp_var)
      ;; beamfile = temporary(temp_bm)
      pixelfile = temporary(temp_pix)
    endif else begin
      ;; no pol identifiers
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
      
      nfile_dims = size(datafile, /dimension)
      if n_elements(nfile_dims) gt 1 then datafile = datafile[0,*]
      nfiles = n_elements(datafile)
      if nfiles gt 2 then message, 'Only 1 or 2 datafiles supported'
      if nfiles eq 2 then if datafile[0] eq datafile[1] then begin
        print, 'datafiles are identical'
        datafile = datafile[0]
        nfiles = 1
      endif
      
      if n_elements(weightfile) gt 0 then if n_elements(weightfile) ne nfiles then message, 'weightfile must have the same number of elements as datafile'
      if n_elements(variancefile) gt 0 then  if n_elements(variancefile) ne nfiles then message, 'variancefile must have the same number of elements as datafile'
      if n_elements(beamfile) gt 0 then if n_elements(beamfile) ne nfiles then message, 'beamfile must have the same number of elements as datafile'
      
      temp = strarr(npol, nfiles)
      temp_wt = strarr(npol, nfiles)
      temp_var = strarr(npol, nfiles)
      ;; temp_bm = strarr(npol, nfiles)
      for i=0, nfiles-1 do begin
        temp[*,i] = datafile[i]
        if n_elements(weightfile) gt 0 then temp_wt[*,i] = weightfile[i] else temp_wt = temp
        if n_elements(variancefile) gt 0 then temp_var[*,i] = variancefile[i] else temp_var = temp
        ;; if n_elements(beamfile) gt 0 then temp_bm[*,i] = beamfile[i] else temp_bm = temp
        if n_elements(pixelfile) gt 0 then temp_pix[*,i] = pixelfile[i] else temp_pix = temp
      endfor
      datafile = temp
      weightfile = temporary(temp_wt)
      variancefile = temporary(temp_var)
      ;; beamfile = temporary(temp_bm)
      pixelfile = temporary(temp_pix)
    endelse
    
    if n_elements(savefilebase_in) gt 1 then message, 'only one savefilebase allowed'
    
    if nfiles eq 1 then infile_label = '' else begin
      data_filebase = [cgRootName(datafile[0,0]), cgRootName(datafile[0,1])]
      
      if max(pol_exist) gt 0 then for i=0, 1 do data_filebase[i] = strjoin(strsplit(data_filebase[i], '_?[xy][xy]_?', /regex, /extract), '_')
      
      fileparts_1 = strsplit(data_filebase[0], '_', /extract)
      fileparts_2 = strsplit(data_filebase[1], '_', /extract)
      match_test = strcmp(fileparts_1, fileparts_2)
      wh_diff = where(match_test eq 0, count_diff, complement = wh_same, ncomplement = count_same)
      
      if count_diff gt 0 then infile_label = [strjoin(fileparts_1[wh_diff]), strjoin(fileparts_2[wh_diff])] $
      else infile_label = strarr(2)
    endelse
    
    if n_elements(savefilebase_in) eq 0 then begin
      if nfiles eq 1 then begin
        data_filebase = cgRootName(datafile[0], directory=datafile_dir)
        
        ;; make this general to all pols
        if max(pol_exist) gt 0 then data_filebase= strjoin(strsplit(data_filebase, '_?[xy][xy]_?', /regex, /extract), '_')
        
        if n_elements(save_path) ne 0 then froot = save_path $
        else froot = datafile_dir
        uvf_froot = froot
        
        general_filebase = data_filebase
        
      endif else begin
        if n_elements(save_path) ne 0 then froot = save_path $
        else temp = cgRootName(datafile[0], directory=froot)
        
        if count_diff eq 0 then general_filebase = data_filebase[0] + '_joint' else begin
          if count_same gt 0 then general_filebase = strjoin(fileparts_1[wh_same], '_') + '__' + strjoin(fileparts_1[wh_diff]) $
            + '_' + strjoin(fileparts_2[wh_diff]) + '_joint' $
          else general_filebase = data_filebase[0] + data_filebase[1] + '_joint'
        endelse
        
      endelse
    endif else begin
      if n_elements(save_path) gt 0 then froot = save_path else begin
        savefilebase_in_base = cgRootName(savefilebase_in[0], directory=froot)
        if froot eq '.' then temp = cgRootName(datafile[0], directory=froot)
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
          file_struct_arr = rts_file_setup(info_file, save_path = save_path, weight_savefilebase = weight_savefilebase_in, $
            uvf_savefilebase = uvf_savefilebase_in, savefilebase = savefilebase_in)
            
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
    type_pol_str = strarr(npol, ntypes)
    for i=0, npol-1 do type_pol_str[i, *] = type_inc + '_' + pol_inc[i]
    
    data_varname = pol_inc + '_data'
    weight_varname = pol_inc + '_weights'
    variance_varname = pol_inc + '_variances'
    
    pixel_varname = strarr(npol) + 'pixel_nums'
    
    if max(pol_exist) gt 1 then begin
      if n_elements(weight_savefilebase_in) gt 0 then begin
        if n_elements(weight_savefilebase_in) ne nfiles then $
          message, 'if weight_savefilebase is specified it must have the same number of elements as datafiles'
          
        pols = stregex(weight_savefilebase_in, '[xy][xy]', /extract, /fold_case)
        temp = strarr(npol, nfiles)
        for i=0, npol-1 do begin
          wh_pol =  where(pols eq pol_inc[i], count_pol)
          if count_pol ne nfiles then message, 'The same number of weight_savefilebase must be included per pol'
          temp[i,*] = weight_savefilebase_in[wh_pol]
        endfor
        weight_savefilebase_in = temp
      endif
      if n_elements(uvf_savefilebase_in) gt 0 then begin
        if n_elements(uvf_savefilebase_in) ne nfiles then $
          message, 'if uvf_savefilebase is specified it must have the same number of elements as datafiles'
          
        pols = stregex(uvf_savefilebase_in, '[xy][xy]', /extract, /fold_case)
        temp = strarr(npol, nfiles)
        for i=0, npol-1 do begin
          wh_pol =  where(pols eq pol_inc[i], count_pol)
          if count_pol ne nfiles then message, 'The same number of uvf_savefilebase must be included per pol'
          temp[i,*] = uvf_savefilebase_in[wh_pol]
        endfor
        weight_savefilebase_in = temp
      endif
    endif else begin
      if n_elements(weight_savefilebase_in) gt 0 then begin
        if n_elements(weight_savefilebase_in) ne nfiles then $
          message, 'if weight_savefilebase is specified it must have the same number of elements as datafiles'
        temp = strarr(npol, nfiles)
        for i=0, nfiles-1 do temp[*,i] = weight_savefilebase_in[i]
        weight_savefilebase_in = temp
      endif
      if n_elements(uvf_savefilebase_in) gt 0 then begin
        if n_elements(uvf_savefilebase_in) ne nfiles then $
          message, 'if uvf_savefilebase is specified it must have the same number of elements as datafiles'
        temp = strarr(npol, nfiles)
        for i=0, nfiles-1 do temp[*,i] = uvf_savefilebase_in[i]
        uvf_savefilebase_in = temp
      endif
    endelse
    
    for j=0, nfiles*npol-1 do begin
      pol_i = j mod npol
      file_i = j / npol
      void = getvar_savefile(datafile[pol_i, file_i], names = varnames)
      
      wh = where(strlowcase(varnames) eq strlowcase(data_varname[pol_i]), count)
      if count eq 0 then message, data_varname[pol_i] + ' is not present in datafile (datafile=' + datafile[pol_i, file_i] + ')'
      
      data_size = getvar_savefile(datafile[pol_i, file_i], data_varname[pol_i], /return_size)
      if data_size[0] gt 0 then this_data_dims = data_size[1:data_size[0]]
      
      if j eq 0 then data_dims = this_data_dims else if total(abs(this_data_dims - data_dims)) ne 0 then message, 'data dimensions in files do not match'
      
      wh = where(strlowcase(varnames) eq strlowcase(weight_varname[pol_i]), count)
      
      if count eq 0 then message, weight_varname[pol_i] + ' is not present in datafile (datafile=' + datafile[pol_i, file_i] + ')'
      
      wt_size = getvar_savefile(datafile[pol_i, file_i], weight_varname[pol_i], /return_size)
      if wt_size[0] gt 0 then wt_dims = wt_size[1:wt_size[0]]
      if total(abs(wt_dims - data_dims)) ne 0 then message, 'weight and data dimensions in files do not match'
      
      wh = where(strlowcase(varnames) eq strlowcase(variance_varname[pol_i]), count)
      
      if count eq 0 then begin
        print, variance_varname[pol_i] + ' is not present in datafile (datafile=' + datafile[pol_i, file_i] + '). Using weights^2 instead of variances'
        no_var = 1
      endif else begin
        no_var = 0
        var_size = getvar_savefile(datafile[pol_i, file_i], variance_varname[pol_i], /return_size)
        if var_size[0] gt 0 then var_dims = var_size[1:var_size[0]]
        if total(abs(var_dims - data_dims)) ne 0 then message, 'variance and data dimensions in files do not match'
      endelse
      
      
      wh = where(strlowcase(varnames) eq strlowcase(pixel_varname[pol_i]), count)
      
      if count eq 0 then message, pixel_varname[pol_i] + ' is not present in datafile (datafile=' + datafile[pol_i, file_i] + ')'
      
      pix_size = getvar_savefile(datafile[pol_i, file_i], pixel_varname[pol_i], /return_size)
      if pix_size[0] gt 0 then pix_dims = pix_size[1:pix_size[0]]
      if total(abs(pix_dims - data_dims[0])) ne 0 then message, 'Number of Healpix pixels does not match data dimension'
      
      if j eq 0 then frequencies = getvar_savefile(datafile[pol_i, file_i], 'frequencies') else $
        if total(abs(frequencies-getvar_savefile(datafile[pol_i, file_i], 'frequencies'))) ne 0 then $
        message, 'frequencies do not match between datafiles'
      n_freq = n_elements(frequencies)
      
      if j eq 0 then nside = getvar_savefile(datafile[pol_i, file_i], 'nside') else if nside ne getvar_savefile(datafile[pol_i, file_i], 'nside') then $
        message, 'nside does not match between datafiles'
      if j eq 0 then time_resolution = getvar_savefile(datafile[pol_i, file_i], 'time_resolution') else $
        if time_resolution ne getvar_savefile(datafile[pol_i, file_i], 'time_resolution') then $
        message, 'time_resolution does not match between datafiles'
        
      if j eq 0 then n_vis_freq = fltarr(npol, nfiles, n_freq)
      n_vis_freq[pol_i, file_i, *] = getvar_savefile(datafile[pol_i, file_i], 'n_vis_arr')
      
      if j eq 0 then kpix_arr = getvar_savefile(datafile[pol_i, file_i], 'kpix_arr') else $
        if total(abs(kpix_arr- getvar_savefile(datafile[pol_i, file_i], 'kpix_arr'))) ne 0 then $
        message, 'kpix_arr does not match between datafiles'
      if j eq 0 then kspan_arr = getvar_savefile(datafile[pol_i, file_i], 'kspan_arr') else $
        if total(abs(kspan_arr - getvar_savefile(datafile[pol_i, file_i], 'kspan_arr'))) ne 0 then $
        message, 'kspan_arr does not match between datafiles'
      if j eq 0 then time_integration = getvar_savefile(datafile[pol_i, file_i], 'time_integration') else $
        if time_integration ne getvar_savefile(datafile[pol_i, file_i], 'time_integration') then $
        message, 'time_integration does not match between datafiles'
      if j eq 0 then max_baseline = getvar_savefile(datafile[pol_i, file_i], 'max_baseline') else $
        if max_baseline ne getvar_savefile(datafile[pol_i, file_i], 'max_baseline') then $
        message, 'max_baseline does not match between datafiles'
        
      if j eq 0 then obs_ra = getvar_savefile(datafile[pol_i, file_i], 'obs_ra') else if obs_ra ne getvar_savefile(datafile[pol_i, file_i], 'obs_ra') then $
        message, 'obs_ra does not match between datafiles'
      if j eq 0 then obs_dec = getvar_savefile(datafile[pol_i, file_i], 'obs_dec') else if obs_dec ne getvar_savefile(datafile[pol_i, file_i], 'obs_dec') then $
        message, 'obs_dec does not match between datafiles'
      if j eq 0 then zen_ra = getvar_savefile(datafile[pol_i, file_i], 'zen_ra') else if zen_ra ne getvar_savefile(datafile[pol_i, file_i], 'zen_ra') then $
        message, 'zen_ra does not match between datafiles'
      if j eq 0 then zen_dec = getvar_savefile(datafile[pol_i, file_i], 'zen_dec') else if zen_dec ne getvar_savefile(datafile[pol_i, file_i], 'zen_dec') then $
        message, 'zen_dec does not match between datafiles'
        
      if j eq 0 then pixels = getvar_savefile(datafile[pol_i, file_i], pixel_varname[pol_i]) else begin
        if total(abs(pixels - getvar_savefile(datafile[pol_i, file_i], pixel_varname[pol_i]))) ne 0 then $
          message, 'pixel nums do not match between datafiles, using common set.'
      endelse
    endfor
    
    if max(kspan_arr) - min(kspan_arr) ne 0 then begin
      print, 'kspan varies with frequency, using mean'
      kspan = mean(kspan_arr)
    endif
    
    if max(kpix_arr) - min(kpix_arr) ne 0 then begin
      print, 'kpix varies with frequency, using mean'
      kpix = mean(kpix_arr)
    endif
    
    ;; convert to MHz if in Hz
    if mean(frequencies) gt 1000. then frequencies = frequencies/1e6
    
    npix = n_elements(temporary(pixels))
    max_baseline_lambda = max_baseline * max(frequencies*1e6) / (3e8)
    
    ;; made up for now
    ;;freq_resolution = 8e3;; native resolution of visibilities in Hz
    ;;time_resolution = 0.5;; native resolution of visibilities in s
    
    freq_diffs = ((frequencies - shift(frequencies, 1))[1:*]) * 1e6 ;in Hz
    if total(abs(freq_diffs - freq_diffs[0])) gt 0 then print, 'inconsistent freq channel differences, using the smallest.'
    freq_resolution = min(freq_diffs)
    
    ;; pointing offset from zenith (for calculating horizon distance for wedge line)
    max_theta = angle_difference(obs_dec, obs_ra, zen_dec, zen_ra, /degree)
    
    ;; degpix
    degpix = (2.*!pi/(kspan*2.))*(180./!pi)
    
    ;; n_obs
    n_obs = lonarr(npol, nfiles) + 1
    
    ;; n_vis
     n_vis = total(n_vis_freq, 3)
    
    metadata_struct = {datafile:datafile, weightfile: datafile, variancefile:datafile, $
      cube_varname:data_varname, weight_varname:weight_varname, variance_varname:variance_varname, $
      frequencies:frequencies, freq_resolution:freq_resolution, time_resolution:time_resolution, $
      n_vis:n_vis, n_vis_freq:n_vis_freq, max_baseline_lambda:max_baseline_lambda, max_theta:max_theta, degpix:degpix, kpix:kpix, kspan:kspan, $
      general_filebase:general_filebase, type_pol_str:type_pol_str, infile_label:infile_label, $
      type_inc:type_inc, n_obs:n_obs, pol_inc:pol_inc, nfiles:nfiles}
      
    metadata_struct = create_struct(metadata_struct, 'pixelfile', datafile, 'pixel_varname', pixel_varname, 'nside', nside)
    
    if no_var then metadata_struct = create_struct(metadata_struct, 'no_var', 1)
    
    if n_elements(vis_noise) gt 0 then metadata_struct = create_struct(metadata_struct, 'vis_noise', vis_noise)
    
    if n_elements(info_file) eq 0 then info_file = froot + general_filebase + '_info.idlsave'
    
    save, filename = info_file, metadata_struct, pol_inc
  endif
  
  npol = n_elements(metadata_struct.pol_inc)
  nfiles = metadata_struct.nfiles
  ntypes = n_elements(metadata_struct.type_inc)
  if tag_exist(metadata_struct, 'nside') then healpix = 1 else healpix = 0
  if tag_exist(metadata_struct, 'no_var') then no_var = 1 else no_var = 0
  
  pol_exist = stregex(metadata_struct.datafile, '[xy][xy]', /boolean, /fold_case)
  
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
  
  if n_elements(delta_uv_lambda) ne 0 then uv_tag = '_deluv' + number_formatter(delta_uv_lambda) else uv_tag = ''
  ;;if n_elements(max_uv_lambda) ne 0 then uv_tag = '_maxuv' + number_formatter(max_uv_lambda)
  uvf_tag = uv_tag + fch_tag
  
  if keyword_set(std_power) then power_tag = power_tag + '_stdp' else power_tag = ''
  power_tag = power_tag + sw_tag
  
  if keyword_set(dft_ian) then dft_label = '_ian' else dft_label = ''
  
  wt_file_label = '_weights_' + strlowcase(metadata_struct.pol_inc)
  file_label = strarr(npol, ntypes)
  for i=0, npol-1 do file_label[i,*] = '_' + strlowcase(metadata_struct.type_inc) + '_' + strlowcase(metadata_struct.pol_inc[i])
  savefilebase = metadata_struct.general_filebase + uvf_tag + file_label
  
  if n_elements(uvf_savefilebase_in) lt nfiles then begin
    if n_elements(save_path) ne 0 then uvf_froot = replicate(save_path, npol, nfiles, ntypes) else begin
      uvf_froot = strarr(npol, nfiles, ntypes)
      ;; test datafile path to see if it exists. if not, use froot
      for pol_i=0, npol-1 do begin
        for file_i=0, nfiles-1 do begin
          datafile_path = file_dirname(metadata_struct.datafile[pol_i, file_i], /mark_directory)
          if file_test(datafile_path, /directory) then uvf_froot[pol_i, file_i, *] = datafile_path else uvf_froot[pol_i, file_i, *] = froot
        endfor
      endfor
    endelse
    uvf_savefilebase = strarr(npol, nfiles, ntypes)
    uvf_label = strarr(npol, nfiles, ntypes)
    for pol_i=0, npol-1 do begin
      for file_i=0, nfiles-1 do begin
        if max(pol_exist) eq 1 then uvf_savefilebase[pol_i, file_i,*] = cgRootName(metadata_struct.datafile[pol_i, file_i]) + uvf_tag + '_' + metadata_struct.type_inc + dft_label $
        else uvf_savefilebase[pol_i, file_i,*] = cgRootName(metadata_struct.datafile[pol_i, file_i]) + uvf_tag + file_label[pol_i, *] + dft_label
        uvf_label[pol_i, file_i,*] = metadata_struct.infile_label[file_i] + '_' + file_label[pol_i, *]
      endfor
    endfor
  endif else begin
    if n_elements(save_path) gt 0 then uvf_froot = replicate(save_path, npol, nfiles, ntypes) else begin
      temp = file_dirname(uvf_savefilebase_in, /mark_directory)
      uvf_froot = strarr(npol, nfiles, ntypes)
      for pol_i=0, npol-1 do begin
        for file_i=0, nfiles-1 do begin
          if temp[i] ne '.' then uvf_froot[i,*] = temp[i] else begin
            datafile_path = file_dirname(metadata_struct.datafile[pol_i, file_i], /mark_directory)
            if file_test(datafile_path, /directory) then uvf_froot[pol_i, file_i, *] = datafile_path else uvf_froot[pol_i, file_i, *] = froot
          endelse
        endfor
      endfor
    endelse
    uvf_savefilebase = file_basename(uvf_savefilebase_in) + uvf_tag + file_label + dft_label
  endelse
  
  ;; add sw tag to general_filebase so that plotfiles have uvf_tag in them
  general_filebase = metadata_struct.general_filebase + uvf_tag
  
  
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
    temp = cgrootname(metadata_struct.weightfile[0, 0], directory = wt_froot)
    if file_test(wt_froot, /directory) eq 0 then wt_froot = froot
    
    if n_elements(save_path) gt 0 then wt_froot = save_path else begin
      wt_froot = strarr(npol, nfiles)
      for pol_i=0, npol-1 do begin
        for file_i=0, nfiles-1 do begin
          wtfile_path = file_dirname(metadata_struct.weightfile[pol_i, file_i], /mark_directory)
          if file_test(wtfile_path, /directory) then wt_froot[pol_i, file_i] = wtfile_path else wt_froot[pol_i, file_i] = froot
        endfor
      endfor
    endelse
    
    weight_savefilebase = strarr(npol, nfiles)
    for pol_i=0, npol-1 do begin
      for file_i=0, nfiles-1 do begin
        if max(pol_exist) eq 1 then weight_savefilebase[pol_i, file_i] = cgRootName(metadata_struct.weightfile[pol_i, file_i]) + uvf_tag + '_weights' + dft_label $
        else weight_savefilebase[pol_i, file_i] = cgRootName(metadata_struct.weightfile[pol_i, file_i]) + uvf_tag + wt_file_label[pol_i] + dft_label
      endfor
    endfor
    
  endif else begin
    if n_elements(save_path) gt 0 then wt_froot = save_path else begin
      temp = file_dirname(weight_savefilebase_in, /mark_directory)
      
      wt_froot = strarr(npol, nfiles)
      for pol_i=0, npol-1 do begin
        for file_i=0, nfiles-1 do begin
          if temp[pol_i, file_i] ne '.' then wt_froot[pol_i, file_i] = temp[pol_i, file_i] else begin
            wtfile_path = file_dirname(metadata_struct.weightfile[pol_i, file_i], /mark_directory)
            if file_test(wtfile_path, /directory) then wt_froot[pol_i, file_i] = wtfile_path else wt_froot[pol_i, file_i] = froot
          endelse
        endfor
      endfor
    endelse
    weight_savefilebase = file_basename(weight_savefilebase_in) + uvf_tag
  endelse
  
  uvf_weight_savefile = wt_froot + weight_savefilebase + '_uvf.idlsave'
  uf_weight_savefile = wt_froot + weight_savefilebase + '_uf_plane.idlsave'
  vf_weight_savefile = wt_froot + weight_savefilebase + '_vf_plane.idlsave'
  uv_weight_savefile = wt_froot + weight_savefilebase + '_uv_plane.idlsave'
  
  
  for pol_i=0, npol-1 do begin
    for type_i=0, ntypes-1 do begin
    
      res_uvf_inputfiles = strarr(nfiles, 2)
      res_uvf_varname = strarr(nfiles, 2)
      
      file_struct = {datafile:reform(metadata_struct.datafile[pol_i,*]), weightfile:reform(metadata_struct.weightfile[pol_i,*]), $
        variancefile:reform(metadata_struct.variancefile[pol_i,*]), $;; beamfile:reform(metadata_struct.beamfile[pol_i,*]), $
        datavar:metadata_struct.cube_varname[pol_i], weightvar:metadata_struct.weight_varname[pol_i], $
        variancevar:metadata_struct.variance_varname[pol_i], $;; beamvar:metadata_struct.beam_varname[pol_i], $
        frequencies:metadata_struct.frequencies, freq_resolution:metadata_struct.freq_resolution, time_resolution:metadata_struct.time_resolution, $
        n_obs:reform(metadata_struct.n_obs[pol_i,*]), n_vis:reform(metadata_struct.n_vis[pol_i,*]), n_vis_freq:reform(metadata_struct.n_vis_freq[pol_i,*,*]), $
        max_baseline_lambda:metadata_struct.max_baseline_lambda, max_theta:metadata_struct.max_theta, $
        degpix:metadata_struct.degpix, kpix:metadata_struct.kpix, kspan:metadata_struct.kspan, $
        uf_savefile:reform(uf_savefile[pol_i,*,type_i]), vf_savefile:reform(vf_savefile[pol_i,*,type_i]), uv_savefile:reform(uv_savefile[pol_i,*,type_i]), $
        uf_raw_savefile:reform(uf_raw_savefile[pol_i,*,type_i]), vf_raw_savefile:reform(vf_raw_savefile[pol_i,*,type_i]), $
        uv_raw_savefile:reform(uv_raw_savefile[pol_i,*,type_i]), $
        uf_sum_savefile:uf_sum_savefile[pol_i, type_i], vf_sum_savefile:vf_sum_savefile[pol_i,type_i], $
        uv_sum_savefile:uv_sum_savefile[pol_i,type_i], uf_diff_savefile:uf_diff_savefile[pol_i,type_i], $
        vf_diff_savefile:vf_diff_savefile[pol_i,type_i], uv_diff_savefile:uv_diff_savefile[pol_i,type_i], $
        uf_weight_savefile:reform(uf_weight_savefile[pol_i, *]), vf_weight_savefile:reform(vf_weight_savefile[pol_i, *]), $
        uv_weight_savefile:reform(uv_weight_savefile[pol_i, *]), $
        ;;beam_savefile:reform(beam_savefile[pol_i, *]), $
        kcube_savefile:kcube_savefile[pol_i,type_i], power_savefile:power_savefile[pol_i,type_i], fits_power_savefile:fits_power_savefile[pol_i,type_i],$
        savefile_froot:froot, savefilebase:savefilebase[pol_i,type_i], general_filebase:general_filebase, $
        weight_savefilebase:reform(weight_savefilebase[pol_i, *]), $
        res_uvf_inputfiles:res_uvf_inputfiles, res_uvf_varname:res_uvf_varname, $
        file_label:file_label[pol_i,type_i], uvf_label:reform(uvf_label[pol_i,*,type_i]), wt_file_label:wt_file_label[pol_i], $
        uvf_tag:uvf_tag, power_tag:power_tag, type_pol_str:metadata_struct.type_pol_str[pol_i,type_i], $
        pol_index:pol_i, type_index:type_i, pol:metadata_struct.pol_inc[pol_i], type:metadata_struct.type_inc[type_i], nfiles:nfiles}
        
      if healpix or not keyword_set(uvf_input) then file_struct = create_struct(file_struct, 'uvf_savefile', reform(uvf_savefile[pol_i,*,type_i]), $
        'uvf_weight_savefile', uvf_weight_savefile[pol_i,*])
        
      if healpix then file_struct = create_struct(file_struct, 'pixelfile', metadata_struct.pixelfile[pol_i,*], 'pixelvar', $
        metadata_struct.pixel_varname[pol_i,*], 'nside', metadata_struct.nside)
        
      if no_var then file_struct = create_struct(file_struct, 'no_var', 1)
      
      if n_elements(freq_mask) gt 0 then file_struct = create_struct(file_struct, 'freq_mask', freq_mask)
      
      if tag_exist(metadata_struct, 'vis_noise') gt 0 then file_struct = create_struct(file_struct, 'vis_noise', reform(metadata_struct.vis_noise[pol_i,*,*]))
      
      if pol_i eq 0 and type_i eq 0 then file_struct_arr = replicate(file_struct, ntypes, npol) else file_struct_arr[type_i, pol_i] = file_struct
      
    endfor
  endfor
  
  return, file_struct_arr
  
end
