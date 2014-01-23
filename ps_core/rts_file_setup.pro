function rts_file_setup, filename, pol_inc, save_path = save_path, $
    weight_savefilebase = weight_savefilebase_in, variance_savefilebase = variance_savefilebase_in, $
    uvf_savefilebase = uvf_savefilebase_in, savefilebase = savefilebase_in, $
    spec_window_type = spec_window_type
    
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
    
      ;; test to see if directories have changed (ie info_file created on different system)
      ;; if so set them to something reasonable
      if n_elements(save_path) gt 0 then begin
        if file_struct_arr.savefile_froot ne save_path then froot = save_path
        if file_dirname(uvf_savefile) ne save_path then uvf_froot = save_path
        if file_dirname(uvf_weight_savefile) ne save_path then wt_froot = save_path
      endif else begin
        if n_elements(savefilebase_in) then begin
          savefilebase_dir = file_dirname(savefilebase_in)
          if file_struct_arr.savefile_froot ne savefilebase_dir then froot = savefilebase_dir
        endif else begin
          info_dir = file_dirname(info_file)
          if file_struct_arr.savefile_froot ne info_dir then froot = info_dir
        endelse
        
        if n_elements(uvf_savefilebase_in) then begin
          uvf_savefilebase_dir = file_dirname(savefilebase_in)
          if file_dirname(uvf_savefile) ne uvf_savefilebase_dir then uvf_froot = uvf_savefilebase_dir
        endif else begin
          info_dir = file_dirname(info_file)
          if file_dirname(uvf_savefile) ne info_dir then uvf_froot = info_dir
        endelse
        
        if n_elements(weight_savefilebase_in) then begin
          wt_savefilebase_dir = file_dirname(savefilebase_in)
          if file_dirname(uvf_weight_savefile) ne wt_savefilebase_dir then wt_froot = wt_savefilebase_dir
        endif else begin
          info_dir = file_dirname(info_file)
          if file_dirname(uvf_weight_savefile) ne info_dir then wt_froot = info_dir
        endelse
        
      endelse
      
      ;; check again (with new directories) to see if datafile(s), uvf, kcube or power files exist
      files_test = file_test([reform(file_struct_arr.datafile, n_elements(file_struct_arr.datafile)), $
        reform(file_struct_arr.uvf_savefile, n_elements(file_struct_arr.uvf_savefile)), $
        reform(file_struct_arr.kcube_savefile, n_elements(file_struct_arr.kcube_savefile)), $
        reform(file_struct_arr.power_savefile, n_elements(file_struct_arr.power_savefile))])
        
    endif
    
    if max(files_test) eq 0 then message, 'Cannot find any datafile(s), uvf, kcube or power files, please specify paths using keywords.'
    
    if n_elements(froot) ne 0 then begin
      file_struct_arr.savefile_froot = froot
      uf_sum_savefile = froot + file_basename(file_struct_arr.uf_sum_savefile)
      vf_sum_savefile = froot + file_basename(file_struct_arr.vf_sum_savefile)
      uv_sum_savefile = froot + file_basename(file_struct_arr.uv_sum_savefile)
      uf_diff_savefile = froot + file_basename(file_struct_arr.uf_diff_savefile)
      vf_diff_savefile = froot + file_basename(file_struct_arr.vf_diff_savefile)
      uv_diff_savefile = froot + file_basename(file_struct_arr.uv_diff_savefile)
      kcube_savefile = froot + file_basename(file_struct_arr.kcube_savefile)
      power_savefile = froot + file_basename(file_struct_arr.power_savefile)
      fits_power_savefile = froot + file_basename(file_struct_arr.fits_power_savefile)
    endif
    
    if n_elements(uvf_froot) ne 0 then begin
      uvf_savefile = uvf_froot + file_basename(file_struct_arr.uvf_savefile)
      uf_savefile = uvf_froot + file_basename(file_struct_arr.uf_savefile)
      vf_savefile = uvf_froot + file_basename(file_struct_arr.vf_savefile)
      uv_savefile = uvf_froot + file_basename(file_struct_arr.uv_savefile)
      uf_raw_savefile = uvf_froot + file_basename(file_struct_arr.uf_raw_savefile)
      vf_raw_savefile = uvf_froot + file_basename(file_struct_arr.vf_raw_savefile)
      uv_raw_savefile = uvf_froot + file_basename(file_struct_arr.uv_raw_savefile)
    endif
    
    if n_elements(wt_froot) ne 0 then uvf_weight_savefile = wt_froot + file_basename(file_struct_arr.uvf_weight_savefile)
    
    return, file_struct_arr
  endelse
  
  nfiles = n_elements(datafile)
  if nfiles eq 1 then datafile = datafile[0]
  
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
  
  type_inc = ['res']
  ntypes = n_elements(type_inc)
  ncubes = npol * ntypes
  type_pol_str = strarr(ncubes)
  for i=0, npol-1 do type_pol_str[ntypes*i:i*ntypes+(ntypes-1)] = type_inc + '_' + pol_inc[i]
  
  wt_file_label = '_weights_' + strlowcase(pol_inc)
  file_label = '_' + strlowcase(type_pol_str)
  
  if n_elements(spec_window_type) ne 0 then begin
    type_list = ['Hann', 'Hamming', 'Blackman', 'Nutall', 'Blackman-Nutall', 'Blackman-Harris']
    sw_tag_list = ['hann', 'ham', 'blm', 'ntl', 'bn', 'bh']
    
    wh_type = where(strlowcase(type_list) eq strlowcase(spec_window_type), count_type)
    if count_type eq 0 then message, 'Spectral window type not recognized.' else begin
      spec_window_type = type_list[wh_type[0]]
      sw_tag = '_' + sw_tag_list[wh_type[0]]
    endelse
  endif else sw_tag = ''
  
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
    uvf_savefilebase = strarr(nfiles, ncubes) + (file_basename(uvf_savefilebase_in) + fch_tag + file_label)
  endif
  
  
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
  
  uvf_weight_savefile = wt_froot + weight_savefilebase + '_uvf.idlsave'
  
  for i=0, nfiles-1 do begin
    if i eq 0 then frequencies = getvar_savefile(datafile[i], 'frequencies') else $
      if total(abs(frequencies-getvar_savefile(datafile[i], 'frequencies'))) ne 0 then $
      message, 'frequencies do not match between datafiles'
    if i eq 0 then nside = getvar_savefile(datafile[i], 'nside') else if nside ne getvar_savefile(datafile[i], 'nside') then $
      message, 'nside does not match between datafiles'
    if i eq 0 then time_resolution = getvar_savefile(datafile[i], 'time_resolution') else $
      if time_resolution ne getvar_savefile(datafile[i], 'time_resolution') then $
      message, 'time_resolution does not match between datafiles'
    ;    if i eq 0 then obs_ra = getvar_savefile(datafile[i], 'obs_ra') else if obs_ra ne getvar_savefile(datafile[i], 'obs_ra') then $
    ;      message, 'obs_ra does not match between datafiles'
    ;    if i eq 0 then obs_dec = getvar_savefile(datafile[i], 'obs_dec') else if obs_dec ne getvar_savefile(datafile[i], 'obs_dec') then $
    ;      message, 'obs_dec does not match between datafiles'
    ;    if i eq 0 then zen_ra = getvar_savefile(datafile[i], 'zen_ra') else if zen_ra ne getvar_savefile(datafile[i], 'zen_ra') then $
    ;      message, 'zen_ra does not match between datafiles'
    ;    if i eq 0 then zen_dec = getvar_savefile(datafile[i], 'zen_dec') else if zen_dec ne getvar_savefile(datafile[i], 'zen_dec') then $
    ;      message, 'zen_dec does not match between datafiles'
    if i eq 0 then pixels = getvar_savefile(datafile[i], 'pixel_nums') else begin
      if total(abs(pixels - getvar_savefile(datafile[i], 'pixel_nums'))) ne 0 then $
        message, 'pixel nums do not match between datafiles, using common set.'
    endelse
  endfor
  
  n_freq = n_elements(frequencies)
  freq_resolution = frequencies[1] - frequencies[0]
  npix = n_elements(temporary(pixels))
  
  data_varname = pol_inc + '_data'
  weight_varname = pol_inc + '_weights'
  variance_varname = pol_inc + '_variances'
  pixel_varname = strarr(nfiles) + 'pixel_nums'
  
  ;; these are totally made up for now
  ;time_resolution = 16
  n_vis = fltarr(nfiles)+(112*111)*n_freq*(6.)
  max_baseline_lambda = 1500 * max(frequencies*1e6) / (3e8)
  kspan = 300.
  kpix = 3.46697 * 1.41555e+08 / (3e8)
  
  ;; pointing offset from zenith (for calculating horizon distance for wedge line)
  max_theta = 0
  
  for i=0, ncubes-1 do begin
    pol_index = i / ntypes
    type_index = i mod ntypes
    
    file_struct = {datafile: datafile, weightfile: datafile, variancefile:datafile, pixelfile:datafile, $
      datavar:data_varname[pol_index], variancevar:variance_varname[pol_index], weightvar:weight_varname[pol_index], $
      pixelvar:pixel_varname, frequencies:frequencies, freq_resolution:freq_resolution, $
      time_resolution:time_resolution, n_vis:n_vis, max_baseline_lambda:max_baseline_lambda, max_theta:max_theta, $
      kpix:kpix, kspan:kspan, nside:nside, $
      uvf_savefile:uvf_savefile[*,i], uvf_weight_savefile:uvf_weight_savefile[*,pol_index], $
      uf_savefile:uf_savefile[*,i], vf_savefile:vf_savefile[*,i], uv_savefile:uv_savefile[*,i], $
      uf_raw_savefile:uf_raw_savefile[*,i], vf_raw_savefile:vf_raw_savefile[*,i], $
      uv_raw_savefile:uv_raw_savefile[*,i], $
      uf_sum_savefile:uf_sum_savefile[i], vf_sum_savefile:vf_sum_savefile[i], $
      uv_sum_savefile:uv_sum_savefile[i], uf_diff_savefile:uf_diff_savefile[i], $
      vf_diff_savefile:vf_diff_savefile[i], uv_diff_savefile:uv_diff_savefile[i], $
      kcube_savefile:kcube_savefile[i], power_savefile:power_savefile[i], fits_power_savefile:fits_power_savefile[i],$
      savefile_froot:froot, savefilebase:savefilebase[i], general_filebase:general_filebase, $
      weight_savefilebase:weight_savefilebase[*,pol_index], $
      file_label:file_label[i], wt_file_label:wt_file_label[pol_index]}
      
    if i eq 0 then file_struct_arr = replicate(file_struct, ncubes) else file_struct_arr[i] = file_struct
  endfor
  
  save, filename = info_file, file_struct_arr, pol_inc
  
  return, file_struct_arr
  
end
