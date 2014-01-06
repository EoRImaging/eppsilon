function rts_file_setup, datafile, pol_inc, save_path = save_path, $
    weight_savefilebase = weight_savefilebase_in, variance_savefilebase = variance_savefilebase_in, $
    uvf_savefilebase = uvf_savefilebase_in, savefilebase = savefilebase_in, $
    spec_window_type = spec_window_type
    
  nfiles = n_elements(datafile)
  if nfiles eq 1 then datafile = datafile[0]
  
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
  
  if n_elements(freq_ch_range) ne 0 then begin
    if min(freq_ch_range) lt 0 or max(freq_ch_range) - min(freq_ch_range) lt 3 then message, 'invalid freq_ch_range'
    fch_tag = '_ch' + number_formatter(min(freq_ch_range)) + '-' + number_formatter(max(freq_ch_range))
  endif else fch_tag = ''
  
  if n_elements(savefilebase_in) eq 0 or n_elements(uvf_savefilebase_in) lt nfiles $
    or n_elements(weight_savefilebase_in) lt nfiles then begin
    if nfiles eq 1 then begin
      if n_elements(save_path) ne 0 then froot = save_path $
      else froot = file_dirname(datafile, /mark_directory)
      uvf_froot = froot
      wt_froot = froot
      infilebase = file_basename(datafile)
      temp2 = strpos(infilebase, '_imagecube.idlsave', /reverse_search)
      
      general_filebase = strmid(infilebase, 0, temp2) + fch_tag
      if n_elements(savefilebase_in) eq 0 then savefilebase = general_filebase + file_label + sw_tag
      
      ;; if we're only dealing with one file and uvf_savefilebase isn't specified then use same base for uvf files
      if n_elements(uvf_savefilebase_in) eq 0 then uvf_savefilebase = strarr(nfiles, ncubes) + (general_filebase + file_label)
      if n_elements(weight_savefilebase_in) eq 0 then weight_savefilebase = strarr(nfiles, ncubes) + (general_filebase + wt_file_label)

    endif else begin
      if n_elements(save_path) ne 0 then froot = save_path $
      else froot = file_dirname(datafile[0], /mark_directory)
      infilebase = file_basename(datafile)
      temp2 = strpos(infilebase, '.', /reverse_search)
      
      if n_elements(savefilebase_in) eq 0 then begin
        fileparts_1 = strsplit(strmid(infilebase[0], 0, temp2[0]), '_', /extract)
        fileparts_2 = strsplit(strmid(infilebase[1], 0, temp2[1]), '_', /extract)
        match_test = strcmp(fileparts_1, fileparts_2)
        wh_diff = where(match_test eq 0, count_diff, complement = wh_same, ncomplement = count_same)
        if count_diff eq 0 then general_filebase = strmid(infilebase[0], 0, temp2[0]) + '_joint' else begin
          if count_same gt 0 then general_filebase = strjoin(fileparts_1[wh_same], '_') + '__' + strjoin(fileparts_1[wh_diff]) $
            + '_' + strjoin(fileparts_2[wh_diff]) + '_joint' $
          else general_filebase = infilebase[0] + infilebase[1] + '_joint'
        endelse
        
        general_filebase = general_filebase + fch_tag
        savefilebase = general_filebase + file_label + sw_tag
      endif
      
      if n_elements(uvf_savefilebase_in) eq 0 then begin
        if n_elements(save_path) ne 0 then uvf_froot = replicate(save_path, nfiles, ncubes) else begin
          uvf_froot = strarr(nfiles, ncubes)
          for i=0, nfiles-1 do uvf_froot[i, *] = file_dirname(datafile[i], /mark_directory)
        endelse

        uvf_savefilebase = strarr(nfiles, ncubes)
        for i=0, nfiles-1 do uvf_savefilebase[i, *] = strmid(infilebase[i], 0, temp2[i]) + fch_tag + file_label
      endif
      
      if n_elements(weight_savefilebase_in) eq 0 then begin
        if n_elements(save_path) ne 0 then wt_froot = replicate(save_path, nfiles, npol) else begin
          wt_froot = strarr(nfiles, npol)
          for i=0, nfiles-1 do wt_froot[i, *] = file_dirname(datafile[i], /mark_directory)
        endelse
        weight_savefilebase = strarr(nfiles, npol)
        for i=0, nfiles-1 do weight_savefilebase[i, *] = strmid(infilebase[i], 0, temp2[i]) + fch_tag + wt_file_label
      endif
    endelse
  endif
  
  
  if n_elements(savefilebase_in) eq 1 then begin
    if n_elements(save_path) gt 0 then froot = save_path else begin
      temp = file_dirname(savefilebase_in, /mark_directory)
      if temp ne '.' then froot = temp else froot = file_dirname(datafile[0], /mark_directory)
    endelse
    savefilebase = file_basename(savefilebase_in) + fch_tag + file_label + sw_tag
    general_filebase = file_basename(savefilebase_in) + fch_tag
  endif
  
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
  
  
  if n_elements(weight_savefilebase_in) eq nfiles then begin
    if n_elements(save_path) gt 0 then wt_froot = replicate(save_path, nfiles, npol) else begin
      temp = file_dirname(uvf_savefilebase_in, /mark_directory)
      wt_froot = strarr(nfiles, npol)
      for i=0, nfiles-1 do if temp[i] ne '.' then wt_froot[i,*] = temp[i] $
      else wt_froot[i,*] = file_dirname(datafile[i], /mark_directory)
    endelse
    weight_savefilebase = strarr(nfiles, ncubes) + (file_basename(weight_savefilebase_in) + fch_tag + wt_file_label)
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
  
  
  return, file_struct_arr
  
end
