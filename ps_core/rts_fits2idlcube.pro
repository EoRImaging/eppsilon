function rts_fits2idlcube, datafiles, weightfiles, variancefiles, pol_inc, save_path = save_path, refresh = refresh

  if n_elements(datafiles) eq 0 then message, 'datafiles must be passed in'
  if n_elements(weightfiles) ne n_elements(datafiles) then message, 'different number of weightfiles and datafiles'
  if n_elements(variancefiles) ne n_elements(datafiles) then message, 'different number of variancefiles and datafiles'
  
  ;; check for even/odd files, separate them if present
  even_mask = stregex(datafiles, 'even', /boolean)
  n_even = total(even_mask)
  odd_mask = stregex(datafiles, 'odd', /boolean)
  n_odd = total(odd_mask)
  
  if n_even ne n_odd then message, 'number of even and odd files are not equal'
  if n_even eq 0 then begin
    nfiles = 1
    n_freq = n_elements(datafiles)
    datafile_arr = datafiles
    weightfile_arr = weightfiles
    variancefile_arr = variancefiles
  endif else begin
    nfiles = 2
    n_freq = n_even
    
    datafile_arr = strarr(n_freq, nfiles)
    datafile_arr[*,0] = datafiles[where(even_mask gt 0)]
    datafile_arr[*,1] = datafiles[where(odd_mask gt 0)]
    
    wt_even_mask = stregex(weightfiles, 'even', /boolean)
    wt_n_even = total(wt_even_mask)
    wt_odd_mask = stregex(weightfiles, 'odd', /boolean)
    wt_n_odd = total(wt_odd_mask)
    if wt_n_even ne n_even or wt_n_odd ne n_odd then message, 'number of even and odd weight files do not match datafiles'
    
    weightfile_arr = strarr(n_freq, nfiles)
    weightfile_arr[*,0] = weightfiles[where(even_mask gt 0)]
    weightfile_arr[*,1] = weightfiles[where(odd_mask gt 0)]
    
    var_even_mask = stregex(variancefiles, 'even', /boolean)
    var_n_even = total(var_even_mask)
    var_odd_mask = stregex(variancefiles, 'odd', /boolean)
    var_n_odd = total(var_odd_mask)
    if var_n_even ne n_even or var_n_odd ne n_odd then message, 'number of even and odd variance files do not match datafiles'
    
    variancefile_arr = strarr(n_freq, nfiles)
    variancefile_arr[*,0] = variancefiles[where(even_mask gt 0)]
    variancefile_arr[*,1] = variancefiles[where(odd_mask gt 0)]
  endelse
  
  
  pol_inc = ['xx', 'yy']
  npol = n_elements(pol_inc)
  
  type_inc = ['res']
  ntypes = n_elements(type_inc)
  ncubes = npol * ntypes
  type_pol_str = strarr(ncubes)
  for i=0, npol-1 do type_pol_str[ntypes*i:i*ntypes+(ntypes-1)] = type_inc + '_' + pol_inc[i]
  
  if n_elements(save_path) ne 0 then froot = save_path $
  else froot = file_dirname(file_dirname(datafile_arr[0]), /mark_directory)
  
  idl_cube_savefile=strarr(nfiles)
  for file_i=0, nfiles-1 do begin
  
    infilebase = file_basename(datafile_arr[*,file_i])
    temp2 = strpos(infilebase, '.', /reverse_search)
    
    for i=0, n_freq-1 do begin
      this_fileparts = strsplit(strmid(infilebase[i], 0, temp2[i]), '_', /extract)
      if i eq 0 then fileparts = this_fileparts else begin
        this_test = long(strcmp(fileparts, this_fileparts))
        if i eq 1 then match_test = this_test else match_test = match_test + this_test
      endelse
    endfor
    
    wh_diff = where(match_test lt n_freq-1, count_diff, complement = wh_same, ncomplement = count_same)
    if count_same gt 0 then general_filebase = strjoin(fileparts[wh_same], '_') else stop
    
    idl_cube_savefile[file_i] = froot + general_filebase + '_cube.idlsave'
    
    test_idlcube = file_test(idl_cube_savefile[file_i]) *  (1 - file_test(idl_cube_savefile[file_i], /zero_length))
    if test_idlcube eq 0 or keyword_set(refresh) then begin
    
      frequencies = fltarr(n_freq)
      ;frequencies2 = fltarr(n_freq)
      n_vis_arr = fltarr(n_freq)
      kpix_arr = fltarr(n_freq)
      kspan_arr = fltarr(n_freq)
      for i=0, n_freq-1 do begin
        ;        temp = strsplit(datafile_arr[i, file_i], '_',/extract)
        ;        freq_loc = where(strpos(strlowcase(strsplit(datafile_arr[i, file_i], '_',/extract)), 'mhz') gt -1, count)
        ;        if count gt 1 then stop else freq_loc = freq_loc[0]
        ;        frequencies[i] = float(temp[freq_loc])
      
        data = mrdfits(datafile_arr[i, file_i], 1, hdr, /silent)
        
        this_nside = fxpar(hdr, 'nside')
        if i eq 0 then nside = this_nside else if this_nside ne nside then message, 'nside values do not agree across datafiles'
        this_dims = fxpar(hdr, 'naxis*')
        if i eq 0 then begin
          dims = this_dims
          if n_elements(dims) eq 2 then npix = dims[1]
        endif else if total(abs(this_dims-dims)) ne 0 then message, 'dimensions do not agree across datafiles'
        
        col_types = fxpar(hdr, 'ttype*')
        
        pix_col = where(strpos(strlowcase(col_types), 'pix') gt -1, count)
        if count ne 1 then stop else pix_col = pix_col[0]
        ordering = strtrim(fxpar(hdr, 'tunit' + number_formatter(pix_col+1)))
        
        pixels = data.(pix_col)
        if i eq 0 then pixel_nums = temporary(pixels) $
        else if total(abs(pixel_nums-pixels)) ne 0 then message, 'pixel nums do not agree across datafiles'
        
        xx_col = where(strpos(strlowcase(col_types), 'weighted pp') gt -1, count)
        if count ne 1 then stop else xx_col = xx_col[0]
        
        if i eq 0 then xx_data = fltarr(npix, n_freq)
        
        xx_data[*,i] = data.(xx_col)
        
        yy_col = where(strpos(strlowcase(col_types), 'weighted qq') gt -1, count)
        if count ne 1 then stop else yy_col = yy_col[0]
        
        if i eq 0 then yy_data = fltarr(npix, n_freq)
        yy_data[*,i] =  data.(yy_col)
        
        undefine, data
        
        data2 = mrdfits(datafile_arr[i, file_i], 2, hdr2, /silent)
        col_types2 = fxpar(hdr2, 'ttype*')
        
        ;        freq_col = where(strpos(strlowcase(col_types2), 'freq') gt -1, count)
        ;        if count ne 1 then stop else freq_col = freq_col[0]
        ;        freq_arr = data2.(freq_col)
        ;        if max(abs(freq_arr - freq_arr[0])) gt 0 then print, 'frequencies in file ' + datafile_arr[i, file_i] + ' vary by ' + number_formatter(max(abs(freq_arr - freq_arr[0])))
        ;        frequencies2[i] = freq_arr[0]
        
        time_col = where(strpos(strlowcase(col_types2), 'fracmjd') gt -1, count)
        if count ne 1 then stop else time_col = time_col[0]
        time_arr = data2.(time_col)
        time_diffs = ((time_arr - shift(time_arr, 1))[1:*])*3600*24 ;fractional days to seconds
        time_diffs = round(time_diffs*2.)/2. ;round to nearest 1/2 second
        if total(abs(time_diffs - time_diffs[0])) gt 0 then print, 'inconsistent time resolutions in file ' + datafile_arr[i, file_i] + ' using the smallest value.'
        if i eq 0 then time_resolution = min(time_diffs) else if abs(time_resolution - min(time_diffs)) ne 0 then $
          message, 'time resolutions are not the same between frequency files'
          
        obs_ra_col = where(strpos(strlowcase(col_types2), 'ha_pt') gt -1, count)
        if count ne 1 then stop else obs_ra_col = obs_ra_col[0]
        obs_ra_arr = data2.(obs_ra_col)
        if total(abs(obs_ra_arr - obs_ra_arr[0])) gt 0 then message, 'inconsistent obs_ra in file ' + datafile_arr[i, file_i]
        if i eq 0 then obs_ra = obs_ra_arr[0] else if abs(obs_ra - obs_ra_arr[0]) ne 0 then $
          message, 'obs_ra are not the same between frequency files'
          
        obs_dec_col = where(strpos(strlowcase(col_types2), 'dec_pt') gt -1, count)
        if count ne 1 then stop else obs_dec_col = obs_dec_col[0]
        obs_dec_arr = data2.(obs_dec_col)
        if total(abs(obs_dec_arr - obs_dec_arr[0])) gt 0 then message, 'inconsistent obs_dec in file ' + datafile_arr[i, file_i]
        if i eq 0 then obs_dec = obs_dec_arr[0] else if abs(obs_dec - obs_dec_arr[0]) ne 0 then $
          message, 'obs_dec are not the same between frequency files'
          
        zen_ra_col = where(strpos(strlowcase(col_types2), 'ra_ph') gt -1, count)
        if count ne 1 then stop else zen_ra_col = zen_ra_col[0]
        zen_ra_arr = data2.(zen_ra_col)
        if total(abs(zen_ra_arr - zen_ra_arr[0])) gt 0 then message, 'inconsistent zen_ra in file ' + datafile_arr[i, file_i]
        if i eq 0 then zen_ra = zen_ra_arr[0] else if abs(zen_ra - zen_ra_arr[0]) ne 0 then $
          message, 'zen_ra are not the same between frequency files'
          
        zen_dec_col = where(strpos(strlowcase(col_types2), 'dec_ph') gt -1, count)
        if count ne 1 then stop else zen_dec_col = zen_dec_col[0]
        zen_dec_arr = data2.(zen_dec_col)
        if total(abs(zen_dec_arr - zen_dec_arr[0])) gt 0 then message, 'inconsistent zen_dec in file ' + datafile_arr[i, file_i]
        if i eq 0 then zen_dec = zen_dec_arr[0] else if abs(zen_dec - zen_dec_arr[0]) ne 0 then $
          message, 'zen_dec are not the same between frequency files'
          
        undefine, data2
        
        data3 = mrdfits(datafile_arr[i, file_i], 3, hdr3, /silent)
        frequencies[i] = fxpar(hdr3, 'freq')
        n_vis_arr[i] = fxpar(hdr3, 'n_bsnl')
        kpix_arr[i] = fxpar(hdr3, 'uv_pixsz')
        kspan_arr[i] = fxpar(hdr3, 'uv_span')
        
        if i eq 0 then time_integration = fxpar(hdr3, 'time_int') else if (time_integration-fxpar(hdr3, 'time_int')) ne 0 then $
          message, 'time integrations are not the same between frequency files'
        if i eq 0 then max_baseline = fxpar(hdr3, 'max_bsnl') else if (max_baseline-fxpar(hdr3, 'max_bsnl')) ne 0 then $
          message, 'max baseline are not the same between frequency files'
          
      endfor
      
      for i=0, n_freq-1 do begin
        ;        temp = strsplit(weightfile_arr[i, file_i], '_',/extract)
        ;        freq_loc = where(strpos(strlowcase(strsplit(weightfile_arr[i, file_i], '_',/extract)), 'mhz') gt -1, count)
        ;        if count gt 1 then stop else freq_loc = freq_loc[0]
        ;        if float(temp[freq_loc]) ne frequencies[i] then message, 'data and weights frequencies do not match'
      
        data = mrdfits(weightfile_arr[i, file_i], 1, hdr, /silent)
        
        this_nside = fxpar(hdr, 'nside')
        if this_nside ne nside then message, 'weights nside value does not agree with datafiles'
        this_dims = fxpar(hdr, 'naxis*')
        if total(abs(this_dims-dims)) ne 0 then message, 'weights dimensions do not agree with datafiles'
        
        col_types = fxpar(hdr, 'ttype*')
        
        pix_col = where(strpos(strlowcase(col_types), 'pix') gt -1, count)
        if count ne 1 then stop else pix_col = pix_col[0]
        ordering = strtrim(fxpar(hdr, 'tunit' + number_formatter(pix_col+1)))
        
        pixels = data.(pix_col)
        if total(abs(pixel_nums-pixels)) ne 0 then message, 'weights pixel nums do not agree with datafiles'
        
        xx_col = where(strpos(strlowcase(col_types), 'weighted pp') gt -1, count)
        if count ne 1 then stop else xx_col = xx_col[0]
        
        if i eq 0 then xx_weights = fltarr(npix, n_freq)
        
        xx_weights[*,i] = data.(xx_col)
        
        yy_col = where(strpos(strlowcase(col_types), 'weighted qq') gt -1, count)
        if count ne 1 then stop else yy_col = yy_col[0]
        
        if i eq 0 then yy_weights = fltarr(npix, n_freq)
        yy_weights[*,i] =  data.(yy_col)
        
        undefine, data
        
        data2 = mrdfits(weightfile_arr[i, file_i], 2, hdr2, /silent)
        col_types2 = fxpar(hdr2, 'ttype*')
        
        ;        freq_col = where(strpos(strlowcase(col_types2), 'freq') gt -1, count)
        ;        if count ne 1 then stop else freq_col = freq_col[0]
        ;        freq_arr = data.(freq_col)
        ;        if max(abs(freq_arr - freq_arr[0])) gt 0 then print, 'frequencies in file ' + weightfile_arr[i, file_i] + ' vary by ' + number_formatter(max(abs(freq_arr - freq_arr[0])))
        ;        if abs(frequencies2[i] - freq_arr[0]) gt 1e-3 then message, 'frequencies are not the same in weights and data files.'
        
        time_col = where(strpos(strlowcase(col_types2), 'fracmjd') gt -1, count)
        if count ne 1 then stop else time_col = time_col[0]
        time_arr = data2.(time_col)
        time_diffs = ((time_arr - shift(time_arr, 1))[1:*])*3600*24 ;fractional days to seconds
        time_diffs = round(time_diffs*2.)/2. ;round to nearest 1/2 second
        if total(abs(time_diffs - time_diffs[0])) gt 0 then print, 'inconsistent time resolutions in file ' + weightfile_arr[i, file_i] + ' using the smallest value.'
        if abs(time_resolution - min(time_diffs)) ne 0 then message, 'time resolutions are not the same in weights and data files'
        
        obs_ra_col = where(strpos(strlowcase(col_types2), 'ha_pt') gt -1, count)
        if count ne 1 then stop else obs_ra_col = obs_ra_col[0]
        obs_ra_arr = data2.(obs_ra_col)
        if total(abs(obs_ra_arr - obs_ra_arr[0])) gt 0 then message, 'inconsistent obs_ra in file ' + weightfile_arr[i, file_i]
        if abs(obs_ra - obs_ra_arr[0]) ne 0 then message, 'obs_ra is not the same in weights and data files'
        
        obs_dec_col = where(strpos(strlowcase(col_types2), 'dec_pt') gt -1, count)
        if count ne 1 then stop else obs_dec_col = obs_dec_col[0]
        obs_dec_arr = data2.(obs_dec_col)
        if total(abs(obs_dec_arr - obs_dec_arr[0])) gt 0 then message, 'inconsistent obs_dec in file ' + weightfile_arr[i, file_i]
        if abs(obs_dec - obs_dec_arr[0]) ne 0 then message, 'obs_dec is not the same in weights and data files'
        
        zen_ra_col = where(strpos(strlowcase(col_types2), 'ra_ph') gt -1, count)
        if count ne 1 then stop else zen_ra_col = zen_ra_col[0]
        zen_ra_arr = data2.(zen_ra_col)
        if total(abs(zen_ra_arr - zen_ra_arr[0])) gt 0 then message, 'inconsistent zen_ra in file ' + weightfile_arr[i, file_i]
        if abs(zen_ra - zen_ra_arr[0]) ne 0 then message, 'zen_ra is not the same in weights and data files'
        
        zen_dec_col = where(strpos(strlowcase(col_types2), 'dec_ph') gt -1, count)
        if count ne 1 then stop else zen_dec_col = zen_dec_col[0]
        zen_dec_arr = data2.(zen_dec_col)
        if total(abs(zen_dec_arr - zen_dec_arr[0])) gt 0 then message, 'inconsistent zen_dec in file ' + weightfile_arr[i, file_i]
        if abs(zen_dec - zen_dec_arr[0]) ne 0 then message, 'zen_dec is not the same in weights and data files'
        
        undefine, data2
        
        data3 = mrdfits(weightfile_arr[i, file_i], 3, hdr3, /silent)
        if fxpar(hdr3, 'freq') ne frequencies[i] then begin
          if abs(fxpar(hdr3, 'freq') - frequencies[i])/frequencies[i] lt 0.001 then begin
            print, 'data and weights frequencies different by ' + $
              number_formatter(abs(fxpar(hdr3, 'freq') - frequencies[i])*100/frequencies[i], format = '(e10.2)') + ' %, using data value'
          endif else message, 'data and weights frequencies do not match'
        endif
        
        if fxpar(hdr3, 'n_bsnl') ne n_vis_arr[i] then message, 'data and weights n_vis do not match'
        if fxpar(hdr3, 'uv_pixsz') ne kpix_arr[i] then begin
          if abs(fxpar(hdr3, 'uv_pixsz') - kpix_arr[i])/kpix_arr[i] lt 0.001 then begin
            print, 'data and weights uv pixel sizes different by ' + $
              number_formatter(abs(fxpar(hdr3, 'uv_pixsz') - kpix_arr[i])*100/kpix_arr[i], format = '(e10.2)') + ' %, using data value'
          endif else message, 'data and weights uv pixel sizes do not match'
        endif
        
        if fxpar(hdr3, 'uv_span') ne kspan_arr[i] then begin
          if abs(fxpar(hdr3, 'uv_span') - kspan_arr[i])/kspan_arr[i] lt 0.001 then begin
            print, 'data and weights uv spans different by ' + $
              number_formatter(abs(fxpar(hdr3, 'uv_span') - kspan_arr[i])*100/kspan_arr[i], format = '(e10.2)') + ' %, using data value'
          endif else message, 'data and weights uv spans do not match'
        endif
        if (time_integration-fxpar(hdr3, 'time_int')) ne 0 then message, 'data and weights time integrations do not match'
        if (max_baseline-fxpar(hdr3, 'max_bsnl')) ne 0 then message, 'data and weights max baseline do not match'
        
        
      endfor
      
      for i=0, n_freq-1 do begin
        ;        temp = strsplit(variancefile_arr[i, file_i], '_',/extract)
        ;        freq_loc = where(strpos(strlowcase(strsplit(variancefile_arr[i, file_i], '_',/extract)), 'mhz') gt -1, count)
        ;        if count gt 1 then stop else freq_loc = freq_loc[0]
        ;        if float(temp[freq_loc]) ne frequencies[i] then message, 'data and variances frequencies do not match'
      
        data = mrdfits(variancefile_arr[i, file_i], 1, hdr, /silent)
        
        this_nside = fxpar(hdr, 'nside')
        if this_nside ne nside then message, 'variances nside value does not agree with datafiles'
        this_dims = fxpar(hdr, 'naxis*')
        if total(abs(this_dims-dims)) ne 0 then message, 'variances dimensions do not agree with datafiles'
        
        col_types = fxpar(hdr, 'ttype*')
        
        pix_col = where(strpos(strlowcase(col_types), 'pix') gt -1, count)
        if count ne 1 then stop else pix_col = pix_col[0]
        ordering = strtrim(fxpar(hdr, 'tunit' + number_formatter(pix_col+1)))
        
        pixels = data.(pix_col)
        if total(abs(pixel_nums-pixels)) ne 0 then message, 'variances pixel nums do not agree with datafiles'
        
        xx_col = where(strpos(strlowcase(col_types), 'weighted pp') gt -1, count)
        if count ne 1 then stop else xx_col = xx_col[0]
        
        if i eq 0 then xx_variances = fltarr(npix, n_freq)
        
        xx_variances[*,i] = data.(xx_col)
        
        yy_col = where(strpos(strlowcase(col_types), 'weighted qq') gt -1, count)
        if count ne 1 then stop else yy_col = yy_col[0]
        
        if i eq 0 then yy_variances = fltarr(npix, n_freq)
        yy_variances[*,i] =  data.(yy_col)
        
        undefine, data
        
        data2 = mrdfits(variancefile_arr[i, file_i], 2, hdr2, /silent)
        col_types2 = fxpar(hdr2, 'ttype*')
        
        ;        freq_col = where(strpos(strlowcase(col_types2), 'freq') gt -1, count)
        ;        if count ne 1 then stop else freq_col = freq_col[0]
        ;        freq_arr = data.(freq_col)
        ;        if max(abs(freq_arr - freq_arr[0])) gt 0 then print, 'frequencies in file ' + variancefile_arr[i, file_i] + ' vary by ' + number_formatter(max(abs(freq_arr - freq_arr[0])))
        ;        if abs(frequencies2[i] - freq_arr[0]) gt 1e-3 then message, 'frequencies are not the same in variance and data files.'
        
        time_col = where(strpos(strlowcase(col_types2), 'fracmjd') gt -1, count)
        if count ne 1 then stop else time_col = time_col[0]
        time_arr = data2.(time_col)
        time_diffs = ((time_arr - shift(time_arr, 1))[1:*])*3600*24 ;fractional days to seconds
        time_diffs = round(time_diffs*2.)/2. ;round to nearest 1/2 second
        if total(abs(time_diffs - time_diffs[0])) gt 0 then print, 'inconsistent time resolutions in file ' + variancefile_arr[i, file_i] + ' using the smallest value.'
        if abs(time_resolution - min(time_diffs)) ne 0 then message, 'time resolutions are not the same in variance and data files'
        
        obs_ra_col = where(strpos(strlowcase(col_types2), 'ha_pt') gt -1, count)
        if count ne 1 then stop else obs_ra_col = obs_ra_col[0]
        obs_ra_arr = data2.(obs_ra_col)
        if total(abs(obs_ra_arr - obs_ra_arr[0])) gt 0 then message, 'inconsistent obs_ra in file ' + variancefile_arr[i, file_i]
        if abs(obs_ra - obs_ra_arr[0]) ne 0 then message, 'obs_ra is not the same in variance and data files'
        
        obs_dec_col = where(strpos(strlowcase(col_types2), 'dec_pt') gt -1, count)
        if count ne 1 then stop else obs_dec_col = obs_dec_col[0]
        obs_dec_arr = data2.(obs_dec_col)
        if total(abs(obs_dec_arr - obs_dec_arr[0])) gt 0 then message, 'inconsistent obs_dec in file ' + variancefile_arr[i, file_i]
        if abs(obs_dec - obs_dec_arr[0]) ne 0 then message, 'obs_dec is not the same in variance and data files'
        
        zen_ra_col = where(strpos(strlowcase(col_types2), 'ra_ph') gt -1, count)
        if count ne 1 then stop else zen_ra_col = zen_ra_col[0]
        zen_ra_arr = data2.(zen_ra_col)
        if total(abs(zen_ra_arr - zen_ra_arr[0])) gt 0 then message, 'inconsistent zen_ra in file ' + variancefile_arr[i, file_i]
        if abs(zen_ra - zen_ra_arr[0]) ne 0 then message, 'zen_ra is not the same in variance and data files'
        
        zen_dec_col = where(strpos(strlowcase(col_types2), 'dec_ph') gt -1, count)
        if count ne 1 then stop else zen_dec_col = zen_dec_col[0]
        zen_dec_arr = data2.(zen_dec_col)
        if total(abs(zen_dec_arr - zen_dec_arr[0])) gt 0 then message, 'inconsistent zen_dec in file ' + variancefile_arr[i, file_i]
        if abs(zen_dec - zen_dec_arr[0]) ne 0 then message, 'zen_dec is not the same in variance and data files'
        
        undefine, data2
        
        data3 = mrdfits(variancefile_arr[i, file_i], 3, hdr3, /silent)
        if fxpar(hdr3, 'freq') ne frequencies[i] then begin
          if abs(fxpar(hdr3, 'freq') - frequencies[i])/frequencies[i] lt 0.001 then begin
            print, 'data and variance frequencies different by ' + $
              number_formatter(abs(fxpar(hdr3, 'freq') - frequencies[i])*100/frequencies[i], format = '(e10.2)') + ' %, using data value'
          endif else message, 'data and variance frequencies do not match'
        endif
        
        if fxpar(hdr3, 'n_bsnl') ne n_vis_arr[i] then message, 'data and variance n_vis do not match'

        if fxpar(hdr3, 'uv_pixsz') ne kpix_arr[i] then begin
          if abs(fxpar(hdr3, 'uv_pixsz') - kpix_arr[i])/kpix_arr[i] lt 0.001 then begin
            print, 'data and variance uv pixel sizes different by ' + $
              number_formatter(abs(fxpar(hdr3, 'uv_pixsz') - kpix_arr[i])*100/kpix_arr[i], format = '(e10.2)') + ' %, using data value'
          endif else message, 'data and variance uv pixel sizes do not match'
        endif
        
        if fxpar(hdr3, 'uv_span') ne kspan_arr[i] then begin
          if abs(fxpar(hdr3, 'uv_span') - kspan_arr[i])/kspan_arr[i] lt 0.001 then begin
            print, 'data and variance uv spans different by ' + $
              number_formatter(abs(fxpar(hdr3, 'uv_span') - kspan_arr[i])*100/kspan_arr[i], format = '(e10.2)') + ' %, using data value'
          endif else message, 'data and variance uv spans do not match'
        endif

        if (time_integration-fxpar(hdr3, 'time_int')) ne 0 then message, 'data and variance time integrations do not match'
        if (max_baseline-fxpar(hdr3, 'max_bsnl')) ne 0 then message, 'data and variance max baseline do not match'
        
        
      endfor
      
      save, file = idl_cube_savefile[file_i], nside, frequencies, time_resolution, n_vis_arr, kpix_arr, kspan_arr, time_integration, max_baseline, $
        obs_ra, obs_dec, zen_ra, zen_dec, $
        pixel_nums, xx_data, yy_data, xx_weights, yy_weights, xx_variances, yy_variances
        
    endif
  endfor
  
  
  if nfiles eq 2 then begin
    pixel_nums1 = getvar_savefile(idl_cube_savefile[0], 'pixel_nums')
    pixel_nums2 = getvar_savefile(idl_cube_savefile[1], 'pixel_nums')
    
    if total(abs(pixel_nums2-pixel_nums1)) ne 0 then begin
      print, 'pixel nums do not match between datafiles, using common set.'
      match2, pixel_nums1, pixel_nums2, sub1, sub2
      
      for i=0, 1 do begin
        if i eq 0 then pix_use = where(sub1 gt -1, count_pix, ncomplement = count_drop_pix) $
        else pix_use = where(sub2 gt -1, count_pix, ncomplement = count_drop_pix)
        
        if count_pix eq 0 then message, 'no common pixels between datafiles'
        if count_drop_pix gt 0 then begin
          restore, idl_cube_savefile[i]
          
          pixel_nums = pixel_nums[pix_use]
          xx_data = xx_data[pix_use, *]
          yy_data = yy_data[pix_use, *]
          xx_weights = xx_weights[pix_use, *]
          yy_weights = yy_weights[pix_use, *]
          xx_variances = xx_variances[pix_use, *]
          yy_variances = yy_variances[pix_use, *]
          
          save, file = idl_cube_savefile[file_i], nside, frequencies, time_resolution, n_vis_arr, kpix_arr, kspan_arr, time_integration, max_baseline, $
            obs_ra, obs_dec, zen_ra, zen_dec, $
            pixel_nums, xx_data, yy_data, xx_weights, yy_weights, xx_variances, yy_variances
            
        endif
      endfor
    endif
  endif
  return, idl_cube_savefile
  
end
