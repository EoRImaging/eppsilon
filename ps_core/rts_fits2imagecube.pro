function rts_fits2imagecube, datafiles, weightfiles, variancefiles, pol_inc, save_path = save_path

  file_dims = size(datafiles, /dimension)
  if n_elements(file_dims) gt 1 then nfiles = file_dims[1] else nfiles = 1
  n_freq = file_dims[0]
  
  if n_freq eq 0 then message, 'no files in specified directory: ' + data_dir
  if n_elements(weightfiles) ne n_freq*nfiles then message, 'different number of weightfiles and datafiles'
  if n_elements(variancefiles) ne n_freq*nfiles then message, 'different number of variancefiles and datafiles'
  
  pol_inc = ['xx', 'yy']
  npol = n_elements(pol_inc)
  
  type_inc = ['res']
  ntypes = n_elements(type_inc)
  ncubes = npol * ntypes
  type_pol_str = strarr(ncubes)
  for i=0, npol-1 do type_pol_str[ntypes*i:i*ntypes+(ntypes-1)] = type_inc + '_' + pol_inc[i]
  
  if n_elements(save_path) ne 0 then froot = save_path $
  else froot = file_dirname(file_dirname(datafiles[0]), /mark_directory)
  
  image_cube_savefile=strarr(nfiles)
  for file_i=0, nfiles-1 do begin
  
    infilebase = file_basename(datafiles[*,file_i])
    temp2 = strpos(infilebase, '.', /reverse_search)
    
    for i=0, n_freq-1 do begin
      this_fileparts = strsplit(strmid(infilebase[i], 0, temp2[i]), '_', /extract)
      if i eq 0 then fileparts = this_fileparts else begin
        this_test = strcmp(fileparts, this_fileparts)
        if i eq 1 then match_test = this_test else match_test = match_test + this_test
      endelse
    endfor
    
    wh_diff = where(match_test lt n_freq-1, count_diff, complement = wh_same, ncomplement = count_same)
    if count_same gt 0 then general_filebase = strjoin(fileparts[wh_same], '_') else stop
    
    image_cube_savefile[file_i] = froot + general_filebase + '_imagecube.idlsave'
    
    test_imagecube = file_test(image_cube_savefile[file_i]) *  (1 - file_test(image_cube_savefile[file_i], /zero_length))
    if test_imagecube eq 0 then begin
    
      frequencies = fltarr(n_freq)
      frequencies2 = fltarr(n_freq)
      for i=0, n_freq-1 do begin
        temp = strsplit(datafiles[i, file_i], '_',/extract)
        freq_loc = where(strpos(strlowcase(strsplit(datafiles[i, file_i], '_',/extract)), 'mhz') gt -1, count)
        if count gt 1 then stop else freq_loc = freq_loc[0]
        frequencies[i] = float(temp[freq_loc])
        
        data = mrdfits(datafiles[i, file_i], 1, hdr, /silent)
        
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
        
        data = mrdfits(datafiles[i, file_i], 2, hdr2, /silent)
        col_types2 = fxpar(hdr2, 'ttype*')
        
        freq_col = where(strpos(strlowcase(col_types2), 'freq') gt -1, count)
        if count ne 1 then stop else freq_col = freq_col[0]
        freq_arr = data.(freq_col)
        if total(abs(freq_arr - freq_arr[0])) gt 1e-3 then message, 'inconsistent frequencies in file ' + datafiles[i, file_i]
        frequencies2[i] = freq_arr[0]
        
        time_col = where(strpos(strlowcase(col_types2), 'fracmjd') gt -1, count)
        if count ne 1 then stop else time_col = time_col[0]
        time_arr = data.(time_col)
        time_diffs = ((time_arr - shift(time_arr, 1))[1:*])*3600*24 ;fractional days to seconds
        time_diffs = round(time_diffs*2.)/2. ;round to nearest 1/2 second
        if total(abs(time_diffs - time_diffs[0])) gt 0 then message, 'inconsistent time resolutions in file ' + datafiles[i, file_i]
        if i eq 0 then time_resolution = time_diffs[0] else if abs(time_resolution - time_diffs[0]) ne 0 then $
          message, 'time resolutions are not the same between frequency files'
          
          
;        obs_ra_col = where(strpos(strlowcase(col_types2), 'ha_pt') gt -1, count)
;        if count ne 1 then stop else obs_ra_col = obs_ra_col[0]
;        obs_ra_arr = data.(obs_ra_col)
;        if total(abs(obs_ra_arr - obs_ra_arr[0])) gt 0 then message, 'inconsistent obs_ra in file ' + datafiles[i, file_i]
;        if i eq 0 then obs_ra = obs_ra_arr[0] else if abs(obs_ra - obs_ra_arr[0]) ne 0 then $
;          message, 'obs_ra are not the same between frequency files'
;          
;        obs_dec_col = where(strpos(strlowcase(col_types2), 'dec_pt') gt -1, count)
;        if count ne 1 then stop else obs_dec_col = obs_dec_col[0]
;        obs_dec_arr = data.(obs_dec_col)
;        if total(abs(obs_dec_arr - obs_dec_arr[0])) gt 0 then message, 'inconsistent obs_dec in file ' + datafiles[i, file_i]
;        if i eq 0 then obs_dec = obs_dec_arr[0] else if abs(obs_dec - obs_dec_arr[0]) ne 0 then $
;          message, 'obs_dec are not the same between frequency files'
;          
;        zen_ra_col = where(strpos(strlowcase(col_types2), 'ha_ph') gt -1, count)
;        if count ne 1 then stop else zen_ra_col = zen_ra_col[0]
;        zen_ra_arr = data.(zen_ra_col)
;        if total(abs(zen_ra_arr - zen_ra_arr[0])) gt 0 then message, 'inconsistent zen_ra in file ' + datafiles[i, file_i]
;        if i eq 0 then zen_ra = zen_ra_arr[0] else if abs(zen_ra - zen_ra_arr[0]) ne 0 then $
;          message, 'zen_ra are not the same between frequency files'
;          
;        zen_dec_col = where(strpos(strlowcase(col_types2), 'dec_ph') gt -1, count)
;        if count ne 1 then stop else zen_dec_col = zen_dec_col[0]
;        zen_dec_arr = data.(zen_dec_col)
;        if total(abs(zen_dec_arr - zen_dec_arr[0])) gt 0 then message, 'inconsistent zen_dec in file ' + datafiles[i, file_i]
;        if i eq 0 then zen_dec = zen_dec_arr[0] else if abs(zen_dec - zen_dec_arr[0]) ne 0 then $
;          message, 'zen_dec are not the same between frequency files'
          
      endfor
      
      for i=0, n_freq-1 do begin
        temp = strsplit(weightfiles[i, file_i], '_',/extract)
        freq_loc = where(strpos(strlowcase(strsplit(weightfiles[i, file_i], '_',/extract)), 'mhz') gt -1, count)
        if count gt 1 then stop else freq_loc = freq_loc[0]
        if float(temp[freq_loc]) ne frequencies[i] then message, 'data and weights frequencies do not match'
        
        data = mrdfits(weightfiles[i, file_i], 1, hdr, /silent)
        
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
        
        data = mrdfits(weightfiles[i, file_i], 2, hdr2, /silent)
        col_types2 = fxpar(hdr2, 'ttype*')
        
        freq_col = where(strpos(strlowcase(col_types2), 'freq') gt -1, count)
        if count ne 1 then stop else freq_col = freq_col[0]
        freq_arr = data.(freq_col)
        if total(abs(freq_arr - freq_arr[0])) gt 1e-3 then message, 'inconsistent frequencies in file ' + weightfiles[i, file_i]
        if abs(frequencies2[i] - freq_arr[0]) gt 1e-3 then message, 'frequencies are not the same in weights and data files.'
        
        time_col = where(strpos(strlowcase(col_types2), 'fracmjd') gt -1, count)
        if count ne 1 then stop else time_col = time_col[0]
        time_arr = data.(time_col)
        time_diffs = ((time_arr - shift(time_arr, 1))[1:*])*3600*24 ;fractional days to seconds
        time_diffs = round(time_diffs*2.)/2. ;round to nearest 1/2 second
        if total(abs(time_diffs - time_diffs[0])) gt 0 then message, 'inconsistent time resolutions in file ' + weightfiles[i, file_i]
        if abs(time_resolution - time_diffs[0]) ne 0 then message, 'time resolutions are not the same in weights and data files'
        
        
;        obs_ra_col = where(strpos(strlowcase(col_types2), 'ha_pt') gt -1, count)
;        if count ne 1 then stop else obs_ra_col = obs_ra_col[0]
;        obs_ra_arr = data.(obs_ra_col)
;        if total(abs(obs_ra_arr - obs_ra_arr[0])) gt 0 then message, 'inconsistent obs_ra in file ' + weightfiles[i, file_i]
;        if abs(obs_ra - obs_ra_arr[0]) ne 0 then message, 'obs_ra is not the same in weights and data files'
;        
;        obs_dec_col = where(strpos(strlowcase(col_types2), 'dec_pt') gt -1, count)
;        if count ne 1 then stop else obs_dec_col = obs_dec_col[0]
;        obs_dec_arr = data.(obs_dec_col)
;        if total(abs(obs_dec_arr - obs_dec_arr[0])) gt 0 then message, 'inconsistent obs_dec in file ' + weightfiles[i, file_i]
;        if abs(obs_dec - obs_dec_arr[0]) ne 0 then message, 'obs_dec is not the same in weights and data files'
;        
;        zen_ra_col = where(strpos(strlowcase(col_types2), 'ha_ph') gt -1, count)
;        if count ne 1 then stop else zen_ra_col = zen_ra_col[0]
;        zen_ra_arr = data.(zen_ra_col)
;        if total(abs(zen_ra_arr - zen_ra_arr[0])) gt 0 then message, 'inconsistent zen_ra in file ' + weightfiles[i, file_i]
;        if abs(zen_ra - zen_ra_arr[0]) ne 0 then message, 'zen_ra is not the same in weights and data files'
;        
;        zen_dec_col = where(strpos(strlowcase(col_types2), 'dec_ph') gt -1, count)
;        if count ne 1 then stop else zen_dec_col = zen_dec_col[0]
;        zen_dec_arr = data.(zen_dec_col)
;        if total(abs(zen_dec_arr - zen_dec_arr[0])) gt 0 then message, 'inconsistent zen_dec in file ' + weightfiles[i, file_i]
;        if abs(zen_dec - zen_dec_arr[0]) ne 0 then message, 'zen_dec is not the same in weights and data files'
        
      endfor
      
      for i=0, n_freq-1 do begin
        temp = strsplit(variancefiles[i, file_i], '_',/extract)
        freq_loc = where(strpos(strlowcase(strsplit(variancefiles[i, file_i], '_',/extract)), 'mhz') gt -1, count)
        if count gt 1 then stop else freq_loc = freq_loc[0]
        if float(temp[freq_loc]) ne frequencies[i] then message, 'data and variances frequencies do not match'
        
        data = mrdfits(variancefiles[i, file_i], 1, hdr, /silent)
        
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
        
        data = mrdfits(variancefiles[i, file_i], 2, hdr2, /silent)
        col_types2 = fxpar(hdr2, 'ttype*')
        
        freq_col = where(strpos(strlowcase(col_types2), 'freq') gt -1, count)
        if count ne 1 then stop else freq_col = freq_col[0]
        freq_arr = data.(freq_col)
        if total(abs(freq_arr - freq_arr[0])) gt 1e-3 then message, 'inconsistent frequencies in file ' + variancefiles[i, file_i]
        if abs(frequencies2[i] - freq_arr[0]) gt 1e-3 then message, 'frequencies are not the same in variance and data files.'
        
        time_col = where(strpos(strlowcase(col_types2), 'fracmjd') gt -1, count)
        if count ne 1 then stop else time_col = time_col[0]
        time_arr = data.(time_col)
        time_diffs = ((time_arr - shift(time_arr, 1))[1:*])*3600*24 ;fractional days to seconds
        time_diffs = round(time_diffs*2.)/2. ;round to nearest 1/2 second
        if total(abs(time_diffs - time_diffs[0])) gt 0 then message, 'inconsistent time resolutions in file ' + variancefiles[i, file_i]
        if abs(time_resolution - time_diffs[0]) ne 0 then message, 'time resolutions are not the same in variance and data files'
        
        
;        obs_ra_col = where(strpos(strlowcase(col_types2), 'ha_pt') gt -1, count)
;        if count ne 1 then stop else obs_ra_col = obs_ra_col[0]
;        obs_ra_arr = data.(obs_ra_col)
;        if total(abs(obs_ra_arr - obs_ra_arr[0])) gt 0 then message, 'inconsistent obs_ra in file ' + variancefiles[i, file_i]
;        if abs(obs_ra - obs_ra_arr[0]) ne 0 then message, 'obs_ra is not the same in variance and data files'
;        
;        obs_dec_col = where(strpos(strlowcase(col_types2), 'dec_pt') gt -1, count)
;        if count ne 1 then stop else obs_dec_col = obs_dec_col[0]
;        obs_dec_arr = data.(obs_dec_col)
;        if total(abs(obs_dec_arr - obs_dec_arr[0])) gt 0 then message, 'inconsistent obs_dec in file ' + variancefiles[i, file_i]
;        if abs(obs_dec - obs_dec_arr[0]) ne 0 then message, 'obs_dec is not the same in variance and data files'
;        
;        zen_ra_col = where(strpos(strlowcase(col_types2), 'ha_ph') gt -1, count)
;        if count ne 1 then stop else zen_ra_col = zen_ra_col[0]
;        zen_ra_arr = data.(zen_ra_col)
;        if total(abs(zen_ra_arr - zen_ra_arr[0])) gt 0 then message, 'inconsistent zen_ra in file ' + variancefiles[i, file_i]
;        if abs(zen_ra - zen_ra_arr[0]) ne 0 then message, 'zen_ra is not the same in variance and data files'
;        
;        zen_dec_col = where(strpos(strlowcase(col_types2), 'dec_ph') gt -1, count)
;        if count ne 1 then stop else zen_dec_col = zen_dec_col[0]
;        zen_dec_arr = data.(zen_dec_col)
;        if total(abs(zen_dec_arr - zen_dec_arr[0])) gt 0 then message, 'inconsistent zen_dec in file ' + variancefiles[i, file_i]
;        if abs(zen_dec - zen_dec_arr[0]) ne 0 then message, 'zen_dec is not the same in variance and data files'
        
      endfor
      
      save, file = image_cube_savefile[file_i], nside, frequencies, frequencies2, time_resolution, $;obs_ra, obs_dec, zen_ra, zen_dec, $
        pixel_nums, xx_data, yy_data, xx_weights, yy_weights, xx_variances, yy_variances
        
    endif
  endfor
  
  
  if nfiles eq 2 then begin  
    pixel_nums1 = getvar_savefile(image_cube_savefile[0], 'pixel_nums')
    pixel_nums2 = getvar_savefile(image_cube_savefile[1], 'pixel_nums')
    
    if total(abs(pixel_nums2-pixel_nums1)) ne 0 then begin
      print, 'pixel nums do not match between datafiles, using common set.'
      match2, pixel_nums1, pixel_nums2, sub1, sub2
      
      for i=0, 1 do begin
        if i eq 0 then pix_use = where(sub1 gt -1, count_pix, ncomplement = count_drop_pix) $
        else pix_use = where(sub2 gt -1, count_pix, ncomplement = count_drop_pix)
        
        if count_pix eq 0 then message, 'no common pixels between datafiles'
        if count_drop_pix gt 0 then begin
          restore, image_cube_savefile[i]
          
          pixel_nums = pixel_nums[pix_use]
          xx_data = xx_data[pix_use, *]
          yy_data = yy_data[pix_use, *]
          xx_weights = xx_weights[pix_use, *]
          yy_weights = yy_weights[pix_use, *]
          xx_variances = xx_variances[pix_use, *]
          yy_variances = yy_variances[pix_use, *]
          
          save, file = image_cube_savefile[i], nside, frequencies, pixel_nums, xx_data, yy_data, xx_weights, yy_weights, $
            xx_variances, yy_variances
            
        endif
      endfor
    endif
  endif
  return, image_cube_savefile
  
end
