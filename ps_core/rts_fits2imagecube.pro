function rts_fits2imagecube, datafiles, weightfiles, variancefiles, pol_inc, save_path = save_path
                  
  n_freq = n_elements(datafiles)
  if n_freq eq 0 then message, 'no files in specified directory: ' + data_dir
  if n_elements(weightfiles) ne n_freq then message, 'different number of weightfiles and datafiles'
  if n_elements(variancefiles) ne n_freq then message, 'different number of variancefiles and datafiles'
  
  pol_inc = ['xx', 'yy']
  npol = n_elements(pol_inc)

  type_inc = ['res']
  ntypes = n_elements(type_inc)
  ncubes = npol * ntypes
  type_pol_str = strarr(ncubes)
  for i=0, npol-1 do type_pol_str[ntypes*i:i*ntypes+(ntypes-1)] = type_inc + '_' + pol_inc[i]

  if n_elements(save_path) ne 0 then froot = save_path $
  else froot = file_dirname(file_dirname(datafiles[0]), /mark_directory)
  
  infilebase = file_basename(datafiles)
  temp2 = strpos(infilebase, '.', /reverse_search)
  
  for i=0, n_freq-1 do begin
     this_fileparts = strsplit(strmid(infilebase[i], 0, temp2[i]), '_', /extract)
     if i eq 0 then fileparts = this_fileparts $
     else begin
        this_test = strcmp(fileparts, this_fileparts)
        if i eq 1 then match_test = this_test else match_test = match_test + this_test  
     endelse
  endfor 
  
  wh_diff = where(match_test lt n_freq-1, count_diff, complement = wh_same, ncomplement = count_same)
  if count_same gt 0 then general_filebase = strjoin(fileparts[wh_same], '_') else stop

  image_cube_savefile = froot + general_filebase + '_imagecube.idlsave'

  test_imagecube = file_test(image_cube_savefile) *  (1 - file_test(image_cube_savefile, /zero_length))
  if test_imagecube eq 0 then begin

     frequencies = fltarr(n_freq)
     for i=0, n_freq-1 do begin
        temp = strsplit(datafiles[i], '_',/extract)
        freq_loc = where(strpos(strlowcase(strsplit(datafiles[i], '_',/extract)), 'mhz') gt -1, count)
        if count gt 1 then stop else freq_loc = freq_loc[0]
        frequencies[i] = float(temp[freq_loc])
        
        data = mrdfits(datafiles[0], 1, hdr)
        
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
        
        xx_col = where(strpos(strlowcase(col_types), 'pp') gt -1, count)
        if count ne 1 then stop else xx_col = xx_col[0]
        
        if i eq 0 then xx_data = fltarr(npix, n_freq)
        
        xx_data[*,i] = data.(xx_col)
        
        yy_col = where(strpos(strlowcase(col_types), 'qq') gt -1, count)
        if count ne 1 then stop else yy_col = yy_col[0]
        
        if i eq 0 then yy_data = fltarr(npix, n_freq)
        yy_data[*,i] =  data.(yy_col)
        
        undefine, data
     endfor
     
     for i=0, n_freq-1 do begin
        temp = strsplit(weightfiles[i], '_',/extract)
        freq_loc = where(strpos(strlowcase(strsplit(weightfiles[i], '_',/extract)), 'mhz') gt -1, count)
        if count gt 1 then stop else freq_loc = freq_loc[0]
        if float(temp[freq_loc]) ne frequencies[i] then message, 'data and weights frequencies do not match'
        
        data = mrdfits(weightfiles[0], 1, hdr)
        
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
        
        xx_col = where(strpos(strlowcase(col_types), 'pp') gt -1, count)
        if count ne 1 then stop else xx_col = xx_col[0]
        
        if i eq 0 then xx_weights = fltarr(npix, n_freq)
        
        xx_weights[*,i] = data.(xx_col)
        
        yy_col = where(strpos(strlowcase(col_types), 'qq') gt -1, count)
        if count ne 1 then stop else yy_col = yy_col[0]
        
        if i eq 0 then yy_weights = fltarr(npix, n_freq)
        yy_weights[*,i] =  data.(yy_col)
        
        undefine, data
     endfor
     
     for i=0, n_freq-1 do begin
        temp = strsplit(variancefiles[i], '_',/extract)
        freq_loc = where(strpos(strlowcase(strsplit(variancefiles[i], '_',/extract)), 'mhz') gt -1, count)
        if count gt 1 then stop else freq_loc = freq_loc[0]
        if float(temp[freq_loc]) ne frequencies[i] then message, 'data and variances frequencies do not match'
        
        data = mrdfits(variancefiles[0], 1, hdr)
        
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
        
        xx_col = where(strpos(strlowcase(col_types), 'pp') gt -1, count)
        if count ne 1 then stop else xx_col = xx_col[0]
        
        if i eq 0 then xx_variances = fltarr(npix, n_freq)
        
        xx_variances[*,i] = data.(xx_col)
        
        yy_col = where(strpos(strlowcase(col_types), 'qq') gt -1, count)
        if count ne 1 then stop else yy_col = yy_col[0]
        
        if i eq 0 then yy_variances = fltarr(npix, n_freq)
        yy_variances[*,i] =  data.(yy_col)
        
        undefine, data
     endfor

     save, file = image_cube_savefile, nside, frequencies, pixel_nums, xx_data, yy_data, xx_weights, yy_weights, $
           xx_variances, yy_variances

  endif

  return, image_cube_savefile

end
