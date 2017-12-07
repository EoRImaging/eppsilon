function casa_fits2idlcube, datafiles, weightfiles, variancefiles, pol_inc, save_path = save_path, refresh = refresh

  if n_elements(datafiles) eq 0 then message, 'datafiles must be passed in'
  if n_elements(weightfiles) ne n_elements(datafiles) then message, 'different number of weightfiles and datafiles'
  if n_elements(variancefiles) ne n_elements(datafiles) then message, 'different number of variancefiles and datafiles'
  
  pol_inc = ['xx', 'yy']
  npol = n_elements(pol_inc)
  
  ;; check for even/odd files,
  even_mask = stregex(datafiles, '_0_', /boolean)
  n_even = total(even_mask)
  odd_mask = stregex(datafiles, '_1_', /boolean)
  n_odd = total(odd_mask)
  
  ;; check for different polarization files
  pol_mask = intarr(n_elements(datafiles), npol)
  wt_pol_mask = intarr(n_elements(weightfiles), npol)
  var_pol_mask = intarr(n_elements(variancefiles), npol)
  n_pol_mask = intarr(npol)
  wt_n_pol_mask = intarr(npol)
  var_n_pol_mask = intarr(npol)
  for i=0, npol-1 do begin
    pol_mask[*,i] = stregex(datafiles, '_' + pol_inc[i] + '_', /boolean)
    n_pol_mask[i] = total(pol_mask[*,i])
    wt_pol_mask[*,i] = stregex(weightfiles, '_' + pol_inc[i] + '_', /boolean)
    wt_n_pol_mask[i] = total(wt_pol_mask[*,i])
    var_pol_mask[*,i] = stregex(variancefiles, '_' + pol_inc[i] + '_', /boolean)
    var_n_pol_mask[i] = total(var_pol_mask[*,i])
  endfor
  
  if n_even ne n_odd then message, 'number of even and odd files are not equal'
  if max(abs(n_pol_mask - n_pol_mask[0])) gt 0 then message, 'number of files for each polarization are not equal'
  if max(n_pol_mask) gt 2 or min(n_pol_mask) lt 1 then message, 'Only 1 or 2 datafiles per polarization is allowed'
  
  if n_even eq 0 then begin
    nfiles = 1
    for i=0, npol do begin
      datafile_arr[i] = datafiles[where(pol_mask[*,i] gt 0)]
      weightsfile_arr[i,1] = weightfiles[where(pol_mask[*,i] gt 0)]
      variancefile_arr[i,1] = variancefiles[where(pol_mask[*,i] gt 0)]
    endfor
  endif else begin
    nfiles = 2
    
    wt_even_mask = stregex(weightfiles, '_0_', /boolean)
    wt_n_even = total(wt_even_mask)
    wt_odd_mask = stregex(weightfiles, '_1_', /boolean)
    wt_n_odd = total(wt_odd_mask)
    if wt_n_even ne n_even or wt_n_odd ne n_odd then message, 'number of even and odd weight files do not match datafiles'
    
    var_even_mask = stregex(variancefiles, '_0_', /boolean)
    var_n_even = total(var_even_mask)
    var_odd_mask = stregex(variancefiles, '_1_', /boolean)
    var_n_odd = total(var_odd_mask)
    if var_n_even ne n_even or var_n_odd ne n_odd then message, 'number of even and odd variance files do not match datafiles'
    
    
    datafile_arr = strarr(npol, nfiles)
    weightsfile_arr = strarr(npol, nfiles)
    variancefile_arr = strarr(npol, nfiles)
    for i=0, npol-1 do begin
      datafile_arr[i,0] = datafiles[where(even_mask*pol_mask[*,i] gt 0)]
      datafile_arr[i,1] = datafiles[where(odd_mask*pol_mask[*,i] gt 0)]
      weightsfile_arr[i,0] = weightfiles[where(wt_even_mask*wt_pol_mask[*,i] gt 0)]
      weightsfile_arr[i,1] = weightfiles[where(wt_odd_mask*wt_pol_mask[*,i] gt 0)]
      variancefile_arr[i,0] = variancefiles[where(var_even_mask*var_pol_mask[*,i] gt 0)]
      variancefile_arr[i,1] = variancefiles[where(var_odd_mask*var_pol_mask[*,i] gt 0)]
    endfor
    
  endelse
  
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
    
    for i=0, npol-1 do begin
      this_fileparts = strsplit(strmid(infilebase[i], 0, temp2[i]), '_', /extract)
      if i eq 0 then fileparts = this_fileparts else begin
        this_test = long(strcmp(fileparts, this_fileparts))
        if i eq 1 then match_test = this_test else match_test = match_test + this_test
      endelse
    endfor
    
    wh_diff = where(match_test lt npol-1, count_diff, complement = wh_same, ncomplement = count_same)
    if count_same gt 0 then general_filebase = strjoin(fileparts[wh_same], '_') else stop
    
    idl_cube_savefile[file_i] = froot + general_filebase + '_cube.idlsave'
    
    test_idlcube = file_test(idl_cube_savefile[file_i]) *  (1 - file_test(idl_cube_savefile[file_i], /zero_length))
    if test_idlcube eq 0 or keyword_set(refresh) then begin
    
      ;n_vis_arr = fltarr(npol)
      ;kpix_arr = fltarr(npol)
      ;kspan_arr = fltarr(npol)
      for i=0, npol-1 do begin
        ;        temp = strsplit(datafile_arr[i, file_i], '_',/extract)
        ;        freq_loc = where(strpos(strlowcase(strsplit(datafile_arr[i, file_i], '_',/extract)), 'mhz') gt -1, count)
        ;        if count gt 1 then stop else freq_loc = freq_loc[0]
        ;        frequencies[i] = float(temp[freq_loc])
      
        data = mrdfits(datafile_arr[i, file_i], 0, hdr, /silent)
        
        this_dims = fxpar(hdr, 'naxis*')
        if i eq 0 then begin
          dims = this_dims
        endif else if total(abs(this_dims-dims)) ne 0 then message, 'dimensions do not agree across datafiles'
        
        col_types = fxpar(hdr, 'ctype*')
        
        
        freq_col = where(strpos(strlowcase(col_types), 'freq') gt -1, count)
        if count ne 1 then stop else freq_col = freq_col[0]
        if i eq 0 then begin
          n_freq = dims[freq_col]
          delta_f = fxpar(hdr, 'cdelt'+number_formatter(freq_col+1)) ;;fits columns are 1-indexed
          frequencies = findgen(n_freq)*delta_f + fxpar(hdr, 'crval'+number_formatter(freq_col+1)) - (1-fxpar(hdr, 'crpix'+number_formatter(freq_col+1)))*delta_f
        endif else begin
          if dims[freq_col] ne n_freq then message, 'number of frequencies does not agree across datafiles'
          if abs(delta_f-fxpar(hdr, 'cdelt'+number_formatter(freq_col+1))) gt 0 then message, 'frequency channel width does not agree across datafiles'
          temp = findgen(n_freq)*delta_f + fxpar(hdr, 'crval'+number_formatter(freq_col+1)) - (1-fxpar(hdr, 'crpix'+number_formatter(freq_col+1)))*delta_f
          if max(abs(temp-frequencies)) gt 0 then message, 'frequencies do not agree across datafiles'
        endelse
        
        ra_col = where(strpos(strlowcase(col_types), 'ra') gt -1, count)
        if count ne 1 then stop else ra_col = ra_col[0]
        if i eq 0 then begin
          n_ra = dims[ra_col]
          delta_ra = fxpar(hdr, 'cdelt'+number_formatter(ra_col+1)) ;;fits columns are 1-indexed
          ra_vals = findgen(n_ra)*delta_ra + fxpar(hdr, 'crval'+number_formatter(ra_col+1)) - (1-fxpar(hdr, 'crpix'+number_formatter(ra_col+1)))*delta_ra
        endif else begin
          if dims[ra_col] ne n_ra then message, 'ra dimesion does not agree across datafiles'
          if abs(delta_ra-fxpar(hdr, 'cdelt'+number_formatter(ra_col+1))) gt 0 then message, 'ra pixel width does not agree across datafiles'
          temp = findgen(n_ra)*delta_ra + fxpar(hdr, 'crval'+number_formatter(ra_col+1)) - (1-fxpar(hdr, 'crpix'+number_formatter(ra_col+1)))*delta_ra
          if max(abs(temp-ra_vals)) gt 0 then message, 'ra values do not agree across datafiles'
        endelse
        
        dec_col = where(strpos(strlowcase(col_types), 'dec') gt -1, count)
        if count ne 1 then stop else dec_col = dec_col[0]
        if i eq 0 then begin
          n_dec = dims[dec_col]
          delta_dec = fxpar(hdr, 'cdelt'+number_formatter(dec_col+1)) ;;fits columns are 1-indexed
          dec_vals = findgen(n_dec)*delta_dec + fxpar(hdr, 'crval'+number_formatter(dec_col+1)) - (1-fxpar(hdr, 'crpix'+number_formatter(dec_col+1)))*delta_dec
        endif else begin
          if dims[dec_col] ne n_dec then message, 'dec dimesion does not agree across datafiles'
          if abs(delta_dec-fxpar(hdr, 'cdelt'+number_formatter(dec_col+1))) gt 0 then message, 'dec pixel width does not agree across datafiles'
          temp = findgen(n_dec)*delta_dec + fxpar(hdr, 'crval'+number_formatter(dec_col+1)) - (1-fxpar(hdr, 'crpix'+number_formatter(dec_col+1)))*delta_dec
          if max(abs(temp-dec_vals)) gt 0 then message, 'dec values do not agree across datafiles'
        endelse
        
        case pol_inc[i] of
          'xx': xx_data = temporary(data)
          'yy': yy_data = temporary(data)
        endcase
        
      endfor ;; loop over pol. for data
      
      for i=0, npol-1 do begin
      
        data = mrdfits(weightsfile_arr[i, file_i], 0, hdr, /silent)
        
        this_dims = fxpar(hdr, 'naxis*')
        if total(abs(this_dims-dims)) ne 0 then message, 'weights dimensions do not agree with datafiles'
        
        col_types = fxpar(hdr, 'ctype*')
        
        freq_col = where(strpos(strlowcase(col_types), 'freq') gt -1, count)
        if count ne 1 then stop else freq_col = freq_col[0]
        if dims[freq_col] ne n_freq then message, 'weights number of frequencies does not agree with datafiles'
        if abs(delta_f-fxpar(hdr, 'cdelt'+number_formatter(freq_col+1))) gt 0 then message, 'weights frequency channel width does not agree with datafiles'
        temp = findgen(n_freq)*delta_f + fxpar(hdr, 'crval'+number_formatter(freq_col+1)) - (1-fxpar(hdr, 'crpix'+number_formatter(freq_col+1)))*delta_f
        if max(abs(temp-frequencies)) gt 0 then message, 'weights frequencies do not agree with datafiles'
        
        ra_col = where(strpos(strlowcase(col_types), 'ra') gt -1, count)
        if count ne 1 then stop else ra_col = ra_col[0]
        if dims[ra_col] ne n_ra then message, 'weights ra dimesion does not agree with datafiles'
        if abs(delta_ra-fxpar(hdr, 'cdelt'+number_formatter(ra_col+1))) gt 0 then message, 'weights ra pixel width does not agree with datafiles'
        temp = findgen(n_ra)*delta_ra + fxpar(hdr, 'crval'+number_formatter(ra_col+1)) - (1-fxpar(hdr, 'crpix'+number_formatter(ra_col+1)))*delta_ra
        if max(abs(temp-ra_vals)) gt 0 then message, 'weights ra values do not agree with datafiles'
        
        dec_col = where(strpos(strlowcase(col_types), 'dec') gt -1, count)
        if count ne 1 then stop else dec_col = dec_col[0]
        if dims[dec_col] ne n_dec then message, 'weights dec dimesion does not agree with datafiles'
        if abs(delta_dec-fxpar(hdr, 'cdelt'+number_formatter(dec_col+1))) gt 0 then message, 'weights dec pixel width does not agree with datafiles'
        temp = findgen(n_dec)*delta_dec + fxpar(hdr, 'crval'+number_formatter(dec_col+1)) - (1-fxpar(hdr, 'crpix'+number_formatter(dec_col+1)))*delta_dec
        if max(abs(temp-dec_vals)) gt 0 then message, 'weights dec values do not agree with datafiles'
        
        
        case pol_inc[i] of
          'xx': xx_weights = temporary(data)
          'yy': yy_weights = temporary(data)
        endcase
        
      endfor ;; loop over pol. for weights
      
      for i=0, npol-1 do begin
      
        data = mrdfits(variancefile_arr[i, file_i], 0, hdr, /silent)
        
        this_dims = fxpar(hdr, 'naxis*')
        if total(abs(this_dims-dims)) ne 0 then message, 'variance dimensions do not agree with datafiles'
        
        col_types = fxpar(hdr, 'ctype*')
        
        freq_col = where(strpos(strlowcase(col_types), 'freq') gt -1, count)
        if count ne 1 then stop else freq_col = freq_col[0]
        if dims[freq_col] ne n_freq then message, 'variance number of frequencies does not agree with datafiles'
        if abs(delta_f-fxpar(hdr, 'cdelt'+number_formatter(freq_col+1))) gt 0 then message, 'variance frequency channel width does not agree with datafiles'
        temp = findgen(n_freq)*delta_f + fxpar(hdr, 'crval'+number_formatter(freq_col+1)) - (1-fxpar(hdr, 'crpix'+number_formatter(freq_col+1)))*delta_f
        if max(abs(temp-frequencies)) gt 0 then message, 'variance frequencies do not agree with datafiles'
        
        ra_col = where(strpos(strlowcase(col_types), 'ra') gt -1, count)
        if count ne 1 then stop else ra_col = ra_col[0]
        if dims[ra_col] ne n_ra then message, 'variance ra dimesion does not agree with datafiles'
        if abs(delta_ra-fxpar(hdr, 'cdelt'+number_formatter(ra_col+1))) gt 0 then message, 'variance ra pixel width does not agree with datafiles'
        temp = findgen(n_ra)*delta_ra + fxpar(hdr, 'crval'+number_formatter(ra_col+1)) - (1-fxpar(hdr, 'crpix'+number_formatter(ra_col+1)))*delta_ra
        if max(abs(temp-ra_vals)) gt 0 then message, 'variance ra values do not agree with datafiles'
        
        dec_col = where(strpos(strlowcase(col_types), 'dec') gt -1, count)
        if count ne 1 then stop else dec_col = dec_col[0]
        if dims[dec_col] ne n_dec then message, 'variance dec dimesion does not agree with datafiles'
        if abs(delta_dec-fxpar(hdr, 'cdelt'+number_formatter(dec_col+1))) gt 0 then message, 'variance dec pixel width does not agree with datafiles'
        temp = findgen(n_dec)*delta_dec + fxpar(hdr, 'crval'+number_formatter(dec_col+1)) - (1-fxpar(hdr, 'crpix'+number_formatter(dec_col+1)))*delta_dec
        if max(abs(temp-dec_vals)) gt 0 then message, 'variance dec values do not agree with datafiles'
        
        
        case pol_inc[i] of
          'xx': xx_variances = temporary(data)
          'yy': yy_variances = temporary(data)
        endcase
        
      endfor ;;loop over pol. for variances
      
      save, file = idl_cube_savefile[file_i], frequencies, xx_data, yy_data, xx_weights, yy_weights, xx_variances, yy_variances, ra_vals, dec_vals
      
    endif ;; end refresh
  endfor ;; loop over odd/even
  
  return, idl_cube_savefile
  
end
