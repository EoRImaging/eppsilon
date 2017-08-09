function prep_uvf_cube, filename, varname, z_mpc_mean, healpix, uvf_input, git_hash, kx_mpc, ky_mpc, $
    pix_area_rad = pix_area_rad, pol_index = pol_index, uv_avg = uv_avg, uv_img_clip = uv_img_clip
    
  cube = getvar_savefile(filename, varname)
  
  void = getvar_savefile(filename, names = names)
  wh_hash = where(strmatch(names, '*git_hash', /fold_case) eq 1, count_hash)
  if count_hash gt 0 then begin
    git_name = names[wh_hash[0]]
    git_hash = getvar_savefile(filename, git_name)
  endif else git_hash = ''
  
  if healpix or not uvf_input then begin
  
    if max(abs(cube)) eq 0 then message, varname + ' in file ' + filename + ' is entirely zero.'
    
    kx_rad_vals = getvar_savefile(filename, 'kx_rad_vals')
    ky_rad_vals = getvar_savefile(filename, 'ky_rad_vals')
    
    kx_rad_delta = kx_rad_vals[1] - kx_rad_vals[0]
    this_kx_mpc = temporary(kx_rad_vals) / z_mpc_mean
    n_kx = n_elements(this_kx_mpc)
    
    ky_rad_delta = ky_rad_vals[1] - ky_rad_vals[0]
    this_ky_mpc = temporary(ky_rad_vals) / z_mpc_mean
    n_ky = n_elements(this_ky_mpc)
    
    if min(this_ky_mpc) lt 0 then begin
      ;; negative ky values haven't been cut yet
      ;; need to cut uvf cubes in half because image is real -- we'll cut negative ky
      cube = cube[*, n_ky/2:n_ky-1,*]
      
      this_ky_mpc = this_ky_mpc[n_ky/2:n_ky-1]
    endif
    
    ;; multiply data, weights & variance cubes by pixel_area_rad to get proper units from DFT
    ;; (not squared for variance because they weren't treated as units squared in FHD code)
    cube = cube * pix_area_rad
    
  endif else begin
  
    cube_size = size(cube)
    if cube_size[n_elements(cube_size)-2] eq 10 then begin
      ;; cube is a pointer
      dims2 = size(*cube[0], /dimension)
      n_freq = cube_size[2]
      iter=1
      while min(dims2) eq 0 and iter lt n_freq do begin
        print, 'warning: some frequency slices have null pointers'
        dims2 = size(*cube[iter], /dimension)
        iter = iter+1
      endwhile
      
      temp = complex(fltarr([dims2, n_freq]))
      if nfiles eq 2 then temp2 = complex(fltarr([dims2, n_freq]))
      for i = 0, n_freq-1 do begin
        foo = *cube[pol_index, i]
        if foo ne !null then begin
          temp[*,*,i] = temporary(foo)
        endif
      endfor
      undefine_fhd, cube
      
      cube = temporary(temp)
      
    endif else begin
      dims2 = size(cube, /dimension)
      n_freq = dims2[2]
    endelse
    
    if max(abs(cube)) eq 0 then message, varname + ' in file ' + filename + ' is entirely zero.'
    
    n_kx = dims2[0]
    if abs(file_struct.kpix-1/(n_kx[0] * (abs(file_struct.degpix) * !pi / 180d)))/file_struct.kpix gt 1e-4 then $
      message, 'Something has gone wrong with calculating uv pixel size'
    kx_mpc_delta = (2.*!pi)*file_struct.kpix / z_mpc_mean
    this_kx_mpc = (dindgen(n_kx)-n_kx/2) * kx_mpc_delta
    
    n_ky = dims2[1]
    ky_mpc_delta = (2.*!pi)*file_struct.kpix / z_mpc_mean
    this_ky_mpc = (dindgen(n_ky)-n_ky/2) * ky_mpc_delta
    
    if keyword_set(uv_avg) then begin
      nkx_new = floor(n_kx / uv_avg)
      temp = complex(fltarr(nkx_new, n_ky, n_freq))
      temp_kx = fltarr(nkx_new)
      for i=0, nkx_new-1 do begin
        temp[i,*,*] = total(cube[i*uv_avg:(i+1)*uv_avg-1,*,*], 1) / uv_avg
        temp_kx[i] = total(this_kx_mpc[i*uv_avg:(i+1)*uv_avg-1]) / uv_avg
      endfor
      
      nky_new = floor(n_ky / uv_avg)
      temp3 = complex(fltarr(nkx_new, nky_new, n_freq))
      temp_ky = fltarr(nky_new)
      for i=0, nky_new-1 do begin
        temp3[*,i,*] = total(temp[*,i*uv_avg:(i+1)*uv_avg-1,*], 2) / uv_avg
        temp_ky[i] = total(this_ky_mpc[i*uv_avg:(i+1)*uv_avg-1]) / uv_avg
      endfor
      undefine, temp
      
      cube = temporary(temp3)
      
      if max(abs(cube)) eq 0 then message, varname + ' in file ' + filename + ' is entirely zero after averaging.'
      
      this_kx_mpc = temporary(temp_kx)
      this_ky_mpc = temporary(temp_ky)
      n_kx = nkx_new
      n_ky = nky_new
    endif
    
    if keyword_set(uv_img_clip) then begin
      kx_mpc_delta_old = kx_mpc_delta
      ky_mpc_delta_old = ky_mpc_delta
      
      temp = shift(fft(fft(shift(cube,dims2[0]/2,dims2[1]/2,0), dimension=1),dimension=2),dims2[0]/2,dims2[1]/2,0)
      temp = temp[(dims2[0]/2)-(dims2[0]/uv_img_clip)/2:(dims2[0]/2)+(dims2[0]/uv_img_clip)/2-1, *, *]
      temp = temp[*, (dims2[1]/2)-(dims2[1]/uv_img_clip)/2:(dims2[1]/2)+(dims2[1]/uv_img_clip)/2-1, *]
      
      temp_dims = size(temp, /dimension)
      n_kx = temp_dims[0]
      n_ky = temp_dims[1]
      kx_mpc_delta = kx_mpc_delta * uv_img_clip
      this_kx_mpc = (dindgen(n_kx)-n_kx/2) * kx_mpc_delta
      ky_mpc_delta = ky_mpc_delta * uv_img_clip
      this_ky_mpc = (dindgen(n_kx)-n_kx/2) * ky_mpc_delta
      
      temp = shift(fft(fft(shift(temp,n_kx/2,n_ky/2,0), dimension=1, /inverse), dimension=2, /inverse),n_kx/2,n_ky/2,0)
      
      cube = temp
      
      if max(abs(cube)) eq 0 then message, varname + ' in file ' + filename + ' is entirely zero after image clip.'
    endif
    
    if n_elements(freq_ch_range) ne 0 then begin
      cube = cube[*, *, min(freq_ch_range):max(freq_ch_range)]
    endif
    if n_elements(freq_flags) ne 0 then begin
      flag_arr = rebin(reform(freq_mask, 1, 1, n_elements(file_struct.frequencies)), size(weights_cube1,/dimension), /sample)
      cube = cube * flag_arr
    endif
    
    ;; need to cut uvf cubes in half because image is real -- we'll cut negative ky
    cube = cube[*, n_ky/2:n_ky-1,*]
    this_ky_mpc = this_ky_mpc[n_ky/2:n_ky-1]
    
  endelse
  
  if n_elements(kx_mpc) eq 0 then begin
    kx_mpc = this_kx_mpc
  endif else begin
    if max(abs(this_kx_mpc - kx_mpc)) gt 0 then message, 'kx values in file ' + filename + ' do not match other files'
  endelse
  
  if n_elements(ky_mpc) eq 0 then begin
    ky_mpc = this_ky_mpc
  endif else begin
    if max(abs(this_ky_mpc - ky_mpc)) gt 0 then message, 'ky values in file ' + filename + ' do not match other files'
  endelse
  
  ;; Also need to drop 1/2 of ky=0 line -- drop ky=0, kx<0
  cube[0:n_kx/2-1, 0, *] = 0
  
  return, cube
  
end
