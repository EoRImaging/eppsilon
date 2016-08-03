function image_mask_eppsilon, x_rot, y_rot, wh_close, filter_name=filter_name,old_code=old_code

  if keyword_set(old_code) then begin
    restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_decon_March2016_small/pix_center_vec.sav'
    x_minmax = minmax(pix_center_vec[*,0])
    y_minmax = minmax(pix_center_vec[*,1])
    
    x_smooth_boundary = FLTARR(2)
    y_smooth_boundary = FLTARR(2)
    
    x_extent = x_minmax[1] - x_minmax[0]
    y_extent = y_minmax[1] - y_minmax[0]
    x_smooth_boundary[1] = x_minmax[1] - x_extent/8
    x_smooth_boundary[0] = x_minmax[0] + x_extent/8
    y_smooth_boundary[1] = y_minmax[1] - y_extent/14
    y_smooth_boundary[0] = y_minmax[0] + y_extent/14
    
    
    elements=1000
    dimension=1000
    x_arr=FINDGEN(elements)*(x_extent/elements) + x_minmax[0]
    y_arr=FINDGEN(dimension)*(y_extent/dimension) + y_minmax[0]
    
    width_smooth=x_extent/4.
    rarray=Sqrt((meshgrid(elements,1)*(x_extent/elements)-x_extent/2.)^2.+(meshgrid(dimension,2)*(y_extent/dimension)-y_extent/2.)^2.)
    cut_i0=where(x_arr GT x_smooth_boundary[1],n_cut)
    cut_i1=where(x_arr LT x_smooth_boundary[0],n_cut)
    cut_i2=where(y_arr GT y_smooth_boundary[1],n_cut)
    cut_i3=where(y_arr LT y_smooth_boundary[0],n_cut)
    cut_ix=[cut_i0,cut_i1]
    cut_iy=[cut_i2,cut_i3]
    mask_bt=fltarr(dimension,elements)+1.
    
    mask_bt[cut_ix,*]=0
    mask_bt[*,cut_iy]=0
    IF Keyword_Set(width_smooth) THEN mask_bt=Smooth(mask_bt,200, /Edge_truncate)
    ;cut_i=where(dirty_image_uv*mask_bt EQ 0,n_cut,comp=keep_i,ncomp=n_keep)
    ;di_uv_use*=mask_bt
    
    pix_center_vec_expand = pix_center_vec
    pix_center_vec_expand[*,0] = ( pix_center_vec[*,0] - x_minmax[0] ) * 1000/x_extent
    pix_center_vec_expand[*,1] = ( pix_center_vec[*,1] - y_minmax[0] ) * 1000/y_extent
    
    pix_mask=FLTARR(N_elements(pix_center_vec[*,0]),192)
    
    pix_mask_per_freq=interpolate(mask_bt, pix_center_vec_expand[*,0],pix_center_vec_expand[*,1])
  endif
  
  x_minmax = minmax(x_rot[wh_close])
  y_minmax = minmax(y_rot[wh_close])
  x_extent = x_minmax[1] - x_minmax[0]
  y_extent = y_minmax[1] - y_minmax[0]
  
  if ~keyword_set(filter_name) then filter_name='Blackman-Harris'
  mask_bh = spectral_window(1000, type = filter_name, periodic = periodic)
  xy_mask_bh = mask_bh # transpose(mask_bh)

  pix_center_x = ( x_rot[wh_close] - x_minmax[0] ) * 1000/x_extent
  pix_center_y = ( y_rot[wh_close] - y_minmax[0] ) * 1000/y_extent
  
  pix_mask_per_freq=interpolate(xy_mask_bh, pix_center_x, pix_center_y)
  
  pix_mask=FLTARR(N_elements(x_rot[wh_close]),192)
  for freq_i=0,191 do pix_mask[*,freq_i]=pix_mask_per_freq
  
  return, pix_mask
end