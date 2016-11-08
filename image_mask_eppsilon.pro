function image_mask_eppsilon, x_rot, y_rot, wh_close, n_freq=n_freq, filter_name=filter_name, alpha=alpha
	
	x_minmax = minmax(x_rot[wh_close])
	y_minmax = minmax(y_rot[wh_close])
	x_extent = x_minmax[1] - x_minmax[0]
	y_extent = y_minmax[1] - y_minmax[0]
	
	;Make a 1000 element mask given the filter type.
	if filter_name NE 'none' then mask_bh = spectral_window(1000, type = filter_name, periodic = periodic, alpha=alpha) else  begin
		mask_bh = FLTARR(1000)
		mask_bh[*]=1.
	endelse
	;Make a 1000x1000 mask
	xy_mask_bh = mask_bh # transpose(mask_bh)
	
	;Find locations of pixel centers
	pix_center_x = ( x_rot[wh_close] - x_minmax[0] ) * N_elements(mask_bh)/x_extent
	pix_center_y = ( y_rot[wh_close] - y_minmax[0] ) * N_elements(mask_bh)/y_extent
	
	;Interpolate the 1000x1000 mask to the pixel centeres
	pix_mask_per_freq=interpolate(xy_mask_bh, pix_center_x, pix_center_y)
	
	;Build up the mask over frequency
	pix_mask=FLTARR(N_elements(x_rot[wh_close]),n_freq)
	for freq_i=0,n_freq-1 do pix_mask[*,freq_i]=pix_mask_per_freq
	
	;Normalization
	norm_factor = FLTARR(n_freq)
	for freq_i=0, n_freq-1 do begin
		norm_factor[freq_i] = sqrt(n_freq/total(pix_mask[*,freq_i]^2.))
		pix_mask[*,freq_i] = pix_mask[*,freq_i] * norm_factor[freq_i]
	endfor

	return, pix_mask
end