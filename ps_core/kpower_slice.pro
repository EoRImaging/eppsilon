

function kpower_slice, power_3d, kx_mpc, ky_mpc, kz_mpc, kperp_lambda_conv, delay_params, hubble_param, noise_3d = noise_3d, $
                       noise_expval_3d = noise_expval_3d, weights_3d = weights_3d, slice_axis = slice_axis, slice_inds = slice_inds, $
                       slice_weights = weights, slice_noise_expval = noise_expval, slice_savefile = slice_savefile
 
  dims = size(power_3d, /dimensions)
  wt_dims = size(weights_3d, /dimensions)
  nev_dims = size(noise_expval_3d, /dimensions)
  noise_dims = size(noise_3d, /dimensions)

  if n_elements(dims) ne 3 then message, 'power_3d must be an array with 3 dimensions and ordered (x,y,z)'

  n_kx = dims[0]
  n_ky = dims[1]
  n_kz = dims[2]

  if n_elements(kx_mpc) ne n_kx then $
     message, 'input power file must have a kx_mpc vector with length equal to the first dimension of the power array'
  if n_elements(ky_mpc) ne n_ky then $
     message, 'input power file must have a ky_mpc vector with length equal to the second dimension of the power array'
  if n_elements(kz_mpc) ne n_kz then $
     message, 'input power file must have a kz_mpc vector with length equal to the third dimension of the power array'

  if n_elements(weights_3d) gt 0 then if total(abs(dims-wt_dims)) ne 0 then message, 'weights_3d must have same dimensions as power_3d'
  if n_elements(noise_expval_3d) gt 0 then if total(abs(dims-nev_dims)) ne 0 then $
     message, 'noise_expval_3d must have same dimensions as power_3d'
  if n_elements(noise_3d) gt 0 then if total(abs(dims-noise_dims)) ne 0 then message, 'noise_3d must have same dimensions as power_3d'

  if n_elements(slice_axis) eq 0 then slice_axis = 0
  wh = where(indgen(3) eq slice_axis, count)
  if count eq 0 then message, 'slice_axis not recognized'

  if n_elements(slice_inds) eq 0 then begin
     if slice_axis eq 0 then slice_inds = dims[slice_axis]/2 $
     else slice_inds = 0
  endif

  case slice_axis of
     0: begin
        power = power_3d[slice_inds,*,*]
        if n_elements(weights_3d) gt 0 then weights = weights_3d[slice_inds,*,*]
        if n_elements(noise_expval_3d) gt 0 then noise_expval = noise_expval_3d[slice_inds,*,*]
        if n_elements(noise_3d) gt 0 then noise = noise_3d[slice_inds,*,*]
        xarr = ky_mpc
        yarr = kz_mpc
        slice_name = 'x'
        plane_name = 'ky-kz'
        plot_xname = 'y'
        plot_yname = 'z'
     end
     1: begin
        power = power_3d[*,slice_inds,*]
        if n_elements(weights_3d) gt 0 then weights = weights_3d[*,slice_inds,*]
        if n_elements(noise_expval_3d) gt 0 then noise_expval = noise_expval_3d[*,slice_inds,*]
        if n_elements(noise_3d) gt 0 then noise = noise_3d[*,slice_inds,*]
        xarr = kx_mpc
        yarr = kz_mpc
        slice_name = 'y'
        plane_name = 'kx-kz'
        plot_xname = 'x'
        plot_yname = 'z'
     end
     2: begin
        power = power_3d[*,*,slice_inds]
        if n_elements(weights_3d) gt 0 then weights = weights_3d[*,*,slice_inds]
        if n_elements(noise_expval_3d) gt 0 then noise_expval = noise_expval_3d[*,*,slice_inds]
        if n_elements(noise_3d) gt 0 then noise = noise_3d[*,*,slice_inds]
        xarr = kx_mpc
        yarr = ky_mpc
        slice_name = 'z'
        plane_name = 'kx-ky'
        plot_xname = 'x'
        plot_yname = 'y'
     end
  endcase  

  if n_elements(slice_savefile) eq 1 then $
     save, file = slice_savefile, power, noise, noise_expval, weights, kperp_lambda_conv, delay_params, hubble_param, $
           slice_axis, slice_inds, xarr, yarr, slice_name, plane_name, plot_xname, plot_yname

  return, power

end
