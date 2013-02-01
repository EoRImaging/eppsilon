function uvf_slice, uvf_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, delay_params, slice_axis = slice_axis, $
                    slice_inds = slice_inds, slice_savefile = slice_savefile

 
  dims = size(uvf_cube, /dimensions)

  if n_elements(dims) ne 3 then message, 'uvf_cube must be an array with 3 dimensions and ordered (x,y,z)'

  n_kx = dims[0]
  n_ky = dims[1]
  n_f = dims[2]

  if n_elements(kx_mpc) ne n_kx then $
     message, 'input power file must have a kx_mpc vector with length equal to the first dimension of the power array'
  if n_elements(ky_mpc) ne n_ky then $
     message, 'input power file must have a ky_mpc vector with length equal to the second dimension of the power array'
  if n_elements(frequencies) ne n_f then $
     message, 'input power file must have a kz_mpc vector with length equal to the third dimension of the power array'

  if n_elements(slice_axis) eq 0 then slice_axis = 0
  wh = where(indgen(3) eq slice_axis, count)
  if count eq 0 then message, 'slice_axis not recognized'

  if n_elements(slice_inds) eq 0 then begin
     if slice_axis lt 2 then slice_inds = dims[slice_axis]/2 $
     else slice_inds = 0
  endif

  case slice_axis of
     0: begin
        uvf_slice = uvf_cube[slice_inds,*,*]
        xarr = ky_mpc
        yarr = frequencies
        slice_name = 'u'
        plane_name = 'vf'
        plot_xname = 'v'
        plot_yname = 'f'
     end
     1: begin
        uvf_slice = uvf_cube[*,slice_inds,*]
        xarr = kx_mpc
        yarr = frequencies
        slice_name = 'v'
        plane_name = 'uf'
        plot_xname = 'u'
        plot_yname = 'f'
     end
     2: begin
        uvf_slice = uvf_cube[*,*,slice_inds]
        xarr = kx_mpc
        yarr = ky_mpc
        slice_name = 'f'
        plane_name = 'uv'
        plot_xname = 'u'
        plot_yname = 'v'
     end
  endcase  

  if n_elements(slice_savefile) eq 1 then $
     save, file = slice_savefile, uvf_slice, kperp_lambda_conv, delay_params, slice_axis, slice_inds, xarr, yarr, slice_name, $
           plane_name, plot_xname, plot_yname

  return, uvf_slice
end
