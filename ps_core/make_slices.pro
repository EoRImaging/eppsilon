pro make_slices, file_struct, type=type, file_ind = file_ind, data_cube = data_cube, $
  weights_cube = weights_cube, var_cube = var_cube, $
  noise_3d = noise_3d, noise_expval_3d = noise_expval_3d, $
  kx_mpc = kx_mpc, ky_mpc = ky_mpc, frequencies = frequencies, kz_mpc = kz_mpc, $
  kperp_lambda_conv = kperp_lambda_conv, delay_params = delay_params, $
  hubble_param = hubble_param

  case type of
    'raw': begin
      uf_tot = total(total(abs(weights_cube),3),1)
      wh_uf_n0 = where(uf_tot gt 0, count_uf_n0)
      if count_uf_n0 eq 0 then message, 'uvf weights appear to be entirely zero'
      min_dist_uf_n0 = min(wh_uf_n0, min_loc)
      uf_slice_ind = wh_uf_n0[min_loc]

      uf_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, $
        delay_params, hubble_param, slice_axis = 1, slice_inds = uf_slice_ind, $
          slice_savefile = file_struct.uf_raw_savefile[file_ind])

      uf_weight_slice = uvf_slice(weights_cube, kx_mpc, ky_mpc, frequencies, $
        kperp_lambda_conv, delay_params, hubble_param, slice_axis = 1, $
        slice_inds = uf_slice_ind, slice_savefile = file_struct.uf_weight_savefile[file_ind])

      uf_var_slice = uvf_slice(var_cube, kx_mpc, ky_mpc, frequencies, $
        kperp_lambda_conv, delay_params, hubble_param, slice_axis = 1, $
        slice_inds = uf_slice_ind, slice_savefile = file_struct.uf_var_savefile[file_ind])
      undefine, uf_slice, uf_weight_slice, uf_var_slice

      vf_tot = total(total(abs(weights_cube),3),2)
      wh_vf_n0 = where(vf_tot gt 0, count_vf_n0)
      if count_vf_n0 eq 0 then message, 'uvf weights appear to be entirely zero'
      n_kx = n_elements(kx_mpc)
      min_dist_vf_n0 = min(abs(n_kx/2-wh_vf_n0), min_loc)
      vf_slice_ind = wh_vf_n0[min_loc]

      vf_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, $
        delay_params, hubble_param, slice_axis = 0, slice_inds = vf_slice_ind, $
        slice_savefile = file_struct.vf_raw_savefile[file_ind])

      vf_weight_slice = uvf_slice(weights_cube, kx_mpc, ky_mpc, frequencies, $
        kperp_lambda_conv, delay_params, hubble_param, slice_axis = 0, $
        slice_inds = vf_slice_ind, slice_savefile = file_struct.vf_weight_savefile[file_ind])

      vf_var_slice = uvf_slice(var_cube, kx_mpc, ky_mpc, frequencies, $
        kperp_lambda_conv, delay_params, hubble_param, slice_axis = 0, $
        slice_inds = vf_slice_ind, slice_savefile = file_struct.vf_var_savefile[file_ind])

      if max(abs(vf_slice)) eq 0 then message, 'vf data slice is entirely zero'
      undefine, vf_slice, vf_weight_slice, vf_var_slice

      uv_tot = total(total(abs(weights_cube),2),1)
      wh_uv_n0 = where(uv_tot gt 0, count_uv_n0)
      if count_uv_n0 eq 0 then message, 'uvf weights appear to be entirely zero'
      min_dist_uv_n0 = min(wh_uv_n0, min_loc)
      uv_slice_ind = wh_uv_n0[min_loc]

      uv_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, $
        delay_params, hubble_param, slice_axis = 2, slice_inds = uv_slice_ind, $
        slice_savefile = file_struct.uv_raw_savefile[file_ind])

      uv_weight_slice = uvf_slice(weights_cube, kx_mpc, ky_mpc, frequencies, $
        kperp_lambda_conv, delay_params, hubble_param, slice_axis = 2, $
        slice_inds = uv_slice_ind, slice_savefile = file_struct.uv_weight_savefile[file_ind])

      uv_var_slice = uvf_slice(var_cube, kx_mpc, ky_mpc, frequencies, $
        kperp_lambda_conv, delay_params, hubble_param, slice_axis = 2, $
        slice_inds = uv_slice_ind, slice_savefile = file_struct.uv_var_savefile[file_ind])

      if max(abs(uv_slice)) eq 0 then message, 'uv data slice is entirely zero'
      undefine, uv_slice, uv_weight_slice, uv_var_slice
    end
    'divided': begin
      uf_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, $
        delay_params, hubble_param, slice_axis = 1, $
        slice_inds = uf_slice_ind, slice_savefile = file_struct.uf_savefile[file_ind])
      undefine, uf_slice

      vf_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, $
        delay_params, hubble_param, slice_axis = 0, $
        slice_inds = vf_slice_ind, slice_savefile = file_struct.vf_savefile[file_ind])
      undefine, vf_slice

      uv_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, $
        delay_params, hubble_param, slice_axis = 2, $
        slice_inds = uv_slice_ind, slice_savefile = file_struct.uv_savefile[file_ind])
      undefine, uv_slice
    end
    'sum': begin
      uf_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, $
        delay_params, hubble_param, slice_axis = 1, $
        slice_inds = uf_slice_ind, slice_savefile = file_struct.uf_sum_savefile)
      undefine, uf_slice

      vf_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, $
        delay_params, hubble_param, slice_axis = 0, $
        slice_inds = vf_slice_ind, slice_savefile = file_struct.vf_sum_savefile)
      undefine, vf_slice

      uv_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, $
        delay_params, hubble_param, slice_axis = 2, $
        slice_inds = uv_slice_ind, slice_savefile = file_struct.uv_sum_savefile)
      undefine, uv_slice
    end
    'diff': begin
      uf_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, $
        delay_params, hubble_param, slice_axis = 1, $
        slice_inds = uf_slice_ind, slice_savefile = file_struct.uf_diff_savefile)
      undefine, uf_slice

      vf_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, $
        delay_params, hubble_param, slice_axis = 0, $
        slice_inds = vf_slice_ind, slice_savefile = file_struct.vf_diff_savefile)
      undefine, vf_slice

      uv_slice = uvf_slice(data_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, $
        delay_params, hubble_param, slice_axis = 2, $
        slice_inds = uv_slice_ind, slice_savefile = file_struct.uv_diff_savefile)
      undefine, uv_slice
    end
    'kspace': begin
      y_tot = total(total(abs(data_cube),3),1)
      wh_y_n0 = where(y_tot gt 0, count_y_n0)
      min_dist_y_n0 = min(wh_y_n0, min_loc)
      y_slice_ind = wh_y_n0[min_loc]

      yslice_savefile = file_struct.xz_savefile
      yslice_power = kpower_slice(data_cube, kx_mpc, ky_mpc, kz_mpc, kperp_lambda_conv, $
        delay_params, hubble_param, noise_3d = noise_3d, noise_expval_3d = noise_expval_3d, $
        weights_3d = weights_cube, slice_axis = 1, slice_inds = y_slice_ind, $
        slice_savefile = yslice_savefile)

      x_tot = total(total(abs(data_cube),3),2)
      wh_x_n0 = where(x_tot gt 0, count_x_n0)
      n_kx = n_elements(kx_mpc)
      min_dist_x_n0 = min(abs(n_kx/2-wh_x_n0), min_loc)
      x_slice_ind = wh_x_n0[min_loc]

      xslice_savefile = file_struct.yz_savefile
      xslice_power = kpower_slice(data_cube, kx_mpc, ky_mpc, kz_mpc, kperp_lambda_conv, $
        delay_params, hubble_param, noise_3d = noise_3d, noise_expval_3d = noise_expval_3d, $
        weights_3d = weights_cube, slice_axis = 0, slice_inds = x_slice_ind, $
        slice_savefile = xslice_savefile)

      z_tot = total(total(abs(data_cube),2),1)
      wh_z_n0 = where(z_tot gt 0, count_z_n0)
      min_dist_z_n0 = min(wh_z_n0, min_loc)
      z_slice_ind = wh_z_n0[min_loc]

      zslice_savefile = file_struct.xy_savefile
      zslice_power = kpower_slice(data_cube, kx_mpc, ky_mpc, kz_mpc, kperp_lambda_conv, $
        delay_params, hubble_param, noise_3d = noise_3d, noise_expval_3d = noise_expval_3d, $
        weights_3d = weights_cube, slice_axis = 2, slice_inds = z_slice_ind, $
        slice_savefile = zslice_savefile)
    end
  endcase
end
