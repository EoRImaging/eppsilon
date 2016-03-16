pro ps_slice_differences, slice_file1, slice_file2, savefile_diff = savefile_diff 


  slice_axis = getvar_savefile(slice_file1, 'slice_axis')
  slice_axis2 = getvar_savefile(slice_file2, 'slice_axis')
  if n_elements(slice_axis) ne n_elements(slice_axis2) or max(abs(slice_axis - slice_axis2)) gt 0 then message, 'slice_axis do not match between slices'
  slice_inds = getvar_savefile(slice_file1, 'slice_inds')
  slice_inds2 = getvar_savefile(slice_file2, 'slice_inds')
  if n_elements(slice_inds) ne n_elements(slice_inds2) or max(abs(slice_inds - slice_inds2)) gt 0 then message, 'slice_inds do not match between slices'
  xarr = getvar_savefile(slice_file1, 'xarr')
  xarr2 = getvar_savefile(slice_file2, 'xarr')
  if n_elements(xarr) ne n_elements(xarr2) or max(abs(xarr - xarr2)) gt 1.05e-3 then message, 'xarr does not match between slices'
  yarr = getvar_savefile(slice_file1, 'yarr')
  yarr2 = getvar_savefile(slice_file2, 'yarr')
  if n_elements(yarr) ne n_elements(yarr2) or max(abs(yarr - yarr2)) gt 1.05e-3 then message, 'yarr does not match between slices'
  
  kperp_lambda_conv = getvar_savefile(slice_file1, 'kperp_lambda_conv')
  kperp_lambda_conv2 = getvar_savefile(slice_file2, 'kperp_lambda_conv')
  if n_elements(kperp_lambda_conv) ne n_elements(kperp_lambda_conv2) or max(abs(kperp_lambda_conv - kperp_lambda_conv2)) gt 1.05e-3 then $
    message, 'kperp_lambda_conv does not match between slices'
  delay_params = getvar_savefile(slice_file1, 'delay_params')
  delay_params2 = getvar_savefile(slice_file2, 'delay_params')
  if n_elements(delay_params) ne n_elements(delay_params2) or max(abs(delay_params - delay_params2)) gt 1.05e-3 then $
    message, 'delay_params does not match between slices'
  hubble_param = getvar_savefile(slice_file1, 'hubble_param')
  hubble_param2 = getvar_savefile(slice_file2, 'hubble_param')
  if n_elements(hubble_param) ne n_elements(hubble_param2) or max(abs(hubble_param - hubble_param2)) gt 1.05e-3 then $
    message, 'hubble_param does not match between slices'
    
    
  ;; if the above are all the same these should be too
  slice_name = getvar_savefile(slice_file1, 'slice_name')
  plane_name = getvar_savefile(slice_file1, 'plane_name')
  plot_xname = getvar_savefile(slice_file1, 'plot_xname')
  plot_yname = getvar_savefile(slice_file1, 'plot_yname')
  
  power1 = getvar_savefile(slice_file1, 'power')
  power2 = getvar_savefile(slice_file2, 'power')
  power_diff = power1 - power2
  
  if max(abs(power_diff)) eq 0 then begin
    print, 'The cubes are identical -- power difference is zero everywhere'
  ;continue
  endif
  undefine, power1, power2
  
  weights1 = getvar_savefile(slice_file1, 'weights')
  weights2 = getvar_savefile(slice_file2, 'weights')
  
  ;; variance = 1/weights
  var1 = 1./weights1
  wh_wt1_0 = where(weights1 eq 0, count_wt1_0)
  if count_wt1_0 gt 0 then var1[wh_wt1_0] = 0
  var2 = 1./weights2
  wh_wt2_0 = where(weights2 eq 0, count_wt2_0)
  if count_wt2_0 gt 0 then var2[wh_wt2_0] = 0
  undefine, weights1, weights2
  
  var_diff = var1 + var2
  weight_diff = 1/var_diff
  if count_wt1_0 gt 0 then weight_diff[wh_wt1_0] = 0
  if count_wt2_0 gt 0 then weight_diff[wh_wt2_0] = 0
  undefine, var1, var2, var_diff
  
  
  power = power_diff
  weights = weight_diff
  save, file = savefile_diff, power, weights, kperp_lambda_conv, delay_params, hubble_param, $
    slice_axis, slice_inds, xarr, yarr, slice_name, plane_name, plot_xname, plot_yname
    
end