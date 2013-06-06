
function limits_fom, limits, n_inside = n_inside, n_out_inside = n_out_inside

  common func_data, n_pts, x_in, y_in, x_out, y_out

  x_above_min = (x_in  - limits[0]) gt 0
  y_above_min = (y_in - limits[1]) gt 0
  x_below_max = (limits[2] - x_in) gt 0
  y_below_max = (limits[3] - y_in) gt 0
  
  xy_inside = x_above_min * y_above_min * x_below_max * y_below_max
  n_inside = total(xy_inside)

  x_out_above_min = (x_out  - limits[0]) gt 0
  y_out_above_min = (y_out - limits[1]) gt 0
  x_out_below_max = (limits[2] - x_out) gt 0
  y_out_below_max = (limits[3] - y_out) gt 0
  
  xy_out_inside = x_out_above_min * y_out_above_min * x_out_below_max * y_out_below_max
  n_out_inside = total(xy_out_inside)

  xlen = limits[2] - limits[0]
  ylen = limits[3] - limits[1]

  max_len_ratio = max([ylen, xlen]) / min([ylen, xlen])

  ;;fom = n_pts/4d * double(n_out_inside) + double((n_pts - n_inside)) 
  fom = n_pts/4d * double(n_out_inside) + double((n_pts - n_inside)) + (max_len_ratio-1)*n_pts/100.

  if limits[2] lt limits[0] or limits[3] lt limits[1] then fom = double(n_pts)*double(n_elements(x_out))

  return, fom

end

function healpix_limits, x_pix, y_pix, x_out_pix, y_out_pix

  common func_data, n_pts, x_in, y_in, x_out, y_out

  x_in = x_pix
  y_in = y_pix
  x_out = x_out_pix
  y_out = y_out_pix

  n_pts = n_elements(x_in)
  if n_elements(y_in) ne n_pts then message, 'x_in and y_in must have the same number of elements'
  n_pts_out = n_elements(x_out)
  if n_elements(y_out) ne n_pts_out then message, 'x_out and y_out must have the same number of elements'

  init_lims = [min(x_in), min(y_in), max(x_in), max(y_in)]
  length_scales = [abs(min(x_in) - min(x_out)), abs(min(y_in) - min(y_out)), abs(max(x_out) - max(x_in)), abs(max(y_out) - max(y_in))]

  init_fom = limits_fom(init_lims)

  ;; running amoeba twice seems to be less sensitive to initial guess
  lims = amoeba(1.0e-5, function_name = 'limits_fom', scale=length_scales, p0 = init_lims, function_value=fval, ncalls = ncalls)
  lims2 = amoeba(1.0e-5, function_name = 'limits_fom', scale=length_scales, p0 = lims, function_value=fval, ncalls = ncalls)

  final_fom =  limits_fom(lims2)

  limits = lims2
  return, limits


end
