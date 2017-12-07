function healpix_rot_calc, x_locs, y_locs, ang_range, ang_step

  n_angles = (ang_range[1] - ang_range[0]) / ang_step
  angles = (dindgen(n_angles) * (ang_range[1]-ang_range[0])*!dpi/180.) / n_angles + ang_range[0]*!dpi/180.

  cos_arr = cos(angles)
  sin_arr = sin(angles)

  x_rot = matrix_multiply(x_locs, cos_arr) - matrix_multiply(y_locs, sin_arr)
  y_rot = matrix_multiply(x_locs, sin_arr) + matrix_multiply(y_locs, cos_arr)
  undefine, cos_arr, sin_arr

  area = dblarr(n_angles)
  limits = dblarr(n_angles,4)
  limits[*,0] = min(x_rot, dimension=1)
  limits[*,1] = min(y_rot, dimension=1)
  limits[*,2] = max(x_rot, dimension=1)
  limits[*,3] = max(y_rot, dimension=1)
  area = (limits[*,2] - limits[*,0]) * (limits[*,3] - limits[*,1])

  wh = where(area eq min(area), count)

  if count eq 1 then rot_angle = angles[wh[0]] $
  else if count eq 2 and abs(wh[1] - wh[0]) eq 1 then rot_angle = mean(angles[wh]) $
  else stop

  return, rot_angle
end

function healpix_rot, x_locs, y_locs

  n_pts = n_elements(x_locs)
  if n_elements(y_locs) ne n_pts then message, 'x_locs and y_locs must have the same number of elements'

  ang_range = [0, 90]
  ang_step = 1
  rot_ang1 = healpix_rot_calc(x_locs, y_locs, ang_range, ang_step)

  ang_range = [rot_ang1*180./!dpi - 2*ang_step, rot_ang1*180./!dpi + 2*ang_step]
  ang_step = 0.1
  rot_angle = healpix_rot_calc(x_locs, y_locs, ang_range, ang_step)

  return, rot_angle

end
