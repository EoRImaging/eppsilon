
function get_rot_matrix, theta, phi, inverse = inverse

  rot_matrix = [[cos(theta) * cos(phi)^2d + sin(phi)^2d, (cos(theta)-1d) * cos(phi)*sin(phi), sin(theta) * cos(phi)], $
                [(cos(theta)-1d) * cos(phi)*sin(phi), cos(theta) * sin(phi)^2d + cos(phi)^2d, sin(theta) * sin(phi)], $
                [(-1) * sin(theta) * cos(phi), (-1) * sin(theta) * sin(phi), cos(theta)]]

  if keyword_set(inverse) then return, invert(rot_matrix) else return, rot_matrix

end


pro healpix_setup_ft, pixel_nums, nside, new_pix_vec, limits, kx_rad_vals, ky_rad_vals, quiet = quiet

  pix2vec_ring, nside, pixel_nums, pix_center_vec
  ;; find mid point (work in x/y because of possible jumps in phi)
  vec_mid = [mean(pix_center_vec[*,0]), mean(pix_center_vec[*,1]), mean(pix_center_vec[*,2])]
  theta0 = acos(vec_mid[2])
  phi0 = atan(vec_mid[1], vec_mid[0])
  
  dists = sqrt((pix_center_vec[*,0]-vec_mid[0])^2d + (pix_center_vec[*,1]-vec_mid[1])^2d + (pix_center_vec[*,2]-vec_mid[2])^2d)
  radius = max(dists)
  
  disc_covers = 0
  nloops = 0
  while disc_covers lt 1 do begin
     query_disc, nside, vec_mid, radius, listpix, nlist, /inc
     min_pix = min([pixel_nums, listpix])
     wh2 = where(histogram(listpix, min=min_pix) eq 0 and histogram(pixel_nums, min=min_pix) gt 0, count2)
     if count2 gt 0 then radius = radius * (1 + (nloops+1)*0.1d) else disc_covers = 1
     nloops = nloops+1
  endwhile
  
  ;; remove pixels from listpix that are in my image -- only want nearby pixels not in my image
  min_pix = min([pixel_nums, listpix])
  max_pix = max([pixel_nums, listpix])
  wh = where(histogram(listpix, min = min_pix, max = max_pix) gt 0 and histogram(pixel_nums, min = min_pix, max = max_pix) eq 0, count)
  if count gt 0 then outside_pix = wh + min_pix else stop
  listpix=0
  wh=0
  
  pix2vec_ring, nside, outside_pix, out_center_vec
  
  ;; define new coordinate system
  ;; To get to current location, need to first rotate around z by
  ;; phi, then around y by -theta, then around z by -phi
  ;; use inverse to rotate back to vertical
  rot_matrix = get_rot_matrix(theta0, phi0, /inverse)
  
  new_pix_vec = rot_matrix ## pix_center_vec
  new_out_vec = rot_matrix ## out_center_vec
  
  ;; get current colortable so it can be restored later
  tvlct, r, g, b, /get
  loadct,39

  if not keyword_set(quiet) then begin
     if windowavailable(1) then wset, 1 else window, 1
     surface, dist(5), /nodata, /save, xrange = [-1, 1], yrange = [-1, 1], zrange = [-1, 1], xtitle = 'x', ytitle = 'y'
     plots, out_center_vec[*,0], out_center_vec[*,1], out_center_vec[*,2], psym = 4, color = 200, /T3D
     plots, pix_center_vec[*,0], pix_center_vec[*,1], pix_center_vec[*,2], psym = 4, color = 254, /T3D
     
     plots, new_out_vec[*,0], new_out_vec[*,1], new_out_vec[*,2], psym = 4, color = 100, /T3D
     plots, new_pix_vec[*,0], new_pix_vec[*,1], new_pix_vec[*,2], psym = 4, color = 75, /T3D
  endif

  pred_angle = healpix_rot(new_pix_vec[*,0], new_pix_vec[*,1])
  
  x_rot = new_pix_vec[*,0] * cos(pred_angle) - new_pix_vec[*,1] * sin(pred_angle)
  y_rot = new_pix_vec[*,0] * sin(pred_angle) + new_pix_vec[*,1] * cos(pred_angle)
  x_out_rot = new_out_vec[*,0] * cos(pred_angle) - new_out_vec[*,1] * sin(pred_angle)
  y_out_rot = new_out_vec[*,0] * sin(pred_angle) + new_out_vec[*,1] * cos(pred_angle)
  
  lims = healpix_limits(x_rot, y_rot, x_out_rot, y_out_rot)
  
  wh_inside = where(x_rot ge lims[0] and x_rot le lims[2] and y_rot ge lims[1] and y_rot le lims[3], count_in)
  wh_out_inside = where(x_out_rot ge lims[0] and x_out_rot le lims[2] and y_out_rot ge lims[1] and y_out_rot le lims[3], $
                        count_out)
  
  if count_out ne 0 then begin
     print, 'Amoeba failed to find limits that exclude all outside pixels. Set limits by hand.'
     stop
  endif
  
  if not keyword_set(quiet) then begin
     if windowavailable(2) then wset, 2 else window, 2
     plot, x_rot, y_rot, /nodata, xrange = minmax(x_out_rot), yrange = minmax(y_out_rot), color=0
     oplot, x_out_rot, y_out_rot, psym = 3, color = 100
     oplot, x_rot, y_rot, psym = 3, color = 75
  
     x_range_plot = [lims[0], replicate(lims[2], 2), replicate(lims[0], 2)]
     y_range_plot = [replicate(lims[1], 2), replicate(lims[3], 2), lims[1]]
     oplot, x_range_plot, y_range_plot, psym = -3, color = 0
  endif
  ;; print, n_elements(pixels)
  ;; print, minmax(pixels)
  ;; print, theta0, phi0
  ;; print, vec_mid
  ;; print, pred_angle*180/!pi
  ;; print, lims
  ;; print, count_in
  
  rot_angle = pred_angle
  limits = lims

  new_pix_vec = [[x_rot], [y_rot], [new_pix_vec[*,2]]]
  
  ;; Calculate k step size and range, which are given by spatial resolution & size of field
  ;; Angular resolution is given in Healpix paper in units of arcminutes, need to convert to radians
  ang_resolution = sqrt(3d/!pi) * 3600d/nside * (1d/60d) * (!pi/180d)
  degpix = ang_resolution * 180d / !dpi
  
  x_rad_length = limits[2] - limits[0] + ang_resolution
  y_rad_length = limits[3] - limits[1] + ang_resolution
  
  kxy_rad_range = (2d*!pi) / ang_resolution
  kx_rad_delta = (2d*!pi) / x_rad_length
  ky_rad_delta = (2d*!pi) / y_rad_length
  
  ;; define locations (in k) to take FT
  n_kx = round(kxy_rad_range/kx_rad_delta) + 1
  n_ky = round(kxy_rad_range/ky_rad_delta) + 1
  if (ceil(n_kx/2d)-floor(n_kx/2d)) gt 0 then kx_rad_vals = (dindgen(n_kx)-floor(n_kx/2d)) * kx_rad_delta $
  else kx_rad_vals = (dindgen(n_kx)-n_kx/2+1) * kx_rad_delta
  if (ceil(n_ky/2d)-floor(n_ky/2d)) gt 0 then ky_rad_vals = (dindgen(n_ky)-floor(n_ky/2d)) * ky_rad_delta $
  else ky_rad_vals = (dindgen(n_ky)-n_ky/2+1) * ky_rad_delta
  
 tvlct, r, g, b

end
