function choose_pix_ft, file_struct, pixel_nums = pixel_nums, data_dims = data_dims, $
    uvf_options = uvf_options

  if tag_exist(file_struct, 'nside') ne 0 then healpix = 1 else healpix = 0

  if healpix then begin
    if n_elements(pixel_nums) eq 0 then message, 'If cubes are Healpix pixel_nums must be passed'

    ;; get pixel vectors
    pix2vec_ring, file_struct.nside, pixel_nums, pix_center_vec
    ;; find mid point (work in x/y because of possible jumps in phi)
    vec_mid = [mean(pix_center_vec[*,0]), mean(pix_center_vec[*,1]), mean(pix_center_vec[*,2])]
    theta0 = acos(vec_mid[2])
    phi0 = atan(vec_mid[1], vec_mid[0])

    ;; To go to flat sky, rotate patch to zenith and flatten.
    ;; To get to current location, need to first rotate around z by
    ;; phi, then around y by -theta, then around z by -phi
    ;; use inverse to rotate back to zenith
    rot_matrix = get_rot_matrix(theta0, phi0, /inverse)
    new_pix_vec = rot_matrix ## pix_center_vec

    ;; then rotate to make as rectangular as possible
    pred_angle = healpix_rot(new_pix_vec[*,0], new_pix_vec[*,1])

    x_rot = new_pix_vec[*,0] * cos(pred_angle) - new_pix_vec[*,1] * sin(pred_angle)
    y_rot = new_pix_vec[*,0] * sin(pred_angle) + new_pix_vec[*,1] * cos(pred_angle)

  endif else begin
    if n_elements(data_dims) eq 0 then message, 'If cubes are not Healpix data_dims must be passed'

    ;; gridded image to dft to parallel Healpix computation
    pix_size_rad = abs(file_struct.degpix) * !pi / 180d

    x_vec = (findgen(data_dims[0]) - data_dims[0]/2.) * pix_size_rad
    y_vec = (findgen(data_dims[1]) - data_dims[1]/2.) * pix_size_rad

    x_rot = fltarr(data_dims[0]*data_dims[1])
    y_rot = fltarr(data_dims[0]*data_dims[1])
    x_rot = reform(rebin(x_vec, data_dims[0], data_dims[1], /sample), data_dims[0]*data_dims[1])
    y_rot = reform(rebin(reform(y_vec, 1, data_dims[1]), data_dims[0], data_dims[1], /sample), data_dims[0]*data_dims[1])

  endelse

  ;; get size of a square region that fits in the Healpix image
  if healpix then begin
    ;; get surrounding pixels
    dists = sqrt((pix_center_vec[*,0]-vec_mid[0])^2d + (pix_center_vec[*,1]-vec_mid[1])^2d + (pix_center_vec[*,2]-vec_mid[2])^2d)
    radius = max(dists)*1.1
    query_disc, file_struct.nside, vec_mid, radius, listpix, nlist, /inc
    min_pix = min([pixel_nums, listpix])
    wh2 = where(histogram(listpix, min=min_pix) eq 0 and histogram(pixel_nums, min=min_pix) gt 0, count2)
    while count2 gt 0 do begin
      radius = radius*1.1
      query_disc, file_struct.nside, vec_mid, radius, listpix, nlist, /inc
      min_pix = min([pixel_nums, listpix])
      wh2 = where(histogram(listpix, min=min_pix) eq 0 and histogram(pixel_nums, min=min_pix) gt 0, count2)
    endwhile

    ;; remove pixels from listpix that are in my image -- only want nearby pixels not in my image
    min_pix = min([pixel_nums, listpix])
    max_pix = max([pixel_nums, listpix])
    wh = where(histogram(listpix, min = min_pix, max = max_pix) gt 0 and histogram(pixel_nums, min = min_pix, max = max_pix) eq 0, count_out)
    if count_out gt 0 then outside_pix = wh + min_pix else $
      message, 'Something has gone wrong with finding excluded Healpix pixels in region of interest'
    pix2vec_ring, file_struct.nside, outside_pix, out_center_vec
    new_out_vec = rot_matrix ## out_center_vec
    x_out_rot = new_out_vec[*,0] * cos(pred_angle) - new_out_vec[*,1] * sin(pred_angle)
    y_out_rot = new_out_vec[*,0] * sin(pred_angle) + new_out_vec[*,1] * cos(pred_angle)
    limits = healpix_limits(x_rot, y_rot, x_out_rot, y_out_rot)

    ;; limits should not extend beyond horizon (1/sqrt(2))
    if limits[0] lt -1/sqrt(2) then limits[0] = -1/sqrt(2)
    if limits[1] lt -1/sqrt(2) then limits[1] = -1/sqrt(2)
    if limits[2] gt 1/sqrt(2) then limits[2] = 1/sqrt(2)
    if limits[3] gt 1/sqrt(2) then limits[3] = 1/sqrt(2)

  endif else begin
    limits = [min(x_vec), min(y_vec), max(x_vec), max(y_vec)]
  endelse
  image_len = min([limits[2]-limits[0],limits[3]-limits[1]])

  frequencies = file_struct.frequencies
  if n_elements(freq_ch_range) ne 0 then begin
    frequencies = frequencies[min(freq_ch_range):max(freq_ch_range)]
  endif

  ;; figure out k values to calculate dft
  uv_cellsize_m = 5 ;; based on calculations of beam FWHM by Aaron
  if tag_exist(uvf_options, 'image_window_name') then begin
    ;; if we have a image window we want to go out wider than normal to accommodate the window.
    min_uv_cellsize_m = max([file_struct.kpix * (3e8) / mean(frequencies*1e6), (3e8) / (image_len * mean(frequencies*1e6))])
    std_fraction = min_uv_cellsize_m / uv_cellsize_m

    if std_fraction gt 0.9 then print, 'The Healpix image cubes do not extend much beyond the standard window size.'

    if not tag_exist(uvf_options, 'image_window_frac_size') then begin
      uvf_options = create_uvf_options(uvf_options = uvf_options, $
        image_window_frac_size = std_fraction)
      if std_fraction gt 0.9 then print, 'The calculated window fractional size is ' + number_formatter(std_fraction)
    endif
    uv_cellsize_m = min_uv_cellsize_m
  endif


  if tag_exist(uvf_options, 'delta_uv_lambda') then begin
    delta_kperp_rad = uvf_options.delta_uv_lambda * (2.*!pi)
  endif else begin
    delta_kperp_rad = uv_cellsize_m * mean(file_struct.frequencies*1e6) * (2.*!pi) / (3e8)
  endelse
  ;; go a little beyond max_baseline to account for expansion due to w projection
  ;; max_kperp_rad = (file_struct.max_baseline_lambda/kperp_lambda_conv) * z_mpc_mean * 1.1
  ;; use kspan of Ian's cubes
  if tag_exist(file_struct, 'kspan') then begin
    max_kperp_rad = min([file_struct.kspan/2.,file_struct.max_baseline_lambda])* (2.*!pi)
  endif else max_kperp_rad = min([file_struct.max_baseline_lambda])* (2.*!pi)

  if tag_exist(uvf_options, 'max_uv_lambda') then begin
    max_kperp_rad = min([max_kperp_rad, uvf_options.max_uv_lambda * (2.*!pi)])
  endif

  ;; limit field of view to match calculated k-modes
  xy_len = 2*!pi/delta_kperp_rad

  ;; image may be smaller than expected, may need to adjust delta_kperp_rad
  if image_len lt xy_len then begin
    print, 'Image FoV is smaller than expected, increasing delta kperp to match image FoV'
    delta_kperp_rad = 2*!pi/image_len

    x_range = [-1,1]*image_len/2. + mean(x_rot)
    y_range = [-1,1]*image_len/2. + mean(y_rot)

  endif else begin
    x_range = [-1,1]*xy_len/2. + mean(x_rot)
    y_range = [-1,1]*xy_len/2. + mean(y_rot)
  endelse

  wh_close = where(x_rot le x_range[1] and x_rot ge x_range[0] and $
    y_rot le y_range[1] and y_rot ge y_range[0], count_close, $
    ncomplement = count_far)

  if n_elements(wh_close) eq 0 and count_far eq 0 then wh_close = lindgen(n_elements(x_rot))

  n_kperp = round(max_kperp_rad / delta_kperp_rad) * 2 + 1
  kx_rad_vals = (findgen(n_kperp) - (n_kperp-1)/2) * delta_kperp_rad

  ;; need to cut uvf cubes in half because image is real -- we'll cut in v
  ;; drop the unused half before the DFT to save time
  ky_rad_vals = kx_rad_vals[n_kperp/2:n_kperp-1]

  ret_struct = {wh_close: wh_close, x_use: x_rot[wh_close], y_use: y_rot[wh_close], $
    kx_rad_vals: kx_rad_vals, ky_rad_vals: ky_rad_vals, delta_kperp_rad: delta_kperp_rad, n_kperp: n_kperp}

  return, ret_struct

end
