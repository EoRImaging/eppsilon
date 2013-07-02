pro write_ps_fits, fits_savefile, power_3d, weights_3d, noise_expval, kx_mpc, ky_mpc, kz_mpc, kperp_lambda_conv, delay_params, $
                   hubble_param, noise_3d = noise_3d

  ;; make a basic primary header & 2 extension headers (for weights & noise expectation value)
  mkhdr, header, power_3d, /extend
  mkhdr, header2, weights_3d, /image
  mkhdr, header3, noise_expval, /image
  if n_elements(noise_3d) ne 0 then mkhdr, header4, noise_3d, /image

  ;; put them all into a structure to loop over adding parameters
  header_struct = {power:header, weights: header2, noise_expval:header3}
  if n_elements(noise_3d) ne 0 then begin
     header_struct = create_struct(header_struct, 'noise', header4)
     n_hdr = 4
  endif else n_hdr = 3

  ;; add on units parameter
  fxaddpar, header_struct.power, 'BUNIT', 'mK^2 Mpc^3', 'Power. For [mK^2 h^-3 Mpc^3], multiply by ' + number_formatter(hubble_param) + '^3'
  fxaddpar, header_struct.weights, 'BUNIT', 'mK^-4 Mpc^-6', 'Weights (1/variance). For [mK^-4 h^6 Mpc^-6], multiply by ' + $
            number_formatter(hubble_param) + '^-6'
  fxaddpar, header_struct.noise_expval, 'BUNIT', 'mK^2 Mpc^3', 'Noise expectation value. For [mK^2 h^-3 Mpc^3], multiply by ' + $
            number_formatter(hubble_param) + '^3'
  if n_elements(noise_3d) ne 0 then fxaddpar, header_struct.noise, 'BUNIT', 'mK^2 Mpc^3', 'Observed noise. For [mK^2 h^-3 Mpc^3], ' + $
     'multiply by ' + number_formatter(hubble_param) + '^3'

  ;; setup primary WCS coordinates
  axis_names = ['kx', 'ky', 'kz']
  ref_pixs = intarr(3)
  ref_vals = dblarr(3)

  wh_kx0 = where(kx_mpc eq 0, count_kx0)
  if count_kx0 eq 1 then ref_pixs[0] = wh_kx0[0] else if count_kx0 gt 1 then message, 'multiple 0s in kx_mpc'

  wh_ky0 = where(ky_mpc eq 0, count_ky0)
  if count_ky0 eq 1 then ref_pixs[1] = wh_ky0[0] else if count_ky0 gt 1 then message, 'multiple 0s in ky_mpc'

  wh_kz0 = where(kz_mpc eq 0, count_kz0)
  if count_kz0 eq 1 then ref_pixs[2] = wh_kz0[0] else if count_kz0 gt 1 then message, 'multiple 0s in kz_mpc'

  ref_vals = [kx_mpc[ref_pixs[0]], ky_mpc[ref_pixs[1]], kz_mpc[ref_pixs[2]]]
  ref_pixs = ref_pixs+1 ;; 1-indexed

  deltas = [kx_mpc[1]-kx_mpc[0], ky_mpc[1]-ky_mpc[0], kz_mpc[1]-kz_mpc[0]]

  ;; add header parameters
  for j=0, n_hdr-1 do begin
     fxaddpar, header_struct.(j), 'WCSNAME', 'kspace'
     for i=0, 2 do begin
        fxaddpar, header_struct.(j), 'CNAME'+number_formatter(i+1), axis_names[i]
        fxaddpar, header_struct.(j), 'CUNIT'+number_formatter(i+1), 'Mpc^-1', 'For [h Mpc^-1], divide by ' + $
                  number_formatter(hubble_param)
        fxaddpar, header_struct.(j), 'CRPIX'+number_formatter(i+1), ref_pixs[i]
        fxaddpar, header_struct.(j), 'CRVAL'+number_formatter(i+1), ref_vals[i]
        fxaddpar, header_struct.(j), 'CDELT'+number_formatter(i+1), deltas[i]
     endfor
  endfor


  ;; setup secondary WCS coordinates
  axis_names2 = ['u', 'v', 'delay']
  axis_units2 = ['lambda', 'lambda', 's']
  
  ref_pix2 = ref_pixs
  ref_vals2 = ref_vals
  deltas2 = deltas

  ref_vals2[0:1] = ref_vals2[0:1] * kperp_lambda_conv
  deltas2[0:1] = deltas[0:1] * kperp_lambda_conv
  
  if ref_vals[2] ne 0 then stop
  deltas2[2] = delay_params[0]

  ;; add header parameters
  for j=0, n_hdr-1 do begin
     fxaddpar, header_struct.(j), 'WCSNAMEA', 'wavelength_delay'
     for i=0, 2 do begin
        fxaddpar, header_struct.(j), 'CNAME'+number_formatter(i+1)+'A', axis_names2[i]
        fxaddpar, header_struct.(j), 'CUNIT'+number_formatter(i+1)+'A', axis_units2[i]
        fxaddpar, header_struct.(j), 'CRPIX'+number_formatter(i+1)+'A', ref_pix2[i]
        fxaddpar, header_struct.(j), 'CRVAL'+number_formatter(i+1)+'A', ref_vals2[i]
        fxaddpar, header_struct.(j), 'CDELT'+number_formatter(i+1)+'A', deltas2[i]
     endfor
  endfor

  ;; write primary HDU
  writefits, fits_savefile, power_3d, header_struct.power

  ;; now add extension for weights
  writefits, fits_savefile, weights_3d, header_struct.weights, /append
  
  ;; now add extension for noise expectation value
  writefits, fits_savefile, noise_expval, header_struct.noise_expval, /append

  ;; now add extension for noise
  if n_hdr eq 4 then writefits, fits_savefile, noise_3d, header_struct.noise, /append


end
