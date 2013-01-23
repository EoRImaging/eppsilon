; Partially vectorized 2D discrete fourier transform of unevenly spaced data
;    A hybrid of the regular and vector approaches
; locations1 & 2 are x/y values of data points
; k1 & k2 are kx/ky values to test at


function discrete_ft_2d_fast, locations1, locations2, data, k1, k2, timing = timing

  time0 = systime(1)

  data_dims = size(data, /dimensions)
  n_dims = n_elements(data_dims)
  if n_dims gt 2 then message,  'Data can have up to 2 dimensions -- 1 for x & y combined and 1 for z'

  n_k1 = n_elements(k1)
  n_k2 = n_elements(k2)
  n_pts = data_dims[0]
  if n_dims eq 1 then n_slices = 1 else n_slices = data_dims[1]

  if n_elements(locations1) ne n_pts or n_elements(locations2) ne n_pts then message, $
     'locations1 & 2 must have same number of elements as first dimension of data.'

  ft = dcomplex(dblarr(n_k1, n_k2, n_slices))
  y_exp = exp(-1*dcomplex(0,1)*rebin(reform(k2, 1, n_k2), n_pts, n_k2)*rebin(locations2, n_pts, n_k2))
  x_exp = exp(-1*dcomplex(0,1)*rebin(reform(k1, 1, n_k1), n_pts, n_k1)*rebin(locations1, n_pts, n_k1))
  ;; y_exp = exp(-1*dcomplex(0,1)*rebin(rebin(reform(k2, 1, n_k2), n_pts, n_k2)*rebin(locations2, n_pts, n_k2), n_pts, n_k2, n_slices))
  ;; x_exp = exp(-1*dcomplex(0,1)*rebin(rebin(reform(k1, 1, n_k1), n_pts, n_k1)*rebin(locations1, n_pts, n_k1), n_pts, n_k1, n_slices))

  ;; want progress reports every so often + first 5 steps
  nprogsteps = 20
  progress_steps = [indgen(5)+1, round(n_k1 * findgen(nprogsteps) / double(nprogsteps))]
  times = dblarr(n_k1)

  time_preloop = systime(1) - time0
  print, 'pre-loop time: ' + strsplit(string(time_preloop), /extract)

  for i=0, n_k1-1 do begin
     ;;generate a vector of x_exp values and a 2d array of y_exp values
     wh = where(progress_steps eq i, count)
     if count gt 0 then begin
        ave_t = mean(times[0:i-1])
        print, 'progress: on loop ' + strsplit(string(i), /extract) + ' of ' + strsplit(string(n_k1), /extract) + $
               ' (~ ' + strsplit(string(round(100d*i/n_k1)), /extract) + '% done) approx. time remaining: ' + $
               strsplit(string(ave_t*(n_k1-i)))
        if i gt 0 then print, 'average time per loop: ' + strsplit(string(mean(ave_t)), /extract)
     endif

     temp=systime(1)

     for k=0, n_slices-1 do begin
        ft[i,*,k] = (data[*,k]*x_exp[*,i]) # y_exp
     endfor

     ;;ft[i,*,*] = (data*x_exp[*,i,*]) # y_exp

     times[i] = systime(1) - temp
     ;; print, 'average time per loop: ' + strsplit(string(mean(times[0:i])), /extract)
     ;; time_complete = mean(times[0:i])*n_k1 + time_preloop
     ;; print, 'projected completion time (hours): ' + strsplit(string(time_complete/3600d), /extract)
     ;; stop
  endfor

  time1= systime(1)
  timing = time1-time0
  
  return, ft
end

