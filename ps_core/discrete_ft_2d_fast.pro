; Partially vectorized 2D discrete fourier transform of unevenly spaced data
;    A hybrid of the regular and vector approaches
; locations1 & 2 are x/y values of data points
; k1 & k2 are kx/ky values to test at


function discrete_ft_2d_fast, locations1, locations2, data, k1, k2, timing = timing, mem_param = mem_param

  ;; mem_param is trade-off between speed and memory usage. It
  ;; essentially specifies how many loops to do (0-2). Defaut is 2
  ;; (minimum memory usage)
  if n_elements(mem_param) eq 0 then mem_param = 2
  if (mem_param lt 0) or (mem_param) gt 2 then $
     message, 'mem_param specifies how many loops to use, allowed values are 0-2'

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

  x_exp = exp(-1*dcomplex(0,1)*matrix_multiply(locations1, k1, /btranspose))
  y_exp = exp(-1*dcomplex(0,1)*matrix_multiply(locations2, k2, /btranspose))

  if mem_param lt 2 then begin
     y_exp = rebin_complex(reform(y_exp, n_pts, n_k2, 1), n_pts, n_k2, n_slices)
  endif

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
        t_left = ave_t*(n_k1-i)
        if t_left lt 60 then t_left_str = number_formatter(t_left, format='(d8.2)') + ' sec' $
        else if t_left lt 3600 then t_left_str = number_formatter(t_left/60d, format='(d8.2)') + ' min' $
        else t_left_str = number_formatter(t_left/360d, format='(d8.2)') + ' hours'

        print, 'progress: on loop ' + number_formatter(i) + ' of ' + number_formatter(n_k1) + $
               ' (~ ' + number_formatter(round(100d*i/n_k1)) + '% done)'
        if i gt 0 then print, 'average time per loop: ' + number_formatter(mean(ave_t)) + '; approx. time remaining: ' + t_left_str
     endif

     temp=systime(1)

     if mem_param eq 2 then begin
        for k=0, n_slices-1 do begin
           ft[i,*,k] = matrix_multiply(data[*,k]*x_exp[*,i], y_exp)
        endfor
     endif else begin
        x_exp_loop = rebin_complex(reform(x_exp[*,i], n_pts, 1), n_pts, n_slices)
        ft[i,*,*] = transpose(matrix_multiply(data*x_exp_loop, y_exp, /atranspose))
     endelse

     times[i] = systime(1) - temp
 
  endfor

  time1= systime(1)
  timing = time1-time0
  
  return, ft
end

