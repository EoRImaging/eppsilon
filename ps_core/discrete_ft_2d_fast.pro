; Partially vectorized 2D discrete fourier transform of unevenly spaced data
;    A hybrid of the regular and vector approaches
; locations1 & 2 are x/y values of data points
; k1 & k2 are kx/ky values to test at


function discrete_ft_2d_fast, locations1, locations2, data, k1, k2, timing = timing, fchunk = fchunk

  print, 'Beginning discrete 2D FT'

  time0 = systime(1)

  data_dims = size(data, /dimensions)
  n_dims = n_elements(data_dims)
  if n_dims gt 2 then message,  'Data can have up to 2 dimensions -- 1 for x & y combined and 1 for z'

  n_k1 = n_elements(k1)
  n_k2 = n_elements(k2)
  n_pts = data_dims[0]
  if n_dims eq 1 then n_slices = 1 else n_slices = data_dims[1]

  ;; fchunk is how many frequencies to process simultaneously.
  ;; It allows trade-off between speed and memory usage.
  ;; Defaut is 1 (minimum memory usage), max is number of f channels
  if n_elements(fchunk) eq 0 then fchunk = 1 else fchunk = round(fchunk)
  if (fchunk lt 0) or (fchunk) gt n_slices then $
     message, 'fchunk specifies how many frequencies to process simultaneously. Allowed values are 0-' + number_formatter(n_slices)

  if n_elements(locations1) ne n_pts or n_elements(locations2) ne n_pts then message, $
     'locations1 & 2 must have same number of elements as first dimension of data.'

  ft = complex(fltarr(n_k1, n_k2, n_slices))

  x_loc_k = float(matrix_multiply(locations1, k1, /btranspose))
  y_loc_k = float(matrix_multiply(locations2, k2, /btranspose))

  if fchunk gt 1 then begin
     n_chunks = ceil(n_slices/fchunk)
     fchunk_sizes = intarr(n_chunks) + fchunk
     if n_chunks gt 1 then fchunk_sizes[n_chunks-1] = n_slices - total(fchunk_sizes[0:n_chunks-2])
     fchunk_edges = [0, total(fchunk_sizes, /cumulative)-1]
  endif

  ;;for now require fchunk =1 or nfreq
  ;; if fchunk eq 1 then mem_param=2 else if fchunk eq n_slices then mem_param=1 else stop

  ;; if mem_param eq 2 then begin
  ;;   x_exp = exp(-1*dcomplex(0,1)*x_loc_k)
  ;;   y_exp = exp(-1*dcomplex(0,1)*y_loc_k)
  ;; endif else begin
  ;;    y_exp = exp(-1*dcomplex(0,1)*rebin(y_loc_k, n_pts, n_k2, n_slices))
  ;; endelse

  ;; want progress reports every so often + first 5 steps
  nsteps = n_k1;;*n_chunks
  nprogsteps = 20
  progress_steps = [indgen(5)+1, round(nsteps * findgen(nprogsteps) / double(nprogsteps))]
  times = fltarr(nsteps)

  time_preloop = systime(1) - time0
  print, 'pre-loop time: ' + strsplit(string(time_preloop), /extract)

  ;; for i=0, n_k1-1 do begin
  ;;    ;;generate a vector of x_exp values and a 2d array of y_exp values
  ;;    wh = where(progress_steps eq i, count)
  ;;    if count gt 0 then begin
  ;;       ave_t = mean(times[0:i-1])
  ;;       t_left = ave_t*(n_k1-i)
  ;;       if t_left lt 60 then t_left_str = number_formatter(t_left) + 's' $
  ;;       else if t_left lt 3600 then t_left_str = number_formatter(t_left/60d) + 'm' $
  ;;       else t_left_str = number_formatter(t_left/3600d) + 'h'

  ;;       print, 'progress: on loop ' + number_formatter(i) + ' of ' + number_formatter(n_k1) + $
  ;;              ' (~ ' + number_formatter(round(100d*i/n_k1)) + '% done)'
  ;;       if i gt 0 then print, 'memory used: ' + number_formatter(memory(/current)/1.e9, format='(d8.1)') + ' GB; mean loop time: ' + $
  ;;                                number_formatter(mean(ave_t), format='(d8.2)') + '; approx. time remaining: ' + t_left_str
  ;;    endif

  ;;    temp=systime(1)

  ;;    if mem_param eq 2 then begin
  ;;       for k=0, n_slices-1 do begin
  ;;          ft[i,*,k] = matrix_multiply(data[*,k]*x_exp[*,i], y_exp)
  ;;       endfor
  ;;    endif else begin
  ;;       x_exp_loop = exp(-1*dcomplex(0,1)*rebin(x_loc_k[*,i], n_pts, n_slices))
  ;;       ft[i,*,*] = transpose(matrix_multiply(data*x_exp_loop, y_exp, /atranspose))
  ;;    endelse

  ;;    times[i] = systime(1) - temp
  ;; endfor

  for j=0, n_chunks-1 do begin
     y_exp = exp(-1*dcomplex(0,1)*rebin(y_loc_k, n_pts, n_k2, fchunk_sizes[j]))

     for i=0, n_k1-1 do begin
        ;;generate a vector of x_exp values and a 2d array of y_exp values
        this_step = (i+1)*(j+1)-1
        wh = where(progress_steps eq this_step, count)
        if count gt 0 then begin
           ave_t = mean(times[0:this_step])
           t_left = ave_t*(nsteps-this_step)
           if t_left lt 60 then t_left_str = number_formatter(t_left, format='(d8.2)') + ' sec' $
           else if t_left lt 3600 then t_left_str = number_formatter(t_left/60d, format='(d8.2)') + ' min' $
           else t_left_str = number_formatter(t_left/3600d, format='(d8.2)') + ' hours'
           
           print, 'progress: on step ' + number_formatter(this_step) + ' of ' + number_formatter(nsteps) + $
                  ' (~ ' + number_formatter(round(100d*this_step/(nsteps))) + '% done)'
           if i gt 0 then print, 'memory used: ' + number_formatter(memory(/current)/1.e9, format='(d8.1)') + $
                                 ' GB; approx. time remaining: ' + t_left_str
        endif

        temp=systime(1)

        x_exp_loop = exp(-1*dcomplex(0,1)*rebin(x_loc_k[*,i], n_pts, fchunk_sizes[j]))
        ft[i,*,fchunk_edges[j]:fchunk_edges[j+1]] = transpose(matrix_multiply(data*x_exp_loop, y_exp, /atranspose))

        times[this_step] = systime(1) - temp
 
     endfor
  endfor

  time1 = systime(1)
  timing = time1-time0
  
  if timing lt 60 then timing_str = number_formatter(timing, format='(d8.2)') + ' sec' $
  else if timing lt 3600 then timing_str = number_formatter(timing/60, format='(d8.2)') + ' min' $
  else timing_str= number_formatter(timing/3600, format='(d8.2)') + ' hours'

  print, 'discrete 2D FT time: ' + timing_str

  return, ft
end

