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

  n_chunks = ceil(n_slices/float(fchunk))
  fchunk_sizes = intarr(n_chunks) + fchunk
  if n_chunks gt 1 then fchunk_sizes[n_chunks-1] = n_slices - total(fchunk_sizes[0:n_chunks-2])
  fchunk_edges = [0, total(fchunk_sizes, /cumulative)]

  ;; want progress reports every so often + on 3rd step
  nsteps = n_k1*n_chunks
  nprogsteps = 20
  progress_steps = [3, round(nsteps * findgen(nprogsteps) / double(nprogsteps))]
  inner_times = fltarr(nsteps)
  outer_times = fltarr(n_chunks)

  if fchunk_sizes[0] eq 1 then y_exp = reform(exp(-1.*complex(0,1)*y_loc_k), n_pts, n_k2, fchunk_sizes[0]) $
  else y_exp = exp(-1.*complex(0,1)*rebin(y_loc_k, n_pts, n_k2, fchunk_sizes[0]))

  time_preloop = systime(1) - time0
  print, 'pre-loop time: ' + strsplit(string(time_preloop), /extract)
   
  for j=0, n_chunks-1 do begin
     temp=systime(1)
 
     if fchunk_sizes[j] ne fchunk_sizes[0] then $
        if fchunk_sizes[j] eq 1 then y_exp = reform(exp(-1.*complex(0,1)*y_loc_k), n_pts, n_k2, fchunk_sizes[j]) $
        else y_exp = exp(-1.*complex(0,1)*rebin(y_loc_k, n_pts, n_k2, fchunk_sizes[j]))
   
     for i=0, n_k1-1 do begin
        this_step = j*n_k1 + i
        wh = where(progress_steps eq this_step, count)
        if count gt 0 then begin
           print, 'progress: on step ' + number_formatter(this_step) + ' of ' + number_formatter(nsteps) + $
                  ' (~ ' + number_formatter(round(100d*this_step/(nsteps))) + '% done)'
           if this_step gt 0 then begin
              ave_t = mean(inner_times[0:this_step-1])
              t_left = ave_t*(nsteps-this_step)
              if t_left lt 60 then t_left_str = number_formatter(t_left, format='(d8.2)') + ' sec' $
              else if t_left lt 3600 then t_left_str = number_formatter(t_left/60d, format='(d8.2)') + ' min' $
              else t_left_str = number_formatter(t_left/3600d, format='(d8.2)') + ' hours'
              
              print, 'memory used: ' + number_formatter(memory(/current)/1.e9, format='(d8.1)') + ' GB; ave step time: ' + $
                     number_formatter(ave_t, format='(d8.2)') + '; approx. time remaining: ' + t_left_str

           endif
        endif

        temp=systime(1)

        if fchunk_sizes[j] eq 1 then begin
           x_inds = findgen(n_pts) + n_pts*i
           data_inds = findgen(n_pts) + n_pts*fchunk_edges[j]
           term1 = reform(data[data_inds]*exp(-1.*complex(0,1)*x_loc_k[x_inds]), n_pts, fchunk_sizes[j])
           term1_old = reform(data[*,fchunk_edges[j]]*exp(-1.*complex(0,1)*x_loc_k[*,i]), n_pts, fchunk_sizes[j])
        endif else begin
           x_inds = matrix_multiply(findgen(n_pts) + n_pts*i, fltarr(fchunk_sizes[j])+1)
           data_inds = matrix_multiply(findgen(n_pts), fltarr(fchunk_sizes[j])+1) + $
                       n_pts*transpose(matrix_multiply(findgen(fchunk_sizes[j]) + fchunk_edges[j], fltarr(n_pts)+1))
           term1 = data[data_inds] * exp(-1.*complex(0,1)*x_loc_k[x_inds])
           term1_old = data[*,fchunk_edges[j]:fchunk_edges[j+1]-1] * exp(-1.*complex(0,1)*rebin(x_loc_k[*,i], n_pts, fchunk_sizes[j]))
        endelse
        undefine, data_inds, x_inds
        inds = reform(i + n_k1*matrix_multiply(findgen(n_k2), fltarr(fchunk_sizes[j])+1) + n_k1 * n_k2 * $
                      transpose(matrix_multiply(findgen(fchunk_sizes[j]) + fchunk_edges[j], fltarr(n_k2)+1)), n_k2*fchunk_sizes[j])
        ft[inds] = transpose(matrix_multiply(term1, y_exp, /atranspose))
        undefine, inds

        inner_times[this_step] = systime(1) - temp
 
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

