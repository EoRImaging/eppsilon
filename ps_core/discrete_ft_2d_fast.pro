; Partially vectorized 2D discrete fourier transform of unevenly spaced data
;    A hybrid of the regular and vector approaches
; locations1 & 2 are x/y values of data points
; k1 & k2 are kx/ky values to test at


function discrete_ft_2d_fast, locations1, locations2, data, k1, k2, max_k_mag = max_k_mag, timing = timing, fchunk = fchunk

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
  nsteps = n_chunks
  nprogsteps = 20
  progress_steps = [3, round(nsteps * findgen(nprogsteps) / double(nprogsteps))]
  inner_times = fltarr(nsteps)
  step1_times = fltarr(nsteps)
  step2_times = fltarr(nsteps)
  step3_times = fltarr(nsteps)

  y_exp = exp(-1.*complex(0,1)*y_loc_k)
  x_exp = exp(-1.*complex(0,1)*x_loc_k)

  x_exp_real = real_part(x_exp)
  x_exp_imag = imaginary(x_exp)
 
  undefine, x_exp

  time_preloop = systime(1) - time0
  print, 'pre-loop time: ' + strsplit(string(time_preloop), /extract)
  
  ;; check memory to reset the highwater mark
  temp = memory(/current)/1.e9

  for j=0, n_chunks-1 do begin
    
     this_step = j
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

           print, 'peak memory used: ' + number_formatter(memory(/highwater)/1.e9, format='(d8.1)') + ' GB; ave step time: ' + $
                  number_formatter(ave_t, format='(d8.2)') + '; approx. time remaining: ' + t_left_str
        endif
     endif

     temp=systime(1)

     if fchunk_sizes[j] eq 1 then begin
        data_inds = dindgen(n_pts) + n_pts*fchunk_edges[j]

        data_expand = rebin(data[temporary(data_inds)], n_pts, n_k1, /sample)

        term1_real = data_expand * x_exp_real
        term1_imag = temporary(data_expand) * x_exp_imag
     endif else begin
        data_inds = rebin(dindgen(n_pts), n_pts, fchunk_sizes[j], /sample) + n_pts * $
                    rebin(reform(findgen(fchunk_sizes[j]) + fchunk_edges[j], 1, fchunk_sizes[j]), n_pts, fchunk_sizes[j], /sample)
     
        data_expand = rebin(reform(data[temporary(data_inds)], n_pts, 1, fchunk_sizes[j]), n_pts, n_k1, fchunk_sizes[j], /sample)

        term1_real = reform(data_expand * rebin(x_exp_real, n_pts, n_k1,fchunk_sizes[j], /sample), n_pts, n_k1*fchunk_sizes[j])
        
        term1_imag = reform(temporary(data_expand) * rebin(x_exp_imag, n_pts, n_k1,fchunk_sizes[j], /sample), $
                            n_pts, n_k1*fchunk_sizes[j])

     endelse
     term1 = temporary(term1_real)+complex(0,1)*temporary(term1_imag)

     temp2 = systime(1)

        ;; get ky vals that are inside max_k_mag
        ;; if n_elements(max_k_mag) gt 0 then begin
        ;;    hist = histogram(sqrt(abs(k1[i]^2. + k2^2.))/max_k_mag, min=0, max=1.5, reverse_indices=ri)
        ;;    count_inside = hist[0]
        ;;    if count_inside gt 0 then y_inds = ri[ri[0]:ri[1]-1] else continue
                   

        ;;    ;; inds_temp = matrix_multiply(dindgen(n_pts), fltarr(count_inside)) + n_pts*matrix_multiply(fltarr(n_pts), double(y_inds))
        ;;    ;; yexp_inds = reform(matrix_multiply(reform(temporary(inds_temp), n_pts*count_inside), $
        ;;    ;;                                    n_pts*count_inside*(dindgen(fchunk_sizes[j]) + fchunk_edges[j])), $
        ;;    ;;                    n_pts,count_inside, fchunk_sizes[j])
        ;;    ;; y_exp_use = y_exp[yexp_inds]

        ;;    y_exp_use = y_exp[*, y_inds, *]
        ;;    inds = i + n_k1*matrix_multiply(double(y_inds), fltarr(fchunk_sizes[j])+1) + n_k1 * n_k2 * $
        ;;           transpose(matrix_multiply(dindgen(fchunk_sizes[j]) + fchunk_edges[j], fltarr(count_inside)+1))
        ;;    ft[inds] = transpose(matrix_multiply(term1, y_exp_use, /atranspose))


        ;; endif else begin

     if fchunk_sizes[j] eq 1 then begin
        inds = dindgen(n_k1*n_k2) + n_k1 * n_k2 * fchunk_edges[j]

        temp3 = systime(1)
        
        ft[temporary(inds)] = matrix_multiply(temporary(term1), y_exp,/atranspose)

     endif else begin
        inds = reform(rebin(dindgen(n_k1*n_k2), n_k1*n_k2, fchunk_sizes[j], /sample) + $
                      n_k1 * n_k2 * rebin(reform(dindgen(fchunk_sizes[j]) + fchunk_edges[j], 1, fchunk_sizes[j]), $
                                          n_k1*n_k2, fchunk_sizes[j], /sample), n_k1, n_k2, fchunk_sizes[j])

        temp3 = systime(1)
        ft[temporary(inds)] = transpose(reform(matrix_multiply(temporary(term1), y_exp,/atranspose), $
                                               n_k1, fchunk_sizes[j], n_k2), [0,2,1])
     endelse

     temp4 = systime(1)

     inner_times[this_step] = temp4 - temp
     step1_times[this_step] = temp2-temp
     step2_times[this_step] = temp3-temp2
     step3_times[this_step] = temp4-temp3
     
  endfor
  
  time1 = systime(1)
  timing = time1-time0
  
  if timing lt 60 then timing_str = number_formatter(timing, format='(d8.2)') + ' sec' $
  else if timing lt 3600 then timing_str = number_formatter(timing/60, format='(d8.2)') + ' min' $
  else timing_str= number_formatter(timing/3600, format='(d8.2)') + ' hours'

  print, 'discrete 2D FT time: ' + timing_str

  return, ft
end

