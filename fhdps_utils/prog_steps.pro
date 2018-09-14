;; this_step: 0-based iteration
;; nsteps: number of iterations to be completed
;; loop_times: vector of length nsteps with values filled in up to this_step
;; nprogsteps: optional, number of evenly spaced progress steps to report (default = 5)
;; progress_steps: optional, list of step indices to report progress on
;;   (default = [nprogsteps evenly spaced values starting at 0])
pro prog_steps, this_step, nsteps, loop_times, nprogsteps=nprogsteps, $
  progress_steps=progress_steps

  if n_elements(progress_steps) eq 0 then begin
    ;; get progress reports periodically + on 3rd step (for early estimate for time to finish)
    if n_elements(nprogsteps) eq 0 then begin
      ;; default to 20%
      nprogsteps = 5
    endif
    progress_steps = round(nsteps * findgen(nprogsteps) / double(nprogsteps))
  endif

  wh = where(progress_steps eq this_step, count)
  if count gt 0 then begin
     pct_done = round(100d*(this_step+1)/(nsteps))

     print, 'progress: finished step ' + number_formatter(this_step+1) + ' of ' + number_formatter(nsteps) + $
            ' (~ ' + number_formatter(pct_done) + '% done)'

      ave_t = mean(loop_times[0:this_step])
      t_left = ave_t*(nsteps-this_step-1)
      if t_left lt 60 then t_left_str = number_formatter(t_left, format='(d8.2)') + ' sec' $
      else if t_left lt 3600 then t_left_str = number_formatter(t_left/60d, format='(d8.2)') + ' min' $
      else t_left_str = number_formatter(t_left/3600d, format='(d8.2)') + ' hours'

      print, 'peak memory used: ' + number_formatter(memory(/highwater)/1.e9, format='(d8.1)') + ' GB; ave step time: ' + $
             number_formatter(ave_t, format='(d8.2)') + '; approx. time remaining: ' + t_left_str
  endif

end
