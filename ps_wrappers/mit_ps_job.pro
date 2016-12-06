pro mit_ps_job
  ;; wrapper for mit_wrapper to take in shell parameters

  compile_opt strictarr
  args = Command_Line_Args(count=nargs)
  folder_name=args[0]
  obs_range=args[1]
  filter_name=args[2]
  ;delta_uv_lambda =args[3]
  ;if (nargs eq 3) then n_obs=args[2]
  print,'folder_name = '+folder_name
  print,'obs_range = '+obs_range
  print,'filter name = '+filter_name
  ;print,'delta_uv_lambda = '+delta_uv_lambda
  
  ;if (nargs eq 3) then begin
  ;   n_obs=args[2]
  ;   ps_wrapper,folder_name,obs_range,n_obs=n_obs,/plot_kpar_power,/plot_kperp_power,/png,/plot_k0_power,/exact_obsnames,loc_name='mit'
  ;endif else begin
  ps_wrapper,folder_name,obs_range,/plot_kpar_power,/plot_kperp_power,/png,/plot_k0_power,/exact_obsnames,loc_name='mit', filter_name=filter_name;, delta_uv_lambda=delta_uv_lambda
;endelse
  
end
