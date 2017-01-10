pro mit_ps_job
  ;; wrapper for mit_wrapper to take in shell parameters

  compile_opt strictarr
  args = Command_Line_Args(count=nargs)
  folder_name=args[0]
  obs_range=args[1]
  
  print,'folder_name = '+folder_name
  print,'obs_range = '+obs_range
  
  if (nargs eq 5) then begin
    cube_type=args[2]
    pol=args[3]
    evenodd=args[4]
    print,'cube_type = '+cube_type
    print,'pol = '+pol
    print,'evenodd = '+evenodd
    single_cube_dft,folder_name,obs_range,/exact_obsnames,loc_name='mit',$
      cube_type=cube_type,pol=pol,evenodd=evenodd
  endif else begin
    ps_wrapper,folder_name,obs_range,/plot_kpar_power,/plot_kperp_power,/png,/plot_k0_power,/exact_obsnames,loc_name='mit'
  endelse
  
end
