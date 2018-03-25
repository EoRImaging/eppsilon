pro mit_ps_job
  ;wrapper for ps_wrapper to take in shell parameters from a grid engine bash script

  ;***Read-in typical parameters
  compile_opt strictarr
  args = Command_Line_Args(count=nargs)
  folder_name=args[0]
  obs_range=args[1]
  
  print,'folder_name = '+folder_name
  print,'obs_range = '+obs_range
  ;***
  
  ;***If only the typical parameters are set, run an evenodd/pol continuous lightweight wrapper. Will contain
  ;   significant aliasing due to a harsh window cut, but will be the fastest continuous option
  if (nargs EQ 2) then ps_wrapper,folder_name,obs_range,/plot_kpar_power,/plot_kperp_power,/png,/plot_k0_power,$
    /exact_obsnames,loc_name='aws'
  ;***
    
  ;***If the typical parameters plus the window filter are set, run an evenodd/pol continuous heavyweight wrapper
  ;   Will not contain significant aliasing due to a harsh window cut, but will take over 40 hrs
  if (nargs eq 3) then begin
    image_window_name=args[2]
    print,'image_window_name = '+image_window_name
    ps_wrapper,folder_name,obs_range,/plot_kpar_power,/plot_kperp_power,/png,/plot_k0_power,/exact_obsnames, $
      loc_name='aws',image_window_name=image_window_name
  endif
  ;***
  
  ;***If the typical parameters plus parallelization parameters are set, run a paralleled lightweight
  ;   wrapper. Significant aliasing, but by far the fastest option (~30 minutes).
  if (nargs eq 5) then begin
    cube_type=args[2]
    pol=args[3]
    evenodd=args[4]
    print,'cube_type = '+cube_type
    print,'pol = '+pol
    print,'evenodd = '+evenodd
    single_cube_dft,folder_name,obs_range,/exact_obsnames,loc_name='aws',$
      cube_type=cube_type,pol=pol,evenodd=evenodd
  endif
  ;***
  
  ;***If the full suite of parameters are set, run a paralleled heavyweight wrapper. This is the fastest
  ;   way to run without significant aliasing.
  if (nargs eq 6) then begin
    cube_type=args[2]
    pol=args[3]
    evenodd=args[4]
    image_window_name=args[5]
    print,'cube_type = '+cube_type
    print,'pol = '+pol
    print,'evenodd = '+evenodd
    print,'image_window_name = '+image_window_name
    single_cube_dft,folder_name,obs_range,/exact_obsnames,loc_name='aws',$
      cube_type=cube_type,pol=pol,evenodd=evenodd,image_window_name=image_window_name
  endif
  ;***
  
end
