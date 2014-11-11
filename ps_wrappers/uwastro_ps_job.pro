pro uwastro_ps_job
  ;; wrapper for uwastro_wrapper to take in shell parameters

  compile_opt strictarr
  args = Command_Line_Args(count=nargs)
  folder_name=args[0]
  obs_range=args[1]
  uwastro_wrapper,folder_name,obs_range,/png
 
  
end
