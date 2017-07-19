pro slurm_ps_job
  ;; wrapper for ps_wrapper to take in shell parameters

  compile_opt strictarr
  args = Command_Line_Args(count=nargs)
  folder_name=args[0]
  obs_range=args[1]
  ;if (nargs eq 3) then n_obs=args[2]
  print,'folder_name = '+folder_name
  print,'obs_range = '+obs_range

  ;;;;;; FT Parameters
  ;freq_ch_range=[ 161, 199 ]   ; 10MHz, z ~ 8.36 to 8.99
;  freq_ch_range=[ 141, 179 ]   ; 10MHz, z ~ 9.83 to 10.37
  freq_ch_range=[ 161, 238 ]   ; 20MHz, z ~ 8.36 to ??
  no_spec_window=1

  ;;;;;; Which plots?:
  plot_stdset=1
  plot_slices=0
;  slice_type='raw'
;  uvf_plot_type='abs'	; These two only used if plot_slices=1
  plot_1to2d=0
  plot_2d_masked=0
  plot_kpar_power=0
  plot_kperp_power=0
  plot_k0_power = 0    ; Causes issues when running with other 1d plots
  plot_noise_1d = 0
  plot_sim_noise = 0
  plot_binning_hist= 0   ;for debugging

  ;;;;;; 2D plotting options
  plot_wedge_line=1
  kperp_linear_axis=0
  kpar_linear_axis=0
;  kperp_plot_range=[5,10]
;  kperp_lambda_plot_range=
;  kpar_plot_range=
  baseline_axis=1
  delay_axis=1
  cable_length_axis=0
;  ref_eor_file="/users/alanman/data/alanman/BubbleCube/TiledHpxCubes/pspecs_21cmFAST/ps_no_halos_z008.24_nf0.550294_useTs1_zetaX2.0e+56_alphaX1.2_TvirminX3.0e+04_aveTb009.92_Pop2_256_300Mpc_v3"

  ;;;;;; 1D plotting options
  set_krange_1dave=0
;  range_1d=
  plot_1d_delta=0
  plot_1d_error_bars=0
  plot_1d_nsigma=2
;  plot_eor_1d=1
  plot_flat_1d=0
  no_text_1d=0
  sim=1
  png=1
;  data_range=[1,3]

  hinv=0

  extra = var_bundle()


  ps_wrapper,folder_name,obs_range,/exact_obsnames,loc_name='oscar', _Extra=extra


;  if (nargs eq 3) then begin
;     n_obs=args[2]
;     ps_wrapper,folder_name,obs_range,n_obs=n_obs,/plot_kpar_power,/plot_kperp_power,/png,/plot_k0_power,/plot_eor_1d,/exact_obsnames,loc_name='oscar'
;  endif else begin
;     ps_wrapper,folder_name,obs_range,/plot_kpar_power,/plot_kperp_power,/png,/plot_k0_power,/plot_eor_1d,/exact_obsnames,loc_name='oscar'
;  endelse
 
end
