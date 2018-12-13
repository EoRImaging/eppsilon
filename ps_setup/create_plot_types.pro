function create_plot_types, plot_types = plot_types, $
    plot_stdset = plot_stdset, plot_slices = plot_slices, $
    slice_type = slice_type, uvf_plot_type = uvf_plot_type, plot_1to2d = plot_1to2d, $
    plot_2d_masked = plot_2d_masked, plot_kpar_power = plot_kpar_power, $
    plot_kperp_power = plot_kperp_power, plot_k0_power = plot_k0_power, $
    plot_noise_1d = plot_noise_1d, plot_sim_noise = plot_sim_noise, $
    plot_binning_hist = plot_binning_hist, return_new = return_new

  update_tags = list()
  update_values = list()
  if n_elements(plot_types) eq 0 then begin
    ;; default to not plotting slices
    if n_elements(plot_slices) eq 0 then plot_slices=0

    ;; default to making standard plot set if [plot_slices, plot_1to2d, plot_2d_masked] isn't set
    if keyword_set(plot_slices) or keyword_set(plot_1to2d) or keyword_set(plot_2d_masked) then begin
      if n_elements(plot_stdset) eq 0 then plot_stdset = 0
    endif else begin
      if n_elements(plot_stdset) eq 0 then plot_stdset = 1
    endelse

    ;; default to not plotting 1 to 2d plots, 2d mask plots, kpar power,
    ;; kperp power, kpar=0 power, 1D thermal noise, 1D variance sim noise,
    ;; or binning histogram plots
    if n_elements(plot_1to2d) eq 0 then plot_1to2d=0
    if n_elements(plot_2d_masked) eq 0 then plot_2d_masked=0
    if n_elements(plot_kpar_power) eq 0 then plot_kpar_power=0
    if n_elements(plot_kperp_power) eq 0 then plot_kperp_power=0
    if n_elements(plot_k0_power) eq 0 then plot_k0_power=0
    if n_elements(plot_noise_1d) eq 0 then plot_noise_1d=0
    if n_elements(plot_sim_noise) eq 0 then plot_sim_noise=0
    if n_elements(plot_binning_hist) eq 0 then plot_binning_hist=0

    plot_types = {plot_stdset: plot_stdset, plot_slices: plot_slices, $
      plot_1to2d: plot_1to2d, plot_2d_masked: plot_2d_masked, $
      plot_kpar_power:plot_kpar_power, plot_kperp_power: plot_kperp_power, $
      plot_k0_power: plot_k0_power, plot_noise_1d: plot_noise_1d, $
      plot_sim_noise: plot_sim_noise, plot_binning_hist: plot_binning_hist}
  endif else begin
    if n_elements(plot_stdset) gt 0 then begin
      update_tags.add, 'plot_stdset'
      update_values.add, plot_stdset
    endif
    if n_elements(plot_slices) gt 0 then begin
      update_tags.add, 'plot_slices'
      update_values.add, plot_slices
    endif
    if n_elements(plot_1to2d) gt 0 then begin
      update_tags.add, 'plot_1to2d'
      update_values.add, plot_1to2d
    endif
    if n_elements(plot_2d_masked) gt 0 then begin
      update_tags.add, 'plot_2d_masked'
      update_values.add, plot_2d_masked
    endif
    if n_elements(plot_kpar_power) gt 0 then begin
      update_tags.add, 'plot_kpar_power'
      update_values.add, plot_kpar_power
    endif
    if n_elements(plot_kperp_power) gt 0 then begin
      update_tags.add, 'plot_kperp_power'
      update_values.add, plot_kperp_power
    endif
    if n_elements(plot_k0_power) gt 0 then begin
      update_tags.add, 'plot_k0_power'
      update_values.add, plot_k0_power
    endif
    if n_elements(plot_noise_1d) gt 0 then begin
      update_tags.add, 'plot_noise_1d'
      update_values.add, plot_noise_1d
    endif
    if n_elements(plot_sim_noise) gt 0 then begin
      update_tags.add, 'plot_sim_noise'
      update_values.add, plot_sim_noise
    endif
    if n_elements(plot_binning_hist) gt 0 then begin
      update_tags.add, 'plot_binning_hist'
      update_values.add, plot_binning_hist
    endif
  endelse

  if n_elements(slice_type) gt 0 then begin
    update_tags.add, 'slice_type'
    update_values.add, slice_type
  endif
  if n_elements(uvf_plot_type) gt 0 then begin
    update_tags.add, 'uvf_plot_type'
    update_values.add, uvf_plot_type
  endif

  if n_elements(update_tags) gt 0 or keyword_set(return_new) then begin
    update_tags = update_tags.toarray()

    new_plot_types = update_option_struct(plot_types, update_tags, update_values, $
      return_new = return_new)

    return, new_plot_types
  endif else begin
    return, plot_types
  endelse
end
