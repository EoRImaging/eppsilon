pro ps_ratio_wrapper, folder_names_in, obs_names_in, ps_foldernames=ps_foldernames, $
    exact_obsnames = exact_obsnames, version_test = version_test, $
    cube_types = cube_types,  pols = pols, $
    all_pol_diff_ratio = all_pol_diff_ratio, freq_ch_range = freq_ch_range, $
    freq_flags = freq_flags, freq_flag_name = freq_flag_name, $
    freq_flag_repeat = freq_flag_repeat, $
    freq_avg_factor = freq_avg_factor, force_even_freqs = force_even_freqs, $
    spec_window_types = spec_window_types, delta_uv_lambda = delta_uv_lambda, $
    full_image = full_image, image_clip = image_clip, $
    ave_removal = ave_removal, std_power = std_power, diff_ratio = diff_ratio, $
    diff_range = diff_range, png = png, eps = eps, pdf = pdf, $
    data_range = data_range, $
    color_type = color_type, invert_colorbar = invert_colorbar, $
    kperp_linear_axis = kperp_linear_axis, $
    kpar_linear_axis = kpar_linear_axis, sim = sim, wt_cutoffs = wt_cutoffs, $
    wt_measures = wt_measures, window_num = window_num, $
    uvf_input = uvf_input, diff_save_path = diff_save_path, plot_path = diff_plot_path

  compare_setup_structures, folder_names_in, obs_names_in, $
    ps_foldernames = ps_foldernames, version_test = version_test, $
    spec_window_types = spec_window_types, delta_uv_lambda = delta_uv_lambda, $
    max_uv_lambda = max_uv_lambda, full_image = full_image, $
    image_clip = image_clip, ave_removal = ave_removal, $
    freq_dft = freq_dft, dft_z_use = dft_z_use, std_power = std_power, $
    all_type_pol = all_type_pol, freq_ch_range = freq_ch_range, $
    freq_flags = freq_flags, freq_flag_name = freq_flag_name, $
    freq_flag_repeat = freq_flag_repeat, $
    freq_avg_factor = freq_avg_factor, force_even_freqs = force_even_freqs, $
    diff_plot_path = diff_plot_path, diff_save_path = diff_save_path, $
    folder_names = folder_names, obs_info = obs_info, uvf_options = uvf_options, $
    freq_options = freq_options, ps_options = ps_options

  plot_options = create_plot_options(plot_path = diff_plot_path, $
    png = png, eps = eps, pdf = pdf)

  plot_2d_options = create_plot_2d_options(kperp_linear_axis = kperp_linear_axis, $
    kpar_linear_axis = kpar_linear_axis, data_range = data_range, color_type = color_type)

  ps_ratio_plots, folder_names, obs_info, cube_types, ps_foldernames=ps_foldernames, $
    pols, all_pol_diff_ratio = all_pol_diff_ratio, $
    uvf_options = uvf_options, freq_options = freq_options, ps_options = ps_options, $
    plot_options = plot_options, plot_2d_options = plot_2d_options, $
    save_path = diff_save_path, plot_filebase = plot_filebase, $
    diff_ratio = diff_ratio, diff_range = diff_range, invert_colorbar = invert_colorbar, $
    quiet = quiet, window_num = window_num
end
