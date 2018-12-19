# eppsilon Keyword Dictionary
eppsilon uses keywords to create unique run-specific settings. This dictionary describes the purpose of each keyword, as well as their logic or applicable ranges.

This is a work in progress; please add keywords as you find them with their corresponding definition. If there's a question about a definition or keyword, label it with !Q.

# ps_wrapper keywords
Many of these keywords are re-used in other wrappers, but since this is the primary eppsilon wrapper they are defined here.

## input data definitions:

**folder_name**: **Required** This defines what folder the data live in. For FHD, it should be the top level folder for that run (which contains folders like Metadata and Healpix). On machines where there is support for the standard folder paths (i.e. there is a relevant *_folder_locs.pro file), this can just be the short folder name. It can also be a full path.

**loc_name**: This is a keyword to indicate which location the code is running on to help with defining the standard folder paths. It is detected from the hostname by default, but that apparently fails in some cases, so this keyword allows the user to specify it. Options are 'enterprise' (enterprise really means any of the ASU machines.)

**obs_name**: This defines which sets of cubes to use within the specified folder_name. It is not required, but is useful if there are more than one run present in the folder_name. It needs to be a string that uniquely identifies the run from all other runs in the folder (For FHD, this matches the obs_name defined in FHD/fhd_core/HEALPix/integrate_healpix_cubes.pro)

**data_subdirs**: This defines the subdirectory (within the folder_name) where the cubes are located. For FHD, if uvf_input is not set, this defaults to 'Healpix' + path_sep(), otherwise it defaults to ''

**exact_obsnames**: This is a flag (valid values are 0/1, default=0) indicating that the obs_name that was provided was the exact obs_name, rather than a unique string that should be in the obs_name.

**beamfiles**: This is a keyword giving the location and names of the beam files. For FHD inputs, it is only used by the code if the beamfiles are not found in their standard location (it is actually overwritten by those locations if they are found.)

**pol_inc**: A list of which polarizations to calculate power spectra for (e.g. 'xx', 'yy'). For FHD the default is all available polarizations, for RTS it is ['xx', 'yy'].

**type_inc**: A list of which cube types to calculate power spectra for (e.g. 'dirty', 'model', 'res'). For FHD the default is all available types, for RTS it is ['res'].

**freq_ch_range**: Specifies what range of frequency channels to calculate the power spectra for. Default is all available channels (this is modified if the coarse_harm_width keyword is used).

**freq_flags**: A list of frequency channels to flag in calcuating the power spectrum. Errors will be generated if there are fewer than 3 unflagged frequencies.

**freq_flag_name**: Only used if freq_flags is used. String to use in the output file names to identify the files with the frequency flagging applied.

**uvf_input**: This is a flag (valid values are 0/1, default=0) indicating that the input cubes are uvf cubes, rather than image space cubes. This is only supported for FHD inputs and only works for single obsid cubes. This is used most often for simulation testing & exploration.

**no_evenodd**: This is a flag (valid values are 0/1, default=0) indicating that only one
set of files (rather than both evens and odds) are expected to be present. If not set,
eppsilon will error fairly quickly rather than proceeding with just the one set of files.


## RTS specific keywords:

**rts**: This is a flag (valid values are 0/1, default=0) indicating that we are using RTS inputs.

**no_wtvar_rts**: This is a flag (valid values are 0/1, default=0) indicating that RTS weight cubes are not present.

**norm_rts_with_fhd**: This is a flag (valid values are 0/1, default=0) indicating that the vis_noise value should be set to the value usually calculated by FHD rather than the RTS value. This usually brings the normalization of the power spectra from RTS runs in line with the normalization of FHD runs


## CASA specific keywords:

**casa**: This is a flag (valid values are 0/1, default=0) indicating that we are using CASA inputs. This mode has not been used in a long time and it may be out of date.


## Simulation specific keywords:

**sim**: This is a flag (valid values are 0/1, default=0) indicating that the inputs are simulations, which results in minor plotting changes that are better suited to simulations.

**sim_use_weight_cutoff**: This is a flag (valid values are 0/1, default=1) that is only used if the sim keyword is set and it causes the code to use the weights cutoff with density correction in the simulation (it does both with and without). Set it to 0 to causes the code to do no cutoff or density correction.

**fix_sim_input**: This is a backwards compatability flag (valid values are 0/1, default=0) to fix some errors with very old simulations where the normalization was off by a factor of 2 because of a mistake in the simulation code. It should only be used with care if working with old simulations. It will soon be deprecated.


## Refresh keywords:
These flags (all default=0) tell the code to redo parts of the analysis that would normally be skipped because the outputs from that part of the code already exist. They cascade, so if you set an early one it will cause all affected downstream code to also re-run.

**refresh_dft**: Recalculate the DFT. Should only be used if the input cube files have changed. This is the longest step computationally.

**refresh_beam**: Redo the beam related calculations. This causes the power spectrum to be recalculated.

**refresh_ps**: Recalculate the power spectrum. This redoes everything downstream of the DFT and beam calculations.

**refresh_binning**: Redo the 1D & 2D binning.

**refresh_info**: Redo getting all the metadata together. This generally does not cause any major recalculations and is only really useful if folder structures have changed.


## DFT and UV grid related keywords:

**delta_uv_lambda**: The default size of the uv pixels to calculate in wavelengths. If full_image is not set, this is defaulted to (5 meters)*freq/c, otherwise it is decreased to match the size of the input image cube. **Note**: the actual uv pixel size used by the code may be larger than this value if the size of the input image cube is smaller than would be needed to support delta_uv_lambda. In that case, the following warning will be printed: 'Image FoV is smaller than expected, increasing delta kperp to match image FoV'.

**max_uv_lambda**: The maximum distance from (0,0) in uv space that we should calculate in the DFT. **Note**: The actual uv extent used may be smaller than this value if the maximum uv extent of the input image cubes is smaller or if the longest gridded baseline is smaller. So the actual uv extent will be the minimum of max_uv_lambda, the maximum uv extent of the input image cubes and the longest gridded baseline.

**full_image**: This is a flag (valid values are 0/1, default=0) indicating that the full image should be used and that the uv spacing should be set based on the image size. This keyword is only needed if no image window function is applied. This keyword cannot be set if **delta_uv_lambda** is set.

**uv_avg**: A rarely used testing/exploration keyword that is only used if the uvf_input flag is set. This is a factor by which to average up pixels in the uv plane. A value of 2 results in pixels that are larger by a factor of 2 in each direction and 4x fewer total pixels.

**uv_img_clip**: A rarely used testing/exploration keyword that is only used if the uvf_input flag is set. This is a factor by which to clip the cubes in image space (i.e. 2D FFT to image space, clip, then 2D FFT back to uv). A value of 2 results in images that are smaller by a factor of 2 in each direction, 4x smaller overall, and uv pixels that are larger by a factor of 2 in each direction.

**dft_fchunk**: A very rarely used keyword, this tells the code how many frequencies to calculate at once. It used to be useful for memory management, but rewrites of the DFT code have made this largely obsolete.

**no_dft_progress**: This is a flag (valid values are 0/1, default=0) to suppress the printing of progress report statements during the DFT.


## Power spectrum calculation keywords:

**ave_removal**: This is a flag (valid values are 0/1, default=1) indicating that the average along the frequency axis of the cubes should be removed before the frequency Fourier Transform and then added back in to the k_parallel=0 mode. This has been well tested and substantially reduces the bleed of the flat spectrum foregrounds into the window.

**wt_cutoffs**: This is a keyword to support a very subtle power spectrum normalization issue identified by Adrian Liu and quantified for MWA with simulations. It turns out that when baselines substantially overlap there is a suppression of the power of stochastic fields (like the EoR). We measure this overlap using the weights cube and we found through simulations that for MWA antennas, above a certain weight value the power suppression asymptotes to a factor of 0.5. A conservative weight value to use to be sure we are in that regime is a weight value of 1, and if we throw away bins below that value we lose very little sensitivity with the MWA. The wt_cutoffs keyword is the weights value to use for the cutoff value (it is calculated in uvf space, so it applies to all kz values for a given kx,ky pixel). It is defaulted to 1.0 and is only used during the binnning to 1D or 2D power spectra (for 2D, it is used to normalize areas of the ps properly but lower values are not thrown out, for 1D, pixels with lower values are excluded and the power spectrum is normalized by a factor of 2). To turn this normalization correction off, set wt_cutoffs=0.

**wt_measures**: See the description for wt_cutoffs to get the full context. This keyword defines how the weight value for the power normalization is quantified in the uv plane given a uvf cube. Options are 'min' (minimum value along frequency axis) and 'ave' (average value along frequency axis), default is 'min'. We don't see much change between 'min' and 'ave', but 'min' is more conservative.

**spec_window_type**: Name of spectral window to apply in the frequency direction before Fourier Transforming. Default is 'Blackman-Harris'. Full set of options is: ['Hann', 'Hamming', 'Blackman', 'Nutall', 'Blackman-Nutall', 'Blackman-Harris', 'Blackman-Harris^2'].

**no_spec_window**: This is a flag (valid values are 0/1, default=0) indicating that no spectral window should be applied (see spec_window_type keyword).

**allow_beam_approx**: This is a flag (valid values are 0/1, default=0) indicating that a less good beam integral approximation should be used if the beam_integral is not present or is zero in the obs structure. If this keyword is not set and the beam_integral is not present or is zero, eppsilon will fail fairly quickly. If this keyword is set and the beam_integral is not present or is zero, eppsilon will print a warning but continue using a less good approximation.

**std_power**: A rarely used testing/exploration flag (valid values are 0/1, default=0) indicating that the power should be calculated using the sine and cosine terms of the FFT along the frequency axis, rather than the Lomb-Scargle components.

**no_wtd_avg**: A rarely used testing/exploration flag (valid values are 0/1, default=0) indicating that the power should be calculated as the sum of the square of the two Fourier Transform components (either Lomb-Scargle components or sine/cosine components if the std_power flag is set), rather than as a variance weighted sum of the square of the components.

**inverse_covar_weight**: *This is very much under development and is not yet trustworthy.* A rarely used testing/exploration flag (valid values are 0/1, default=0) that constructs a frequency covariance matrix assuming flat spectrum foregrounds, transforms it to k_parallel and uses it to weight the power spectra.


## Power spectrum binning keywords:

**no_kzero**: This is a flag (valid values are 0/1, default=0) indicating that the k_parallel=0 mode should be removed before binning. This is a rarely used keyword because we typically use other keywords to control it.

### 2D binning keywords:

**log_kpar**: This is a flag (valid values are 0/1, default=0) indicating that logarithmic bins should be used in the k_parallel direction for 2D power spectra.

**log_kperp**: This is a flag (valid values are 0/1, default=0) indicating that logarithmic bins should be used in the k_perpendicular direction for 2D power spectra.

**kpar_bin**: This sets the size of the bins for 2D power spectra in the k_parallel direction. For linear binning, this is the actual binsize, for log binning, this gives the bins per decade as 1/kpar_bin. The default is the kz pixel size for linear binning and 0.1 for log binning (giving 10 bins per decade).

**kperp_bin**: This sets the size of the bins for 2D power spectra in the k_perpendicular direction. For linear binning, this is the actual binsize, for log binning, this gives the bins per decade as 1/kperp_bin. The default is the min of the kx and ky pixel size (which are usually the same) for linear binning and 0.1 for log binning (giving 10 bins per decade).


### 1D binning keywords

**wedge_angles**: This is both a binning and a plotting keyword. List of angles on the sky to use for relevant binning cuts and plotted as lines on 2d power spectra. These are specified in degrees, the default list is for the MWA primary beam width and the horizon: [20, max_theta+90d] (where max_theta is the maximum zenith pointing angle since the MWA can point away from the zenith).

**coarse_harm_width**: Controls the number of kz bins around the MWA coarse band lines to exclude from 1D binning. The actual number of bins excluded = coarse_harm_width*2 - 1, so a value of 1 leads to 1 bin excluded, a value of 2 leads to 3 bins excluded, etc.

**log_k1d**: This is a flag (valid values are 0/1, default=0) indicating that logarithmic bins should be used for 1D power spectra.

**k1d_bin**: This sets the size of the bins for 1D power spectra. For linear binning, this is the actual binsize, for log binning, this gives the bins per decade as 1/kperp_bin. The default is the min of the (kx, ky, kz) pixel size for linear binning and 0.1 for log binning (giving 10 bins per decade).

**kpar_range_1dave**: This is a 2D vector giving the range of k_parallel values to include in the 1D binning.

**kperp_range_1dave**: This is a 2D vector giving the range of k_perpendicular values to include in the 1D binning. This keyword and kperp_range_lambda_1dave cannot both be set.

**kperp_range_lambda_1dave**: This is a 2D vector giving the range of k_perpendicular values *specified in wavelengths* to include in the 1D binning. This keyword and kperp_range_1dave cannot both be set.

**kx_range_1dave**: This is a 2D vector giving the range of k_x values to include in the 1D binning. This keyword and kx_range_lambda_1dave cannot both be set.

**kx_range_lambda_1dave**: This is a 2D vector giving the range of k_x values *specified in wavelengths* to include in the 1D binning. This keyword and kx_range_1dave cannot both be set.

**ky_range_1dave**: This is a 2D vector giving the range of k_y values to include in the 1D binning. This keyword and ky_range_lambda_1dave cannot both be set.

**ky_range_lambda_1dave**: This is a 2D vector giving the range of k_y values *specified in wavelengths* to include in the 1D binning. This keyword and ky_range_1dave cannot both be set.

**kperp_range_lambda_kparpower**: This is a 2D vector giving the range of k_perpendicular values *specified in wavelengths* to include in the 1D k_parallel power binning (resulting in a 1D power as a function k_parallel rather than k).

**kpar_range_kperppower**: This is a 2D vector giving the range of k_parallel values to include in the 1D k_perpendicular power binning (resulting in a 1D power as a function k_perpendicular rather than k).


## Outputs keywords and options:

**cube_power_info**: An output structure that contains some information from the run that we use for quantifying simulations. The fields are: [ave_power, wt_ave_power, uv_pix_area, uv_area, ave_weights, ave_weights_freq, wt_ave_power_freq, ave_power_freq, wt_ave_power_uvf, ave_power_uvf, nbsl_lambda2, nbsl_lambda2_freq] and optionally flat_power.

**ps_foldername**: This is the subdirectory under the folder_name where the eppsilon outputs are saved. Defaults to 'ps' + path_sep().

**save_path**: Path to save eppsilon outputs to. Default is folder_names + path_sep() + ps_foldername + path_sep().

**savefilebase**: Base name for eppsilon output files. Default is taken from the cube file names before the polarization tags.

**no_binning**: This is a flag (valid values are 0/1, default=0) indicating that binning should not be performed. This may be desirable if the large number of resulting files poses a problem (an issue for some clusters). In general binning is fast, so can be done later on another machine quickly. If this keyword is set no 1D or 2D binned plots will be displayed.

**save_slices**: This is a flag (valid values are 0/1, default=1) indicating that slices of the 3D cubes should not be saved at various stages. It may be desirable to turn this off if the number of resulting files poses a problem (an issue for some clusters).


## Plotting options:

**hinv**: This is a flag (valid values are 0/1, default=1) to use h^-1 Mpc rather than physical Mpc in plot units.

**plot_path**: Path to save eppsilon plot files to. Default is save_paths + path_sep() + 'plots' + path_sep()

**plot_filebase**: Base name for eppsilon plot files. Default is to savefilebase.

**note**: An output string giving some limited information about the eppsilon run that is useful for labelling plots.

**individual_plots**: This is a flag (valid values are 0/1, default=0) indicating that separate plot files should be made for each plot, rather than combining different polarization and cube type plots into a single file.

**png**: This is a flag (valid values are 0/1, default=0) to have eppsilon plots saved as png files rather than just being displayed on the screen. Only one of png, eps and pdf can be set.

**eps**: This is a flag (valid values are 0/1, default=0) to have eppsilon plots saved as eps files rather than just being displayed on the screen. Only one of png, eps and pdf can be set.

**pdf**: This is a flag (valid values are 0/1, default=0) to have eppsilon plots saved as pdf files rather than just being displayed on the screen. Only one of png, eps and pdf can be set.

### Flags & keywords controlling which plots to make:

**plot_stdset**: This is a flag (valid values are 0/1, default=1) to have eppsilon make the canonical set of plots. It can be usful to turn this off if you only want to look at a few unusual plots.

**plot_slices**: This is a flag (valid values are 0/1, default=0) to have eppsilon make plots of 2d slices through one of the 3d cubes.

**slice_type**: Only used if plot_slices is set. Controls which type of cube to plot slices of. Options are: ['raw', 'divided', 'sumdiff', 'weights', 'variance', 'power', 'var_power'], default is 'sumdiff'. 'power' and 'var_power' are slices in k-space, the others are in uvf space.

**uvf_plot_type**: Only used if plot_slices is set and slice_type is not 'kspace': Controls what to plot for slices of uvf cubes. Options are: ['abs', 'phase', 'real', 'imaginary', 'normalized'], default is 'abs'. 'normalized' means peak normalized (for a real plane only) -- useful for weights or variances.

**plot_1to2d**: This is a flag (valid values are 0/1, default=0) to have eppsilon make plots showing where the 1D bins come from up in 2D plots.

**plot_2d_masked**: This is a flag (valid values are 0/1, default=0) to have eppsilon make 2D plots using the 1D masking scenarios.

**plot_kpar_power**: This is a flag (valid values are 0/1, default=0) to have eppsilon make a 1D plot of the power as a function of k_parallel.

**plot_kperp_power**: This is a flag (valid values are 0/1, default=0) to have eppsilon make a 1D plot of the power as a function of k_perpendicular.

**plot_k0_power**: This is a flag (valid values are 0/1, default=0) to have eppsilon make a 1D plot showing the power in the k_parallel=0 modes of the 2D power spectra.

**plot_noise_1d**: This is a flag (valid values are 0/1, default=0) to have eppsilon make 1D plots showing the observed noise in the 1D power spectrum (power of the difference cubes).

**plot_sim_noise**: This is a flag (valid values are 0/1, default=0) to have eppsilon make 1D plots showing the simulated noise in the 1D power spectrum (from simulations based on the input variance cubes).

**plot_binning_hist**:  This is a flag (valid values are 0/1, default=0) to have eppsilon make histogram plots of the values going into each bin. These are *very* non-standard and are used for testing/debugging to understand how binning affects power spectrum levels and noise. These are the plots we used to understand why we could have 2D power spectra with noise-like pixels in the window but end up with 1D power spectra with no noise-like bins (it's just because you have more voxels going into the average, so the noise drops, but we were worried about biases and these plots showed there weren't any obvious ones.)

### 2D plotting options:

**plot_wedge_line**: This is a flag (valid values are 0/1, default=1) controlling whether or not the wedge lines (corresponding to the wedge_angles) are drawn on 2D plots.

#### axes options:
**kperp_linear_axis**: This is a flag (valid values are 0/1, default=0) controlling whether or not the k_perpendicular axis is linear (as opposed to logarithmic).

**kpar_linear_axis**: This is a flag (valid values are 0/1, default=0) controlling whether or not the k_parallel axis is linear (as opposed to logarithmic).

**kperp_plot_range**: A 2D vector giving the range of k_perpendicular values to show on 2D plots. Default is [5. wavelengths equivilent, maximum uv extent]

**kperp_lambda_plot_range**: A 2D vector *specified in wavelengths* giving the range of k_perpendicular values to show on 2D plots. Default is [5., maximum uv extent].

**kpar_plot_range**: A 2D vector giving the range of k_parallel values to show on 2D plots. Default is the full range.

**baseline_axis**: This is a flag (valid values are 0/1, default=1) controlling whether or not the top axis is labelled in wavelengths units.

**delay_axis**: This is a flag (valid values are 0/1, default=1) controlling whether or not the right axis is labelled in delay units.

**cable_length_axis**: This is a flag (valid values are 0/1, default=0) controlling whether or not the right axis is labelled in cable length units showing where reflections for different cable lengths would appear.


#### colorbar options:

**color_type**: Can be set to 'log', 'linear' or 'integer' to control the color mapping. The default in most cases is 'log', but it is 'linear' and 'integer' for some plots that map 1d bin characteristics to 2d.

#### colorbar range options:

**set_data_ranges**: This is a flag (valid values are 0/1, default=1) controlling whether the default colorbar ranges is set in ps_wrapper.pro.

**data_range**: A 2D vector giving the range to use for the colorbar for the 2D power spectrum plots. Defaults to the maximum range of the data or to the value in ps_wrapper if set_data_ranges is set.

**sigma_range**: A 2D vector giving the range to use for the colorbar for the 2D power spectrum error plots. Defaults to the maximum range of the data or to the value in ps_wrapper if set_data_ranges is set.

**nev_range**: A 2D vector giving the range to use for the colorbar for the 2D power spectrum expected noise plots. Defaults to the maximum range of the data or to the value in ps_wrapper if set_data_ranges is set.

**snr_range**: A 2D vector giving the range to use for the colorbar for the 2D power spectrum signal to noise plots. Defaults to the maximum range of the data or to the value in ps_wrapper if set_data_ranges is set.

**noise_range**: A 2D vector giving the range to use for the colorbar for the 2D power spectrum measured noise plots. Defaults to the maximum range of the data or to the value in ps_wrapper if set_data_ranges is set.

**nnr_range**: A 2D vector giving the range to use for the colorbar for the 2D power spectrum noise ratio plots. Defaults to the maximum range of the data or to the value in ps_wrapper if set_data_ranges is set.

**slice_range**: Only used if plot_slices is set. A 2D vector giving the range to use for the colorbar for the specified slice plots. Defaults to the maximum range of the data.

### 1D plot options:

**set_krange_1dave**: This is a flag (valid values are 0/1, default=1) controlling whether the default 1D plot range is set in ps_wrapper.pro.

**range_1d**: A 2D vector giving the range to use for the 1D power axis. Defaults to the maximum range of the data or to the value in ps_wrapper if set_krange_1dave is set.

**plot_1d_delta**: This is a flag (valid values are 0/1, default=0) controlling whether the 1D power axis is scaled by k^3 (theoretician units).

**plot_1d_error_bars**: This is a flag (valid values are 0/1, default=0) controlling whether the 1D plots are plotted with error bars. If this keyword is set, the thermal noise line is *not* drawn.

**plot_1d_nsigma**: This sets what sigma level the is shown for the thermal noise line (or error bars if plot_1d_error_bars is set) on 1D plots. Default is 1 (meaning 1 sigma), typical values to show for upper limits are 2 sigma so a value of 2 for this keyword.

**plot_eor_1d**: This is a flag (valid values are 0/1, default=0) controlling whether or not a theoretical eor 1D power line is drawn. Default is 1 if the sim keyword is set.

**plot_flat_1d**: This is a flag (valid values are 0/1, default=0) controlling whether or not a flat 1D power line is drawn. This is useful if flat power simulations are being plotted.

**no_text_1d**: This is a flag (valid values are 0/1, default=0) to prevent printing the legend text on the plot.
