# eppsilon Keyword Dictionary
eppsilon uses keywords to create unique run-specific settings. This dictionary describes the purpose of each keyword, as well as their logic or applicable ranges.

This is a work in progress; please add keywords as you find them in alphabetical order with their corresponding definition. If there's a question about a definition or keyword, label it with !Q.

## input data definitions:

**folder_name**: **Required** This defines what folder the data live in. For FHD, it should be the top level folder for that run (which contains folders like Metadata and Healpix). On machines where there is support for the standard folder paths (i.e. there is a relevant *_folder_locs.pro file), this can just be the short folder name. It can also be a full path.

**loc_name**: This is a keyword to indicate which location the code is running on to help with defining the standard folder paths. It is detected from the hostname by default, but that apparently fails in some cases, so this keyword allows the user to specify it. Options are 'mit' or 'enterprise' (enterprise really means any of the ASU machines.)

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


## RTS specific keywords:

**rts**: This is a flag (valid values are 0/1, default=0) indicating that we are using RTS inputs.

**no_wtvar_rts**: This is a flag (valid values are 0/1, default=0) indicating that RTS weight cubes are not present.

**norm_rts_with_fhd**: This is a flag (valid values are 0/1, default=0) indicating that the vis_noise value should be set to the value usually calculated by FHD rather than the RTS value. This usually brings the normalization of the power spectra from RTS runs in line with the normalization of FHD runs


## CASA specific keywords:

**casa**: This is a flag (valid values are 0/1, default=0) indicating that we are using CASA inputs. This mode has not been used in a long time and it may be out of date.


## Simulation specific keywords:

**sim**: This is a flag (valid values are 0/1, default=0) indicating that the inputs are simulations, which results in minor plotting changes that are better suited to simulations.

**fix_sim_input**: This is a backwards compatability flag (valid values are 0/1, default=0) to fix some errors with very old simulations where the normalization was off by a factor of 2 because of a mistake in the simulation code. It should only be used with care if working with old simulations. It will soon be deprecated.


## Refresh keywords:
These flags (all default=0) tell the code to redo parts of the analysis that would normally be skipped because the outputs from that part of the code already exist. They cascade, so if you set an early one it will cause all affected downstream code to also re-run.

**refresh_dft**: Recalculate the DFT. Should only be used if the input cube files have changed. This is the longest step computationally.

**refresh_beam**: Redo the beam related calculations. This causes the power spectrum to be recalculated.

**refresh_ps**: Recalculate the power spectrum. This redoes everything downstream of the DFT and beam calculations.

**refresh_binning**: Redo the 1D & 2D binning.

**refresh_info**: Redo getting all the metadata together. This generally does not cause any major recalculations and is only really useful if folder structures have changed.


## DFT and UV grid related keywords:

**image_window_name**: Name of an image space window to apply to the input cubes before the DFT as an anti-aliasing filter. This is not supplied by default, but evidence is growing for a Tukey style window and the default may change to 'Tukey' soon. Any of the spec_window_type windows are also supported, but have never been tested.

**image_window_frac_size**: Only applies if image_window_name is 'Tukey'. The fractional size of the flat part of the Tukey function. Defaults to the ratio of the standard image size (given by a uv pixel size of 5 meters) to the input cube image size (so that the flat section matches the size of the standard image size).

**delta_uv_lambda**: The default size of the uv pixels to calculate in wavelengths. If there's no image_window_name set, this is defaulted to (5 meters)*freq/c, otherwise it is decreased to match the size of the input image cube. **Note**: the actual uv pixel size used by the code may be larger than this value if the size of the input image cube is smaller than would be needed to support delta_uv_lambda. In that case, the following warning will be printed: 'Image FoV is smaller than expected, increasing delta kperp to match image FoV'.

**max_uv_lambda**: The maximum distance from (0,0) in uv space that we should calculate in the DFT. **Note**: The actual uv extent used may be smaller than this value if the maximum uv extent of the input image cubes is smaller or if the longest gridded baseline is smaller. So the actual uv extent will be the minimum of max_uv_lambda, the maximum uv extent of the input image cubes and the longest gridded baseline.

**uv_avg**: A rarely used testing/exploration keyword that is only used if the uvf_input flag is set. This is a factor by which to average up pixels in the uv plane. A value of 2 results in pixels that are larger by a factor of 2 in each direction and 4x fewer total pixels.

**uv_img_clip**: A rarely used testing/exploration keyword that is only used if the uvf_input flag is set. This is a factor by which to clip the cubes in image space (i.e. 2D FFT to image space, clip, then 2D FFT back to uv). A value of 2 results in images that are smaller by a factor of 2 in each direction, 4x smaller overall, and uv pixels that are larger by a factor of 2 in each direction.

**dft_fchunk**: A very rarely used keyword, this tells the code how many frequencies to calculate at once. It used to be useful for memory management, but rewrites of the DFT code have made this largely obsolete.

**no_dft_progress**: This is a flag (valid values are 0/1, default=0) to suppress the printing of progress report statements during the DFT.


## Power spectrum calculation keywords:

**ave_removal**: This is a flag (valid values are 0/1, default=1) indicating that the average along the frequency axis of the cubes should be removed before the frequency Fourier Transform and then added back in to the k_parallel=0 mode. This has been well tested and substantially reduces the bleed of the flat spectrum foregrounds into the window.

**wt_cutoffs**: This is a keyword to support a very subtle power spectrum normalization issue identified by Adrian Liu and quantified for MWA with simulations. It turns out that when baselines substantially overlap there is a suppression of the power of stochastic fields (like the EoR). We measure this overlap using the weights cube and we found through simulations that for MWA antennas, above a certain weight value the power suppression asymptotes to a factor of 0.5. A conservative weight value to use to be sure we are in that regime is a weight value of 1, and if we throw away bins below that value we lose very little sensitivity with the MWA. The wt_cutoffs keyword is the weights value to use for the cutoff value (it is calculated in uvf space, so it applies to all kz values for a given kx,ky pixel). It is defaulted to 1.0 and is only used during the binnning to 1D or 2D power spectra (for 2D, it is used to normalize areas of the ps properly but lower values are not thrown out, for 1D, pixels with lower values are excluded and the power spectrum is normalized by a factor of 2). To turn this normalization correction off, set wt_cutoffs=0.

**wt_measures**: See the description for wt_cutoffs to get the full context. This keyword defines how the weight value for the power normalization is quantified in the uv plane given a uvf cube. Options are 'min' (minimum value along frequency axis) and 'ave' (average value along frequency axis), default is 'min'. We don't see much change between 'min' and 'ave', but 'min' is more conservative.

**spec_window_type**: Name of spectral window to apply in the frequency direction before Fourier Transforming. Default is 'Blackman-Harris'. Full set of options is: ['Hann', 'Hamming', 'Blackman', 'Nutall', 'Blackman-Nutall', 'Blackman-Harris', 'Blackman-Harris^2'].

**no_spec_window**: This is a flag (valid values are 0/1, default=1) indicating that no spectral window should be applied (see spec_window_type keyword).

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


## 1D binning keywords

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


## Output options:

**cube_power_info**: An output structure that contains some information from the run that we use for quantifying simulations. The fields are: [ave_power, wt_ave_power, uv_pix_area, uv_area, ave_weights, ave_weights_freq, wt_ave_power_freq, ave_power_freq, wt_ave_power_uvf, ave_power_uvf, nbsl_lambda2, nbsl_lambda2_freq] and optionally flat_power.

**ps_foldername**: This is the subdirectory under the folder_name where the eppsilon outputs are saved. Defaults to 'ps' + path_sep().

**save_path**: Path to save eppsilon outputs to. Default is folder_names + path_sep() + ps_foldername + path_sep().

**savefilebase**: Base name for eppsilon output files. Default is taken from the cube file names before the polarization tags.

**plot_path**: Path to save eppsilon plot files to. Default is save_paths + path_sep() + 'plots' + path_sep()

**plot_filebase**: Base name for eppsilon plot files. Default is to savefilebase.

**individual_plots**: This is a flag (valid values are 0/1, default=0) indicating that separate plot files should be made for each plot, rather than combining different polarization and cube type plots into a single file.

**png**: This is a flag (valid values are 0/1, default=0) indicating that eppsilon plots should be saved as png files rather than just being displayed on the screen. Only one of png, eps and pdf can be set.

**eps**: This is a flag (valid values are 0/1, default=0) indicating that eppsilon plots should be saved as eps files rather than just being displayed on the screen. Only one of png, eps and pdf can be set.

**pdf**: This is a flag (valid values are 0/1, default=0) indicating that eppsilon plots should be saved as pdf files rather than just being displayed on the screen. Only one of png, eps and pdf can be set.


## Plotting options:

**plot_wedge_line**

**hinv**

**note**

### Which plots to make:

**plot_slices**

**slice_type**

**slice_range**

**uvf_plot_type**

**plot_stdset**


**plot_1to2d**

**plot_2d_masked**


**plot_kpar_power**

**plot_kperp_power**

**plot_k0_power**

**plot_noise_1d**

**plot_sim_noise**

**plot_binning_hist**

**plot_2d_masked**: use the 1D masking scenarios on all 2D plots. <br />
  -*Turn on/off*: 1/0 <br />


### axes options:
**kperp_linear_axis**

**kpar_linear_axis**

**kperp_plot_range**

**kperp_lambda_plot_range**

**kpar_plot_range**

**baseline_axis**

**delay_axis**

**cable_length_axis**


### colorbar range options:

**set_data_ranges**: This is a flag (valid values are 0/1, default=1) indicating that the default colorbar ranges should be set in ps_wrapper.pro.

**data_range**: A 2D vector giving the range to use for the colorbar for the 2D power spectrum plots. Defaults to the maximum range of the data, overridden by set_data_ranges.

**sigma_range**: A 2D vector giving the range to use for the colorbar for the 2D power spectrum error plots. Defaults to the maximum range of the data, overridden by set_data_ranges.

**nev_range**: A 2D vector giving the range to use for the colorbar for the 2D power spectrum expected noise plots. Defaults to the maximum range of the data, overridden by set_data_ranges.

**snr_range**: A 2D vector giving the range to use for the colorbar for the 2D power spectrum signal to noise plots. Defaults to the maximum range of the data, overridden by set_data_ranges.

**noise_range**: A 2D vector giving the range to use for the colorbar for the 2D power spectrum measured noise plots. Defaults to the maximum range of the data, overridden by set_data_ranges.

**nnr_range**: A 2D vector giving the range to use for the colorbar for the 2D power spectrum noise ratio plots. Defaults to the maximum range of the data, overridden by set_data_ranges.

### 1D plot options:

**set_krange_1dave**

**range_1d**

**plot_1d_delta**

**plot_1d_error_bars**

**plot_1d_nsigma**

**plot_eor_1d**

**plot_flat_1d**

**no_text_1d**: This is a flag (valid values are 0/1, default=0) to prevent printing the legend text on the plot.
