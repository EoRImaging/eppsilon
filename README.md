# eppsilon
error propagated power spectrum with interleaved observed noise

# Introduction
eppsilon is a package to calculate 21 cm power spectra and associated
error bars and metrics from 3D image cubes. eppsilon was developed in close
collaboration with FHD ([Fast Holographic Deconvolution](https://github.com/EoRImaging/FHD))
and it takes in FHD outputs seamlessly, but it can also create power spectra
from image cubes produced by other packages, including the MWA RTS
(Murchison Widefield Array Real Time System).

eppsilon takes a very straight-forward error propagation approach that is quite
simple but powerful. To work best, eppsilon needs to be provided with
data, weights and variance cubes which have been constructed in the following way:
  - data: visibilities gridded with the primary beam in uv-space
  - weights: ones gridded with the primary beam in uv-space
  - variances: ones gridded with the square of the primary beam in uv-space
The cubes can be provided either in _uvf_ space or in image space
(healpix or regularly gridded). eppsilon can provide the best metrics if the
cubes are split into two separate cubes (usually by gridding the even time
samples into one cube and the odd time samples into another to produced
pairs of interleaved cubes). FHD can generate all these eppsilon inputs natively.

## Language Details
eppsilon is written in IDL (Interactive Data Language), a common language used
in astronomy and geoscience. IDL is not an open source language, but most
university astronomy departments have IDL licenses. There is an open source
alternative, [GDL](https://github.com/gnudatalanguage/gdl), but we have not
tested how well eppsilon works in GDL. We are supportive of open-source code,
however, so if you would like to run eppsilon using GDL and have any trouble
please file an issue and we will do our best to help.

# Installation:
Installing IDL libraries is fairly simple, you just need to download the files
and then ensure that they are on your IDL path. See David Fanning's excellent
[Coyote's guide to IDL Programming](http://www.idlcoyote.com/code_tips/installcoyote.php)
for a good guide to installing IDL libraries. In fact, we cannot recommend
that site highly enough for all your IDL questions.

## Dependencies
eppsilon requires two standard IDL libraries and the `fhdps_utils` library,
a small utility library shared by FHD and eppsilon.
If you want to start from HEALPix image cubes (this is standard if you're using
FHD with eppsilon), the HEALPix IDL library is also required.
ImageMagick is an optional (non-IDL) dependency, it's only needed if you want to
produce PDF or PNG plots (postscript or encapsulated postscript don't require it).
 - [IDL astronomy library](https://idlastro.gsfc.nasa.gov/)
 - [Coyote Graphics Library](http://www.idlcoyote.com/documents/programs.php)
 - [fhdps_utils](https://github.com/EoRImaging/fhdps_utils)
 - Optional: [HEALPix IDL Library](https://healpix.sourceforge.io/)
 - Optional: [ImageMagick](https://www.imagemagick.org/).

## Check Installations
Open a fresh terminal and start IDL to test the installation.
Suggested test commands:

- `print, cgHasImageMagick()`
 - prints an error if the Coyote library is not installed, prints 0 if the
 Coyote Library is installed but Imagemagick not installed correctly, and prints
 1 if both are installed correctly

- `astrolib`
  - prints an error if the IDL Astronomy library is not installed correctly,
  prints message “ASTROLIB: Astronomy Library system variables have been added”
  if it is installed correctly

- `init_healpix`
 - prints an error if the HEALPix IDL library is not installed correctly

 ## Community Guidelines
Contributions to this package to add new features or address any of the
issues in the [issue log](https://github.com/EoRImaging/eppsilon/issues)
are very welcome, as are bug reports and feature requests.

Please submit improvements as pull requests against the repo. We commit to
reviewing pull requests promptly. Bug reports or feature requests are also very
welcome, please add them to the issue log after verifying that the issue does
not already exist. Comments on existing issues are also welcome.

# API
The primary user interface is through `ps_wrapper`. There are a fair number
of keyword options, for which we strive to have sensible defaulting but which the
power user may want to adjust. The keywords are described in full in the
[dictionary](dictionary.md). Two other wrappers, `ps_diff_wrapper` and
`ps_ratio_wrapper` are useful for jackknife comparisons of different data sets
or different runs using the same data. Finally, `ps_cube_images` and
`ps_cube_movie` can be useful for examining HEALPix image cubes in detail.

# Maintainers
eppsilon is maintained by Bryna Hazelton (University of Washington). Please use
the issue log for code-related discussions. You can contact me privately if needed at brynah (at) phys.washington.edu.
