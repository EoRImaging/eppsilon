# eppsilon installation instructions


## Requirements:

1. this repository (https://github.com/EoRImaging/eppsilon)

2. the astro IDL library (http://idlastro.gsfc.nasa.gov/)

3. the coyote library (http://www.idlcoyote.com/documents/programs.php)

4. the HEALPix library (http://sourceforge.net/projects/healpix/)

## Set up your IDL path
Add the above libraries to your IDL path in the order listed

To add a library to your path, on Windows machines type (on the IDL command line)
!PATH = !PATH + ';' + Expand_Path('+path\to\library\')

And on Unix-based machines type
!PATH = Expand_Path('+path/to/library/') + ':' + !PATH

Be sure to use the correct separator (; or : ) and include the ‘+’ sign at the start of +path/to/library/

Alternatively you can create an IDL startup file to set the IDL path.
For example, see here: http://slugidl.pbworks.com/w/page/28913708/Adding%20Programs%20to%20Your%20IDL%20Path

## Install HEALPix
following instructions at: http://healpix.jpl.nasa.gov/html/install.htm

## Install Imagemagick
if not already present, following instructions at: http://www.imagemagick.org

## Check Installations
Open a fresh terminal and start IDL to test the installation.

Suggested test commands:

 - print,cgHasImageMagick()
  - prints error if coyote library not installed, returns 0 if coyote library installed but Imagemagick not installed correctly, and returns 1 if both are installed correctly

 - astrolib
    - prints error if astro IDL library not installed correctly, prints message “ASTROLIB: Astronomy Library system variables have been added” if installed correctly]

 - init_healpix
  - prints error if HEALPix not installed correctly

If the above commands all work, it’s time to try it out on some data
