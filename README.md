# PSF-based-spectral-extraction-in-IDL
PSF-based spectral extraction code written in IDL

Files in this release:

from ~/procs/imspex: ----------------------------------------------------------

imfspex.pro      - wrapper for imspex, handles file input/output and keywords
imspex.pro       - extracts spectra from one image
psfheader.pro    - parses the PSF file header, returns the spectral trace
nominalpos.pro   - finds nominal source position given a spec. trace and FOVID
pos_faint.pro    - algorithm to find a faint source in an image 
pos_bright.pro   - algorithm to find a bright source in an image
sum_profile.pro  - used by pos_faint to determine error
build_psf.pro    - builds a PSF image given a 2-D PSF, pos'n, and spec. trace
build_psfrow.pro - builds one row of the PSF image
exspecim.pro     - extracts spectra from an image
exspecrow.pro    - extracts one wavelength element of the spectra from an image
opt_mregress.pro - multiple linear regression (orig. written by Ph. Prugniel)
mk_pos.pro       - determines RA and DEC from FITS header and pos'n in slit

from ~/procs/irs: -------------------------------------------------------------

mk_map.pro       - converts wavesamp file to order, wavelength, pos'n maps
mk_script.pro    - generates batch scripts for imspex and other pipeline steps
                   (NOTE - code for other pipeline steps NOT in this release)

from ~/procs/fitstran: --------------------------------------------------------

rd_ipac_tbl.pro  - old, rude, crude code to read wavesamp file

from ~/procs/psf: -------------------------------------------------------------

mask_image.pro   - zeros out-of-order data

from ~/procs/sp: --------------------------------------------------------------

wr_spfits.pro    - writes a spectral FITS file

calibration files: ------------------------------------------------------------

b0_wavsamp_19.tbl - SSC provided wavesamp file
sl.lam.fits       - SL wavelength grid (for resampling)
psf_0_18v1.fits   - SL super-sampled PSF file

test image: -------------------------------------------------------------------

hd173511.61a2.0301.ic.fits

Target:  HD 173511, Campaign 61, 1st pointing (AOR 24373760)
exposure:  SL2, nod 2 (EXPID 0003), 2nd of 3 ramps (DCE 0001)

this image has been flatfielded, nod-differenced to remove the background
  and cleaned to remove rogue pixels

IDL command to extract SL2 and SL-bonus spectra from the image: ---------------

imfspex,'hd173511.61a2.0301.ic.fits','hd173511.61a2.0301.sp.fits',$
[2,3],'psf_0_18v1.fits','sl.lam.fits',wavsamp='b0_wavsamp_18.tbl',fovrange=2

The spectrum is saved in a spectral FITS format
col 0 = wavelength (um)
col 1 = flux density (Jy)
col 2 = uncertainty in col 1 (Jy)
col 3 = spectral order/segment number (2 = SL2, 3 = SL-bonus)
