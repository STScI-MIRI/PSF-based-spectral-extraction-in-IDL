FUNCTION exspecim,image,error,bmask,psf,lmap,order,degree=degree,residual=residual,sky=sky,_extra=e

; 17 Feb 12 passing degree to exspecrow explicitly, passing sky image to imspex
; 24 Aug 09 receiving sky polynomial from exspecrow as skycc
; 22 Jul 09 modifications to run properly for multiple sources
; 10 Jul 09 added residual keyword to return the residual image
; 29 Jun 09 modified to receive and use error and bmask planes
; 22 Jun 09 modified substantially to handle multiple PSFs and background
; 26 Apr 07 implemented _extra keyword to pass keywords to exspecrow
; 21 Apr 07 created as imspec1
;
; exspecim extracts one or more spectra from an image, with/without background
; it calls exspecrow, one row at a time
; only the columns between the first and last non-zero elements are passed
; every row is called, whether valid or not
; when done, it filters the resulting spectral data array for NaN data
;
; INPUT
;   image     - 2D spectral image, out-of-order information is masked
;   error     - 2D image with uncertainties, out-of-order information is masked
;   bmask     - 2D image with bmask data, out-of-order information is masked
;   psf       - 2-3D image with PSF in each row shifted to position of source
;             - this image has one plane per PSF to be extracted
;   lmap      - 2D image containing wavelength of each pixel
;   order     - the order being extracted
;   degree    - the degree of the polynomial fitted to the background (-1=none)
;
; OUTPUT  - a spectral data cube, one plane per source
;           col 0 - wavelength (um), col 1 - flux (DN), col 2 - error (DN)
;           col 3 - order, col 4 - mask information
;           for now moment, mask information is set to zero by exspecrow
;   residual  - the residual image generated to pass upwards
;
; algorithm
; create spectral data array to match number of rows in image
; iterate through rows to load columns
; parse result to filter rows with NaN data in wavelength or flux columns

; check image dimensions and create spectral data array

ncol=n_elements(image[*,0])
nrow=n_elements(image[0,*])

; check psf image/cube to determine dimensionality and number of sources

sz=size(psf)
case sz[0] of
  2    : npsf=1      ; psf is an image, one PSF for same source per row
  3    : npsf=sz[3]  ; psf is a cube, one PSF image per source, multiple PSFs
  else : npsf=0
endcase

; define output columns - definitions fixed in this software version

outcol=5 & lcol=0 & fcol=1 & ecol=2 & ocol=3 & mcol=4

; set up spectral, residual, and sky arrays

specdata = fltarr(outcol,nrow,npsf)
residual = image
sky      = image

; load spectral data array with repeated calls to exspecrow

for j=0,nrow-1 do begin

; find first and last nonzero column in the row

  idx=where(image[*,j] ne 0)
  if (max(idx) gt -1) then begin
    i0=min(idx) & i1=max(idx)
  endif else begin
    i0=0 & i1=20              ; dummy values don't matter since all data=0
  endelse

; load image, error, bmask, and wavelength data

  imrow  = reform(double(image[i0:i1,j]))
  errow  = reform(double(error[i0:i1,j]))
  bmrow  = reform(double(bmask[i0:i1,j]))
  lamrow = reform(double(lmap[i0:i1,j]))

; load psf data - have to account for varying dimensionalities

  case npsf of
    0    : psfrow = imrow-imrow                    ; load with zeroes
    1    : psfrow = reform(double(psf[i0:i1,j]))   ; pass 1-D vector
    else : psfrow = reform(double(psf[i0:i1,j,*])) ; pass all PSFs in given row
  endcase

; call exspecrow and load the results (which will now always be 3D)

  specrow = exspecrow(imrow,errow,bmrow,psfrow,lamrow,order,$
    row=j,sky=skycc,degree=degree,_extra=e)
  specdata[*,j,*] = specrow

; compute residual row 

  model  = imrow-imrow
  for n=0,npsf-1 do model = model + specrow[fcol,n]*psfrow[*,n]

  skyrow = imrow-imrow
  if (degree gt -1) then begin
    x = i0 + findgen(i1-i0+1)
    for d=0,degree-1 do skyrow = skyrow + skycc[d]*x^d
  endif
  sky[i0:i1,j]      = skyrow
  residual[i0:i1,j] = imrow - model - skyrow
endfor

; check for NaNed data - if no data valid, then return scalar "0"

count=0
for n=0,npsf-1 do begin
  idx=where(finite(specdata[lcol,*,n]) ne 0 and finite(specdata[fcol,*,n]) ne 0)
  if (max(idx) gt -1) then begin
    plane=specdata[*,idx,n]
    if (count eq 0) then outdata=plane else outdata=[[[outdata]],[[plane]]]
    count=count+1
  endif
endfor
if (count eq 0) then outdata=0

RETURN,outdata
END
