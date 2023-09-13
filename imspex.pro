FUNCTION imspex,image,order,omap,lmap,tracecoeff,inpsf,subpix,centcol,FINDMETHOD=findmethod,INPOS=inpos,OUTPOS=outpos,NFIND=nfind,COLRANGE=colrange,DEGREE=degree,RESIDUAL=residual,SKY=sky,DIAG=diag,_extra=e

; 19 Feb 12 passing degree explicitly to pos_faint, now always set in imfspex
;           separating input and output position as inpos and outpos
; 16 Jul 09 cleaning up
; 10 Jul 09 default position finder is now pos_faint.pro
;  7 Jul 09 updated to call pos_bright.pro to find position
; 29 Jun 09 updated to receive a 3D spectral image
; 22 Jun 09 updated to handle new kernel, which allows more than one PSF
; 26 Apr 07 added method, diag, and _extra keywords
; 22 Apr 07 ensured that psfs are 2D in pass to build_psf with reform function
; 21 Apr 07 added position keyword
; 20 Apr 07 created
;
; imspex extracts a spectrum (or spectra) from one order of an IRS lo-res image
; it sets up input for a call to exspecim, calls it
; and arranges the output to pass back to the calling function (imfspex)
;
; this test version makes several assumptions
; - the image is cleaned
; - the noise in each pixel is the same (error planes are not yet passed)
; - the source is easily traced, allowing a simple find algorithm
; - or better yet, the position(s) is/are given by the position keyword

; INPUT
;   image      - a 3D spectral image (plane 0 = signal, 1 = error, 2 = bmask)
;   order      - the order number, 1, 2, or 3 for bonus
;   omap       - a 2D image containing the order for each pixel
;   lmap       - a 2D image containing the wavelength in um for each pixel
;   tracecoeff - the polynomial coefficients defining the spectral trace
;   inpsf      - an image containg the PSF for each row
;   subpix     - the number of subpix per pixel in the inpsf image
;   centcol    - the central column of the inpsf image
;   findmethod - keyword, sets method of finding source in slit, DEFAULT=0
;                0 = pos_faint (default), 1 = pos_bright, 
;                2 = old faint-source algorithm, 3 = old bright-source algorithm
;   inpos      - fixed position to extract
;   degree     - degree of polynomial fitted to background - for exspecrow
;   diag       - keyword for diagnostic mode
;                passed to exspecim, pos_faint via _extra
;   _extra     - keywords passed to exspecim
;
; OUTPUT - returns a spectral data array with the extracted spectrum
;              - col 0 = wavelength (um), col 1 = flux (DN)
;              - col 2 = error (zeroed),  col 3 = order number
;              - col 4 = mask data (zeroed)
;  outpos      - array of all positions extracted
;  residual    - residual image passed up from exspecim
;  sky         - sky image passed up from exspecim

; calls the following subroutines:
;   mask_image    
;   im_centroid   old BRIGHT-SOURCE algorithm only 
;   find_trace    old BRIGHT-SOURCE algorithm only 
;   find_faint    old FAINT-SOURCE algorithm
;   pos_faint     DEFAULT algorithm
;   build_psf
;   exspecim      (calls exspecrow)
;
; algorithm
; 1.  setup - find array sizes, prepare the image by calling mask_image
; 2.  find the position of the spectrum in the image 
;     findmethod=0 - new FAINT_SOURCE algorithm - pos_faint - DEFAULT
;     findmethod=1 - new BRIGHT-SOURCE algorithm - pos_bright
;     findmethod=2 - old FAINT-SOURCE algorithm - find_faint
;       use trace coefficients to align each row to center in row 0
;       coadd most rows, fit a gaussian, add 0.5 pix to register zero point
;     findmethod=3 - BRIGHT-SOURCE algorithm 
;       call im_centroid to find the centroids vector (pos'n in each row)
;       pass this and trace coefficients to find_trace to locate the position
; 3.  call build_psf to build a pixel-scale PSF image for the actual data image
; 4.  call exspecim to extract the spectrum (calls exspecrow)
; 5.  arrange output from exspecim

; STEP 1 - setup
; check keywords and set defaults (inpos checked in step 2)

if (keyword_set(diag) eq 0) then diag=0 
if (keyword_set(findmethod) eq 0) then findmethod=0 else findmethod=1
if (keyword_set(nfind) eq 0) then nfind=0
if (n_elements(colrange) eq 0) then colrange=0

; find array sizes
; prepare the image - mask_image zeroes out-of-order data

ncol=n_elements(image[*,0,0])
nrow=n_elements(image[0,*,0])

data  = mask_image(reform(image[*,*,0]),omap,order)
error = mask_image(reform(image[*,*,1]),omap,order)
bmask = mask_image(reform(image[*,*,2]),omap,order)

; STEP 2 - find the source if position not given as a parameter
; default is FAINT-SOURCE algorithm
; note that the other algorithms ignore all fixed positions

findflag=0
if (keyword_set(inpos) eq 0) then findflag=1 else foundpos=inpos
if (nfind gt 0 or n_elements(colrange) gt 1) then findflag=1

if (findflag eq 1) then begin
  case findmethod of 
    1 : begin                     ; new bright-source algorithm
      foundpos=pos_bright(data,inpsf,subpix,centcol,tracecoeff,$
               omap,order,posvector=spectrace,sigpos=sigpos)
    end
    2 : begin                     ; old faint-source algorithm
      foundpos=find_faint(data,tracecoeff) 
    end
    3 : begin                     ; old bright-source algorithm
      centroids=im_centroid(data) ; find centroids vector for our image
      foundpos=find_trace(centroids,tracecoeff) ; find scalar position
    end
    else : begin                  ; new faint-source algorithm - DEFAULT
      foundpos=pos_faint(data,error,bmask,inpsf,subpix,centcol,$
               tracecoeff,omap,order,lmap,position=position,$
               nfind=nfind,colrange=colrange,diag=diag,degree=degree,_extra=e)
    end
  endcase
endif 

; STEP 3 - build the psf image using build_psf
;          if pos contains more than one element, then psf will be a cube

npsf=n_elements(foundpos)
psf=fltarr(ncol,nrow,npsf)
for n=0,npsf-1 do begin
  psfplane=build_psf(reform(inpsf),subpix,centcol,omap,order,$
           foundpos[n],tracecoeff,_extra=e)
  psf[*,*,n]=psfplane
  if (diag eq 9) then $
    writefits,'psf.'+string(n,format='(I1)')+'.im.fits',psfplane
endfor
outpos=foundpos

; STEP 4 - extract the spectrum using exspecim, which calls exspecrow

spectrum = exspecim(data,error,bmask,psf,lmap,order,$
           diag=diag,degree=degree,residual=residual,sky=sky,_extra=e)

; STEP 5 - return the spectrum array (which may contain more than one spectrum)
; THIS PORTION OF THE CODE WILL NEED TO BETTER HANDLE WHAT exspecim RETURNS

RETURN,spectrum
END
