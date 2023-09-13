FUNCTION build_psf,psf,subpix,centcol,omap,order,position,tracecoeff,_extra=e

; 16 Jul 09 added _extra to pass pillow function or fast option to build_psfrow
; 22 Jun 09 modified to allow PSF to shift smoothly vs. integral subpix shifts
;           now implementing Vianney's shift algorithm
; 28 Apr 07 made some minor improvements
; 20 Apr 07 created
;
; build_psf contructs a psf image matching the pixel scale and position of
;   the image defined by the input values position and tracecoeff
; usually called by exspecrow, calls build_psfrow
;
; INPUT
;   psf        - a PSF image at subpixel resolution
;   subpix     - the number of subpixels per pixel
;   centcol    - scalar - the central column of the PSF image
;   omap       - order map
;   order      - order to build PSFs for
;   position   - the position of the real trace in row 0
;   tracecoeff - the polynomical coefficients defining the spectral trace
;
; OUTPUT - returns a 2D image with the PSF mapped to full pixel resolution
;          and shifting up the array as defined by position and tracecoeff

; algorithm
; 1.  build the trace by calling poly
;     first add position to tracecoeff[0] to shift y_intercept to actual trace
; 2.  build regridded PSF image one row at a time
;     in each row, find the subpixel position of the trace
;     then bin input psf subpixels with correct shift
;     resulting PSF will be at actual pixel resolution

; STEP 1 - build trace

psfcol=n_elements(psf[*,0])
nrow=n_elements(psf[0,*])
row=findgen(nrow)

newtrace=tracecoeff
newtrace[0]=position         ; assuming that y intercept should be zero
spectrace=poly(row,newtrace)

; STEP 2 - build the PSF image one row at a time by calling build_psfrow.pro
; output array will have same dimensions as order map

outarray=float(omap-omap) 

for j=0,nrow-1 do begin
  outarray[*,j] = build_psfrow( reform(psf[*,j]),subpix,centcol, $
                  reform(omap[*,j]),order,spectrace[j],_extra=e)
endfor

RETURN,outarray
END
