FUNCTION build_psfrow,psfrow,subpix,centcol,omaprow,order,position,pillow=pillow,fast=fast

; 24 Aug 09 returned tsum to original subpix range
; 22 Aug 09 moved pillow application earlier in the code
;           replaced rebin, shifted tsum range, put 0.5 subpix shift back
; 20 Aug 09 added one column to tsum totals, removed out_x 0.5 subpix shift
; 17 Aug 09 shifted out_x by 0.5 subpixels
; 16 Jul 09 now including the pillow function (intrapixel responsivity)
; 29 Jun 09 slightly modified to follow Vianney's changes to his version
; 22 Jun 09 created
;
; contains the kernel of build_psf
; constructs a spatial profile for a row using the PSF for that row
;
; INPUT
;   psfrow     - a PSF at subpixel resolution for the row in question
;   subpix     - the number of subpixels per pixel
;   centcol    - scalar - the central column of the PSF image
;   omaprow    - the relevant row of the order map
;   order      - order number
;   position   - the fractional center for the reconstructed PSF in the row 
;
; OUTPUT - returns a vector - PSF shifted and mapped at full pixel resolution
;
; ALGORITHM:
; follows Vianney and Don's PSF-shifting algorithm
;   build regridded PSF image one row at a time 
;   shift PSF (in subpixel units) to account for fractional shift of position
;   then rebin back to full pixel resolution and renormalize

; normalize pillow if defined, else define to be 1 for all subpix

psfcol=n_elements(psfrow)
ncol=n_elements(omaprow) ; or else hardcode it to 128

; find range of columns in order

idx=where(omaprow eq order, complement=nidx)
if (max(idx) gt -1) then begin

; stretch and shift the data pixel scale (in subpixel units) - out_x

  out_x = subpix * (findgen(ncol)-position) + centcol + 0.5
  temp  = out_x
  for i=1,subpix-1 do out_x=[out_x,temp+i]
  out_x = out_x[sort(out_x)]
  len=n_elements(out_x)

; shift the PSF to the new data pixel scale (out_x)
; and apply the pillow f'n to every set of subpix pixels

  psf = interpol( psfrow, findgen(psfcol), out_x ) 
  for i=0,ncol-1 do begin
    i0 = i*subpix & i1 = i0 + subpix-1
    psf[i0:i1] = psf[i0:i1] * pillow
  endfor

; rebin the shifted (and pillowed) PSF back to full pixels
; fast mode uses total instead of rebin (as of 22 Aug)
; default (slow) mode uses tsum 

  outpsf = fltarr(ncol)
  if (keyword_set(fast) ne 0) then begin
    for i=0,ncol-1 do begin
      i0 = i*subpix & i1 = i0 + subpix-1
      outpsf[i] = total(psf[i0:i1])
    endfor
  endif else begin
    for i=0,ncol-1 do begin
      i0 = i*subpix & i1 = i0 + subpix-1
      outpsf[i] = tsum(out_x[i0:i1],psf[i0:i1])
    endfor
  endelse

  outpsf[nidx] = 0.0                     ; zero out-of-order data
  outpsf = outpsf / total(psfrow,/nan)   ; normalization (usually close to 1)

endif else outpsf=fltarr(ncol)           ; return zeroes if no valid data in row

RETURN,outpsf
END
