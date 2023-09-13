FUNCTION mask_image,image,omask,order,SMOOTHBOX=smoothbox,NORMALIZE=normalize
;
; 14 Apr 06 added normalize keyword
; 10 Apr 06 modified by G. Sloan, mostly cosmetic, renamed to mask_image.pro
; 16 Feb 06 delivered by M. Devost to G. Sloan
; 15 Nov 06 modified by M. Devost
; c  Sep 06 created by D. Levitan as process_image.pro
;
; mask_image sets all data in an image outside the given order to zero
; if the smoothbox keyword is set, the data are smoothed
;
; INPUT
;   image     - an IRS lores input image
;   omask     - the corresponding order mask
;   order     - 1,2,3 for 1st, 2nd, bonus orders
;   smoothbox - keyword to smooth the image (= smoothing box size)
;   normalize - keyword to normalize each row in the image (to total=1)

; OUTPUT      - returns the processed image

nrow=n_elements(image[0,*])

; load the data into outimage and mask out-of-order pixels

outimage=image
idx=where(omask ne order,zerocount)
if (zerocount gt 0) then outimage[idx]=0.0

; if the norm keyword is set, then normalize each row in the image so total=1

if (keyword_set(normalize)) then for j=0,nrow-1 do $
  outimage[*,j]=outimage[*,j]/total(outimage[*,j],/nan)

; if smoothbox is set, then smooth the image

if (keyword_set(smoothbox)) then outimage=smooth(outimage,smoothbox,/nan)

RETURN,outimage
END
