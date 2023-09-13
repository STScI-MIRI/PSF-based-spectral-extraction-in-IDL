;
; NAME:
;              OPT_MREGRESS
;
; PURPOSE:
;              Perform a multiple linear regression fit
;
; USAGE:
;              Result = MREGRESS(X, Y,                           $ 
;                                 MEASURE_ERRORS=measure_errors, $
;                                 SIGMA=sigma,                   $
;                                 INVERSE=inv,                   $
;                                 STATUS=status)
;
; ARGUMENTS:
;       X:     (Input) The array of independent variable data.  X must
;              be dimensioned as an array of Npoints by Nterms, where
;              there are Nterms coefficients (independent variables) to be
;              found and Npoints of samples.
;              Note that this array is transposed with respect to the
;              input of the IDL routine REGRESS.
;
;       Y:     (Input) The vector of dependent variable points.  Y must have 
;              Npoints elements.
;
;
; KEYWORDS:
;   MEASURE_ERRORS : Set this keyword to a vector containing standard
;       measurement errors for each point Y[i].  This vector must be the same
;       length as X and Y.
;
;       For Gaussian errors (e.g. instrumental uncertainties), MEASURE_ERRORS 
;       should be set to the standard deviations of each point in Y. 
;       For Poisson or statistical weighting MEASURE_ERRORS should be
;       set to sqrt(Y).
;
;   SIGMA : Set this keyword to a named variable to recieve the errors
;       on the returned coefficients
;
;   INVERSE : Set this keyword to a named value to recieve the structure
;       containing the covariance matrix and other terms that can be re-used
;       to solve the system with the same X and other Y.
;
;   STATUS : Set this keyword to a named variable to receive the status
;       of the operation. Possible status values are:
;         - 0 for successful completion, 
;         - 1 for a singular array (which indicates that the inversion
;           is invalid), and 
;         - 2 which is a warning that a small pivot element was used
;           and that significant accuracy was probably lost.
;       If STATUS is not specified then any error messages will be output
;       to the screen.
;
; RETURN:
;       MREGRESS returns a column vector of coefficients that has Nterms
;       elements.
;
; DESCRIPTION:
;       Adapted from the IDL program REGRESS
;
;       Perform a multiple linear regression fit:
;               Y[i] = A[0]*X[0,i] + A[1]*X[1,i] + ... +
;                      A[Nterms-1]*X[Nterms-1,i]
;
;       This is a variant of the IDL funtion REGRESS, where the
;       covariance matrix and other intermediate arrays can be saved
;       to solve faster the system with successive values of Y and
;       the same X.
;       An noticeable difference is also that MREGRESS does not
;       include a constant term (if necessary, the constant may be
;       included as one of the terms of X).
;
; HISTORY:
;       17 Jan 08, Written, by Ph. Prugniel from the IDL REGRESS)
;       spring 2009, modified slightly by Vianney Lebouteiller - weights
;       22 Jun 09, renamed from sm_mregress to opt_mregress by G. Sloan
;-
; CATEGORY:    ULY_UTIL
;------------------------------------------------------------------------------
FUNCTION opt_mregress,x,y, $
                  weights=weights, $
                  SIGMA=sigma, $
                  INVERSE=inv, $
                  STATUS=status

COMPILE_OPT idl2

ON_ERROR,2              ;Return to caller if an error occurs
sy = SIZE(Y)            ;Get dimensions of x and y.
sx = SIZE(X)

ndimX = sx[0]
nterm = (ndimX EQ 1) ? 1 : sx[ndimX]   ; # of terms (coefficients)
nptsX = sx[1]       ;# of observations (samples)
npts = sy[1]            ;# of observations (samples)

IF (nptsX NE npts) THEN MESSAGE, $
        'X and Y have incompatible dimensions.'
IF (ndimX EQ 1) THEN x = REFORM(x, npts, 1, /OVER)   ; change X to a 2D array

nfree = npts-1          ;degs of freedom

invert = 0b
if n_tags(inv) gt 0 then if (size(inv.wx))[1] ne (size(y))[1] then invert = 1b

if n_tags(inv) le 0 or invert then begin
    sw = TOTAL(weights)/npts
    weights = weights/sw
    wgt = REBIN(weights,npts,nterm)

; This commented block is the version with a constant term
;    xmean = TOTAL(wgt*x,1, DOUBLE=double)
;    if (nterm gt 1) then $
;      xx = x - REBIN(transpose(xmean),npts,nterm) $ ;x(j,i) - xmean(i)
;    else xx = x - xmean
;    wx = TEMPORARY(wgt) * xx    ;weights(i)*(x(j,i)-xmean(i))
;    sigmax = SQRT(TOTAL(xx*wx,1)/nfree) ;weights(i)*(x(j,i)-xm)*(x(k,i)-xm)
;    ar = matrix_multiply(wx, xx, /ATRANSPOSE)/(nfree * sigmax #sigmax)

    wx = TEMPORARY(wgt) * x    ;weights(i)*x(j,i)
    
    sigmax = SQRT(TOTAL(x*wx,1)/nfree) ;weights(i)*x(j,i)*(x(k,i)-xm)

    ar = matrix_multiply(wx, x, /ATRANSPOSE)/(nfree * sigmax #sigmax)

    ar = INVERT(ar, status)
    IF (status EQ 1L && ~ARG_PRESENT(status)) THEN BEGIN
        IF (ndimX EQ 1) THEN x = REFORM(x, /OVER) ; change X back to a vector
        MESSAGE, /INFO, "Inversion failed due to singular array."
    END
    IF (status EQ 2L && ~ARG_PRESENT(status)) THEN BEGIN
        MESSAGE, /INFO, $
          "Inversion used a small pivot element. Results may be inaccurate."
    endif

    ; error terms
    sigma = SQRT(ar[LINDGEN(nterm)*(nterm+1)]/(sw*nfree*sigmax^2))

    inv = {a:ar, ww:weights, wx:wx, sx:sigmax, sigma:sigma}
endif

; this commented block is the version with constant term
;ymean = TOTAL(y*inv.ww)   ;y mean
;yy = y-ymean
;sigmay = SQRT(TOTAL(inv.ww * yy^2)*npts/nfree) ;weights*(y(i)-ymean)
;correlation = matrix_multiply(inv.wx, yy, /ATRANSPOSE) / (inv.sx * sigmay * nfree)

sigmay = SQRT(TOTAL(inv.ww * y^2)/nfree)    ;weights*y(i)

correlation = matrix_multiply(inv.wx, y, /ATRANSPOSE) / (inv.sx*sigmay*nfree)
a = (correlation # inv.a)*(sigmay/inv.sx)   ;get coefficients

;yfit = matrix_multiply(a,x, /BTRANSPOSE)      ;compute fit
;const = ymean - TOTAL(a*xmean)                ;constant term
;yfit = yfit + const                           ;add it in

IF (ndimX EQ 1) THEN x = REFORM(x, /OVER)     ; change X back to a vector

return, a
END
