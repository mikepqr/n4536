function myfunc, X, P

;p[0] = sigma_y of the 1st gaussian,
;p[1] = sigma_z of the 1st gaussian,
;p[2] = sigma_y of the 2nd gaussian,
;p[3] = sigma_z of the 2nd gaussian,
;p[4] = weight1 ==> see psfcomb
;p[5] = sigma_y of the 3rd gaussian,
;p[6] = sigma_z of the 3rd gaussian,
;p[7] = weight2  ==> see psfcomb
;p[8] = y-center
;p[9] = z-center

y = X[*,0]-p[8]
z = X[*,1]-p[9]
ysq = y*y
zsq = z*z
ny = n_elements(y)
nz = n_elements(z)
psf1 = fltarr(ny,nz)
psf2 = fltarr(ny,nz)
psf3 = fltarr(ny,nz)
psfcomb = fltarr(ny,nz)
sigsq_y1 = p[0]*p[0]
sigsq_z1 = p[1]*p[1]
sigsq_y2 = p[2]*p[2]
sigsq_z2 = p[3]*p[3]
sigsq_y3 = p[5]*p[5]
sigsq_z3 = p[6]*p[6]
norm1 = 1.0/(2*!dpi*p[0]*p[1])
norm2 = 1.0/(2*!dpi*p[2]*p[3])
norm3 = 1.0/(2*!dpi*p[5]*p[6])

for i = 0,ny-1 do begin
  for j = 0,nz-1 do begin
    psf1[i,j] = norm1*exp(-0.5*y[i]*y[i]/sigsq_y1-0.5*z[j]*z[j]/sigsq_z1)
    psf2[i,j] = norm2*exp(-0.5*y[i]*y[i]/sigsq_y2-0.5*z[j]*z[j]/sigsq_z2)
    psf3[i,j] = norm3*exp(-0.5*y[i]*y[i]/sigsq_y3-0.5*z[j]*z[j]/sigsq_z3)
  endfor
endfor

psfcomb = p[4]*psf1 + p[7]*psf2 + (1-p[4]-p[7])*psf3

return, psfcomb

end



function myfunc1, X, P

;p[0] = sigma_y of the 1st gaussian,
;p[1] = sigma_z of the 1st gaussian,
;p[2] = sigma_y of the 2nd gaussian,
;p[3] = sigma_z of the 2nd gaussian,
;p[4] = weight1
;p[5] = sigma_y of the 3rd gaussian,
;p[6] = sigma_z of the 3rd gaussian,
;p[7] = weight2
;p[8] = y-center
;p[9] = z-center

y = X[*,0]-p[8]
z = X[*,1]-p[9]
ysq = y*y
zsq = z*z
ny = n_elements(y)
nz = n_elements(z)
psf1 = fltarr(ny,nz)
psf2 = fltarr(ny,nz)
psf3 = fltarr(ny,nz)
psfcomb = fltarr(ny,nz)
sigsq_y1 = p[0]*p[0]
sigsq_z1 = p[1]*p[1]
sigsq_y2 = p[2]*p[2]
sigsq_z2 = p[3]*p[3]
sigsq_y3 = p[5]*p[5]
sigsq_z3 = p[6]*p[6]
norm1 = 1.0/(2*!dpi*p[0]*p[1])
norm2 = 1.0/(2*!dpi*p[2]*p[3])
norm3 = 1.0/(2*!dpi*p[5]*p[6])

for i = 0,ny-1 do begin
  for j = 0,nz-1 do begin
    psf1[i,j] = norm1*exp(-0.5*y[i]*y[i]/sigsq_y1-0.5*z[j]*z[j]/sigsq_z1)
    psf2[i,j] = norm2*exp(-0.5*y[i]*y[i]/sigsq_y2-0.5*z[j]*z[j]/sigsq_z2)
    psf3[i,j] = norm3*exp(-0.5*y[i]*y[i]/sigsq_y3-0.5*z[j]*z[j]/sigsq_z3)
  endfor
endfor

psfcomb = p[4]*psf1 + p[7]*psf2 + (1-p[4]-p[7])*psf3

return, psf1

end




function myfunc2, X, P

;p[0] = sigma_y of the 1st gaussian,
;p[1] = sigma_z of the 1st gaussian,
;p[2] = sigma_y of the 2nd gaussian,
;p[3] = sigma_z of the 2nd gaussian,
;p[4] = weight1
;p[5] = sigma_y of the 3rd gaussian,
;p[6] = sigma_z of the 3rd gaussian,
;p[7] = weight2
;p[8] = y-center
;p[9] = z-center

y = X[*,0]-p[8]
z = X[*,1]-p[9]
ysq = y*y
zsq = z*z
ny = n_elements(y)
nz = n_elements(z)
psf1 = fltarr(ny,nz)
psf2 = fltarr(ny,nz)
psf3 = fltarr(ny,nz)
psfcomb = fltarr(ny,nz)
sigsq_y1 = p[0]*p[0]
sigsq_z1 = p[1]*p[1]
sigsq_y2 = p[2]*p[2]
sigsq_z2 = p[3]*p[3]
sigsq_y3 = p[5]*p[5]
sigsq_z3 = p[6]*p[6]
norm1 = 1.0/(2*!dpi*p[0]*p[1])
norm2 = 1.0/(2*!dpi*p[2]*p[3])
norm3 = 1.0/(2*!dpi*p[5]*p[6])

for i = 0,ny-1 do begin
  for j = 0,nz-1 do begin
    psf1[i,j] = norm1*exp(-0.5*y[i]*y[i]/sigsq_y1-0.5*z[j]*z[j]/sigsq_z1)
    psf2[i,j] = norm2*exp(-0.5*y[i]*y[i]/sigsq_y2-0.5*z[j]*z[j]/sigsq_z2)
    psf3[i,j] = norm3*exp(-0.5*y[i]*y[i]/sigsq_y3-0.5*z[j]*z[j]/sigsq_z3)
  endfor
endfor

psfcomb = p[4]*psf1 + p[7]*psf2 + (1-p[4]-p[7])*psf3

return, psf2

end



function myfunc3, X, P

;p[0] = sigma_y of the 1st gaussian,
;p[1] = sigma_z of the 1st gaussian,
;p[2] = sigma_y of the 2nd gaussian,
;p[3] = sigma_z of the 2nd gaussian,
;p[4] = weight1
;p[5] = sigma_y of the 3rd gaussian,
;p[6] = sigma_z of the 3rd gaussian,
;p[7] = weight2
;p[8] = y-center
;p[9] = z-center

y = X[*,0]-p[8]
z = X[*,1]-p[9]
ysq = y*y
zsq = z*z
ny = n_elements(y)
nz = n_elements(z)
psf1 = fltarr(ny,nz)
psf2 = fltarr(ny,nz)
psf3 = fltarr(ny,nz)
psfcomb = fltarr(ny,nz)
sigsq_y1 = p[0]*p[0]
sigsq_z1 = p[1]*p[1]
sigsq_y2 = p[2]*p[2]
sigsq_z2 = p[3]*p[3]
sigsq_y3 = p[5]*p[5]
sigsq_z3 = p[6]*p[6]
norm1 = 1.0/(2*!dpi*p[0]*p[1])
norm2 = 1.0/(2*!dpi*p[2]*p[3])
norm3 = 1.0/(2*!dpi*p[5]*p[6])

for i = 0,ny-1 do begin
  for j = 0,nz-1 do begin
    psf1[i,j] = norm1*exp(-0.5*y[i]*y[i]/sigsq_y1-0.5*z[j]*z[j]/sigsq_z1)
    psf2[i,j] = norm2*exp(-0.5*y[i]*y[i]/sigsq_y2-0.5*z[j]*z[j]/sigsq_z2)
    psf3[i,j] = norm3*exp(-0.5*y[i]*y[i]/sigsq_y3-0.5*z[j]*z[j]/sigsq_z3)
  endfor
endfor

psfcomb = p[4]*psf1 + p[7]*psf2 + (1-p[4]-p[7])*psf3

return, psf3

end


;====================================================

pro fitpsf3, fitsfile, pixscale

spec = readfits (fitsfile,0,h)
s = size(spec)
nyy = s[1]
nzz = s[2]

arr_ind = array_indices(spec,where(spec eq max(spec)))
yymax = arr_ind[0]
zzmax = arr_ind[1]
yy = findgen(nyy)
zz = findgen(nzz)

dyy2=(max(yy)-min(yy)+1)/340.
dzz2=(max(zz)-min(zz)+1)/340.
yy2=fltarr(340)
zz2=fltarr(340)
yy0=min(yy)
zz0=min(zz)
for i=0, 339 do begin
  yy2[i]=yy0+dyy2
  zz2[i]=zz0+dzz2
  yy0=yy2[i]
  zz0=zz2[i]
endfor


X = fltarr(nyy,2)
X2 = fltarr(340,2)
X[*,0] = yy
X[*,1] = zz
X2[*,0] = yy2
X2[*,1] = zz2


specnorm = spec/total(spec)  

errnorm = spec*0.0+1.0

;############################################

;guess values:
P=[3.2, 3.7, 1.2, 1.9, 0.6, 0.3, 0.4, 0.2, yymax, zzmax]

result = MPFITFUN('myfunc',X,specnorm,errnorm,P)

print, "RESULTS = ", result

;#### CALCULATE STANDARD DEVIATION ######
;P = result
;gauss = myfunc(X,P)
;mu = mean(gauss)
; 
;std = 0
;for i = 0,nyy-1 do begin
;  for j = 0,nzz-1 do begin
;    aaa = (gauss[i,j]-mu)^2
;    std = std+aaa
;  endfor
;endfor
;std=std/nyy/nzz
;stddev=sqrt(std)
; 
;chi = 0
;for i = 0,nyy-1 do begin
;  for j = 0,nzz-1 do begin
;    aaa = (gauss[i,j]-mu)^2/std
;    chi = chi+aaa
;  endfor
;endfor
;
;print, stddev, chi

loadct, 13

set_plot,"PS"
device,file = 'fitpsf3.ps',/color,/portrait, xsize=18.0 , ysize=10.8
multiplot, [2,1]

psfcomb=myfunc(X2,result)
psf1=myfunc1(X2,result)
psf2=myfunc2(X2,result)
psf3=myfunc3(X2,result)

arr_ind = array_indices(psfcomb,where(psfcomb eq max(psfcomb)))
yy2max = arr_ind[0]
zz2max = arr_ind[1]

;a cut thru yy = 0
plot, (zz-result[9])*pixscale, specnorm[yymax,*],xrange=[-0.8,0.8], xstyle=1, color = 0,psym=4, thick=2.5, xtitle='x [arcsec]', ytitle='normalised count', charthick=4.0, charsize=1.3
oplot, (zz2-result[9])*pixscale, psfcomb[yy2max,*], color = 254, thick = 2.5
oplot, (zz2-result[9])*pixscale, result[4]*psf1[yy2max,*], color = 0, thick = 2.5, linestyle=2
oplot, (zz2-result[9])*pixscale, result[7]*psf2[yy2max,*], color = 0, thick = 2.5, linestyle=2
oplot, (zz2-result[9])*pixscale, (1-result[4]-result[7])*psf3[yy2max,*], color = 0, thick = 2.5, linestyle=2

multiplot

;a cut thru zz = 0

plot, (yy-result[8])*pixscale, specnorm[*,zzmax], xrange=[-0.8,0.8], xstyle=1, color = 0, psym=4, xtitle ='y [arcsec]', charthick=4.0,charsize=1.3, thick=2.5
oplot, (yy2-result[8])*pixscale, psfcomb[*,zz2max], color = 254, thick = 2.5
oplot, (yy2-result[8])*pixscale, result[4]*psf1[*,zz2max], color = 0, thick = 2.5, linestyle=2
oplot, (yy2-result[8])*pixscale, result[7]*psf2[*,zz2max], color = 0, thick = 2.5, linestyle=2
oplot, (yy2-result[8])*pixscale, (1-result[4]-result[7])*psf3[*,zz2max], color = 0, thick = 2.5, linestyle=2

multiplot, /reset
device,/close

end



