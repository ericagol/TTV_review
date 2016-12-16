readcol,'fits_list.txt',fname,format='a',comment='#'
nfile=n_elements(fname)
skipcbv=1
for iname=0,nfile-1 do begin &$
  fits_read,fname[iname],tab,htab,exten=1 &$
  channel=fxpar(htab,'CHANNEL')&$
  quarter=fxpar(htab,'QUARTER')&$
  if(skipcbv ne 1) then begin&$
;  print,fxpar(htab,'MODULE'),fxpar(htab,'OUTPUT')&$
  cbv_read,quarter,channel,cbv &$
  print,quarter,channel &$
;  if(cbv ne 0) then begin &$
    if(iname eq 0) then vec1=cbv.vector_1 else vec1=[vec1,cbv.vector_1]&$
    if(iname eq 0) then vec2=cbv.vector_2 else vec2=[vec2,cbv.vector_2]&$
    if(iname eq 0) then vec3=cbv.vector_3 else vec3=[vec3,cbv.vector_3]&$
    if(iname eq 0) then vec4=cbv.vector_4 else vec4=[vec4,cbv.vector_4]&$
    if(iname eq 0) then vec5=cbv.vector_5 else vec5=[vec5,cbv.vector_5]&$
    if(iname eq 0) then vec6=cbv.vector_6 else vec6=[vec6,cbv.vector_6]&$
    if(iname eq 0) then vec7=cbv.vector_7 else vec7=[vec7,cbv.vector_7]&$
    if(iname eq 0) then vec8=cbv.vector_8 else vec8=[vec8,cbv.vector_8]&$
  endif &$
;  endif &$
  if(iname eq 0) then time=tbget(htab,tab,'TIME') else time=[time,tbget(htab,tab,'TIME')] &$
;  ftmp=tbget(htab,tab,'SAP_FLUX') &$
  ftmp=tbget(htab,tab,'PDCSAP_FLUX') &$
  if(iname eq 0) then fcor=ftmp/median(ftmp) else fcor=[fcor,ftmp/median(ftmp)] &$
  if(iname eq 0) then qtr=intarr(n_elements(ftmp))+quarter else qtr=[qtr,intarr(n_elements(ftmp))+quarter]&$
  stmp=tbget(htab,tab,'SAP_FLUX_ERR')&$
  inan=where(finite(stmp) eq 0) &$
  if(inan[0] eq 0) then stmp[inan]=median(stmp[where(finite(stmp) eq 1)]) &$
  if(iname eq 0) then sig=stmp/median(ftmp) else sig=[sig,stmp/median(ftmp)] &$
  cadencetmp=tbget(htab,tab,'CADENCENO')&$
  if(iname eq 0) then cadence=cadencetmp else cadence=[cadence,cadencetmp] &$
endfor
;help,fcor
indx=where(finite(fcor) eq 1)
cadence=cadence[indx]
time=time[indx]+54833d0
fcor=fcor[indx]
sig=sig[indx]
if(skipcbv ne 1) then begin  &$
 vec1=vec1[indx] &$
vec2=vec2[indx] &$
vec3=vec3[indx] &$
vec4=vec4[indx] &$
vec5=vec5[indx] &$
vec6=vec6[indx] &$
vec7=vec7[indx] &$
vec8=vec8[indx] &$
qtr = qtr[indx] &$
endif
weight=1d0/sig^2
ntime=n_elements(time)
fflat=fcor
for iqtr=1,12 do begin &$
  jqtr=where(qtr eq iqtr) &$
  if(jqtr[0] ge 0) then begin &$
    fqtr=fcor[jqtr] &$
    medqtr=median(fqtr) &$
    fqtr=fqtr/medqtr-1d0 &$
    if(skipcbv ne 1) then begin &$
;     coeff=regress(transpose([[vec1[jqtr]],[vec2[jqtr]],[vec3[jqtr]],[vec4[jqtr]],[vec5[jqtr]]]),fqtr)&$
     coeff=regress(transpose([[vec1[jqtr]],[vec2[jqtr]],[vec3[jqtr]],[vec4[jqtr]],[vec5[jqtr]],[vec6[jqtr]]]),fqtr)&$
;     coeff=regress(transpose([[vec1[jqtr]]]),fqtr)&$
;     correction1=medqtr*(coeff[0]*vec1[jqtr]) &$
     correction1=medqtr*(coeff[0]*vec1[jqtr]+coeff[1]*vec2[jqtr]+coeff[2]*vec3[jqtr]+coeff[3]*vec4[jqtr]+coeff[4]*vec5[jqtr]+coeff[5]*vec6[jqtr]) &$
     fqtr=fcor[jqtr]-correction1+median(correction1) &$
     plot,fcor[jqtr],psym=4,ys=1&$
     oplot,fqtr,psym=4,col=255 &$
;     char=get_kbrd(1)&$
;     if(char eq '0') then 
;     fqtr=fcor[jqtr]&$
    endif else begin &$
      fqtr=fcor[jqtr] &$
    endelse &$
    fflat[jqtr]=fqtr/median(fqtr) &$
;  plot,time[jqtr]-55000d0,fqtr/median(fqtr),psym=3,ys=1,yr=[.998,1.002] &$
;  oplot,time[jqtr]-55000d0,fqtr/median(fqtr),psym=3,col=255 & char=get_kbrd(1) &$
  endif&$
endfor
;indx=where(fflat gt 1.05d0)
;if(indx[0] gt fflat[indx]=1d0
;plot,time-54900d0,fflat,ys=1,psym=3
