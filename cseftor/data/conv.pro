;
 logi0 = fltarr(137280,5)
 z1i0  = fltarr(137280,5)
 z2i0  = fltarr(137280,5)
 ti0   = fltarr(137280,5)
 sb    = fltarr(2288,5)
;
 dndxlogi0 = fltarr(1647360)
 dndxz1i0  = fltarr(1647360)
 dndxz2i0  = fltarr(1647360)
 dndxti0   = fltarr(1647360)
 dndxsb    = fltarr(27456)
;
 openr,4,'TABLE.OMPS',/f77_unformatted
 readu,4, logi0, z1i0, z2i0, ti0, sb
 close,4
 save,logi0,z1i0,z2i0,ti0,sb,filename='TABLE.SAV'
 openr,4,'DNDX.OMPS',/f77_unformatted
 readu,4, dndxlogi0, dndxz1i0, dndxz2i0, dndxti0, dndxsb
 close,4
 save,dndxlogi0,dndxz1i0,dndxz2i0,dndxti0,dndxsb,filename='DNDX.SAV'
 end
