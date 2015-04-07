  bcw = [308.418,310.518,312.199,312.619,313.879,315.139,315.979,317.239, $
         318.079,320.180,321.020,321.860,325.220,328.160,329.000,331.101, $
         331.941,336.141,363.859,367.219,371.838,376.877]
  logi0 = fltarr(6,10,26,22,4,5)
  z1i0  = fltarr(6,10,26,22,4,5)
  z2i0  = fltarr(6,10,26,22,4,5)
  ti0   = fltarr(6,10,26,22,4,5)
  sb    = fltarr(26,22,4,5)
  openr,10,'TABLE.OMPS',/f77_unformatted
  readu,10,logi0,z1i0,z2i0,ti0,sb
  close,10
  ologi0 = fltarr(6,10,26,22,4,5)
  oz1i0  = fltarr(6,10,26,22,4,5)
  oz2i0  = fltarr(6,10,26,22,4,5)
  oti0   = fltarr(6,10,26,22,4,5)
  osb    = fltarr(26,22,4,5)
  openr,20,'TABLE.ORIG',/f77_unformatted
  readu,20,ologi0,oz1i0,oz2i0,oti0,osb
  close,20
  for i = 0,9 do begin
    for j = 0,5 do begin
      print,i,j
      plot,bcw,100.*(logi0(j,i,2,*,0,0)-ologi0(j,i,2,*,0,0))/ologi0(j,i,2,*,0,0)
      dd = ''
      read,dd
    endfor
  endfor
  end
