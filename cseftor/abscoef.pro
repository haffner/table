  bcw    = [308.418,310.518,312.199,312.619,313.879,315.139,315.979,317.239, $
            318.079,320.180,321.020,321.860,325.220,328.160,329.000,331.101, $
            331.941,336.141,363.859,363.859,363.859,363.859]

  tmp = fltarr(6,1521)
  openr,10,'BASS_2.DAT'
  readf,10,tmp
  close,10

  wavl = tmp(0,*)/10.
  c0   = tmp(1,*)
  c1   = tmp(2,*)
  c2   = tmp(3,*)

  cc0 = interpol(c0,wavl,bcw)
  cc1 = interpol(c1,wavl,bcw)
  cc2 = interpol(c2,wavl,bcw)

  openw,10,'abs_coef.cjs'

  for i = 0,21 do printf,10,bcw(i),cc0(i),cc1(i),cc2(i)

  close,10

  end
