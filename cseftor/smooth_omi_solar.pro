  tmp = fltarr(3,27001)

  openr,10,'SolarRefSpec_OPFv9_v04_TEST_ref2.dat'
  readf,10,tmp
  close,10

  wav = tmp(0,*)
  f1  = tmp(1,*)
  f2  = tmp(2,*)

  sf1 = smooth(f1,4.8)
  sf2 = smooth(f2,4.8)

  openw,20,'Smooth_flux_for_omi.dat'
  for i = 0,27000 do printf,20,wav(i),sf1(i),sf2(i),format='(f8.2,f14.5,e16.9)'
  close,20

  end
