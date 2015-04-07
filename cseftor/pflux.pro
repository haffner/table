  device,true_color=24
  device,decomposed=0

  rr = [255,000,000,200]
  gg = [255,000,000,000]
  bb = [255,000,200,000]

  tvlct,rr,gg,bb

  t1 = fltarr(3,27001)
  t2 = fltarr(3,27001)

  openr,10,'SolarRefSpec_OPFv9_v04_TEST_ref2.dat'
  readf,10,t1
  close,10

  openr,10,'Smooth_flux.dat'
  readf,10,t2
  close,10

  plot,t1(0,*),t1(1,*),/nodata,color=1,xrange=[300,301]

  oplot,t1(0,*),t1(1,*),color=2

  oplot,t2(0,*),t2(1,*),color=3,thick=3

  end
