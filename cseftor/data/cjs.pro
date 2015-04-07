  fid          = H5F_OPEN('/home/dhaffner/dev/OMTO3-2.1.0/OMTO3/tbl/LUT/NVAL.SO2_w12_RRS.h5')

  index        = H5D_OPEN(fid,'RingAdded')
  RingAdded    = H5D_READ(index)
  index        = H5D_OPEN(fid,'SensorZenith')
  SensorZenith = H5D_READ(index)
  index        = H5D_OPEN(fid,'SolarZenith')
  SolarZenith  = H5D_READ(index)
  index        = H5D_OPEN(fid,'c0')
  c0           = H5D_READ(index)
  index        = H5D_OPEN(fid,'c1')
  c1           = H5D_READ(index)
  index        = H5D_OPEN(fid,'c2')
  c2           = H5D_READ(index)
  index        = H5D_OPEN(fid,'f0flux')
  f0flux       = H5D_READ(index)
  index        = H5D_OPEN(fid,'kn0')
  kn0          = H5D_READ(index)
  index        = H5D_OPEN(fid,'kn1')
  kn1          = H5D_READ(index)
  index        = H5D_OPEN(fid,'kn2')
  kn2          = H5D_READ(index)
  index        = H5D_OPEN(fid,'knr')
  knr          = H5D_READ(index)
  index        = H5D_OPEN(fid,'knr2')
  knr2         = H5D_READ(index)
  index        = H5D_OPEN(fid,'pressure')
  pressure     = H5D_READ(index)
 
  nlgi0         = fltarr(6,10,21,12,4)
  nz1           = fltarr(6,10,21,12,4)
  nz2           = fltarr(6,10,21,12,4)
  ntr           = fltarr(6,10,21,12,4)
  nsb           = fltarr(21,12,4)

  openr,10,'TABLE_OMI.DAT',/f77_unformatted
  readu,10,nlgi0,nz1,nz2,ntr,nsb
  close,10

  print,'nlgi0 = ',nlgi0(*,0,0,0,0)

  outfn = 'TABLE_OMI.h5'

  H5_RingAdded                      = {_NAME:'RingAdded'        ,_TYPE:'Dataset', _DATA:RingAdded}
  H5_SensorZenith                   = {_NAME:'SensorZenith'     ,_TYPE:'Dataset', _DATA:SensorZenith}
  H5_SolarZenith                    = {_NAME:'SolarZenith'      ,_TYPE:'Dataset', _DATA:SolarZenith}
  H5_c0                             = {_NAME:'c0'               ,_TYPE:'Dataset', _DATA:c0}
  H5_c1                             = {_NAME:'c1'               ,_TYPE:'Dataset', _DATA:c1}
  H5_c2                             = {_NAME:'c2'               ,_TYPE:'Dataset', _DATA:c2}
  H5_f0flux                         = {_NAME:'f0flux'           ,_TYPE:'Dataset', _DATA:f0flux}
  H5_kn0                            = {_NAME:'kn0'              ,_TYPE:'Dataset', _DATA:kn0}
  H5_kn1                            = {_NAME:'kn1'              ,_TYPE:'Dataset', _DATA:kn1}
  H5_kn2                            = {_NAME:'kn2'              ,_TYPE:'Dataset', _DATA:kn2}
  H5_knr                            = {_NAME:'knr'              ,_TYPE:'Dataset', _DATA:knr}
  H5_knr2                           = {_NAME:'knr2'             ,_TYPE:'Dataset', _DATA:knr2}
  H5_knb                            = {_NAME:'knb'              ,_TYPE:'Dataset', _DATA:knb}
  H5_lgi0                           = {_NAME:'lgi0'             ,_TYPE:'Dataset', _DATA:nlgi0}
  H5_pressure                       = {_NAME:'pressure'         ,_TYPE:'Dataset', _DATA:pressure}
  H5_sb                             = {_NAME:'sb'               ,_TYPE:'Dataset', _DATA:nsb}
  H5_tr                             = {_NAME:'tr'               ,_TYPE:'Dataset', _DATA:ntr}
  H5_wlen                           = {_NAME:'wlen'             ,_TYPE:'Dataset', _DATA:wlen}
  H5_z1                             = {_NAME:'z1'               ,_TYPE:'Dataset', _DATA:nz1}
  H5_z2                             = {_NAME:'z2'               ,_TYPE:'Dataset', _DATA:nz2}

  APID  = {_TYPE:                   'Group',H5_RingAdded:          H5_RingAdded         , $
                                            H5_SensorZenith:       H5_SensorZenith      , $
                                            H5_SolarZenith:        H5_SolarZenith       , $
                                            H5_c0:                 H5_c0                , $
                                            H5_c1:                 H5_c1                , $
                                            H5_c2:                 H5_c2                , $
                                            H5_f0flux:             H5_f0flux            , $
                                            H5_kn0:                H5_kn0               , $
                                            H5_kn1:                H5_kn1               , $
                                            H5_kn2:                H5_kn2               , $
                                            H5_knr:                H5_knr               , $
                                            H5_knr2:               H5_knr2              , $
                                            H5_knb:                H5_knb               , $
                                            H5_lgi0:               H5_lgi0              , $
                                            H5_pressure:           H5_pressure          , $
                                            H5_sb:                 H5_sb                , $
                                            H5_tr:                 H5_tr                , $
                                            H5_wlen:               H5_wlen              , $
                                            H5_z1:                 H5_z1                , $
                                            H5_z2:                 H5_z2                }

  h5_create, outfn(0),APID

  end
