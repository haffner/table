  fid          = H5F_OPEN('/omps/dev/cseftor/TOMSEPL2/tbl/LUT/TOMS_NVAL.EP.h5')
; fid          = H5F_OPEN('/omps/dev/cseftor/TOMSN7L2/tbl/LUT/TOMS_NVAL.N7.h5')

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
  index        = H5D_OPEN(fid,'kn0')
  kn0          = H5D_READ(index)
  index        = H5D_OPEN(fid,'kn1')
  kn1          = H5D_READ(index)
  index        = H5D_OPEN(fid,'kn2')
  kn2          = H5D_READ(index)
  index        = H5D_OPEN(fid,'knb')
  knb          = H5D_READ(index)
  index        = H5D_OPEN(fid,'knr')
  knr          = H5D_READ(index)
  index        = H5D_OPEN(fid,'knr2')
  knr2         = H5D_READ(index)
  index        = H5D_OPEN(fid,'pressure')
  pressure     = H5D_READ(index)
  index        = H5D_OPEN(fid,'wlen')
  wlen         = H5D_READ(index)
 
  lgi0         = fltarr(6,10,21,6,4)
  z1           = fltarr(6,10,21,6,4)
  z2           = fltarr(6,10,21,6,4)
  tr           = fltarr(6,10,21,6,4)
  sb           = fltarr(21,6,4)

  openr,10,'TABLE.EP',/f77_unformatted
  readu,10,lgi0,z1,z2,tr,sb
  close,10

  outfn = 'TABLE_EP.h5'

  H5_RingAdded                      = {_NAME:'RingAdded'        ,_TYPE:'Dataset', _DATA:RingAdded}
  H5_SensorZenith                   = {_NAME:'SensorZenith'     ,_TYPE:'Dataset', _DATA:SensorZenith}
  H5_SolarZenith                    = {_NAME:'SolarZenith'      ,_TYPE:'Dataset', _DATA:SolarZenith}
  H5_c0                             = {_NAME:'c0'               ,_TYPE:'Dataset', _DATA:c0}
  H5_c1                             = {_NAME:'c1'               ,_TYPE:'Dataset', _DATA:c1}
  H5_c2                             = {_NAME:'c2'               ,_TYPE:'Dataset', _DATA:c2}
  H5_kn0                            = {_NAME:'kn0'              ,_TYPE:'Dataset', _DATA:kn0}
  H5_kn1                            = {_NAME:'kn1'              ,_TYPE:'Dataset', _DATA:kn1}
  H5_kn2                            = {_NAME:'kn2'              ,_TYPE:'Dataset', _DATA:kn2}
  H5_knb                            = {_NAME:'knb'              ,_TYPE:'Dataset', _DATA:knb}
  H5_knr                            = {_NAME:'knr'              ,_TYPE:'Dataset', _DATA:knr}
  H5_knr2                           = {_NAME:'knr2'             ,_TYPE:'Dataset', _DATA:knr2}
  H5_lgi0                           = {_NAME:'lgi0'             ,_TYPE:'Dataset', _DATA:lgi0}
  H5_pressure                       = {_NAME:'pressure'         ,_TYPE:'Dataset', _DATA:pressure}
  H5_sb                             = {_NAME:'sb'               ,_TYPE:'Dataset', _DATA:sb}
  H5_tr                             = {_NAME:'tr'               ,_TYPE:'Dataset', _DATA:tr}
  H5_wlen                           = {_NAME:'wlen'             ,_TYPE:'Dataset', _DATA:wlen}
  H5_z1                             = {_NAME:'z1'               ,_TYPE:'Dataset', _DATA:z1}
  H5_z2                             = {_NAME:'z2'               ,_TYPE:'Dataset', _DATA:z2}

  APID  = {_TYPE:                   'Group',H5_RingAdded:          H5_RingAdded         , $
                                            H5_SensorZenith:       H5_SensorZenith      , $
                                            H5_SolarZenith:        H5_SolarZenith       , $
                                            H5_c0:                 H5_c0                , $
                                            H5_c1:                 H5_c1                , $
                                            H5_c2:                 H5_c2                , $
                                            H5_kn0:                H5_kn0               , $
                                            H5_kn1:                H5_kn1               , $
                                            H5_kn2:                H5_kn2               , $
                                            H5_knb:                H5_knb               , $
                                            H5_knr:                H5_knr               , $
                                            H5_knr2:               H5_knr2              , $
                                            H5_lgi0:               H5_lgi0              , $
                                            H5_pressure:           H5_pressure          , $
                                            H5_sb:                 H5_sb                , $
                                            H5_tr:                 H5_tr                , $
                                            H5_wlen:               H5_wlen              , $
                                            H5_z1:                 H5_z1                , $
                                            H5_z2:                 H5_z2                }

  h5_create, outfn(0),APID

  end
