  fid          = H5F_OPEN('TABLE_OMI.h5')

  index        = H5D_OPEN(fid,'RingAdded')
  RingAdded    = H5D_READ(index)
  RingAdded    = 1l
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
  index        = H5D_OPEN(fid,'lgi0')
  lgi0         = H5D_READ(index)
  index        = H5D_OPEN(fid,'pressure')
  pressure     = H5D_READ(index)
  index        = H5D_OPEN(fid,'sb')
  sb           = H5D_READ(index)
  index        = H5D_OPEN(fid,'tr')
  tr           = H5D_READ(index)
  index        = H5D_OPEN(fid,'wlen')
  wlen         = H5D_READ(index)
  index        = H5D_OPEN(fid,'z1')
  z1           = H5D_READ(index)
  index        = H5D_OPEN(fid,'z2')
  z2           = H5D_READ(index)

  i0     = 10.^lgi0
  z1     = i0 * (z1 + kn1)
  z2     = i0 * (z2 + kn2)
  tr     = i0 * tr * (1.0 + knr)
  i0     = i0 * (1. + kn0)
  z1     = z1/i0
  z2     = z2/i0
  tr     = tr/i0
  lgi0   = alog10(i0)

  H5_RingAdded                      = {_NAME:'RingAdded'        ,_TYPE:'Dataset', _DATA:RingAdded}
  H5_SensorZenith                   = {_NAME:'SensorZenith'     ,_TYPE:'Dataset', _DATA:SensorZenith}
  H5_SolarZenith                    = {_NAME:'SolarZenith'      ,_TYPE:'Dataset', _DATA:SolarZenith}
  H5_c0                             = {_NAME:'c0'               ,_TYPE:'Dataset', _DATA:c0}
  H5_c1                             = {_NAME:'c1'               ,_TYPE:'Dataset', _DATA:c1}
  H5_c2                             = {_NAME:'c2'               ,_TYPE:'Dataset', _DATA:c2}
  H5_kn0                            = {_NAME:'kn0'              ,_TYPE:'Dataset', _DATA:kn0}
  H5_kn1                            = {_NAME:'kn1'              ,_TYPE:'Dataset', _DATA:kn1}
  H5_kn2                            = {_NAME:'kn2'              ,_TYPE:'Dataset', _DATA:kn2}
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
                                            H5_knr2:               H5_knr2              , $
                                            H5_lgi0:               H5_lgi0              , $
                                            H5_pressure:           H5_pressure          , $
                                            H5_sb:                 H5_sb                , $
                                            H5_tr:                 H5_tr                , $
                                            H5_wlen:               H5_wlen              , $
                                            H5_z1:                 H5_z1                , $
                                            H5_z2:                 H5_z2                }

  h5_create, 'TABLE_OMI_raman.h5',APID

  end
