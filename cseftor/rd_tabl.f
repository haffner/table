      program rd_table
      real*8 i0_in(1521,10,6), z1_in(1521,10,6), z2_in(1521,10,6), 
     1 sb_in(1521), t_in(1521,10,6)
      real*8 wavl(1521), scan(6), sza(10)
      integer nwavl, nscan, nsza
      open (10,file='M325.atm10.dat',form='unformatted',status='old')
      call rdtabl(33,i0_in,z1_in,z2_in,t_in,sb_in)
      end
      subroutine read(lun,i0_in,z1_in,z2_in,t_in,sb_in)
c     -- subroutine to read master table
      real*8 i0_in(1521,10,6), z1_in(1521,10,6), z2_in(1521,10,6), 
     1 sb_in(1521), t_in(1521,10,6)
      real*8 wavl(1521), scan(6), sza(10)
      integer nwavl, nscan, nsza
      read (lun) nwavl, nscan, nsza
      read (lun) scan
      read (lun) sza
      read (lun) wavl
      read (lun) i0_in
      read (lun) z1_in
      read (lun) z2_in
      read (lun) t_in
      read (lun) sb_in
      return
      end
