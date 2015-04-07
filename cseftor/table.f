      program table
c
c     purpose: read i0,z1,z2,t and s from the master tables, take slit 
c     and solar flux average for given lamda0 & lamda0+dwl, & write outputs.  
c     there are 26(ozidex)x10(sza)x6(sza0)x2(p)x3047=101016 lines in table.  
c
c     last modified 09/16/96...dave flittner
c     purpose: Change file names to varialbes; allow user to choose satellite;
c     load reflecting channel(s) into logi0,z1i0,z2i0,ti0,sb (these channels)
c     are not averaged like the other absorbing channels).
c
c
c input
      real*8 i0_in(1521,10,6), z1_in(1521,10,6), z2_in(1521,10,6),
     1 sb_in(1521), t_in(1521,10,6)
      real wlflux(2067),flux(2067)
      character*4 prf(26)
      character*5 atm(2)
c
c output
      real*4 logi0(6,10,26,5,2), z1i0(6,10,26,5,2)
      real*4 z2i0(6,10,26,5,2),tti0(6,10,26,5,2),sb(26,5,2)
      real*4 logi0_360(6,10,2), z1i0_360(6,10,2)
      real*4 z2i0_360(6,10,2), sb_360(26,5,2), ti0_360(6,10,2)
      real*4 logi0_ref(6,10,2), z1i0_ref(6,10,2)
      real*4 z2i0_ref(6,10,2), sb_ref(2), ti0_ref(6,10,2)
c
      real i0
      character*37  direc(2)
      character*60 sunfile,outfile,specfile
      character satellite*5
      integer isatellite,lun
      real*4 cntrwls(6,4)
      real*4 dlam, weight
c
      common/array/wl(1521),fitflux(1521),cntrwl(6),output(12,6),nmb(6),
     1 stiwl(1521,6),c(1521,3),csf(12,3),dlam(201),weight(201)
c     -- OMPS1 wavelengths
      data (cntrwls(i,1),i=1,6)
c    &/3087.0,3126.1,3176.2,3224.2,3313.4,3601.5/
     &/3126.1,3126.1,3176.2,3224.2,3313.4,3601.5/
c     -- OMPS2 wavelengths
      data (cntrwls(i,2),i=1,6)
     &/3120.00,3210.00,3225.00,3320.00,3360.00,3640.00/
c     -- OMPS3 wavelengths
      data (cntrwls(i,3),i=1,6)
     &/3123.40,3173.50,3310.60,3396.60,0.,0./
c     -- OMPS4 wavelengths
      data (cntrwls(i,4),i=1,6)
     &/3123.53,3174.00,3311.29,3397.32,0.,0./
c     -- Master table path
      data direc/'/home/cseftor/dev/tomrad/v2.24/atm10/',
     &           '/home/cseftor/dev/tomrad/v2.24/atm04/'/
      data pi/3.1415926/
      data prf/'L225','L275','L325','L375','L425','L475','M125','M175',
     1 'M225','M275','M325','M375','M425','M475','M525','M575',
     2 'H125','H175','H225','H275','H325','H375','H425','H475',
     3 'H525','H575'/
      data atm/'atm10','atm04'/
c
      sunfile='SOLSTICE'
      specfile='BASS_2.DAT'
c
c get the satellite name
c
      open (40,file='new_spec.dat',status='old')
      read (40,*) dlam
      read (40,*)
      read (40,*) weight
      close (40)
      do i = 1,201
        dlam(i) = dlam(i)*10.
      enddo
      isatellite=0
      do while(isatellite.lt.1.or.isatellite.gt.4)
         write(6,2000)
         read(5,*)isatellite
         if(isatellite.eq.1)then
           satellite='OMPS1'
         else if(isatellite.eq.2)then
           satellite='OMPS2'
         else if(isatellite.eq.3)then
           satellite='OMPS3'
         else if(isatellite.eq.4)then
           satellite='OMPS4'
         else
           write(6,2010)
         endif
      enddo
      do i=1,6
         cntrwl(i)=cntrwls(i,isatellite)
      enddo
c     write(outfile,'(a6,a2,52x)')'TABLE.',satellite
c
      open (30,file=sunfile,status='old')
      open (31,file='TABLE.'//satellite,status='unknown',
     1 form='unformatted')
      open (32,file=specfile,status='old')
c     -- parameters needed for calculation of itotal & partial derivatives
c     -- (phi=90 degrees)
      t = -50.
      phi = pi/2.
c     -- read solar flux data
      i = 0
3     continue
         read (30,*,end=400) aa,bb
         i = i+1
         wlflux(i) = aa*10.
         flux(i) = bb
         goto 3
400   continue
      close (30)
c     -- read wl & c coefficients from 'BASS_2.DAT': c(x,1)=c0,
c     -- c(x,2)=c1, c(x,3)=c2
      j = 0
1     continue
         read(32,*,end=800)a,b,x,d,e
         j = j+1
         wl(j) = a
         c(j,1) = b
         c(j,2) = x
         c(j,3) = d
         goto 1
800   continue
      close (32)
      do 10 i=1,1521
         jj = 0
         do 15 k=1,2067
            if (jj.ne.1 .and. wlflux(k).ge.wl(i)) then
               jj = 1
               dl = wlflux(k) - wlflux(k-1)
               fitflux(i) = flux(k-1)+(flux(k)-flux(k-1)) *
     1          (wl(i)-wlflux(k-1))/dl
            endif
15       continue
10    continue
      do 25 k=1,6
         kk=0
         if (cntrwl(k).ge.3000.) then
            do 20 j=1,1521
               if (kk.ne.1 .and. wl(j).ge.cntrwl(k)) then
                  kk=1
                  nmb(k)=j
               endif
20          continue
         endif
         write(6,*)'table. nmb=',nmb(k)
25    continue
c     -- read from the master tables
      do 85 ipres=1,2
         lun=33
         do 75 iozprf=1,26
            open (lun,file=direc(ipres)//prf(iozprf)//'.'//atm(ipres)
     1       //'.dat',status='old',form='unformatted')
            print *,direc(ipres)//prf(iozprf)//'_'//atm(ipres)//'.dat'
c           -- read table values
            call rdtabl(lun,i0_in,z1_in,z2_in,t_in,sb_in)
            close(lun)
            do 65 isza=1,10
               do 55 iscn=1,6
c                 -- keep track of record count
                  irecor=
     1             (ipres-1)*26*10*6+(iozprf-1)*6*10+(isza-1)*6+iscn
                  do 40 iwav=1,1521
                     stiwl(iwav,2) = i0_in(iwav,isza,iscn)
                     stiwl(iwav,3) = z1_in(iwav,isza,iscn)
                     stiwl(iwav,4) = z2_in(iwav,isza,iscn)
                     stiwl(iwav,5) = t_in(iwav,isza,iscn)
                     stiwl(iwav,6) = sb_in(iwav)
40                continue
                  call fsavg(irecor)
                  do 45 icntrwl=1,5
                     if (nmb(icntrwl).ne.0) then
                        logi0(iscn,isza,iozprf,icntrwl,ipres)=
     1                   alog10(output(icntrwl,2)/pi)
                        i0 = output(icntrwl,2)
                        z1i0(iscn,isza,iozprf,icntrwl,ipres)=
     1                   output(icntrwl,3) / i0
                        z2i0(iscn,isza,iozprf,icntrwl,ipres)=
     1                   output(icntrwl,4) / i0 
                        tti0(iscn,isza,iozprf,icntrwl,ipres)=
     1                   output(icntrwl,5) / i0
                        if (isza.eq.1.and.iscn.eq.1) 
     1                   sb(iozprf,icntrwl,ipres)=output(icntrwl,6)
                        totali0 = output(icntrwl,2) + 
     1                   output(icntrwl,3)*cos(phi) + 
     2                   output(icntrwl,4)*cos(2*phi) + resti0
                     endif
45                continue
                  logi0_ref(iscn,isza,ipres) = 
     1             alog10(output(6,2)/pi)
                  i0 = output(icntrwl,2)
                  z1i0_ref(iscn,isza,ipres)=
     1             output(6,3) / i0
                  z2i0_ref(iscn,isza,ipres)=
     1             output(6,4) / i0 
                  ti0_ref(iscn,isza,ipres)=
     1             output(6,5) / i0
                  if (isza.eq.1.and.iscn.eq.1) 
     1             sb_ref(ipres)=output(6,6)
55             continue
65          continue
75       continue
85    continue
c        -- output arrays i0, z1, z2, t, sb; & n-values for conparison with
c        -- current table
      write (31) logi0,z1i0,z2i0,tti0,sb
      write (31) logi0_ref,z1i0_ref,z2i0_ref,ti0_ref,sb_ref
      close (31)
      stop
1000  format(1x,f8.2,3x,f10.2)
1100  format(1x,5e12.6)
c1100 format(1p,6e12.6)
2000  format('Enter the satellite index (integer):',/,
     &'OMPS1=1',/,'OMPS2=2',/,'OMPS3=3',/,'OMPS4=4')
2010  format('Wrong number.  Try again')
3100  format(1x,6e13.6)
      end
      subroutine fsavg(irecor)
c
c     calculate the flux & slit average of i0,z1,z2,t,s; & c0,c1,c2 at
c     lamda0 
c
      common/array/wl(1521),fitflux(1521),cntrwl(6),output(12,6),nmb(6),
     1 stiwl(1521,6),c(1521,3),csf(12,3),dlam(201),weight(201)
      real sumc(12,3), dlam, weight
c     -- initialize array csf(12,3) and sumc(12,3)
c     -- index for lamda:1 to 12; index for parameter:1 to 3
      if (irecor.eq.1) then
         do 20 ii=1,3
            do 10 jj=1,12
               sumc(jj,ii) = 0.0
               csf(jj,ii) = 0.0
10          continue
20       continue
      endif
      do 200 l=1,6
         if (cntrwl(l).ne.0.) then
            k1 = nmb(l)-11
            k2 = nmb(l)+11
            xlamda0 = cntrwl(l)
            output(l,1) = xlamda0
            mm = nmb(l)+2
            do 250 m=2,6
               sum = 0.
               sumw = 0.
               do 100 j=k1,k2
                  if (fitflux(j).ne.0.0) then
                     d  = wl(j)-xlamda0
                     if (abs(d).lt.25.) then
                        do i = 1,201
                           if (d.gt.dlam(i) .and. d.lt.dlam(i+1)) then
                              frac = (d - dlam(i)) / (dlam(i+1)-dlam(i))
c                             if (m.lt.6) then
                                 s = weight(i) + 
     1                               frac*(weight(i+1)-weight(i))
                                 w = s*fitflux(j)
c                             else
c                                s = weight(i) + 
c    1                               frac*(weight(i+1)-weight(i))
c                                w = s
c                             endif
                              goto 1000
                           endif
                        enddo
1000                    continue
                        sum = sum+w * stiwl(j,m)
                        sumw = sumw+w
                        if (m.eq.2.and.irecor.eq.1) then
	                   do 25 ii=1,3
	                      sumc(l,ii) = sumc(l,ii) + w*c(j,ii)
	                      sumc(l+6,ii) = sumc(l+6,ii) + wp*c(j,ii)
25                         continue
                        endif
                     endif
                  endif
100            continue
               if (sumw.ne.0.00000) output(l,m)=sum/sumw
               if (m.eq.2.and.irecor.eq.1) then
	          do 35 ii=1,3
	             if (sumw.ne.0.00000)csf(l,ii) = sumc(l,ii)/sumw
35                continue
               endif
250         continue
         endif
200   continue
      end
      subroutine rdtabl(lun,i0_in,z1_in,z2_in,t_in,sb_in)
c     -- subroutine to read master table
      real*8 i0_in(1521,10,6), z1_in(1521,10,6), z2_in(1521,10,6), 
     1 sb_in(1521), t_in(1521,10,6)
      real*8 wavl(1521)
      real*8 sza(10),satza(6)
      integer nwavl, nscan, nsza
      read (lun) nwavl, nscan, nsza
      read (lun) satza
      read (lun) sza
      read (lun) wavl
      read (lun) i0_in
      read (lun) z1_in
      read (lun) z2_in
      read (lun) t_in
      read (lun) sb_in
      return
      end
