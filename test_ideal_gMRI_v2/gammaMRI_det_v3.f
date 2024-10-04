! Pablo Galve Lahoz @ Grupo de Fisica Nuclear, UCM, Madrid, Spain
! gammaMRI_v0.f: This code reads an input file with the vectors time
! distribution and it returns the detected events in a cylinder.
! v1: I add the time interpolation and inputs
! v2: Input activity image
! v3: We weights the detection prob_map by 1/R^2
! v3.1: We weights the detection prob_map by 1/R 
! v3.2: Let's try with 1/cos_th/r_xyemv/cos_emv_det (differential surface 
! accounting the dz and r_xyemv*d_phiemv projected on the orthogonal 
! direction of the detector cylinder). cos_emv_det is the cosine 
! between the emv and det in the transverse plane.
! VERSIONS 3 AND 3.1 ARE DOING NOTHING DIFFERENT
! v3.3: We add single detector of det_height x det_height/det_radius
!
! det_v1: We add many detectors in ring structure
! det_v2: We set that all the detections are placed from the center of the
! detectors.
! det_v3: I add 2 additional inputs for gyro and B0

      module gammaMRI
      parameter(pi=acos(-1.))
      real*8 amlife !mean life in seconds
      real*8 gyro, time_step, B0, wlarmor, wl_rad, T1, T2, a0, a2
      integer nx,ny,nz,nt
      character*1000 file_magvectors
      real gvec(4) !,totgvec(nx,ny,nt,4),gvec_xyt(4,nx,ny,nt)
      real,allocatable:: totgvec(:,:,:,:),gvec_xyt(:,:,:,:)
      real,allocatable:: act_image(:,:),em_image(:,:),det_image(:,:)
      real*8,allocatable:: phase_rotframe(:)
      real*8 ang_step
      real det_radius, det_height, xfov, yfov, det_rphi, det_phi
      real det_segphi
      real aux !auxiliary
      integer ndetectors
      integer,allocatable:: ndetections_det(:)
      real*8,allocatable:: prob_map_det(:,:,:)

      contains

      subroutine poisson_dist(time,act)
      real*8 act
      real*8 time,dt,rand

      call random_number(rand)
      if(act.gt.0)then
        !dt=-log(1-rand*0.999999)/act !NO USAR: rand+0.00001 (> 1.0)
        dt=-dlog(1-rand*0.999999)/act !NO USAR: rand+0.00001 (> 1.0)
      else
        dt=1
      endif
!      dt=1e-6
!      dt=poisson_time(act0)
      time=time+dt
      !act=act*exp(-dt/amlife)
      act=act*dexp(-dt/amlife)
      !return
      end subroutine !poisson_dist

      subroutine init_gvectors
      !gyro=1.37e6 !Hz/T
      !B0=0.05 !T
      T1=1.0 !s
      T2=0.2 !s
      a0=1.
      !a2=1. !a2=a0 => maximum anysotropy 
      print*,'-------------------'
      print*,'---init_gvectors---'
      print*,'-------------------'
      wlarmor=B0*gyro
      wl_rad=wlarmor*2*pi
      ang_step=wl_rad*time_step
      print*,'time_step (s)=',time_step
      !print*,'gyro (Hz/T)=',gyro
      !print*,'B0 (T)=',B0
      print*,'wlarmor = B0*gyro'
      print*,'wlarmor (Hz)=',wlarmor
      print*,'wl_rad (rad/s)=',wl_rad
      print*,'ang_step (rad)=',ang_step

      print*,'allocating totgvec, gvec_xyt'
      allocate(totgvec(nx,ny,nt,4),gvec_xyt(4,nx,ny,nt))
      allocate(phase_rotframe(nt))

      print*,'Reading ',trim(file_magvectors)
      open(13,file=trim(file_magvectors),
     o form='unformatted',access='stream')

      phase_rotframe(1)=0.
      do 333 ix=1,nx
      do 333 iy=1,ny
      do 333 iz=1,nz
      do 333 it=1,nt
       read(13)gvec
       totgvec(ix,iy,it,:)=gvec(:)
       gvec_xyt(:,ix,iy,it)=gvec(:)
       if(it.gt.1)phase_rotframe(it)=phase_rotframe(it-1)+ang_step
333   continue
      close(13)

      print*,"We write a few files..."

      open(13,file='output_xyt_gvec.raw',
     o  form='unformatted',access='stream')
      write(13)totgvec
      close(13)

      open(13,file='output_xyt_wgvec.raw',
     o  form='unformatted',access='stream')
      write(13)totgvec(:,:,:,1)
      close(13)

      open(13,file='output_xyt_xgvec.raw',
     o  form='unformatted',access='stream')
      write(13)totgvec(:,:,:,2)
      close(13)

      open(13,file='output_xyt_ygvec.raw',
     o  form='unformatted',access='stream')
      write(13)totgvec(:,:,:,3)
      close(13)

      open(13,file='output_xyt_zgvec.raw',
     o  form='unformatted',access='stream')
      write(13)totgvec(:,:,:,4)
      close(13)


      open(12,file='output_xyt_wgvec.txt')
      open(13,file='output_xyt_xgvec.txt')
      open(14,file='output_xyt_ygvec.txt')
      open(15,file='output_xyt_zgvec.txt')
      open(16,file='output_ix4iy4t_gvec.txt')
      open(18,file='output_ix1iy1t_gvec.txt')
      open(17,file='output_xyt_gvec.txt')
      open(19,file='output_phase_rotframe.txt')

      do 334 it=1,nt
       write(12,'(I10,900E15.5)')it,totgvec(:,:,it,1) !4)
       write(13,'(I10,900E15.5)')it,totgvec(:,:,it,2) !1)
       write(14,'(I10,900E15.5)')it,totgvec(:,:,it,3) !2)
       write(15,'(I10,900E15.5)')it,totgvec(:,:,it,4) !3)
       write(17,'(I10,3600E15.5)')it,gvec_xyt(:,:,:,it)
       aux=sqrt(totgvec(4,4,it,1)**2+totgvec(4,4,it,2)**2+
     o  totgvec(4,4,it,3)**2)
       write(16,'(I10,4E15.5)')it,totgvec(4,4,it,:) !,aux
       write(18,'(I10,4E15.5)')it,totgvec(1,1,it,:) !,aux
       write(19,'(I10,E15.5)')it,phase_rotframe(it)
334   continue

      close(12)
      close(13)
      close(14)
      close(15)
      close(16)
      close(17)
      close(18)
      close(19)

      end subroutine !init_gvectors

      subroutine gen_emission(ix,iy)
      integer ix,iy
      amax=maxval(act_image)
      prob=0
      aval_rand=1.
      do while(prob.lt.aval_rand)
       call random_number(rand)
       ix=rand*nx+1
       call random_number(rand)
       iy=rand*ny+1
       prob=act_image(ix,iy)
       call random_number(rand)
       aval_rand=rand*amax
      enddo
      end subroutine !gen_emission

      subroutine calcula_gvector(ix,iy,time,magvec)
      integer ix,iy
      real*8 time,dt,phase,w_tot
      real magvec(3)

      it=int(time/time_step)+1
      dt=time-(it-1)*time_step
      w_grad=totgvec(ix,iy,it,1)*1000. !rad/s
      w_tot=wl_rad+w_grad !rad/s
      phase=dble(phase_rotframe(it))+w_tot*dt !rad
      magvec(1)=totgvec(ix,iy,it,2)*cos(phase)+
     o          totgvec(ix,iy,it,3)*sin(phase)
      magvec(2)=-totgvec(ix,iy,it,2)*sin(phase)+
     o          totgvec(ix,iy,it,3)*cos(phase)
      magvec(3)=totgvec(ix,iy,it,4)

      !magvec=(/cos(aux), sin(aux), 0./)
      !magvec=(/0., 0., 1./)

      end subroutine calcula_gvector

      function aprob_em(cth)
      real cth
      c2th=2*cth**2 - 1.
      aprob_em = (a0 - a2*c2th)/(a0+a2)
      aprob_em = max(aprob_em,0.)
      end function aprob_em

      subroutine gen_direction(magvec,emvec)
      real emvec(3),magvec(3),magvec_norm(3)
      amagnorm=sqrt(magvec(1)**2 + magvec(2)**2 + magvec(3)**2)
      magvec_norm=magvec/amagnorm
      amagnorm=sqrt(magvec_norm(1)**2 + magvec_norm(2)**2 + 
     o  magvec_norm(3)**2)
      if(abs(amagnorm-1).gt.1e-6)print*,'Error amagnorm=',amagnorm
      f2th_m=0.
      rg=1.
      do while(f2th_m.lt.rg)
       call random_number(rand)
       cth=rand*2.-1
       sth=sqrt(1-cth*cth)
       call random_number(rand)
       phi=rand*2*pi
       cph=cos(phi)
       sph=sin(phi)
       emvec=(/sth*cph, sth*sph, cth/)
       vxm=0
       do ik=1,3
        vxm=vxm+emvec(ik)*magvec_norm(ik)
       enddo
       !cth_m=vxm
       !c2th_m=2*cth_m*cth_m - 1.
       !f2th_m= (a0 - a2*c2th_m)/(a0+a2)
       f2th_m=aprob_em(vxm)
       call random_number(rand)
       rg=rand
      enddo 
      end subroutine !gen_direction

      subroutine gen_position(ix,iy,xem,yem)
      integer ix,iy
      real xem,yem
      xem=(ix-0.5)*xfov/nx - xfov/2.
      yem=(iy-0.5)*yfov/ny - yfov/2.
      end subroutine gen_position

      subroutine gen_detposition(idet,xdet,ydet,zdet)
      integer idet
      real xdet,ydet,zdet
      phi=(idet-1)*det_segphi + 0.5*det_phi - pi
      xdet=det_radius*cos(phi)
      ydet=det_radius*sin(phi)
      zdet=0.
      end subroutine gen_detposition

      subroutine gen_detection(emvec,ix,iy,time,xem,yem,
     o xdet,ydet,zdet,idet,prob_map)
      real emvec(3)
      integer ix,iy,idet
      real*8 time
      real auxemvec(3),auxmagvec(3)
      real prob_map(nx,ny),xem,yem
      prob_map=-1.
      idet=-1
      xdet=0.
      ydet=0.
      zdet=0.

      !xem=(ix-0.5)*xfov/nx - xfov/2.
      !yem=(iy-0.5)*yfov/ny - yfov/2.
      call gen_position(ix,iy,xem,yem)
      zem=0.
      if(xem*xem+yem*yem.gt.det_radius*det_radius)return !emission out of scanner

      !Solving crosspoint with the cylinder...
      aa=emvec(1)**2+emvec(2)**2
      if(aa.lt.1e-10)return !no detection
      bb=2*(emvec(1)*xem+emvec(2)*yem)
      cc=xem**2+yem**2-det_radius**2
      square=bb**2-4*aa*cc
      if(square.lt.0)return !no detection
      tt=(-bb+sqrt(square))/2/aa
      if(tt.lt.0)tt=(-bb-sqrt(square))/2/aa
      xdet=xem+tt*emvec(1)
      ydet=yem+tt*emvec(2)
      zdet=zem+tt*emvec(3)
      if(abs(zdet).gt.det_height/2.)return !no detection
      phi=atan2(ydet,xdet)+pi
      if(phi.lt.0.or.phi.gt.2*pi)print*,'Error phi=',phi
      idet=int(phi/det_segphi)
      phi=phi-idet*det_segphi
      idet=idet+1
      if(phi.gt.det_phi)return
      call gen_detposition(idet,xdet2,ydet2,zdet2)
      !rphi=phi*det_radius
      !if(rphi.gt.det_rphi)return

!      r_xydet=sqrt(xdet2*xdet2 + ydet2*ydet2)
      asum=0.
      do iyk=1,ny
      do ixk=1,nx
       !xemk=(ixk-0.5)*xfov/nx - xfov/2.
       !yemk=(iyk-0.5)*yfov/ny - yfov/2.
       call gen_position(ixk,iyk,xemk,yemk)
       zemk=zem
       auxemvec(1)=xdet2-xemk
       auxemvec(2)=ydet2-yemk
       auxemvec(3)=zdet2-zemk
       amodaux=sqrt(auxemvec(1)**2+auxemvec(2)**2+auxemvec(3)**2)
!       r_xyemv=sqrt(auxemvec(1)**2+auxemvec(2)**2)
!       cos_th=abs(auxemvec(3))/amodaux+1e-6
!       cos_emv_det=(auxemvec(1)*xdet+auxemvec(2)*ydet)/r_xydet/r_xyemv
!       if(cos_emv_det.gt.1.001.or.cos_emv_det.lt.0)then
!        print*,
!     o  'Error cos_emv_det ',cos_emv_det
!        print*,'auxemvec(:)=',auxemvec(:)
!        print*,xdet,ydet,zdet
!        print*,'rxydet=',rxydet
!       endif
       auxemvec=auxemvec/amodaux 
       call calcula_gvector(ixk,iyk,time,auxmagvec)
       amodaux2=sqrt(auxmagvec(1)**2+auxmagvec(2)**2+auxmagvec(3)**2)
       auxmagvec=auxmagvec/amodaux2
       cth=auxemvec(1)*auxmagvec(1)+auxemvec(2)*auxmagvec(2)+
     o     auxemvec(3)*auxmagvec(3)
       !diffS=cos_th*r_xyemv/det_radius/cos_emv_det
       prob_map(ixk,iyk)=aprob_em(cth)*prob_map_det(ixk,iyk,idet)!/diffS
       !prob_map(ixk,iyk)=aprob_em(cth)*diffS
       asum=asum+prob_map(ixk,iyk)
       !if(prob_map(ixk,iyk).lt.0.or.prob_map(ixk,iyk).gt.1.)then
       if(prob_map(ixk,iyk).lt.0)then !.or.prob_map(ixk,iyk).gt.1.)then
        print*,'Error prob_map=',prob_map(ixk,iyk),cth
        print*,'cth=',cth
        print*,'aprob_em(cth)=',aprob_em(cth)
        print*,'auxemvec(:)=',auxemvec(:)
        print*,xdet,ydet,zdet
        print*,xdet2,ydet2,zdet2
        print*,xem,yem,zem
        print*,'amodaux = ',amodaux
        print*,'auxmagvec(:)=',auxmagvec(:)
        print*,'amodaux2 = ',amodaux2
!        print*,'r_xyemv=',r_xyemv
!        print*,'cos_emv_det=',cos_emv_det
        stop
       endif
      enddo
      enddo
      !prob_map=prob_map/asum
      end subroutine gen_detection

      end module gammaMRI

!---------------------------------------------------------------------
!-----MAIN------------------------------------------------------------
!---------------------------------------------------------------------
      program main 
      use gammaMRI
      parameter(ntbins=100,nphibins=100,nzbins=100)
      character*1000 file_out,chaux,act_file
      real*8 time,act
      integer*8 ndecays,ndetected
      integer decay_dist(ntbins)
      real magvec(3),emvec(3)
      real,allocatable:: prob_map(:,:)
      real*8,allocatable:: prob_maptot(:,:)
      integer:: zdet_dist(nzbins),phidet_dist(nphibins),
     o phizdet_dist(nphibins,nzbins)

      call random_seed

      print*,"------------------------------------"
      print*,"-------gamma MRI generator----------"
      print*,"------------------------------------"
      nparam=iargc()
      call getarg(1,file_magvectors)
      write(*,*)'1: Input file:',trim(file_magvectors)
      call getarg(2,chaux)
      read(chaux,*)time_step
      write(*,*)'2: time_step (s)=',time_step
      call getarg(3,chaux)
      read(chaux,*)nx
      write(*,*)'3: nx=',nx
      call getarg(4,chaux)
      read(chaux,*)ny
      write(*,*)'4: ny=',ny
      call getarg(5,chaux)
      read(chaux,*)nt
      write(*,*)'5: nt=',nt
      call getarg(6,chaux)
      read(chaux,*)amlife
      write(*,*)'6: amlife (s)=',amlife
      call getarg(7,chaux)
      read(chaux,*)act
      write(*,*)'7: act (Bq)=',act
      call getarg(8,chaux)
      read(chaux,*)det_radius
      write(*,*)'8: det_radius (cm)=',det_radius
      call getarg(9,chaux)
      read(chaux,*)det_height
      write(*,*)'9: det_height (cm)=',det_height
      call getarg(10,chaux)
      read(chaux,*)xfov
      write(*,*)'10: xfov (cm)=',xfov
      call getarg(11,chaux)
      read(chaux,*)yfov
      write(*,*)'11: yfov (cm)=',yfov
      call getarg(12,act_file)
      write(*,*)'12: act_file: ',trim(act_file)
      call getarg(13,chaux)
      read(chaux,*)ndetectors
      write(*,*)'13: ndetectors=',ndetectors
      call getarg(14,chaux)
      read(chaux,*)a2
      write(*,*)'14: a2=',a2
      call getarg(15,chaux)
      read(chaux,*)if_reco
      write(*,*)'15: if_reco=',if_reco
      call getarg(16,chaux)
      read(chaux,*)gyro
      write(*,*)'16: gyro (Hz/T) =',gyro
      call getarg(17,chaux)
      read(chaux,*)B0
      write(*,*)'17: B0 (T) =',B0
     
      nz=1
      det_rphi=det_height
      print*,'det_rphi (cm)=',det_rphi
      det_phi=det_rphi/det_radius
      print*,'det_phi (rad)=',det_phi
      det_segphi=2*pi/ndetectors
      print*,'det_segphi (rad)=',det_segphi
      !amlife=0.173 !s
      time_tot=nt*time_step !s
      tbin=time_tot/(ntbins)
      !act=1500000/time_tot !1/s
      decay_dist=0
      zdet_dist=0
      phidet_dist=0
      phizdet_dist=0
      !time_th=0.03
      !det_radius=4 !cm
      !det_height=16 !cm
      zbin=det_height/nzbins
      phibin=2*pi/nphibins
      !xfov=0.03 !3. !cm
      !yfov=0.03 !3. !cm
      allocate(prob_map(nx,ny),prob_maptot(nx,ny))
      allocate(prob_map_det(nx,ny,ndetectors))
      allocate(ndetections_det(ndetectors))
      print*,"Let's start Poisson list"
      !print*,"mean life (s) = ",amlife
      !print*,"activity (1/s) = ", act
      print*,"Total time (s) = ", time_tot
      !print*,"Detector radius (cm) = ", det_radius
      !print*,"Detector height (cm) = ", det_height
      !print*,"X FOV (cm) = ", xfov
      !print*,"Y FOV (cm) = ", yfov

      call init_gvectors

!-----Define act_image
      allocate(act_image(nx,ny))
      allocate(em_image(nx,ny))
      allocate(det_image(nx,ny))
      !act_image=0.
      !act_image(int(nx/2)+1,int(ny/2)+1)=1
      !act_image(int(nx/2)+2,int(ny/2)+1)=0.5
      !act_image(int(nx/2)+0,int(ny/2)+1)=0.2
      !act_image(:,2)=1.
      !act_image=1.
      open(13,file=trim(act_file),access='stream',form='unformatted')
      read(13)act_image
      close(13)
      print*,'amax=',maxval(act_image)
      em_image=0.
      det_image=0.
!------

!------if_reco
      prob_map_det=0.
      if(if_reco.eq.1)then
       open(13,file="prob_map_det.raw",access='stream',
     o  form='unformatted')
       read(13)prob_map_det
       close(13)
      endif
!------

      aux=0
!      open(13,file='list_decayskk.txt')
!      open(14,file='list_dirkk.txt')
      open(15,file='list_gammaMRI.raw',form='unformatted',
     o access='stream')
      open(16,file='list_detevents.txt')
      time=0.
      ndecays=0
      ndetected=0
      prob_maptot=0.
      ndetections_det=0
      !do i=1,1000
      do while(time.lt.time_tot)
       call gen_emission(ix,iy)
       em_image(ix,iy)=em_image(ix,iy)+1
       call poisson_dist(time,act)
       !print*,i,time,act
       itbin=time/tbin+1
       if(itbin.le.ntbins)decay_dist(itbin)=decay_dist(itbin)+1
       ndecays=ndecays+1
       if(time.gt.time_tot)cycle !poisson_dist estimated after the do while
       it_gvector=int(time/time_step)+1
       if(it_gvector.gt.nt)cycle
       call calcula_gvector(ix,iy,time,magvec)
       aux=aux+2*pi/100.
       call gen_direction(magvec,emvec)
       call gen_detection(emvec,ix,iy,time,xem,yem,
     o   xdet,ydet,zdet,idet,prob_map)
!       write(13,*) ndecays,time,act
!       write(14,*) emvec(:)
       akk=0.
       do iyk=1,ny
       do ixk=1,nx
        akk=akk+prob_map(ixk,iyk)
       enddo
       enddo
       !if(prob_map(1,1).ge.0)then
       if(akk.ge.0)then !detection
        prob_maptot=prob_maptot+prob_map
        write(15) prob_map
        phi=atan2(ydet,xdet)
        !phi=atan2(emvec(2),emvec(1))
        write(16,"(I10,2F17.7,2I10,6F17.7)")ndecays,time,act,ix,iy,
     o   xem,yem,xdet,ydet,zdet,phi
        det_image(ix,iy)=det_image(ix,iy)+1
        izbin=(zdet+det_height/2.)/zbin+1
        if(izbin.ge.1.or.izbin.le.zbin)then
         zdet_dist(izbin)=zdet_dist(izbin)+1
        endif
        iphibin=(phi+pi)/phibin+1
        phidet_dist(iphibin)=phidet_dist(iphibin)+1
        phizdet_dist(iphibin,izbin)=phizdet_dist(iphibin,izbin)+1
        ndetected=ndetected+1
        ndetections_det(idet)=ndetections_det(idet)+1
        if(if_reco.ne.1)prob_map_det(ix,iy,idet)=
     o    prob_map_det(ix,iy,idet)+1
       endif
      enddo
      print*,"ndecays = ",ndecays
      print*,"ndetected = ",ndetected
      print*,"Sensitivity = ",real(ndetected)/real(ndecays)
      print*,"End time (s) = ",time
      print*,"End activity (1/s) = ",act
!      close(13)
!      close(14)
      close(15)
      close(16)

      open(13,file='ndetections_det.txt')
      do i=1,ndetectors
       write(13,*)i,(i-1)*det_segphi+0.5*det_phi,ndetections_det(i),
     o  real(ndetections_det(i))/ndecays
      enddo
      close(13)

      open(13,file='act_dist.txt')
      do i=1,ntbins
       write(13,*)i,(i-0.5)*tbin,decay_dist(i),decay_dist(i)/tbin
      enddo
      close(13)

      open(13,file='zdet_dist.txt')
      do i=1,nzbins
       write(13,*)i,(i-0.5)*zbin-det_height/2.,zdet_dist(i),
     o  zdet_dist(i)/zbin
      enddo
      close(13)

      open(13,file='phidet_dist.txt')
      do i=1,nphibins
       write(13,*)i,(i-0.5)*phibin,phidet_dist(i),
     o  phidet_dist(i)/phibin
      enddo
      close(13)

      open(13,file='phizdet_dist.raw',access='stream',
     o form='unformatted')
      write(13)phizdet_dist
      close(13)

      open(13,file='prob_maptot.raw',access='stream',form='unformatted')
      prob_map=prob_maptot
      write(13)prob_map !tot
      close(13)

      open(13,file='em_image.raw',access='stream',form='unformatted')
      write(13)em_image
      close(13)

      open(13,file='det_image.raw',access='stream',form='unformatted')
      write(13)det_image
      close(13)

      if(if_reco.ne.1)then
      open(13,file='prob_map_det.raw',access='stream',
     o form='unformatted')
!      do idet=1,ndetectors
!       do iy=1,ny
!       do ix=1,nx 
!        prob_map(ix,iy)=dble(prob_map_det(ix,iy,idet))/ndecays
!       enddo
!       enddo
       prob_map_det=prob_map_det/ndecays
       write(13)prob_map_det
!      enddo
      close(13)
      endif
      end !main
