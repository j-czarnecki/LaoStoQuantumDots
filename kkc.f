      program p3d
      use ifport
      implicit double precision(a,b,d-h,o-z)
      implicit double complex(c)
      parameter(nlmx=40000,n=3500000)
      parameter(nbm=2000,nsize=1000 000,milion=10 000 000)
      parameter ( ncc_max = 370 000 000,nmax=1 000 000)
      parameter ( npd=n)
      dimension pl(27)
      dimension r(n,4),nban(n),nbrx(n,0:12),nbs(n,0:6),nbry(n,0:12),
     > nbd(n,0:12),nba(n,0:12),rpo(n,2)
      dimension rs(27,3)
      dimension rp(3)
      dimension ifpm(128)
      dimension xg(5,2),job(8)
      dimension BM(n,1),ipiv(n),cpr2(n,1),cpr(n,1)
      dimension chbm(nbm,nbm),ebm(nbm),cmax(100),cmaxt(100,10)
      dimension ncin(0:4,2)
      parameter(lwork=2*nbm)
      double complex work(lwork),com2e(100,100),cmat2e(100,100)
      dimension rwork(3*nbm)
      integer ncc_max
      dimension nbrzeg(n,2),nbas(n)
      allocatable cnz(:),nzin(:,:),cs(:,:),ch(:,:),cnzb(:),nzinb(:,:),
     >iro(:),ico(:),eband(:,:),cnzz(:),nzz(:,:),enesx(:)
     >,cfeast(:,:),res(:),cband(:,:,:),cih(:,:)
      allocatable
     >chiaja(:),csiaja(:),ia(:),ja(:),
     >iab(:),jab(:),efeast(:)
      dimension dcont(500,8)
      dimension dcont2e(500,8)
      dimension czvecs(npd,9),ptn(2),rdm(n,2)
      integer inf(4)
      dimension ipivz(npd),ipivedsr(100)
      dimension nbaza(nbm,2),mwband(8)
      double complex tra(8,8,3),cmat(100,100),cb(100,1),ct(100),cdt(100)
      double complex transit(70,70,3),com(100,100)
      dimension dlt(-1000:1000),ej(3,3,3,3)
      allocatable cpsi(:,:,:,:)
      dimension komb(nsize,6), iperm(720,0:6)

      allocatable crc(:,:),crc2(:,:),
     >ncinl(:,:),vpt(:),vptp(:),codc(:,:,:,:),cro(:),cpot(:)
      allocate(crc(nmax,860))
      allocate(crc2(nmax,860))
      allocate(ncinl(nmax,0:860))
      allocate(cpsi(nmax/6,3,2,60))
      allocate( cnzz(milion))
      allocate(nzz(milion,2))
      allocate(enesx(milion))
       allocate(cnz(ncc_max))
       allocate(cnzb(ncc_max))
       allocate(nzin(ncc_max,2))
       allocate(iro(ncc_max))
       allocate(ico(ncc_max))
       allocate(nzinb(ncc_max,2))
       ej(1,1,1,1)=0.336
       ej(2,2,2,2)=0.336
       ej(3,3,3,3)=0.336
       ej(1,2,1,2)=0.306
       ej(1,3,1,3)=0.306
       ej(2,1,2,1)=0.306
       ej(2,3,2,3)=0.306
       ej(3,1,3,1)=0.306
       ej(3,2,3,2)=0.306
       ej(1,2,2,1)=0.015
       ej(1,3,3,1)=0.015
       ej(2,1,1,2)=0.015
       ej(2,3,3,2)=0.015
       ej(3,1,1,3)=0.015
       ej(3,2,2,3)=0.015

      dlt(0)=1
      ci=(0.0,1.0)
      pi=4*atan(1.0)
      xm=1

c     as=sqrt(1/xm/wx)*4./nxs
      ay=10/(2*nz+1)/.05292
      ax=10.0/(2*nx+1)/.05292
c     ny=nx*sqrt(wx/wy)
c     nz=nx*sqrt(wx/wz)
      !!WRITE(*,*) ,nx,ny,nz,
     >ax*.05292,az*.05292
      nx=110
      ny=50

      a=0.39/.05292
      ax=a
      ay=a

      in=0
c budujemy sieć węzłów numer jeden
      in2=0
      do iy=-ny,ny
      do ix=-nx,nx
      do id=1,3
      do is=1,-1,-2
      in=in+1
      r(in,1)=ix*a
      r(in,2)=iy*a
      r(in,3)=id
      r(in,4)=is
      num=id*2+(-is+1)/2-1
      nbas(in)=num
c     !!WRITE(*,*) id,is,num
      enddo
      enddo
      in2=in2+1
      rpo(in2,1)=ix*a
      rpo(in2,2)=iy*a
c     write(13,*) r(in,1),r(in,2)
      enddo
      enddo
      !!WRITE(*,*) ((io,isp,io=1,3),isp=1,2)

      do i=1,in
      nbrx(i,0)=0
      nbry(i,0)=0
      nbs(i,0)=0
      nbd(i,0)=0
      nba(i,0)=0
      do j=1,in
      drx=abs(r(i,1)-r(j,1))
      dry=abs(r(i,2)-r(j,2))
      drxb=(r(i,1)-r(j,1))
      dryb=(r(i,2)-r(j,2))
      dr=drx+dry
      if(i.ne.j.and.drx.gt.a*.01.and.drx.lt.a*1.1.and.dry.lt.a*.1) then
      nbrx(i,0)=nbrx(i,0)+1
      nl=nbrx(i,0)
      nbrx(i,nl)=j
      endif

      if(i.ne.j.and.drx.gt.a*.01.and.dry.gt.a*.01
     >.and.drx.lt.a*1.5.and.dry.lt.a*1.5) then
      if(drxb*dryb.gt.1) then
      nbd(i,0)=nbd(i,0)+1
      nl=nbd(i,0)
      nbd(i,nl)=j
      else
      nba(i,0)=nba(i,0)+1
      nl=nba(i,0)
      nba(i,nl)=j
      endif
      endif

      if(i.ne.j.and.dry.gt.a*.01.and.dry.lt.a*1.1.and.drx.lt.a*.1) then
      nbry(i,0)=nbry(i,0)+1
      nl=nbry(i,0)
      nbry(i,nl)=j
      endif
      if(i.ne.j.and.drx.lt.a/10.and.dry.lt.a/10) then
      nbs(i,0)=nbs(i,0)+1
      nl=nbs(i,0)
      nbs(i,nl)=j
c     !!WRITE(*,*) i,j,nbas(i),nbas(j)
      endif
      enddo
      enddo
      tl=875/27211.6
      th=40/27211.6
c     tl=th
      td=th
      delte=47/27211.6
      dso=10/27211.6
      drso=20/27211.6


      do bpz=7.89,8
      ncinl=0
      crc=0
      crc2=0
      bz=2*bpz*.0578838/27211.6
      bpx=0
      bx=2*bpx*.0578838/27211.6
      bpy=0
      by=2*bpy*.0578838/27211.6
      g=3
      ze=g/2


      xma=0.286
      v0=50/27211.6
      v1=54/27211.6
      Rga=12/.05292
      sx=80*a/2

      psz=
     >+v0*exp(-0.5*((+sx)**2+0**2)/Rga**2)
     >+v1*exp(-0.5*((-sx)**2+0**2)/Rga**2)
      inw=in
      do i=1,in
      ig1=i
      ig2=i
      if(nbas(i).eq.1.or.nbas(i).eq.2) celem2=4*tl
      if(nbas(i).eq.3.or.nbas(i).eq.4) celem2=2*tl+2*th+delte
      if(nbas(i).eq.5.or.nbas(i).eq.6) celem2=2*tl+2*th+delte

      celem2=celem2-v1*exp(-0.5*((r(i,1)-sx)**2+r(i,2)**2)/Rga**2)
      celem2=celem2-v0*exp(-0.5*((r(i,1)+sx)**2+r(i,2)**2)/Rga**2)+psz
c     !!WRITE(*,*) in
c     inw=in


      ib=nbas(i)
      kb=ib
      if(ib.eq.1.and.kb.eq.1) celem2=celem2+ze*bz/2
      if(ib.eq.2.and.kb.eq.2) celem2=celem2-ze*bz/2
      if(ib.eq.3.and.kb.eq.3) celem2=celem2+ze*bz/2
      if(ib.eq.4.and.kb.eq.4) celem2=celem2-ze*bz/2
      if(ib.eq.5.and.kb.eq.5) celem2=celem2+ze*bz/2
      if(ib.eq.6.and.kb.eq.6) celem2=celem2-ze*bz/2

      celem=1
      call dopiszl(ig1,ig2,ncinl,celem,celem2,crc,crc2,nmax)
      enddo

c do so loc
      do i=1,in
      do j=1,nbs(i,0)
      k=nbs(i,j)
      ib=nbas(i)
      kb=nbas(k)
      celem=0
      celem2=0
      if(ib.eq.1.and.kb.eq.4) celem2=ci/2
      if(ib.eq.1.and.kb.eq.6) celem2=-1/2.
      if(ib.eq.2.and.kb.eq.3) celem2=ci/2
      if(ib.eq.2.and.kb.eq.5) celem2=1/2.
      if(ib.eq.3.and.kb.eq.2) celem2=-ci/2
      if(ib.eq.3.and.kb.eq.5) celem2=ci/2
      if(ib.eq.4.and.kb.eq.1) celem2=-ci/2
      if(ib.eq.4.and.kb.eq.6) celem2=-ci/2
      if(ib.eq.5.and.kb.eq.2) celem2=1/2.
      if(ib.eq.5.and.kb.eq.3) celem2=-ci/2
      if(ib.eq.6.and.kb.eq.1) celem2=-1/2.
      if(ib.eq.6.and.kb.eq.4) celem2=ci/2
      celem2=celem2*dso/3
      if(1.eq.1) then
      if(ib.eq.1.and.kb.eq.3) celem2=celem2+ci*bx/2.
      if(ib.eq.2.and.kb.eq.4) celem2=celem2+ci*bx/2.
      if(ib.eq.3.and.kb.eq.1) celem2=celem2-ci*bx/2.
      if(ib.eq.4.and.kb.eq.2) celem2=celem2-ci*bx/2.

      if(ib.eq.1.and.kb.eq.5) celem2=celem2-ci*by/2.
      if(ib.eq.2.and.kb.eq.6) celem2=celem2-ci*by/2
      if(ib.eq.5.and.kb.eq.1) celem2=celem2+ci*by/2
      if(ib.eq.6.and.kb.eq.2) celem2=celem2+ci*by/2

      if(ib.eq.3.and.kb.eq.5) celem2=celem2+ci*bz/2
      if(ib.eq.4.and.kb.eq.6) celem2=celem2+ci*bz/2
      if(ib.eq.5.and.kb.eq.3) celem2=celem2-ci*bz/2
      if(ib.eq.6.and.kb.eq.4) celem2=celem2-ci*bz/2
      endif

      if(ib.eq.1.and.kb.eq.2) celem2=celem2+ze*bx/2
      if(ib.eq.2.and.kb.eq.1) celem2=celem2+ze*bx/2
      if(ib.eq.3.and.kb.eq.4) celem2=celem2+ze*bx/2
      if(ib.eq.4.and.kb.eq.3) celem2=celem2+ze*bx/2
      if(ib.eq.5.and.kb.eq.6) celem2=celem2+ze*bx/2
      if(ib.eq.6.and.kb.eq.5) celem2=celem2+ze*bx/2

      if(ib.eq.1.and.kb.eq.2) celem2=celem2-ci*ze*by/2
      if(ib.eq.2.and.kb.eq.1) celem2=celem2+ci*ze*by/2
      if(ib.eq.3.and.kb.eq.4) celem2=celem2-ci*ze*by/2
      if(ib.eq.4.and.kb.eq.3) celem2=celem2+ci*ze*by/2
      if(ib.eq.5.and.kb.eq.6) celem2=celem2-ci*ze*by/2
      if(ib.eq.6.and.kb.eq.5) celem2=celem2+ci*ze*by/2



c     !!WRITE(*,*) i,k,cdabs(celem2),' dso',ib,kb
c     !!WRITE(*,*) celem2*27211.6,r(ib,1),r(kb,1)
      ig1=i
      ig2=k
      call dopiszl(ig1,ig2,ncinl,celem,celem2,crc,crc2,nmax)
      enddo
      enddo



      !!WRITE(*,*) 'him'
      do i=1,in
      do j=1,nbrx(i,0)
      k=nbrx(i,j)
      celem=0
c     !!WRITE(*,*) i,j,k
      celem2=0
      if(abs(r(i,3)-r(k,3)).lt.0.01.and.abs(r(i,4)-r(k,4)).lt.0.01) then
      if(nbas(i).le.4) celem2=-tl
      if(nbas(i).ge.5) celem2=-th
      endif
      sign=-1
      if(r(i,1).gt.r(k,1))sign=1
      if(nbas(i).eq.1.and.nbas(k).eq.5) celem2=celem2-drso/2 *sign
      if(nbas(i).eq.2.and.nbas(k).eq.6) celem2=celem2-drso/2 *sign
      if(nbas(i).eq.5.and.nbas(k).eq.1) celem2=celem2+drso/2 *sign
      if(nbas(i).eq.6.and.nbas(k).eq.2) celem2=celem2+drso/2 *sign
      ig1=i
      ig2=k
      cp=
     >cdexp(-ci*bz*(
     >(R(k,1)-R(i,1))
     >*(R(k,2)+R(i,2))/4
     >-(R(k,2)-R(i,2))
     >*(R(k,1)+R(i,1))/4))
      celem2=celem2*cp
      call dopiszl(ig1,ig2,ncinl,celem,celem2,crc,crc2,nmax)
      enddo
      enddo

      do i=1,in
      do j=1,nbry(i,0)
      k=nbry(i,j)
      celem=0
      celem2=0
      if(abs(r(i,3)-r(k,3)).lt.0.01.and.abs(r(i,4)-r(k,4)).lt.0.01)then
      if(nbas(i).le.2.or.nbas(i).ge.5) celem2=-tl
      if(nbas(i).eq.3.or.nbas(i).eq.4) celem2=-th
      endif
      sign=-1
      if(r(i,2).gt.r(k,2))sign=1
      if(nbas(i).eq.1.and.nbas(k).eq.3) celem2=celem2-drso/2 *sign
      if(nbas(i).eq.2.and.nbas(k).eq.4) celem2=celem2-drso/2 *sign
      if(nbas(i).eq.3.and.nbas(k).eq.1) celem2=celem2+drso/2 *sign
      if(nbas(i).eq.4.and.nbas(k).eq.2) celem2=celem2+drso/2 *sign
      ig1=i
      ig2=k
      cp=
     >cdexp(-ci*bz*(
     >(R(k,1)-R(i,1))
     >*(R(k,2)+R(i,2))/4
     >-(R(k,2)-R(i,2))
     >*(R(k,1)+R(i,1))/4))
      celem2=celem2*cp
      call dopiszl(ig1,ig2,ncinl,celem,celem2,crc,crc2,nmax)
      enddo
      enddo

      do i=1,in
      do j=1,nbd(i,0)
      k=nbd(i,j)
      celem=0
      celem2=0
      if(abs(r(i,4)-r(k,4)).lt.0.01)then
      if(nbas(i).eq.3.and.nbas(k).eq.5) celem2=celem2-td/2
      if(nbas(i).eq.4.and.nbas(k).eq.6) celem2=celem2-td/2
      if(nbas(i).eq.5.and.nbas(k).eq.3) celem2=celem2-td/2
      if(nbas(i).eq.6.and.nbas(k).eq.4) celem2=celem2-td/2
      ig1=i
      ig2=k
      cp=
     >cdexp(-ci*bz*(
     >(R(k,1)-R(i,1))
     >*(R(k,2)+R(i,2))/4
     >-(R(k,2)-R(i,2))
     >*(R(k,1)+R(i,1))/4))
      celem2=celem2*cp
      call dopiszl(ig1,ig2,ncinl,celem,celem2,crc,crc2,nmax)
      endif
      enddo
      enddo
      do i=1,in
      do j=1,nba(i,0)
      k=nba(i,j)
      celem=0
      celem2=0
      if(abs(r(i,4)-r(k,4)).lt.0.01)then
      if(nbas(i).eq.3.and.nbas(k).eq.5) celem2=celem2+td/2
      if(nbas(i).eq.4.and.nbas(k).eq.6) celem2=celem2+td/2
      if(nbas(i).eq.5.and.nbas(k).eq.3) celem2=celem2+td/2
      if(nbas(i).eq.6.and.nbas(k).eq.4) celem2=celem2+td/2
      ig1=i
      ig2=k
      cp=
     >cdexp(-ci*bz*(
     >(R(k,1)-R(i,1))
     >*(R(k,2)+R(i,2))/4
     >-(R(k,2)-R(i,2))
     >*(R(k,1)+R(i,1))/4))
      celem2=celem2*cp
      call dopiszl(ig1,ig2,ncinl,celem,celem2,crc,crc2,nmax)
      endif
      enddo
      enddo
22    format(200psi12.6)

      m0=40

c przepisanie potencjalu na siatke dla szredingera





      call feastinit(ifpm)
      ifpm(4)=30

      jnz=0
      do 761 i=1,inw
      do 762 j=1,ncinl(i,0)
      jnz=jnz+1
      k=ncinl(i,j)
      nzin(jnz,1)=i
      nzin(jnz,2)=k
      cnz(jnz)=crc2(i,j)
      cnzb(jnz)=crc(i,j)
c     write(14,*) i,k,cnz(jnz)
762   continue
761   continue




      write(3,*) 'AAAAAR',jnz,inw
      job(1)=2
      job(2)=1
      job(3)=1
      job(5)=jnz
      job(6)=0
      do i=1,jnz
      iro(i)=nzin(i,1)
      ico(i)=nzin(i,2)
      enddo
      if(.not.allocated(ia)) allocate(ia(((inw+1))))
      if(.not.allocated(ja)) allocate(ja(jnz))
      if(.not.allocated(chiaja)) allocate(chiaja(jnz))
      if(.not.allocated(csiaja)) allocate(csiaja(jnz))
      call mkl_zcsrcoo(job,inw,chiaja,ja,ia,jnz,cnz,iro,ico,info)
      !!WRITE(*,*) inw,nelem,nmax,info,jnz
      if(.not.allocated(iab)) allocate(iab((inw+1)))
      if(.not.allocated(jab)) allocate(jab(jnz))
      call mkl_zcsrcoo(job,inw,csiaja,jab,iab,jnz,cnzb,iro,ico,info)

      do  i=1,jnz
c     !!WRITE(*,*)csiaja(i),jab(i),'s'
c     !!WRITE(*,*)chiaja(i),jab(i),'h'
      enddo
      !!WRITE(*,*) inw,nelem,nmax,info,iab(inw+1),ia(inw+1),jnz

      emin=-0.02/27.2116
      emax=+.02/27.2116


      if(.not.allocated(efeast)) allocate (efeast(m0))
      if(.not.allocated(res))allocate (res(m0))
      if(.not.allocated(cfeast))allocate (cfeast(inw,m0))
      !!WRITE(*,*) 'pfeast'
      call zfeast_hcsrgv('F',inw,chiaja,ia,ja,csiaja,iab,jab,
     >ifpm,epsout,
     > loop, emin, emax, m0,efeast,cfeast,
     >  m, res, info)

      do ip=1,in2
            do isp=1,2
                  do io=1,3
                        do i=1,inw
                              rin1=r(i,1)
                              rin2=r(i,2)
                              iorb=r(i,3)
                              is=r(i,4)
                              if(abs(rin1-rpo(ip,1))+abs(rin2-rpo(ip,2)).lt.a/10)then
                                    if(io .eq. iorb .and. isp .eq. (1-is)/2 + 1)then
                                          rdm(ip,1)=r(i,1)
                                          rdm(ip,2)=r(i,2)
                                          do ist=1,m0
                                                cpsi(ip,io,isp,ist)=cfeast(i,ist)
c     !!WRITE(*,*) ip,io,isp,ist
                                          enddo
c     write(1200+r(i,3)*10+r(i,4),222) rin1,rin2,cdabs(cfeast(i,1))
                                    endif
                              endif
                        enddo
                  enddo
            enddo
      enddo

      is=1
      do ip=1,-in2
            write(13,222) rpo(ip,1),rpo(ip,2),
            >cdabs(cpsi(ip,1,1,is)),
            >cdabs(cpsi(ip,2,1,is)),
            >cdabs(cpsi(ip,3,1,is)),
            >cdabs(cpsi(ip,1,2,is)),
            >cdabs(cpsi(ip,2,2,is)),
            >cdabs(cpsi(ip,3,2,is))
      enddo

      dcont=0
      do i=1,m
            do j=1,inw
                  ib=2*r(j,3)+r(j,4)
                  io=r(j,3)
                  is=r(j,4)
                  ib=2*(io-1)+1+(is+1)/2
                  dcont(i,ib)=dcont(i,ib)+cdabs(cfeast(j,i))**2
            enddo
      enddo
      xn=0
      do i=1,jnz
            k=nzin(i,1)
            l=nzin(i,2)
            xn=xn+dconjg(cfeast(k,1))*cfeast(l,1)*cnzb(i)
      enddo
      !!WRITE(*,*) 'z cnzb ',xn,m

      !!WRITE(*,*) 'zfe info ',info,iband,(efeast(m1)*27211.6,m1=1,m)

      do  i=1,m
            xn=0
            do ib=1,6
                  xn=xn+dcont(i,ib)
            enddo
            !!WRITE(*,*) 'm ',i,xn
            do j=1,inw
c     cfeast(j,i)=cfeast(j,i)/sqrt(xn)
            enddo
            do ib=1,6
                  dcont(i,ib)=dcont(i,ib)/xn
            enddo
      enddo
c     do ix=1,m
      write(820,242) bpx,bpy,bpz,((efeast(ix)-delte)*27211.6,(dcont(ix,ic),ic=1,6),ix=1,m)
c    >(dcont(ix,num),num=1,6)
c     dcont(i,1)=dcont(i,1)+dcont(i,2)
c     dcont(i,2)=dcont(i,3)+dcont(i,6)
c     dcont(i,3)=dcont(i,4)+dcont(i,5)
c     dcont(i,4)=dcont(i,7)+dcont(i,8)

      !!WRITE(*,*) 'zfe info ',info,iband,(efeast(m1)*27211.6,m1=1,m)


      xn=0
      do i=1,jnz
            k=nzin(i,1)
            l=nzin(i,2)
            xn=xn+dconjg(cfeast(k,1))*cfeast(l,1)*cnzb(i)
      enddo
      !!WRITE(*,*) 'z cnzb ',xn

c    >*tra(imo,jmo,1)
c    >*sls(1,mod(i,8)+1,mod(j,8)+1)

      do ix=1,m
            write(520,242) ix,efeast(ix)*27211.6,
            >(dcont(ix,num),num=1,6)
      enddo
      xnc=0
      write(11,*) 'norm ', xnc

      com=0
      dl=100/.05292
      damp=100/27211.6*.01
      do i=1,m
      do j=1,m
      do k=1,in2
      do isp=1,2
      do io=1,3
      com(i,j)=com(i,j)+rdm(k,1)*dconjg(cpsi(k,io,isp,i))
     >*cpsi(k,io,isp,j)*damp/dl
      enddo
      enddo
      enddo
      write(14,*) i,j,cdabs(com(i,j))
      enddo
      enddo
c     stop


      !################ TIME DEPENDENT CALCULATION ???
      dt=400
      w=(efeast(3)-efeast(1))
      do w0=.5,-30,.01
            w=w0/27211.6
            ct=0
            cmax=0
            ct(1)=1
            tmax=.005/2.418e-11/dt
            do tit=1,tmax
                  cmat=0
                  t=(tit-1)*dt
                  do k=1,m
                        cmat(k,k)=1
                        do l=1,m
                              cmat(k,l)=cmat(k,l)-dt/2/ci*sin(w*(t+dt))*com(k,l)
                              >*cdexp(ci*(efeast(k)-efeast(l))*(t+dt))
                        enddo
                  enddo

                  cb=0
                  do k=1,m
                        cb(k,1)=ct(k)
                        do l=1,m
                              cb(k,1)=cb(k,1)+dt/2/ci*sin(w*t)*com(k,l)*ct(l)
                              >*cdexp(ci*(efeast(k)-efeast(l))*t)
                        enddo
                  enddo

                  call zgesv(m,1,cmat,100,ipiv,cb,100,info)
                  do k=1,m
                        ct(k)=cb(k,1)
                  enddo
c     write(14,242)
c    > it*1.0,(cdabs(ct(i))**2,i=1,m),info*1.0
            do i=1,m
                  if(cdabs(ct(i)).gt.cdabs(cmax(i))) cmax(i)=cdabs(ct(i))
            enddo
            do i=1,10
                  if(abs(tit-tmax/i).lt.0.5) then
                        write(100+i,242)w0,t*2.411e-11,(cdabs(cmax(k))**2,k=1,m)
                  endif
            enddo
      enddo


      write(15,242)
     > w0,(cdabs(cmax(i))**2,i=1,m)

      enddo
      enddo
      allocate(codc(m,m,m,m))
      allocate(cro(in2))
      allocate(cpot(in2))
      codc=0
      eps=100

      nost=1

       !!WRITE(*,*) nost,no

      call srand (1233)

      wlos=0
      r12los=0
      iflag=0
      r1los=0
      do 921 ilos=1,100000
      x1=(rand(iflag)-0.5)*ax
      x2=(rand(iflag)-0.5)*ax
      z2=(rand(iflag)-0.5)*az
      y1=(rand(iflag)-0.5)*ax
      y2=(rand(iflag)-0.5)*ax
      z1=(rand(iflag)-0.5)*az
      r12=sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
      r1=sqrt((x1)**2+(y1)**2+z1**2)
      wlos=wlos+1/r12
      r1los=r1los+1/r1
      r12los=r12los+r12
c     write(15,*) wlos/ilos,dx,r12los/ilos
921   continue
      Zupar=wlos/ilos

      if(1.eq.1) then
      do 871 i=1,m
      do 871 k=1,m

      xn=0
      do 8179 ir1=1,in2
      cl1u=0
      cl1d=0
      cr1u=0
      cr1d=0
      do iorb=1,3
      cl1u=cl1u+cpsi(ir1,iorb,1,i)
      cl1d=cl1d+cpsi(ir1,iorb,2,i)
      cr1u=cr1u+cpsi(ir1,iorb,1,k)
      cr1d=cr1d+cpsi(ir1,iorb,2,k)
      enddo
      do 8179 iorb=1,3
      cro(ir1)=dconjg(cl1u)*cr1u+dconjg(cl1d)*cr1d
8179   continue
c     do ir1=1,in2
c     xn=xn+cro(ir1)
c     enddo
c     stop

      cpot=0
      wyklucz=1e-3
      xn=0
      do 878 ir1=1,in2
      do 878 ir2=1,in2

      r12=sqrt((r(ir1,1)-r(ir2,1))**2+
     >    +(r(ir1,2)-r(ir2,2))**2)

      odl=1/eps/(r12+1/(wyklucz)*exp(-r12*10))
c     odl=1
      cpot(ir2)=cpot(ir2)+odl*cro(ir1)
878   continue

      do 872 j=1,m
      do 872 l=1,m
      iw=i
      jw=j
      kw=k
      lw=l
      do 873 ir2=1,in2
      cl2u=0
      cl2d=0
      cr2u=0
      cr2d=0
      do iorb=1,3
      cl2u=cl2u+cpsi(ir2,iorb,1,j)
      cl2d=cl2d+cpsi(ir2,iorb,2,j)
      cr2u=cr2u+cpsi(ir2,iorb,1,l)
      cr2d=cr2d+cpsi(ir2,iorb,2,l)
      enddo
      cro2=dconjg(cl2u)*cr2u+dconjg(cl2d)*cr2d
      codc(iw,jw,kw,lw)=codc(iw,jw,kw,lw)+cpot(ir2)
     >*cro2
873   continue
872   continue

      write(23,*) iw,jw,kw,lw
871   continue
      do i=1,m
      do j=1,m
      write(23,*) i,j
      do k=1,m
      do l=1,m
      do ir2=1,in2
      do iorbi=1,3
      do iorbj=1,3
      do iorbk=1,3
      do iorbl=1,3
      do isp1=1,2
      do isp2=1,2
      codc(i,j,k,l)=codc(i,j,k,l)+
     > dconjg(cpsi(ir2,iorbi,isp1,i))
     >*dconjg(cpsi(ir2,iorbj,isp2,j))
     >*cpsi(ir2,iorbk,isp1,k)
     >*cpsi(ir2,iorbl,isp2,l)
     >*ej(iorbi,iorbj,iorbk,iorbl)/eps
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      do  i=1,m
      do  j=1,m
      do  k=1,m
      do  l=1,m
c     write(920,9121) i,j,k,l,cdabs(codc(i,j,k,l))
      enddo
      enddo
      enddo
      enddo
c     write(920,*) '****'
      endif
c     stop
      write(920,*) cdabs(codc(1,1,1,1))*27211.6

      ne=2
      call perm(ne,iperm)
      lperm=ddgamma(ne)
      call wk(ne,m,komb,numerod,nsize)


      cnzz=0
      jnzz=0
      nzz=0

      do  i=1,numerod
      ene=0
      do 111 ie=1,ne
      ene=ene+efeast(komb(i,ie))
111   continue
      jnzz=jnzz+1
      cnzz(jnzz)=ene
      nzz(jnzz,1)=i
      nzz(jnzz,2)=i
      enddo
      !!WRITE(*,*) 'hiam'

      do 3 i=1,numerod
      do 3 j=i,numerod
c     !!WRITE(*,*) i,j
      ccod=0
      do 4 i1=1,ne
      do 4 i2=i1+1,ne
c petla po permutacjach elementu l
      do 5 ip=1,lperm
c     !!WRITE(*,*) i1,i2,ip
c sprawdzenie delt dla elektronow
c o numerach innych niz i oraz j
      llt=0
      do 6 ie=1,ne
      if(ie.ne.i1.and.ie.ne.i2.and.
     >komb(i,ie).eq.komb(j,iperm(ip,ie)))
     > llt=llt+1
6     continue

      if(llt.eq.ne-2) then
      iaa=komb(i,i1)
      ibb=komb(i,i2)
      icc=komb(j,iperm(ip,i1))
      idd=komb(j,iperm(ip,i2))
      ccod=ccod+codc(iaa,ibb,icc,idd)*(-1)**iperm(ip,0)
      endif

5     continue
4     continue
      if(abs(ccod)*27211.6.gt.1e-6) then
      jnzz=jnzz+1
      cnzz(jnzz)=ccod
      nzz(jnzz,1)=i
      nzz(jnzz,2)=j
c     !!WRITE(*,*) i,j,ccod

      if(j.ne.i) then
      jnzz=jnzz+1
      cnzz(jnzz)=dconjg(ccod)
      nzz(jnzz,1)=j
      nzz(jnzz,2)=i
      endif
      endif

3     continue

      num=numerod
      nstt=20
      if(nstt.gt.numerod) nstt=numerod
      !!WRITE(*,*) numerod,ncc_max,jnzz,nstt,' w2e'
      if(allocated(cih)) deallocate(cih)
      if(.not.allocated(cih)) allocate(cih(numerod,nstt))
      call w2e(num,cnzz,nzz,milion,numerod,jnzz,enesx,nstt,cih)
      dcont2e=0
      do i=1,nstt
      do ik=1,numerod
      do ie=1,ne
      iaa=komb(ik,ie)
      do ic=1,6
      dcont2e(i,ic)=cdabs(cih(ik,i))**2*dcont(iaa,ic)+
     >dcont2e(i,ic)
      enddo
      enddo
      enddo
      enddo

      write(822,221) bpz,(enesx(i)*27211.6-ne*delte*27211.6,
     >(dcont2e(i,ic),ic=1,6),i=1,nstt)

      com2e=0
      do i=1,nstt
      do j=1,nstt
      do ip=1,1
      do jp=1,lperm
      do iwi=1,numerod
      do iwj=1,numerod
      do ie=1,ne

      llt=0
      do ie1=1,ne
      if(ie.ne.ie1.and.
     >komb(iwi,iperm(ip,ie1)).eq.komb(iwj,iperm(jp,ie1)))
     > llt=llt+1
      enddo
      if(llt.eq.ne-1) then
      com2e(i,j)=com2e(i,j)+(-1)**(iperm(ip,0)+iperm(jp,0))
     >*dconjg(cih(iwi,i))*
     >cih(iwj,j)*com(komb(iwi,iperm(ip,ie)),komb(iwj,iperm(jp,ie)))
      endif

      enddo
      enddo
      enddo
      enddo
      enddo
      write(16,*) i,j,cdabs(com2e(i,j))
      enddo
      enddo
c     stop
      dt=400
      do w0=0.001,0.7,.2/2400
      w=w0/27211.6
      ct=0
      cmax=0
      ct(1)=1
      tmax=.005/2.418e-11/dt
      do tit=1,tmax
      cmat=0
      t=(tit-1)*dt
      do k=1,m
      cmat(k,k)=1
      do l=1,m
      cmat(k,l)=cmat(k,l)-dt/2/ci*sin(w*(t+dt))*com2e(k,l)
     >*cdexp(ci*(enesx(k)-enesx(l))*(t+dt))
      enddo
      enddo
      cb=0
      do k=1,m
      cb(k,1)=ct(k)
      do l=1,m
      cb(k,1)=cb(k,1)+dt/2/ci*sin(w*t)*com2e(k,l)*ct(l)
     >*cdexp(ci*(enesx(k)-enesx(l))*t)
      enddo
      enddo
      call zgesv(m,1,cmat,100,ipiv,cb,100,info)
      do k=1,m
      ct(k)=cb(k,1)
      enddo
c     write(14,242)
c    > it*1.0,(cdabs(ct(i))**2,i=1,m),info*1.0
      do i=1,m
      if(cdabs(ct(i)).gt.cdabs(cmax(i))) cmax(i)=cdabs(ct(i))
      enddo
      do i=1,10
      if(abs(tit-tmax/i).lt.0.5) then
      write(100+i,242)w0,t*2.411e-11,(cdabs(cmax(k))**2,k=1,m)
      endif
      enddo
      enddo
      write(17,242)
     > w0,(cdabs(cmax(i))**2,i=1,m)
      enddo
      deallocate(codc)
      deallocate(cpot)
      deallocate(cro)



221   format(1200g16.8)



9121  format(4i10,1g20.12)

      deallocate(ia)
      deallocate(ja)
      deallocate(chiaja)
      deallocate(csiaja)
      deallocate (cfeast)
      deallocate(iab)
      deallocate(jab)
      deallocate (efeast)
      deallocate (res)
4321  continue
      continue

222   format(1200f20.12)
242   format(1200g20.12)
      end program

      function f1(xi)
      implicit double precision(a-h,o-z)
      f1=0.5*xi*(xi-1)
      end

      function f2(xi)
      implicit double precision(a-h,o-z)
      f2=-(xi-1)*(xi+1)
      end

      function f3(xi)
      implicit double precision(a-h,o-z)
      f3=0.5*xi*(xi+1)
      end

      function psi1(xi)
      implicit double precision(a-h,o-z)
      psi1=(1-xi)/2
      end

      function psi2(xi)
      implicit double precision(a-h,o-z)
      psi2=(1+xi)/2
      end

      function phi(i,xi1,xi2,xi3)
      implicit double precision(a-h,o-z)
      if(i.eq.1) then
      phi=f1(xi1)*f1(xi2)*f1(xi3)
      return
      endif
      if(i.eq.2) then
      phi=f3(xi1)*f1(xi2)*f1(xi3)
      return
      endif
      if(i.eq.3) then
      phi=f3(xi1)*f3(xi2)*f1(xi3)
      return
      endif
      if(i.eq.4) then
      phi=f1(xi1)*f3(xi2)*f1(xi3)
      return
      endif
      if(i.eq.5) then
      phi=f1(xi1)*f1(xi2)*f3(xi3)
      return
      endif
      if(i.eq.6) then
      phi=f3(xi1)*f1(xi2)*f3(xi3)
      return
      endif
      if(i.eq.7) then
      phi=f3(xi1)*f3(xi2)*f3(xi3)
      return
      endif
      if(i.eq.8) then
      phi=f1(xi1)*f3(xi2)*f3(xi3)
      return
      endif
      if(i.eq.9) then
      phi=f2(xi1)*f1(xi2)*f1(xi3)
      return
      endif
      if(i.eq.10) then
      phi=f1(xi1)*f2(xi2)*f1(xi3)
      return
      endif
      if(i.eq.11) then
      phi=f2(xi1)*f2(xi2)*f1(xi3)
      return
      endif
      if(i.eq.12) then
      phi=f3(xi1)*f2(xi2)*f1(xi3)
      return
      endif
      if(i.eq.13) then
      phi=f2(xi1)*f3(xi2)*f1(xi3)
      return
      endif
      if(i.eq.14) then
      phi=f1(xi1)*f1(xi2)*f2(xi3)
      return
      endif
      if(i.eq.15) then
      phi=f2(xi1)*f1(xi2)*f2(xi3)
      return
      endif
      if(i.eq.16) then
      phi=f3(xi1)*f1(xi2)*f2(xi3)
      return
      endif
      if(i.eq.17) then
      phi=f1(xi1)*f2(xi2)*f2(xi3)
      return
      endif
      if(i.eq.18) then
      phi=f2(xi1)*f2(xi2)*f2(xi3)
      return
      endif
      if(i.eq.19) then
      phi=f3(xi1)*f2(xi2)*f2(xi3)
      return
      endif
      if(i.eq.20) then
      phi=f1(xi1)*f3(xi2)*f2(xi3)
      return
      endif
      if(i.eq.21) then
      phi=f2(xi1)*f3(xi2)*f2(xi3)
      return
      endif
      if(i.eq.22) then
      phi=f3(xi1)*f3(xi2)*f2(xi3)
      return
      endif
      if(i.eq.23) then
      phi=f2(xi1)*f1(xi2)*f3(xi3)
      return
      endif
      if(i.eq.24) then
      phi=f1(xi1)*f2(xi2)*f3(xi3)
      return
      endif
      if(i.eq.25) then
      phi=f2(xi1)*f2(xi2)*f3(xi3)
      return
      endif
      if(i.eq.26) then
      phi=f3(xi1)*f2(xi2)*f3(xi3)
      return
      endif
      if(i.eq.27) then
      phi=f2(xi1)*f3(xi2)*f3(xi3)
      return
      endif
      end


      function phil(i,xi1,xi2,xi3)
      implicit double precision(a-h,o-z)
      if(i.eq.1) then
      phil=psi1(xi1)*psi1(xi2)*psi1(xi3)
      return
      endif
      if(i.eq.2) then
      phil=psi2(xi1)*psi1(xi2)*psi1(xi3)
      return
      endif
      if(i.eq.3) then
      phil=psi2(xi1)*psi2(xi2)*psi1(xi3)
      return
      endif
      if(i.eq.4) then
      phil=psi1(xi1)*psi2(xi2)*psi1(xi3)
      return
      endif
      if(i.eq.5) then
      phil=psi1(xi1)*psi1(xi2)*psi2(xi3)
      return
      endif
      if(i.eq.6) then
      phil=psi2(xi1)*psi1(xi2)*psi2(xi3)
      return
      endif
      if(i.eq.7) then
      phil=psi2(xi1)*psi2(xi2)*psi2(xi3)
      return
      endif
      if(i.eq.8) then
      phil=psi1(xi1)*psi2(xi2)*psi2(xi3)
      return
      endif
      end


      subroutine wypisz(i,j,jnz,celem,cnz,nz)
      double complex cnz(100 000),celem
      dimension nz(100 000,2)
      do 1 k=1,jnz
      if(nz(k,1).eq.i.and.nz(k,2).eq.j) then
      ibyl=1
      celem=cnz(k)
      return
      endif
1     continue
      end


      subroutine cc_mv ( m, n, ncc, icc, ccc, acc, x, b )

c*********************************************************************72
c
cc CC_MV multiplies a CC matrix by a vector
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 July 2014
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Iain Duff, Roger Grimes, John Lewis,
c    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
c    October 1992
c
c  Parameters:
c
c    Input, integer M, the number of rows.
c
c    Input, integer N, the number of columns.
c
c    Input, integer NCC, the number of CC values.
c
c    Input, integer ICC(NCC), the CC rows.
c
c    Input, integer CCC(N+1), the compressed CC columns
c
c    Input, complex ACC(NCC), the CC values.
c
c    Input, complex X(N), the vector to be multiplied.
c
c    Output, complex B(M), the product A*X.
c
      implicit none

      integer m
      integer n
      integer ncc

      double complex acc(ncc)
      double complex b(m)
      integer ccc(n+1)
      integer i
      integer icc(ncc)
      integer j
      integer k
      double complex x(n)

      do i = 1, m
        b(i) = 0.0E+00
      end do

      do j = 1, n
        do k = ccc(j), ccc(j+1) - 1
          i = icc(k)
          b(i) = b(i) + acc(k) * x(j)
        end do
      end do

      return
      end
      subroutine cc_print ( m, n, ncc, icc, ccc, acc, title )

c*********************************************************************72
c
cc CC_PRINT prints a sparse matrix in CC format.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 July 2014
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows in the matrix.
c
c    Input, integer N, the number of columns in the matrix.
c
c    Input, integer NCC, the number of CC elements.
c
c    Input, integer ICC(NCC), the CC rows.
c
c    Input, integer CCC(N+1), the compressed CC columns.
c
c    Input, complex ACC(NCC), the CC values.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer n
      integer ncc

      double complex acc(ncc)
      integer ccc(n+1)
      integer i
      integer icc(ncc)
      integer j
      integer jnext
      integer k
      integer m
      character * ( * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) '     #     I     J       Ar                Ai'
      write ( *, '(a)' )
     &  '  ----  ----  ----  --------------  --------------'
      write ( *, '(a)' ) ' '

      j = 1
      jnext = ccc(2)

      do k = 1, ncc

        i = icc(k)

10      continue

        if ( jnext .le. k ) then
          j = j + 1
          jnext = ccc(j+1)
          go to 10
        end if

        write ( *, '(2x,i4,2x,i4,2x,i4,2x,g16.8,2x,g16.8)' )
     &    k, i, j, cdabs(acc(k))

      end do

      return
      end

      subroutine cnv(isp,asp,acc,icc,ccc,nmaxcc,ncc,m)
      double complex acc(nmaxcc),asp(nmaxcc)
      integer isp(nmaxcc,2),icc(nmaxcc),ccc(nmaxcc)
      double complex aspc(nmaxcc)
      integer ispc(nmaxcc,2)
      icznik=0
      do 122 io=1,ncc
c     !!WRITE(*,*) io,cdabs(asp(io)),isp(io,1),isp(io,2)
122   continue
c     stop
      do 2 jk=1,m
      do 2 i=1,ncc
      if(isp(i,2).eq.jk) then
      icznik=icznik+1
      ispc(icznik,1)=isp(i,1)
      ispc(icznik,2)=isp(i,2)
      aspc(icznik)=asp(i)
c     !!WRITE(*,*) i,jk,icznik,isp(i,1),isp(i,2),asp(i)
      endif
c3     continue
2     continue
      do 30 i=1,ncc
c     !!WRITE(*,*) i,ispc(i,1),ispc(i,2),aspc(i)
30    continue
c      stop

      do 3 i=1,ncc
      acc(i)=aspc(i)
      icc(i)=ispc(i,1)
3     continue
      ccc(1)=1
      do 4 jk=1,m
      do 5 i=1,ncc
      if(ispc(i,2).gt.jk) then
      ccc(jk+1)=i
      goto 4
      endif
5     continue
4     continue
      ccc(m+1)=ncc+1
      end



C**********************************************************************
C
      SUBROUTINE ZRANDN (N,ZX,SEED)
C
C     Purpose:
C     Fills the vector ZX with random numbers  between 0 and 1.  If the
C     SEED is given, it should be odd and positive.  The generator is a
C     fairly unsophisticated one, from Pearson's  "Numerical methods in
C     engineering and science" book.
C
C     Parameters:
C     N    = the dimension of the vector (input).
C     ZX   = the vector to fill with random numbers (output).
C     SEED = the seed for the generator (input).
C
C     Noel M. Nachtigal
C     April 23, 1993
C
C**********************************************************************
C
      INTRINSIC DBLE, DCMPLX, IABS, MOD
C
      INTEGER N, SEED
      DOUBLE COMPLEX ZX(N)
C
C     Local variables.
C
      INTEGER I, J
      DOUBLE PRECISION IMAGX, REALX
C
C     Local variables that are saved from one call to the next.
C
      DOUBLE PRECISION DMAX
      INTEGER IM, IMAX, IS
      SAVE DMAX, IM, IMAX, IS
      DATA IM/0/
C
C     Initialize the generator data.
C
      IF (IM.EQ.0) THEN
         J  = 0
         IM = 1
         DO 10 I = 1, 31
            J = J + 1
            IF (IM*2.LE.IM) GO TO 20
            IM = IM * 2
 10      CONTINUE
 20      IMAX = (IM-1) * 2 + 1
         DMAX = DBLE(IMAX)
         DO 30 I = 1, MOD(J,3)
            J = J - 1
            IM = IM / 2
 30      CONTINUE
         IM = IM + 5
         IS = IABS(MOD(IM*30107,IMAX))
      END IF
C
C     Check whether there is a new seed.
C
      IF (SEED.GT.0) IS = (SEED / 2) * 2 + 1
C
C     Here goes the rest.
C
      DO 40 I = 1, N
         REALX = DBLE(IS) / DMAX
         IS    = IABS(MOD(IM*IS,IMAX))
         IMAGX = DBLE(IS) / DMAX
         IS    = IABS(MOD(IM*IS,IMAX))
         ZX(I) = DCMPLX(REALX,IMAGX)
 40   CONTINUE
C
      RETURN
      END
C
C**********************************************************************
C**********************************************************************
C
C     Copyright (C) 1992  Roland W. Freund and Noel M. Nachtigal
C     All rights reserved.
C
C     This code is part of a copyrighted package.  For details, see the
C     file "cpyrit.doc" in the top-level directory.
C
C     *****************************************************************
C     ANY USE OF  THIS CODE CONSTITUTES ACCEPTANCE OF  THE TERMS OF THE
C                             COPYRIGHT NOTICE
C     *****************************************************************
C
C**********************************************************************
C
C     This file contains the routine for the TFQMR algorithm.
C
C**********************************************************************
C
      SUBROUTINE ZUTFX (NDIM,NLEN,NLIM,VECS,TOL,INFO)
C
C     Purpose:
C     This subroutine uses the TFQMR algorithm to solve linear systems.
C     It runs the  algorithm to convergence  or until  a user-specified
C     limit on the number of iterations is reached.
C
C     The  code is  set up  to solve  the system  A x = b  with initial
C     guess x_0 = 0.  Here  A x = b  denotes the preconditioned system,
C     and  it  is  connected with the  original system as follows.  Let
C     B y = c be the original unpreconditioned system to be solved, and
C     let y_0 be an arbitrary initial guess for its solution.  Then:
C          A x = b, where  A = M_1^{-1} B M_2^{-1},
C                          x = M_2 (y - y_0), b = M_1^{-1} (c - B y_0).
C     Here M = M_1 M_2 is the preconditioner.
C
C     To recover the final iterate y_n for  the original system B y = c
C     from the final iterate x_n for the preconditioned system A x = b,
C     set
C               y_n = y_0 + M_2^{-1} x_n.
C
C     The algorithm was first described  in the RIACS Technical Report
C     91.18, `A Transpose-Free  Quasi-Minimal  Residual Algorithm  for
C     Non-Hermitian Linear Systems`, by Roland Freund, September 1991,
C     which subsequently appeared in  SIAM J. Sci. Comput., 14 (1993),
C     pp. 470--482.
C
C     Parameters:
C     For a description of  the parameters, see the file "zutfx.doc" in
C     the current directory.
C
C     External routines used:
C     double precision dlamch(ch)
C        LAPACK routine, computes machine-related constants.
C     double precision dznrm2(n,x,incx)
C        BLAS-1 routine, computes the 2-norm of x.
C     subroutine zaxpby(n,z,a,x,b,y)
C        Library routine, computes z = a * x + b * y.
C     double precision zdotu(n,x,incx,y,incy)
C        BLAS-1 routine, computes y^H * x.
C     subroutine zrandn(n,x,seed)
C        Library routine, fills x with random numbers.
C
C     Noel M. Nachtigal
C     April 13, 1993
C
C**********************************************************************
C
      INTRINSIC CDABS, DBLE, DSQRT, MAX0
      EXTERNAL DLAMCH, DZNRM2, ZAXPBY, ZDOTU, ZRANDN
      DOUBLE COMPLEX ZDOTU
      DOUBLE PRECISION DLAMCH, DZNRM2
C
      INTEGER INFO(4), NDIM, NLEN, NLIM
      DOUBLE COMPLEX VECS(NDIM,9)
      DOUBLE PRECISION TOL
C
C     Miscellaneous parameters.
C
      DOUBLE COMPLEX ZONE, ZZERO
      PARAMETER (ZONE = (1.0D0,0.0D0),ZZERO = (0.0D0,0.0D0))
      DOUBLE PRECISION DHUN, DONE, DTEN, DZERO
      PARAMETER (DHUN = 1.0D2,DONE = 1.0D0,DTEN = 1.0D1,DZERO = 0.0D0)
C
C     Local variables, permanent.
C
      INTEGER IERR, N, RETLBL, TF, TRES, VF
      SAVE    IERR, N, RETLBL, TF, TRES, VF
      DOUBLE COMPLEX ALPHA, BETA, ETA, RHO
      SAVE           ALPHA, BETA, ETA, RHO
      DOUBLE PRECISION COS, VAR, R0, RESN, TAU, UCHK, UNRM
      SAVE             COS, VAR, R0, RESN, TAU, UCHK, UNRM
C
C     Local variables, transient.
C
      INTEGER INIT, REVCOM
      DOUBLE COMPLEX ZTMP
      DOUBLE PRECISION DTMP
C
C     Initialize some of the permanent variables.
C
      DATA RETLBL /0/
C
C     Check the reverse communication flag to see where to branch.
C        REVCOM   RETLBL      Comment
C           0        0    first call, go to label 10
C           1       30    returning from AXB, go to label 30
C           1       40    returning from AXB, go to label 40
C           1       60    returning from AXB, go to label 60
C           1       70    returning from AXB, go to label 70
C
      REVCOM  = INFO(2)
      INFO(2) = 0
      IF (REVCOM.EQ.0) THEN
         N = 0
         IF (RETLBL.EQ.0) GO TO 10
      ELSE IF (REVCOM.EQ.1) THEN
         IF (RETLBL.EQ.30) THEN
            GO TO 30
         ELSE IF (RETLBL.EQ.40) THEN
            GO TO 40
         ELSE IF (RETLBL.EQ.60) THEN
            GO TO 60
         ELSE IF (RETLBL.EQ.70) THEN
            GO TO 70
         END IF
      END IF
      IERR = 1
      GO TO 90
C
C     Check whether the inputs are valid.
C
 10   IERR = 0
      IF (NDIM.LT.1)    IERR = 2
      IF (NLEN.LT.1)    IERR = 2
      IF (NLIM.LT.1)    IERR = 2
      IF (NLEN.GT.NDIM) IERR = 2
      IF (IERR.NE.0) GO TO 90
C
C     Extract from INFO the output units TF and VF, the true residual
C     flag TRES, and the left starting vector flag INIT.
C
      VF   = MAX0(INFO(1),0)
      INIT = VF / 100000
      VF   = VF - INIT * 100000
      TRES = VF / 10000
      VF   = VF - TRES * 10000
      TF   = VF / 100
      VF   = VF - TF * 100
C
C     Check the convergence tolerance.
C
      IF (TOL.LE.DZERO) TOL = DSQRT(DLAMCH('E'))
C
C     Start the trace messages and convergence history.
C
      IF (VF.NE.0) WRITE (VF,'(2I8,2E11.4)') 0, 0, DONE, DONE
      IF (TF.NE.0) WRITE (TF,'(2I8,2E11.4)') 0, 0, DONE, DONE
C
C     Set x_0 = 0 and compute the norm of the initial residual.
C
      CALL ZAXPBY (NLEN,VECS(1,5),ZONE,VECS(1,2),ZZERO,VECS(1,5))
      CALL ZAXPBY (NLEN,VECS(1,1),ZZERO,VECS(1,1),ZZERO,VECS(1,1))
      R0 = DZNRM2(NLEN,VECS(1,5),1)
      IF ((TOL.GE.DONE).OR.(R0.EQ.DZERO)) GO TO 90
C
C     Check whether the auxiliary vector must be supplied.
C
      IF (INIT.EQ.0) CALL ZRANDN (NLEN,VECS(1,3),1)
C
C     Initialize the variables.
C
      N    = 1
      RESN = DONE
      RHO  = ZONE
      VAR  = DZERO
      ETA  = ZZERO
      TAU  = R0 * R0
      IERR = 8
      CALL ZAXPBY (NLEN,VECS(1,8),ZZERO,VECS(1,8),ZZERO,VECS(1,8))
      CALL ZAXPBY (NLEN,VECS(1,4),ZZERO,VECS(1,4),ZZERO,VECS(1,4))
      CALL ZAXPBY (NLEN,VECS(1,6),ZZERO,VECS(1,6),ZZERO,VECS(1,6))
C
C     This is one step of the TFQMR algorithm.
C     Compute \beta_{n-1} and \rho_{n-1}.
C
 20   ZTMP = ZDOTU(NLEN,VECS(1,3),1,VECS(1,5),1)
      BETA = ZTMP / RHO
      RHO  = ZTMP
C
C     Compute y_{2n-1}, v_{n-1}, and A y_{2n-1}.
C
      CALL ZAXPBY (NLEN,VECS(1,4),BETA,VECS(1,4),ZONE,VECS(1,8))
      CALL ZAXPBY (NLEN,VECS(1,6),ZONE,VECS(1,5),BETA,VECS(1,6))
C
C     Have the caller carry out AXB, then return here.
C        CALL AXB (VECS(1,6),VECS(1,9))
C
      INFO(2) = 1
      INFO(3) = 6
      INFO(4) = 9
      RETLBL  = 30
      RETURN
 30   CALL ZAXPBY (NLEN,VECS(1,4),BETA,VECS(1,4),ZONE,VECS(1,9))
C
C     Compute \sigma{n-1} and check for breakdowns.
C
      ZTMP = ZDOTU(NLEN,VECS(1,3),1,VECS(1,4),1)
      IF ((CDABS(ZTMP).EQ.DZERO).OR.(CDABS(RHO).EQ.DZERO)) THEN
         IERR = 8
         GO TO 90
      END IF
C
C     Compute \alpha_{n-1}, d_{2n-1} and w_{2n}.
C
      ALPHA = RHO / ZTMP
      ZTMP  = VAR * ETA / ALPHA
      CALL ZAXPBY (NLEN,VECS(1,7),ZONE,VECS(1,6),ZTMP,VECS(1,7))
      CALL ZAXPBY (NLEN,VECS(1,5),ZONE,VECS(1,5),-ALPHA,VECS(1,9))
C
C     Compute \varepsilon_{2n-1}^2, \eta_{2n-1}^2, c_{2n-1}^2, and
C     \tau_{2n-1}^2.
C
      DTMP = DZNRM2(NLEN,VECS(1,5),1)
      DTMP = DTMP * DTMP
      VAR  = DTMP / TAU
      COS  = DONE / ( DONE + VAR )
      TAU  = DTMP * COS
      ETA  = ALPHA * COS
C
C     Compute x_{2n-1} and the upper bound for its residual norm.
C
      CALL ZAXPBY (NLEN,VECS(1,1),ZONE,VECS(1,1),ETA,VECS(1,7))
C
C     Compute the residual norm upper bound.
C     If the scaled upper bound is within one order of magnitude of the
C     target convergence norm, compute the true residual norm.
C
      UNRM = DSQRT(DBLE(2*N) * TAU) / R0
      UCHK = UNRM
      IF ((TRES.EQ.0).AND.(UNRM/TOL.GT.DTEN)) GO TO 50
C
C     Have the caller carry out AXB, then return here.
C        CALL AXB (VECS(1,1),VECS(1,9))
C
      INFO(2) = 1
      INFO(3) = 1
      INFO(4) = 9
      RETLBL  = 40
      RETURN
 40   CALL ZAXPBY (NLEN,VECS(1,9),ZONE,VECS(1,2),-ZONE,VECS(1,9))
      RESN = DZNRM2(NLEN,VECS(1,9),1) / R0
      UCHK = RESN
C
C     Output the trace messages and convergence history.
C
 50   IF (VF.NE.0) WRITE (VF,'(2I8,2E11.4)') N, 2*N-1, UNRM, RESN
      IF (TF.NE.0) WRITE (TF,'(2I8,2E11.4)') N, 2*N-1, UNRM, RESN
C
C     Check for convergence or termination.  Stop if:
C         1. algorithm converged;
C         2. the residual norm upper bound is smaller than the computed
C     residual norm by a factor of at least 100.
C
      IF (RESN.LE.TOL) THEN
         IERR = 0
         GO TO 90
      ELSE IF (UNRM.LT.UCHK/DHUN) THEN
         IERR = 4
         GO TO 90
      END IF
C
C     Compute y_{2n}, A y_{2n}, d_{2n}, and w_{2n+1}.
C
      CALL ZAXPBY (NLEN,VECS(1,6),ZONE,VECS(1,6),-ALPHA,VECS(1,4))
      ZTMP = VAR * COS
      CALL ZAXPBY (NLEN,VECS(1,7),ZONE,VECS(1,6),ZTMP,VECS(1,7))
C
C     Have the caller carry out AXB, then return here.
C        CALL AXB (VECS(1,6),VECS(1,8))
C
      INFO(2) = 1
      INFO(3) = 6
      INFO(4) = 8
      RETLBL  = 60
      RETURN
 60   CALL ZAXPBY (NLEN,VECS(1,5),ZONE,VECS(1,5),-ALPHA,VECS(1,8))
C
C     Compute \varepsilon_{2n}^2, \eta_{2n}^2, c_{2n}^2, and
C     \tau_{2n}^2.
C
      DTMP = DZNRM2(NLEN,VECS(1,5),1)
      DTMP = DTMP * DTMP
      VAR  = DTMP / TAU
      COS  = DONE / ( DONE + VAR )
      TAU  = DTMP * COS
      ETA  = ALPHA * COS
C
C     Compute x_{2n}.
C
      CALL ZAXPBY (NLEN,VECS(1,1),ZONE,VECS(1,1),ETA,VECS(1,7))
C
C     Compute the residual norm upper bound.
C     If the scaled upper bound is within one order of magnitude of the
C     target convergence norm, compute the true residual norm.
C
      UNRM = DSQRT(DBLE(2*N+1) * TAU) / R0
      UCHK = UNRM
      IF ((TRES.EQ.0).AND.(UNRM/TOL.GT.DTEN).AND.(N.LT.NLIM)) GO TO 80
C
C     Have the caller carry out AXB, then return here.
C        CALL AXB (VECS(1,1),VECS(1,9))
C
      INFO(2) = 1
      INFO(3) = 1
      INFO(4) = 9
      RETLBL  = 70
      RETURN
 70   CALL ZAXPBY (NLEN,VECS(1,9),ZONE,VECS(1,2),-ZONE,VECS(1,9))
      RESN = DZNRM2(NLEN,VECS(1,9),1) / R0
      UCHK = UNRM
C
C     Output the trace messages and convergence history.
C
 80   IF (VF.NE.0) WRITE (VF,'(2I8,2E11.4)') N, 2*N, UNRM, RESN
      IF (TF.NE.0) WRITE (TF,'(2I8,2E11.4)') N, 2*N, UNRM, RESN
C
C     Check for convergence or termination.  Stop if:
C         1. algorithm converged;
C         2. the residual norm upper bound is smaller than the computed
C     residual norm by a factor of at least 100;
C         3. algorithm exceeded the iterations limit.
C
      IF (RESN.LE.TOL) THEN
         IERR = 0
         GO TO 90
      ELSE IF (UNRM.LT.UCHK/DHUN) THEN
         IERR = 4
         GO TO 90
      ELSE IF (N.GE.NLIM) THEN
         IERR = 4
         GO TO 90
      END IF
C
C     Update the running counter.
C
      N = N + 1
      GO TO 20
C
C     Done.
C
 90   NLIM    = N
      RETLBL  = 0
      INFO(1) = IERR
C
      RETURN
      END
C
C**********************************************************************
**********************************************************************
C

**********************************************************************
C
C     Copyright (C) 1992  Roland W. Freund and Noel M. Nachtigal
C     All rights reserved.
C
C     This code is part of a copyrighted package.  For details, see the
C     file `cpyrit.doc' in the top-level directory.
C
C     *****************************************************************
C     ANY USE OF  THIS CODE CONSTITUTES ACCEPTANCE OF  THE TERMS OF THE
C                             COPYRIGHT NOTICE
C     *****************************************************************
C
C**********************************************************************
C
      SUBROUTINE ZAXPBY (N,ZZ,ZA,ZX,ZB,ZY)
C
C     Purpose:
C     This subroutine computes ZZ = ZA * ZX + ZB * ZY.  Several special
C     cases are handled separately:
C        ZA =  0.0, ZB =  0.0 => ZZ = 0.0
C        ZA =  0.0, ZB =  1.0 => ZZ = ZY  (this is COPY)
C        ZA =  0.0, ZB = -1.0 => ZZ = -ZY
C        ZA =  0.0, ZB =   ZB => ZZ = ZB * ZY  (this is SCAL)
C        ZA =  1.0, ZB =  0.0 => ZZ = ZX  (this is COPY)
C        ZA =  1.0, ZB =  1.0 => ZZ = ZX + ZY
C        ZA =  1.0, ZB = -1.0 => ZZ = ZX - ZY
C        ZA =  1.0, ZB =   ZB => ZZ = ZX + ZB * ZY (this is AXPY)
C        ZA = -1.0, ZB =  0.0 => ZZ = -ZX
C        ZA = -1.0, ZB =  1.0 => ZZ = -ZX + ZY
C        ZA = -1.0, ZB = -1.0 => ZZ = -ZX - ZY
C        ZA = -1.0, ZB =   ZB => ZZ = -ZX + ZB * ZY
C        ZA =   ZA, ZB =  0.0 => ZZ = ZA * ZX  (this is SCAL)
C        ZA =   ZA, ZB =  1.0 => ZZ = ZA * ZX + ZY  (this is AXPY)
C        ZA =   ZA, ZB = -1.0 => ZZ = ZA * ZX - ZY
C        ZA =   ZA, ZB =   ZB => ZZ = ZA * ZX + ZB * ZY
C     ZZ may be the same as ZX or ZY.
C
C     Parameters:
C     N  = the dimension of the vectors (input).
C     ZZ = the vector result (output).
C     ZA = scalar multiplier for ZX (input).
C     ZX = one of the vectors (input).
C     ZB = scalar multiplier for ZY (input).
C     ZY = the other vector (input).
C
C     Noel M. Nachtigal
C     March 23, 1993
C
C**********************************************************************
C
      INTRINSIC DIMAG, DREAL
C
      INTEGER N
      DOUBLE COMPLEX ZA, ZB, ZX(N), ZY(N), ZZ(N)
C
C     Local variables.
C
      INTEGER I
      DOUBLE PRECISION DAI, DAR, DBI, DBR
C
      IF (N.LE.0) RETURN
C
      DAI = DIMAG(ZA)
      DAR = DREAL(ZA)
      DBI = DIMAG(ZB)
      DBR = DREAL(ZB)
      IF ((DAR.EQ.0.0D0).AND.(DAI.EQ.0.0D0)) THEN
         IF ((DBR.EQ.0.0D0).AND.(DBI.EQ.0.0D0)) THEN
C           ZA = 0.0, ZB = 0.0 => ZZ = 0.0.
            DO 10 I = 1, N
               ZZ(I) = (0.0D0,0.0D0)
 10         CONTINUE
         ELSE IF ((DBR.EQ.1.0D0).AND.(DBI.EQ.0.0D0)) THEN
C           ZA = 0.0, ZB = 1.0 => ZZ = ZY (this is COPY).
            DO 20 I = 1, N
               ZZ(I) = ZY(I)
 20         CONTINUE
         ELSE IF ((DBR.EQ.-1.0D0).AND.(DBI.EQ.0.0D0)) THEN
C           ZA = 0.0, ZB = -1.0 => ZZ = -ZY.
            DO 30 I = 1, N
               ZZ(I) = -ZY(I)
 30         CONTINUE
         ELSE
C           ZA = 0.0, ZB = ZB => ZZ = ZB * ZY (this is SCAL).
            DO 40 I = 1, N
               ZZ(I) = ZB * ZY(I)
 40         CONTINUE
         END IF
      ELSE IF ((DAR.EQ.1.0D0).AND.(DAI.EQ.0.0D0)) THEN
         IF ((DBR.EQ.0.0D0).AND.(DBI.EQ.0.0D0)) THEN
C           ZA = 1.0, ZB = 0.0 => ZZ = ZX (this is COPY).
            DO 50 I = 1, N
               ZZ(I) = ZX(I)
 50         CONTINUE
         ELSE IF ((DBR.EQ.1.0D0).AND.(DBI.EQ.0.0D0)) THEN
C           ZA = 1.0, ZB = 1.0 => ZZ = ZX + ZY.
            DO 60 I = 1, N
               ZZ(I) = ZX(I) + ZY(I)
 60         CONTINUE
         ELSE IF ((DBR.EQ.-1.0D0).AND.(DBI.EQ.0.0D0)) THEN
C           ZA = 1.0, ZB = -1.0 => ZZ = ZX - ZY.
            DO 70 I = 1, N
               ZZ(I) = ZX(I) - ZY(I)
 70         CONTINUE
         ELSE
C           ZA = 1.0, ZB = ZB => ZZ = ZX + ZB * ZY (this is AXPY).
            DO 80 I = 1, N
               ZZ(I) = ZX(I) + ZB * ZY(I)
 80         CONTINUE
         END IF
      ELSE IF ((DAR.EQ.-1.0D0).AND.(DAI.EQ.0.0D0)) THEN
         IF ((DBR.EQ.0.0D0).AND.(DBI.EQ.0.0D0)) THEN
C           ZA = -1.0, ZB = 0.0 => ZZ = -ZX
            DO 90 I = 1, N
               ZZ(I) = -ZX(I)
 90         CONTINUE
         ELSE IF ((DBR.EQ.1.0D0).AND.(DBI.EQ.0.0D0)) THEN
C           ZA = -1.0, ZB = 1.0 => ZZ = -ZX + ZY
            DO 100 I = 1, N
               ZZ(I) = -ZX(I) + ZY(I)
 100        CONTINUE
         ELSE IF ((DBR.EQ.-1.0D0).AND.(DBI.EQ.0.0D0)) THEN
C           ZA = -1.0, ZB = -1.0 => ZZ = -ZX - ZY.
            DO 110 I = 1, N
               ZZ(I) = -ZX(I) - ZY(I)
 110        CONTINUE
         ELSE
C           ZA = -1.0, ZB = ZB => ZZ = -ZX + ZB * ZY
            DO 120 I = 1, N
               ZZ(I) = -ZX(I) + ZB * ZY(I)
 120        CONTINUE
         END IF
      ELSE
         IF ((DBR.EQ.0.0D0).AND.(DBI.EQ.0.0D0)) THEN
C           ZA = ZA, ZB = 0.0 => ZZ = ZA * ZX (this is SCAL).
            DO 130 I = 1, N
               ZZ(I) = ZA * ZX(I)
 130        CONTINUE
         ELSE IF ((DBR.EQ.1.0D0).AND.(DBI.EQ.0.0D0)) THEN
C           ZA = ZA, ZB = 1.0 => ZZ = ZA * ZX + ZY (this is AXPY)
            DO 140 I = 1, N
               ZZ(I) = ZA * ZX(I) + ZY(I)
 140        CONTINUE
         ELSE IF ((DBR.EQ.-1.0D0).AND.(DBI.EQ.0.0D0)) THEN
C           ZA = ZA, ZB = -1.0 => ZZ = ZA * ZX - ZY.
            DO 150 I = 1, N
               ZZ(I) = ZA * ZX(I) - ZY(I)
 150        CONTINUE
         ELSE
C           ZA = ZA, ZB = ZB => ZZ = ZA * ZX + ZB * ZY.
            DO 160 I = 1, N
               ZZ(I) = ZA * ZX(I) + ZB * ZY(I)
 160        CONTINUE
         END IF
      END IF
C
      RETURN
      END
C
C**********************************************************************




      subroutine dopiszl(i,j,nlin,celem,celem2,cnz1,cnz2,nmax)
      double complex cnz1(nmax,860),celem
      double complex cnz2(nmax,860),celem2
      dimension nlin(nmax,0:860)
      if(abs(celem).lt.1e-10.and.abs(celem2).lt.1e-10) return
c     !!WRITE(*,*) i,j,' dl'
      ibyl=0
      ni=nlin(i,0)
c     !!WRITE(*,*) ni,' aa'
      do 1 k=1,ni
      nk=nlin(i,k)
      if(nk.eq.j)then
      cnz1(i,k)=cnz1(i,k)+celem
      cnz2(i,k)=cnz2(i,k)+celem2
      return
      endif
1     continue
      nlin(i,0)=nlin(i,0)+1
      nlin(i,ni+1)=j
      cnz1(i,ni+1)=celem
      cnz2(i,ni+1)=celem2
      end



      subroutine wk(N,M,KNZM,nmb,nsize)
      real *8 xnmb,ddgamma
c w tablicy KNZM zapisuje M po N kombinacji N elementow
c z M-elementowego zbioru
      dimension KNZM(nsize,6)
      if(m.eq.0) then
      nmb=1
      return
      endif
      XNMB=1.d0
c   !!WRITE(*,*) 'komb ',n,m

      if(M-N.gt.N) then
      do 321 i=M-N+1,M
      XNMB=XNMB*i
321   continue
      XNMB=XNMB/ddgamma(N)
      else
      do 322 i=N+1,M
      XNMB=XNMB*i
322   continue
      XNMB=XNMB/ddgamma(M-N)
      endif
      NMB=XNMB


c     !!WRITE(*,*) nmb
c     stop
      if(N.gt.M) then
      !!WRITE(*,*) 'N>M w KNZM'
      stop
      endif
      if(NMB.gt.nsize.
     >or.N.gt.6) then
      !!WRITE(*,*) 'za mala tablica KNZM'
      !!WRITE(*,*) nmb,n
      stop
      endif
c kombinacja startowa
      if(N.eq.1) then
      do 1 i=1,M
1     KNZM(i,1)=i
      endif
      if(N.eq.2) then
      in=0
      do 2 i1=1,M
      do 2 i2=i1+1,M
      in=in+1
      KNZM(in,1)=i1
2     KNZM(in,2)=i2
      endif
      if(N.eq.3) then
      in=0
      do 3 i1=1,M
      do 3 i2=i1+1,M
      do 3 i3=i2+1,M
      in=in+1
      KNZM(in,1)=i1
      KNZM(in,2)=i2
      KNZM(in,3)=i3
3     continue
      endif
      if(N.eq.4) then
      in=0
      do 4 i1=1,M
      do 4 i2=i1+1,M
      do 4 i3=i2+1,M
      do 4 i4=i3+1,M
      in=in+1
      KNZM(in,1)=i1
      KNZM(in,2)=i2
      KNZM(in,3)=i3
      KNZM(in,4)=i4
4     continue
      endif
      if(N.eq.5) then
      in=0
      do 5 i1=1,M
      do 5 i2=i1+1,M
      do 5 i3=i2+1,M
      do 5 i4=i3+1,M
      do 5 i5=i4+1,M

      in=in+1
      KNZM(in,1)=i1
      KNZM(in,2)=i2
      KNZM(in,3)=i3
      KNZM(in,4)=i4
      KNZM(in,5)=i5
5     continue
      endif
      if(N.eq.6) then
      in=0
      do 6 i1=1,M
      do 6 i2=i1+1,M
      do 6 i3=i2+1,M
      do 6 i4=i3+1,M
      do 6 i5=i4+1,M
      do 6 i6=i5+1,M

      in=in+1
      KNZM(in,1)=i1
      KNZM(in,2)=i2
      KNZM(in,3)=i3
      KNZM(in,4)=i4
      KNZM(in,5)=i5
      KNZM(in,6)=i6
6     continue
      endif
      end





      subroutine perm(n,iperm)
c generacja permutacji
      implicit double precision (a-h,o-z)
      dimension iperm(720,0:6),ipe(0:13),ip2(0:13)
      iperm(1,0)=0

c diagonala

      do 10 i=1,n
      iperm(1,i)=i
10    continue
      if(n.eq.1) return
c bierzaca permutowana ip
      ip=1
c znak permutacji
      is=1
c liczba znalezionych permutacji
      np=1

1     continue

c kopiowanie bierzacego

      do 2 i=1,n
      ipe(i)=iperm(ip,i)
2     continue

c petla po probnych zmianach
      do 3 i=2,n
      do 333 ki=0,12
      ip2(ki)=ipe(ki)
333   continue
      ip2(0)=is
      ip2(1)=ipe(i)
      ip2(i)=ipe(1)

c sprawdzenie czy taka permutacja juz byla

      do 4 j=1,np
      ib=0
      do 5 k=1,n
      ib=abs(ip2(k)-iperm(j,k))+ib
5     continue
      if(ib.eq.0) goto 3
4     continue

c dopisanie nowej
      np=np+1
      do 6 j=0,n
      iperm(np,j)=ip2(j)
6     continue
c   write(*,122) (-1)**(iperm(np,0)),(iperm(np,j),j=1,n)

c   write(*,122) (iperm(np,j),j=1,n)

      icontr=ddgamma(n)
c   !!WRITE(*,*) np,icontr
      if(np.eq.icontr) goto 1972
3     continue
      ip=ip+1
      is=mod(iperm(ip,0)+1,2)
      goto1
1972  continue
122   format(20i5)
      end


      function ddgamma(i)
      double precision ddgamma
      ddgamma=1.0
      do 1 j=1,i
      ddgamma=ddgamma*j
1     continue
      end


      subroutine w2e(n,cnz,nz,nsi,numer,jnz,enes,nst,cust)
c n - rozmiar macierzy
c cnz , nz - niezerowe elementy i ich idexy
c nsi- rozmiar fizyczny cnz, nz
c jnz -liczba niezerowych elelmentow macierzowych w cnz
      implicit double precision(a,b,d-h,o-z)
      implicit double complex(c)
      parameter(maxn=2000000,maxnev=50,maxncv=91,ldv=maxn)
      dimension cnz(nsi)
      dimension enes(nsi),nz(nsi,2)
      dimension en1(nsi)
      dimension eval(nsi)

      dimension cust(numer,nst)
      integer           maxn, maxnev, maxncv, ldv
      integer           iparam(11), ipntr(14)
      logical           select(maxncv)
      complex*16 v
      allocatable v(:,:)
      Complex*16
     &                  ax(maxn), d(maxncv),
     &                  workd(3*maxn),
     &                  workev(3*maxncv), resid(maxn),
     &                  workl(3*maxncv*maxncv+5*maxncv)


      Double precision
     &                  rwork(maxncv), rd(maxncv,3)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character         bmat*1, which*2
      integer           ido, n, nx, nev, ncv, lworkl, info, j,
     &                  ierr, nconv, maxitr, ishfts, mode
      Complex*16
     &                  sigma
      Double precision
     &                  tol
      logical           rvec
c
c     %-----------------------------%
c     | BLAS & LAPACK routines used |
c     %-----------------------------%
c
      Double precision
     &                  dznrm2, dlapy2
      external          dznrm2, zaxpy, dlapy2
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c     %--------------------------------------------------%
c     | The number NX is the number of interior points   |
c     | in the discretization of the 2-dimensional       |
c     | convection-diffusion operator on the unit        |
c     | square with zero Dirichlet boundary condition.   |
c     | The number N(=NX*NX) is the dimension of the     |
c     | matrix.  A standard eigenvalue problem is        |
c     | solved (BMAT = 'I').  NEV is the number of       |
c     | eigenvalues to be approximated.  The user can    |
c     | modify NX, NEV, NCV, WHICH to solve problems of  |
c     | different sizes, and to get different parts of   |
c     | the spectrum.  However, The following            |
c     | conditions must be satisfied:                    |
c     |                   N <= MAXN                      |
c     |                 NEV <= MAXNEV                    |
c     |           NEV + 2 <= NCV <= MAXNCV               |
c     %--------------------------------------------------%
c
      allocate(v(ldv,maxncv))

c     !!WRITE(*,*) n,nsi,jnz,nst,allocated(v)
      do 1 i=1,jnz
c     !!WRITE(*,*) i,nz(i,1),nz(i,2),cdabs(cnz(i))
1     continue
c     stop
      nev   = nst
      ncv   = nev*2+1
      if ( n .gt. maxn ) then
         print *, ' ERROR with _NDRV1: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _NDRV1: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _NDRV1: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat  = 'I'
      which = 'SM'
c
c     %---------------------------------------------------%
c     | The work array WORKL is used in ZNAUPD as         |
c     | workspace.  Its dimension LWORKL is set as        |
c     | illustrated below.  The parameter TOL determines  |
c     | the stopping criterion. If TOL<=0, machine        |
c     | precision is used.  The variable IDO is used for  |
c     | reverse communication, and is initially set to 0. |
c     | Setting INFO=0 indicates that a random vector is  |
c     | generated to start the ARNOLDI iteration.         |
c     %---------------------------------------------------%
c
      lworkl  = 3*ncv**2+5*ncv
      tol    = 0.000000001
      ido    = 0
      info   = 0
c
c     %---------------------------------------------------%
c     | This program uses exact shift with respect to     |
c     | the current Hessenberg matrix (IPARAM(1) = 1).    |
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 1 of ZNAUPD is used     |
c     | (IPARAM(7) = 1). All these options can be changed |
c     | by the user. For details see the documentation in |
c     | ZNAUPD.                                           |
c     %---------------------------------------------------%
c
      ishfts = 1
      maxitr = 500000
      mode   = 1
c
      iparam(1) = ishfts
      iparam(3) = maxitr
      iparam(7) = mode
c
c     %-------------------------------------------%
c     | M A I N   L O O P (Reverse communication) |
c     %-------------------------------------------%
c
 10   continue
c
c        %---------------------------------------------%
c        | Repeatedly call the routine ZNAUPD and take |
c        | actions indicated by parameter IDO until    |
c        | either convergence is indicated or maxitr   |
c        | has been exceeded.                          |
c        %---------------------------------------------%
c
         call znaupd ( ido, bmat, n, which, nev, tol, resid, ncv,
     &        v, ldv, iparam, ipntr, workd, workl, lworkl,
     &        rwork,info )
c
         if (ido .eq. -1 .or. ido .eq. 1) then
c
c           %-------------------------------------------%
c           | Perform matrix vector multiplication      |
c           |                y <--- OP*x                |
c           | The user should supply his/her own        |
c           | matrix vector multiplication routine here |
c           | that takes workd(ipntr(1)) as the input   |
c           | vector, and return the matrix vector      |
c           | product to workd(ipntr(2)).               |
c           %-------------------------------------------%
c
           call av2e(n,workd(ipntr(1)),workd(ipntr(2)),cnz,nz,nsi,jnz)


c
c           %-----------------------------------------%
c           | L O O P   B A C K to call ZNAUPD again. |
c           %-----------------------------------------%
c
            go to 10

         end if
c
c     %----------------------------------------%
c     | Either we have convergence or there is |
c     | an error.                              |
c     %----------------------------------------%
c
      if ( info .lt. 0 ) then
c
c        %--------------------------%
c        | Error message, check the |
c        | documentation in ZNAUPD  |
c        %--------------------------%
c
         print *, ' '
         print *, ' Error with _naupd, info = ', info
         print *, ' Check the documentation of _naupd'
         print *, ' '
c
      else
c
c        %-------------------------------------------%
c        | No fatal errors occurred.                 |

ct-Process using ZNEUPD.                |
c        |                                           |
c        | Computed eigenvalues may be extracted.    |
c        |                                           |
c        | Eigenvectors may also be computed now if  |
c        | desired.  (indicated by rvec = .true.)    |
c        %-------------------------------------------%
c
         rvec = .true.
c
         call zneupd (rvec, 'A', select, d, v, ldv, sigma,
     &        workev, bmat, n, which, nev, tol, resid, ncv,
     &        v, ldv, iparam, ipntr, workd, workl, lworkl,
     &        rwork, ierr)

      do 77 ist=1,nst
      wnu=0
      do 78 ix=1,n
      cust(ix,ist)=v(ix,ist)
      wnu=wnu+cdabs(v(ix,ist))**2
78    continue
      do 79 ix=1,n
      cust(ix,ist)=cust(ix,ist)/dsqrt(wnu)
79    continue
      enes(ist)=d(ist)
77    continue

c sortujemy
c     goto 830
      do 334 i=1,nst
      emin=enes(i)
      imin=i
      do 335 j=i+1,nst
      if(enes(j).lt.emin) then
      emin=enes(j)
      imin=j
      endif
335   continue
      e0=enes(i)
      enes(i)=emin
      enes(imin)=e0
      do 400 ix=1,n
      c0=cust(ix,i)
      cust(ix,i)=cust(ix,imin)
      cust(ix,imin)=c0
400   continue
334   continue

      do 829 i=1,nst
       !!WRITE(*,*) i,enes(i)*27211.6
829   continue
830   continue
c      stop

c
c        %----------------------------------------------%
c        | Eigenvalues are returned in the one          |
c        | dimensional array D.  The corresponding      |
c        | eigenvectors are returned in the first NCONV |
c        | (=IPARAM(5)) columns of the two dimensional  |
c        | array V if requested.  Otherwise, an         |
c        | orthogonal basis for the invariant subspace  |
c        | corresponding to the eigenvalues in D is     |
c        | returned in V.                               |
c        %----------------------------------------------%
c
         if ( ierr .ne. 0) then
c
c           %------------------------------------%
c           | Error condition:                   |
c           | Check the documentation of ZNEUPD. |
c           %------------------------------------%
c
             print *, ' '
             print *, ' Error with _neupd, info = ', ierr
             print *, ' Check the documentation of _neupd. '
             print *, ' '
c
         else
c
             nconv = iparam(5)
             do 20 j=1, nconv
c
c               %---------------------------%
c               | Compute the residual norm |
c               |                           |
c               |   ||  A*x - lambda*x ||   |
c               |                           |
c               | for the NCONV accurately  |
c               | computed eigenvalues and  |
c               | eigenvectors.  (iparam(5) |
c               | indicates how many are    |
c               | accurate to the requested |
c               | tolerance)                |
c               %---------------------------%


c
                call av2e(n, v(1,j), ax,cnz,nz,nsi,jnz)
                call zaxpy(n, -d(j), v(1,j), 1, ax, 1)
                rd(j,1) = dble(d(j))
                rd(j,1) = dble(d(j))
                rd(j,2) = dimag(d(j))
                rd(j,3) = dznrm2(n, ax, 1)
                rd(j,3) = rd(j,3) / dlapy2(rd(j,1),rd(j,2))
 20          continue
c
c            %-----------------------------%
c            | Display computed residuals. |
c            %-----------------------------%
c
c            call dmout(6, nconv, 3, rd, maxncv, -6,
c    &            'Ritz values (Real, Imag) and relative residuals')
          end if
c
c        %-------------------------------------------%
c        | Print additional convergence information. |
c        %-------------------------------------------%
c
c        if ( info .eq. 1) then
c            print *, ' '
c            print *, ' Maximum number of iterations reached.'
c            print *, ' '
c        else if ( info .eq. 3) then
c            print *, ' '
c            print *, ' No shifts could be applied during implicit
c    &                  Arnoldi update, try increasing NCV.'
c            print *, ' '
         end if
c
c        print *, ' '
crc
cc
c        print *, ' '
c        print *, ' Size of the matrix is ', n
c        print *, ' The number of Ritz values requested is ', nev
c        print *, ' The number of Arnoldi vectors generated',
c    &            ' (NCV) is ', ncv
c        print *, ' What portion of the spectrum: ', which
c        print *, ' The number of converged Ritz values is ',
c    &              nconv
c        print *, ' The number of Implicit Arnoldi update',
c    &            ' iterations taken is ', iparam(3)
c        print *, ' The number of OP*x is ', iparam(9)
c        print *, ' The convergence criterion is ', tol
c        print *, ' '
c
c     end if
c
c     %---------------------------%
c     | Done with program zndrv1. |
c     %---------------------------%
c
 9000 continue

200   continue
      deallocate(v)
      end


c==========================================================================
c
c     matrix vector subroutine
c
c     The matrix used is the convection-diffusion operator
c     discretized using centered difference.

      subroutine av2e (n, v, w,cnz,nz,nsi,jnz)
      integer           n, j, lo
      Complex*16
     &                  v(n), w(n), one, h2,hkl
      parameter         (one = (1.0D+0, 0.0D+0))
      external          zaxpy, tv
      complex*16 cnz(nsi)
      integer nz(nsi,2)
c
c     Computes w <--- OP*v, where OP is the nx*nx by nx*nx block
c     tridiagonal matrix
c
c                  | T -I          |
c                  |-I  T -I       |
c             OP = |   -I  T       |
c                  |        ...  -I|
c                  |           -I T|
c
c     derived from the standard central difference  discretization
c     of the convection-diffusion operator (Laplacian u) + rho*(du/dx)
c     with zero boundary condition.
c
c      double precision  vpot(ns,ms),bat,dx,c,xk
c      common/cnzz/cnz,
      w=0

      do 77 l=1,jnz
      i=nz(l,1)
      j=nz(l,2)
      w(i)=w(i)+cnz(l)*v(j)
77    continue
       return
      end

