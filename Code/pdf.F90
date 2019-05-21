program iterative_smoothers
   use m_marginalpdf
   use m_set_random_seed2
   use m_moments
   use m_avevar
   use m_func
   use m_steinkern
   use m_cyyreg
   use m_aaprojection
   use m_costf
   use m_getalpha
   use m_tecpdf
   use m_omegafact
   use m_tecmargpdf
   use m_pseudoinv
   use m_pseudoinv2
   use m_tecjointpdf
   use m_tecfunc
   use m_teccostens
   use m_tecsampini
   use m_tecstat
   use m_tecstatIES
   use m_normal
   use m_cov
   use m_integrals
   implicit none
!   integer, parameter :: nrsamp=10000000
   integer nrsamp,esamp,nx,ny
   integer maxiesit   ! 20 
   real gamma_ies    ! step length factor used in IES (gamma_ies=0.5)
   real ave,var

   character(len=7) :: method(0:4)=['ES     ','IES    ','IEnKF  ', 'ESMDA  ', 'EnSTEIN']
   character(len=1) :: variable(1:2)=['x','q']

   real, allocatable, dimension(:) :: avex,adevx,sdevx,varx,skewx,kurtx
   real, allocatable, dimension(:) :: avey,adevy,sdevy,vary,skewy,kurty

   real, allocatable, dimension(:) :: avexIES,varxIES,skewxIES,kurtxIES
   real, allocatable, dimension(:) :: aveyIES,varyIES,skewyIES,kurtyIES

   real, allocatable :: kern(:,:)

   real x0      ! Prior estimate
   real d       ! Measurement
   real rh      ! model parameter
   real beta    ! model parameter
   real siga    ! std dev in prior
   real sigw    ! std dev in model
   real sigo    ! std dev in observation
   real xm      ! half grid domain
   real sft     ! shift in y-direction
   integer funcmode ! which function to use
   integer updatemode ! Update y by ES (0) or evaluate (1)
   integer nmda ! number of mda iterations
   logical :: lies=.true.  ! Run IES or not
   logical :: lmda=.true.  ! Run MDA or not
   logical :: lienks=.true.  ! Run IEnKS or not
   logical :: lenstein=.false.  ! Run EnStein or not
   integer gradient
   integer IESv

!   real x,mode1,mode0,mean1,mean0
   real dx,dy,sump,dq,fac,tmp
   real grad1,grad2(2)
   real cxx,cyy,cqq,cyx,cqy,cqx,alphasum
   real pxx,Pyy,pqq,pyx,pqy,pqx,Pyymat(1,1)
   real cdd,G
   real WoodA,WoodB
   real hessian, cvec(2) , psca
   real Czz(2,2),Pzz(2,2),PIzz(2,2),CIzz(2,2),zi(2),zf(2),Pzy(2,1),Pyz(1,2)
   real, allocatable :: alpha(:)
   real :: alphageo=0.0
   real n1


   real, allocatable :: sf(:),si(:,:),grad(:,:),newgrad(:,:)
   real :: xlength=1.0


   logical lmoderr, testpseudo

   real, allocatable :: x(:)       ! x-locations
   real, allocatable :: y(:)       ! y-locations
   real, allocatable :: q(:)       ! y-locations
   real, allocatable :: margx(:)   ! x-locations
   real, allocatable :: margy(:)   ! y-locations
   real, allocatable :: datum(:)   ! datum pdf
   real, allocatable :: prior(:)   ! prior pdf
   real, allocatable :: pdf(:,:)
   real, allocatable :: cost(:)

   real, allocatable :: xsamp(:),xsampini(:) 
   real, allocatable :: ysamp(:),ysampini(:)
   real, allocatable :: qsamp(:),qsampini(:)
   real, allocatable :: dpert(:)
   integer, allocatable :: iconv(:)
   real pert

   character(len=3) tag3
   character(len=2) tag2
   character(len=1) tag1
   character(len=40) caseid

   integer i,j,k,m,n,loc(1),ival,jval
   integer sumconv

   real jx,djx,djq,ddjx
   real dgx,dgq
   real aveA

   real :: xa,xb,ya,yb,qa,qb ,C(1,1)
   real :: sigq= 0.05

   integer :: nrits=100
   real, allocatable :: costens(:,:)
   real, allocatable :: costite(:,:)
   real, allocatable :: xsampit(:,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! IEnKS variables
   real, allocatable :: YAinv(:,:)
   real, allocatable :: W(:,:)   ! update coefficients
   real, allocatable :: WW(:,:)   ! update coefficients from previous iterate
   real, allocatable :: Yi(:,:)   ! predicted measurements 
   real, allocatable :: YY(:,:)  ! predicted measurement anomalies 
   real, allocatable :: S(:,:)   ! predicted measurement anomalies Y times W^+
   real, allocatable :: S1(:,:)   ! predicted measurement anomalies Y times W^+
   real, allocatable :: S2(:,:)   ! predicted measurement anomalies Y times W^+
   real, allocatable :: S3(:,:)   ! predicted measurement anomalies Y times W^+
   real, allocatable :: YT(:,:),ST(:,:),STO(:,:)
   real, allocatable :: Dens(:,:)   ! Ensemble of perturbed measurements
   real, allocatable :: H(:,:)     ! "Innovation "
   real, allocatable :: E0(:,:)   ! Initial model ensemble 
   real, allocatable :: Ei(:,:)  ! model ensemble iteration i
   real, allocatable :: A0(:,:)   ! Initial model ensemble perturbation
   real, allocatable :: Ai(:,:)   ! model ensemble perturbation at iteration i
   real, allocatable :: Ainv(:,:)   !  Model ensemble perturbation at iteration i
   real, allocatable :: A0inv(:,:)   !  Model ensemble perturbation at iteration i
!   real, allocatable :: AA0(:,:)   !  A^+A
   real, allocatable :: AAi(:,:)   !  Ai^+Ai
   real, allocatable :: aveW(:)   
   real, allocatable :: samples(:,:,:) ! for comparing ES, IES and IEnKF
   real, allocatable :: Omega(:,:)  
   real SS(1,1),meanY,diffW,Y0mean
   real :: truncation=1.0
   real gamma_ienks
   integer maxienksit
   integer innovation
   integer ienks_S
   logical :: lpseudoinv=.false.
   integer, allocatable :: ipiv(:)
   integer info
   logical :: rank1svd=.false.
   logical lcyyreg
   real start,finish,diffx
   real, allocatable :: xx(:),xxold(:),bb(:)
   integer :: ndim
   integer, parameter :: nrobs=1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call set_random_seed2()
   print *,'Normal number', normal()

   open(10,file='infile.in')
      read(10,*)esamp       ; print '(a,i4,i12)',   'number  of samples 10^x    :',esamp,10**esamp
      read(10,*)nx          ; print '(a,i4)',       'x-dimension of grid  (nx)  :',nx
      read(10,*)ny          ; print '(a,i4)',       'y-dimension of grid  (ny)  :',ny
      read(10,*)xa          ; print '(a,f10.3)',    'xa                         :',xa
      read(10,*)xb          ; print '(a,f10.3)',    'xb                         :',xb
      read(10,*)ya          ; print '(a,f10.3)',    'ya                         :',ya
      read(10,*)yb          ; print '(a,f10.3)',    'yb                         :',yb
      read(10,*)x0          ; print '(a,f10.3)',    'prior estimate of x        :',x0
      read(10,*)siga        ; print '(a,f10.3)',    'err std dev of prior       :',siga
      read(10,*)d           ; print '(a,f10.3)',    'measurement of y           :',d
      read(10,*)sigo        ; print '(a,f10.3)',    'measurement std dev        :',sigo
      read(10,*)sigw        ; print '(a,f10.3)',    'model std dev              :',sigw
      read(10,*)lmoderr     ; print '(a,l1)',       'Activate model error       :',lmoderr                  ! Not used
      read(10,*)rh          ; print '(a,f10.3)',    'model parameter rh         :',rh
      read(10,*)beta        ; print '(a,f10.3)',    'model parameter beta       :',beta
      read(10,*)funcmode    ; print '(a,tr7,i3)',   'function to use            :',funcmode
      read(10,*)updatemode  ; print '(a,tr7,i3)',   'update of y(0), func(1)    :',updatemode  
      read(10,*)lcyyreg     ; print '(a,tr10,l1)',  'Regression for Cyy         :',lcyyreg
      read(10,*)lmda        ; print '(a,tr10,l1)',  'Run MDA                    :',lmda
      read(10,*)nmda        ; print '(a,tr7,i3)',   'number of mda iterations   :',nmda
      read(10,*)alphageo    ; print '(a,f10.3)',    'geometrical alpha value    :',alphageo
      read(10,*)lies        ; print '(a,tr10,l1)',  'Run IES                    :',lies
      read(10,*)gradient    ; print '(a,i3)',       'Gradient 0-ana, 1-stochas  :',gradient
      read(10,*)maxiesit    ; print '(a,i5)',       'maxiesit                   :',maxiesit
      read(10,*)gamma_ies   ; print '(a,f10.3)',    'gamma_ies                  :',gamma_ies
      read(10,*)IESv        ; print '(a,i1)',       'IESv                       :',IESv
      read(10,*)lienks      ; print '(a,tr10,l1)',  'Run IEnKS                  :',lienks
      read(10,*)maxienksit  ; print '(a,i5)',       'maxIEnKSit                 :',maxienksit
      read(10,*)gamma_ienks ; print '(a,f10.3)',    'IEnKS_gamma                :',gamma_ienks
      read(10,*)truncation  ; print '(a,f10.3)',    'IEnKS_truncation           :',truncation
      read(10,*)ndim        ; print '(a,i5)',       'ndim                       :',ndim
      read(10,*)IEnKS_S     ; print '(a,i5)',       'IEnKS_S                    :',IEnKS_S
   close(10)

!   if (funcmode == 0 .and. (nx /= 450 .or. ny /= 450)) stop 'check mod_dimensions'
!   if (funcmode == 2 .and. (nx /= 450 .or. ny /= 250)) stop 'check mod_dimensions'
         

!Include model errors in inverse calculation (Tarantola approach)
   cdd=sigo**2
   nrsamp=10**esamp
   nrsamp=1*nrsamp
   allocate(samples(nrsamp,2,0:4))
   allocate(dpert(nrsamp))
   allocate(qsamp(nrsamp),qsampini(nrsamp))
   allocate(xsamp(nrsamp),xsampini(nrsamp)) 
   allocate(ysamp(nrsamp),ysampini(nrsamp))
   allocate(iconv(nrsamp)) 

   allocate(x(nx), y(ny), q(nx))
   allocate(margx(nx), margy(ny))
   allocate(datum(ny))
   allocate(prior(nx))
   allocate(pdf(nx,ny))
   allocate(cost(nx))

   allocate(alpha(nmda))
   allocate(avex(0:nmda),adevx(0:nmda),sdevx(0:nmda),varx(0:nmda),skewx(0:nmda),kurtx(0:nmda))
   allocate(avey(0:nmda),adevy(0:nmda),sdevy(0:nmda),vary(0:nmda),skewy(0:nmda),kurty(0:nmda))

   allocate(avexIES(0:maxiesit),varxIES(0:maxiesit),skewxIES(0:maxiesit),kurtxIES(0:maxiesit))
   allocate(aveyIES(0:maxiesit),varyIES(0:maxiesit),skewyIES(0:maxiesit),kurtyIES(0:maxiesit))

   allocate(xx(nrsamp),bb(nrsamp), xxold(nrsamp))
   allocate (ipiv(nrsamp))
   allocate (Yi(1,nrsamp))
   allocate (YY(1,nrsamp))
   allocate (YAinv(1,ndim))
   allocate (S(1,nrsamp))
   allocate (S1(1,nrsamp))
   allocate (S2(1,nrsamp))
   allocate (S3(1,nrsamp))
   allocate (YT(nrsamp,1),ST(nrsamp,1),STO(nrsamp,1))
   allocate (Dens(1,nrsamp))
   allocate (H(1,nrsamp))
   allocate (E0(ndim,nrsamp))
   allocate (Ei(ndim,nrsamp))
   allocate (A0(ndim,nrsamp))
   allocate (Ai(ndim,nrsamp))
   allocate (Ainv(nrsamp,ndim))
   allocate (A0inv(nrsamp,ndim))
   allocate (aveW(nrsamp))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Theroretical values for statistical moements (trustat.dat)
   call integrals(x0,d,siga,sigo,beta,funcmode,maxiesit)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Grid domain and dx
   qa=-5.0*max(sigq,sigw)
   qb= 5.0*max(sigq,sigw)

   dx=(xb-xa)/real(nx-1); print '(a,f12.5)','dx=',dx
   dy=(yb-ya)/real(ny-1); print '(a,f12.5)','dy=',dy
   dq=(qb-qa)/real(nx-1); print '(a,f12.5)','dq=',dq

   do i=1,nx
      x(i)=xa + real(i-1)*dx
   enddo

   do j=1,ny
      y(j)=ya + real(j-1)*dy
   enddo

   do i=1,nx
      q(i)=qa + real(i-1)*dq
   enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Definition of the strong constraint cost function
   do i=1,nx
      prior(i)=exp(-0.5*(x(i)-x0)**2/siga**2)
      cost(i)=(x(i)-x0)**2/siga**2 + (func(x(i),beta,funcmode)-d)**2/cdd    
   enddo
   sump=sum(prior(:))*dx
   prior=prior/sump
   call tecfunc('costf',cost,x,nx,'x','Cost')
   call tecfunc('prior',prior,x,nx,'x','Prior')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Datum pdf
   do j=1,ny
      datum(j)=exp(-0.5*(y(j)-d)**2/cdd)
   enddo
   sump=sum(datum(:))*dy
   datum=datum/sump
   call tecfunc('datum',datum,y,ny,'y','Datum')




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Definition of the analytical joint unconditional pdf
   do j=1,ny
   do i=1,nx
      pdf(i,j)=exp( -0.5*(x(i)-x0)**2/siga**2            &
                    -0.5*(y(j)-func(x(i),beta,funcmode))**2/max(sigw,sigq)**2 )
   enddo
   enddo
   sump=sum(pdf(:,:))*dx*dy
   pdf=pdf/sump
   call getcaseid(caseid,'PDFJ',alphageo,nmda,esamp,gradient,beta,sigw,0)
   call tecjointpdf(pdf,x,y,nx,ny,caseid)
   call marginalpdf(pdf,margx,margy,nx,ny,x,y,dx,dy)
   call tecfunc('margx',margx,x,nx,'x',caseid)
   call tecfunc('margy',margy,y,ny,'y',caseid)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Definition of the analytical joint conditional pdf
   do j=1,ny
   do i=1,nx
      pdf(i,j)=exp(-0.5*(x(i)-x0)**2/siga**2            &
                   -0.5*(y(j)-func(x(i),beta,funcmode))**2/max(sigw,sigq)**2  &
                   -0.5*(y(j)-d)**2/cdd)    
   enddo
   enddo
   sump=sum(pdf(:,:))*dx*dy
   pdf=pdf/sump


   call getcaseid(caseid,'PDFC',alphageo,nmda,esamp,gradient,beta,sigw,0)
   call tecjointpdf(pdf,x,y,nx,ny,caseid)
   call marginalpdf(pdf,margx,margy,nx,ny,x,y,dx,dy)

   margx(:)=0.0
   do i=1,nx
   do j=1,ny
      if (sigw==0.0) then
         margx(i)=margx(i)+exp(-0.5*(x(i)-x0 )**2/siga**2                   &
                               -0.5*(d-func(x(i),beta,funcmode)-q(j))**2/cdd**2)
      else 
         margx(i)=margx(i)+exp(-0.5*(x(i)-x0 )**2/siga**2                   &
                               -0.5*(q(j)-0.0)**2/max(sigw,0.0001)**2  & 
                               -0.5*(d-func(x(i),beta,funcmode)-q(j))**2/cdd**2)
      endif
   enddo
   enddo
   margx=margx/(sum(margx(:))*dx)
   call tecfunc('margx',margx,x,nx,'x',caseid)
   call tecfunc('margy',margy,y,ny,'y',caseid)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Ensemble generation, prediction, ensemble joint pdf
   do i=1,nrsamp
      xsampini(i)=siga*normal()
   enddo
   call moments(xsampini,nrsamp,avex(0),adevx(0),sdevx(0),varx(0),skewx(0),kurtx(0))
   do i=1,nrsamp
      xsampini(i)=xsampini(i)-avex(0)
   enddo

   do i=1,nrsamp
      qsampini(i)=sigw*normal()
   enddo
   call moments(qsampini,nrsamp,avex(0),adevx(0),sdevx(0),varx(0),skewx(0),kurtx(0))
   do i=1,nrsamp
      qsampini(i)=qsampini(i)-avex(0)
   enddo

   do i=1,nrsamp
      xsampini(i)=x0+xsampini(i)
      ysampini(i)=func(xsampini(i),beta,funcmode) + qsampini(i)
   enddo

   do i=1,nrsamp
      dpert(i)=d+sqrt(cdd)*normal()
   enddo
   print *,'Sampling done'

   call getcaseid(caseid,'INI',alphageo,nmda,esamp,gradient,beta,sigw,0)
   call tecpdf(x,y,nx,ny,xsampini,ysampini,nrsamp,xa,ya,dx,dy,caseid)
   call tecmargpdf('x',xsampini,nrsamp,caseid,xa,xb,nx)
   call tecmargpdf('y',ysampini,nrsamp,caseid,ya,yb,ny)
   call tecmargpdf('q',qsampini,nrsamp,caseid,qa,qb,nx)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! printing the first members cost functions 
   allocate(costens(nx,10))
   allocate(costite(10,nrits))
   allocate(xsampit(10,nrits))
   do n=1,10
      do i=1,nx
         costens(i,n)=(x(i)-xsampini(n))**2/siga**2 + (func(x(i),beta,funcmode)-dpert(n))**2/cdd    
      enddo
      xsampit(n,1)=xsampini(n)
      costite(n,1)= (func(xsampit(n,1),beta,funcmode)-dpert(n))**2/cdd    
   enddo
   call teccostens('costens',costens,x,nx,10,'x')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! E S           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   write(*,'(a)')'++++++++++++++++++++++++++++++++++++++++++++++'
   write(*,'(a)')'ES analysis...'
   do i=1,nrsamp
      xsamp(i)=xsampini(i)
      qsamp(i)=qsampini(i)
      ysamp(i)=func(xsamp(i),beta,funcmode)+qsamp(i)
   enddo

   call cov(Cxx,Cyy,Cqq,Cyx,Cqy,Cqx,xsamp,ysamp,qsamp,nrsamp)

!   write(*,'(a,f10.4)')'Pyy from samples :',Cyy
   if (lcyyreg) cyy=cyyreg(Cxx,Cqq,Cyx,Cqy,Cqx)

   do i=1,nrsamp
      xsamp(i)=xsamp(i) + (cyx/(cyy+cdd))*(dpert(i)-ysamp(i))
      qsamp(i)=qsamp(i) + (cqy/(cyy+cdd))*(dpert(i)-ysamp(i))

      if (updatemode==0) then
         ysamp(i)=ysamp(i)+ (cyy/(cyy+cdd))*(dpert(i)-ysamp(i))
      else
         ysamp(i)=func(xsamp(i),beta,funcmode) + qsamp(i)
      endif

   enddo
   write(*,'(a,10g11.3)')'Es',xsamp(1:10)

   if (sigw < sigq) then
      do i=1,nrsamp
         ysamp(i)=ysamp(i)+sigq*normal()
      enddo
   endif
   call getcaseid(caseid,'ES',alphageo,nmda,esamp,gradient,beta,sigw,0)
   call tecpdf(x,y,nx,ny,xsamp,ysamp,nrsamp,xa,ya,dx,dy,caseid)
   call tecmargpdf('x',xsamp,nrsamp,caseid,xa,xb,nx)
   call tecmargpdf('y',ysamp,nrsamp,caseid,ya,yb,ny)
   call tecmargpdf('q',qsamp,nrsamp,caseid,qa,qb,nx)
   write(*,'(a)')'ES analysis completed'

   samples(:,1,0)=xsamp(:)
   samples(:,2,0)=qsamp(:)

   do n=1,10
      xsampit(n,2)=xsamp(n)
      costite(n,2)= (xsampit(n,2)-xsampini(n))**2/siga**2 + (func(xsampit(n,2),beta,funcmode)-dpert(n))**2/cdd    
   enddo
   call teccostens('costens',costens,x,nx,10,'x')
   call tecsampini('sampiniES',costite,xsampit,2,10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! E S M D A     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (lmda) then
      write(*,'(a)')'++++++++++++++++++++++++++++++++++++++++++++++'
      write(*,'(a)')'MDA analysis...'
      do i=1,nrsamp
         xsamp(i)=xsampini(i)
         qsamp(i)=qsampini(i)
         ysamp(i)=func(xsamp(i),beta,funcmode)+qsamp(i)
      enddo

      alphasum=0.0
      do n=1,nmda
         alpha(n)=getalpha(n,nmda,alphageo)
         alphasum=alphasum+1.0/alpha(n)
         print *,'alpha            :',alpha(n)
         call cov(Cxx,Cyy,Cqq,Cyx,Cqy,Cqx,xsamp,ysamp,qsamp,nrsamp)
         write(*,'(a,f10.4)')'Cyy from samples :',Cyy
         if (lcyyreg) cyy=cyyreg(Cxx,Cqq,Cyx,Cqy,Cqx)
         write(*,'(i3,f10.2,a,2f13.5,e13.5)')n,alpha(n),', cxx= ',cxx,cyx,cqy

         do i=1,nrsamp
            pert=sqrt(alpha(n))*sqrt(cdd)*normal()
            xsamp(i)=xsamp(i) + (cyx/(cyy+alpha(n)*cdd))*(d+pert-ysamp(i))
            qsamp(i)=qsamp(i) + (cqy/(cyy+alpha(n)*cdd))*(d+pert-ysamp(i))
            ysamp(i)=func(xsamp(i),beta,funcmode)+qsamp(i)
         enddo
         call getcaseid(caseid,'MDA',alphageo,nmda,esamp,gradient,beta,sigw,n)
         call tecmargpdf('x',xsamp,nrsamp,caseid,xa,xb,nx)

      enddo


      samples(:,1,3)=xsamp(:)
      samples(:,2,3)=qsamp(:)

      if (sigw < sigq) then
         do i=1,nrsamp
            ysamp(i)=ysamp(i)+sigq*normal()
         enddo
      endif
      call getcaseid(caseid,'MDA',alphageo,nmda,esamp,gradient,beta,sigw,0)
      call tecpdf(x,y,nx,ny,xsamp,ysamp,nrsamp,xa,ya,dx,dy,caseid)
      call tecmargpdf('x',xsamp,nrsamp,caseid,xa,xb,nx)
      call tecmargpdf('y',ysamp,nrsamp,caseid,ya,yb,ny)
      call tecmargpdf('q',qsamp,nrsamp,caseid,qa,qb,nx)
      write(*,'(a,f12.4)')'Alphasum=',alphasum
      write(*,'(a)')'ES-MDA analysis completed'
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! I E S         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (lies) then
   write(*,'(a)')'++++++++++++++++++++++++++++++++++++++++++++++'
!    IES update
      write(*,'(a,2f13.5)',advance='yes')'IES analysis...d,sigo=',d,sigo
      do n=1,nrsamp
         iconv(n)=0
         xsamp(n)=xsampini(n)
         qsamp(n)=qsampini(n)
         ysamp(n)=func(xsamp(n),beta,funcmode)+qsamp(n)
      enddo
      call cov(Cxx,Cyy,Cqq,Cyx,Cqy,Cqx,xsamp,ysamp,qsamp,nrsamp)
      Czz(1,1)=Cxx; Czz(2,2)=Cqq; Czz(1,2)=Cqx; Czz(2,1)=Cqx
      CIzz(1,1)=Cqq; CIzz(2,2)=Cxx; CIzz(1,2)=-Cqx; CIzz(2,1)=-Cqx
      CIzz=CIzz/(Cxx*Cqq-Cqx*Cqx)

      sumconv=0
      do i=1,maxiesit
         if (i > 10) then
            if (mod(i,10) == 0) then
               write(*,'(a,i4,a)',advance='no')'i=',i,'...'
               print *,'converged realizations=', sumconv,nrsamp,real(100*sumconv)/real(nrsamp)
            endif
         endif
         call cov(Pxx,Pyy,Pqq,Pyx,Pqy,Pqx,xsamp,ysamp,qsamp,nrsamp)
         if (lcyyreg) Pyy=cyyreg(Pxx,Pqq,Pyx,Pqy,Pqx)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (IESv == 1) then
!        Old method (Evensen 2018b paper)
            if (sigw > 0.0) then
               print *,'IESv=1 do not support model errors'
               exit
            endif

            dgx=(Pyx)/(Pxx)
            djx=0.0

            do n=1,nrsamp
               if (iconv(n) > 0) cycle ! do nothing for converged realizations

               djx=cdd*(xsamp(n)-xsampini(n))  + dgx*Cxx*(ysamp(n)-dpert(n)) 
               ddjx= Pyy + cdd 

               xsamp(n)=xsamp(n) - gamma_ies*djx/ddjx
               ysamp(n)=func(xsamp(n),beta,funcmode) 

               if ((abs(djx) < 0.000000001) .and. (abs(djq) < 0.000000001)) iconv(n)=1
            enddo

         elseif (IESv == 3 .and. Cqq > 0.0) then
!        Matrix version (takes into account initial correlations between x and q
            Pzz(1,1)=Pxx; Pzz(2,2)=Pqq; Pzz(1,2)=Pqx; Pzz(2,1)=Pqx
            Pzy(1,1)=Pyx; Pzy(2,1)=Pqy
            Pyz(1,1)=Pyx; Pyz(1,2)=Pqy

            do n=1,nrsamp
               if (iconv(n) > 0) cycle ! do nothing for converged realizations

               zf(1)=xsampini(n)
               zf(2)=qsampini(n)

               zi(1)=xsamp(n)
               zi(2)=qsamp(n)

               grad2 = matmul(Pzz,matmul(CIzz,zi-zf)) &
                     - matmul(Pzy , ( matmul(Pyz, matmul(CIzz,zi-zf)) - ysamp(n) + dpert(n) )) / (Pyy + Cdd)

               zi(:) = zi(:) - gamma_ies*grad2(:)

               xsamp(n)=zi(1)
               qsamp(n)=zi(2)
               ysamp(n)=func(xsamp(n),beta,funcmode) + qsamp(n)
      
               if ((abs(grad2(1)) < 0.000000001) .and. (abs(grad2(2)) < 0.000000001)) iconv(n)=1
            enddo

         elseif (IESv == 3 .and. Cqq == 0.0) then
!        Matrix version without model erros
            do n=1,nrsamp
               if (iconv(n) > 0) cycle ! do nothing for converged realizations

               grad1 = (Pxx/Cxx)*(xsamp(n) - xsampini(n)) &
                         -  (Pyx/(Pyy + Cdd))*( (Pyx/Cxx)*(xsamp(n) - xsampini(n)) - ysamp(n) + dpert(n)) 

               xsamp(n) = xsamp(n) - gamma_ies*grad1
               ysamp(n)=func(xsamp(n),beta,funcmode) 
      
               if (abs(grad1) < 0.000000001) iconv(n)=1

            enddo
         endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         sumconv=sum(iconv(1:nrsamp))
         if (sumconv == nrsamp) then
            write(*,'(a,i4,a)',advance='no')'i=',i,'...'
            print *,'converged realizations=', sumconv,nrsamp,real(100*sumconv)/real(nrsamp)
            print *,'Exiting IES iterations'
            exit
         endif

! Dumping costfunction stuff
         if (i < nrits) then
            do j=1,10
               xsampit(j,i+1)=xsamp(j)
               costite(j,i+1)= (xsampit(j,i+1)-xsampini(j))**2/siga**2 + (func(xsampit(j,i+1),beta,funcmode)-dpert(j))**2/cdd    
            enddo
         endif
      enddo
      samples(:,1,1)=xsamp(:)
      samples(:,2,1)=qsamp(:)

! Dumping costfunction stuff
      call tecsampini('sampiniIES',costite,xsampit,min(nrits,i),10)

!    Recomputing ysamp with some noise for nicer plotting
      if (sigw < sigq) then
         do n=1,nrsamp
            ysamp(n)=ysamp(n)+sigq*normal()
         enddo
      endif
      call getcaseid(caseid,'IES',alphageo,nmda,esamp,gradient,beta,sigw,0)
      call tecmargpdf('x',xsamp,nrsamp,caseid,xa,xb,nx)
      call tecmargpdf('y',ysamp,nrsamp,caseid,ya,yb,ny)
      call tecmargpdf('q',qsamp,nrsamp,caseid,qa,qb,nx)
      call tecpdf(x,y,nx,ny,xsamp,ysamp,nrsamp,xa,ya,dx,dy,caseid)
      write(*,'(a)')'IES analysis completed'
   endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! I E n K S     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (lienks) then
!    IEnKS update (IES searching the solution in ensemble subspace)
      allocate (W(nrsamp,nrsamp))
      allocate (WW(nrsamp,nrsamp))
      allocate (Omega(nrsamp,nrsamp))

      write(*,'(a,2f13.5)',advance='yes')'IEnKS analysis...d,sigo=',d,sigo
      do n=1,nrsamp
         Dens(1,n)=dpert(n)
         E0(1,n)=xsampini(n)
         Yi(1,n)=func(E0(1,n),beta,funcmode)
         if (ndim .gt. 1) then
            E0(2,n)=qsampini(n)
            Yi(1,n)=Yi(1,n)+E0(2,n)
         endif
         if (ndim .gt. 2) then
            do i=3,ndim
               E0(i,n)=normal()
            enddo
         endif
      enddo
      n1=sqrt(real(nrsamp-1))
      do i=1,ndim
         A0(i,:)=( E0(i,:) - sum(E0(i,1:nrsamp))/real(nrsamp) )/n1
      enddo

      Ei=E0
      W=0.0

!    Testing AA0
!      allocate (AA0(nrsamp,nrsamp))
!      call cpu_time(start)
!      Ai=A0
!      call aaprojection(Ai,AA0,ndim,nrsamp,1.0)
!      write(*,'(a)')'Checking pseudo inversion A0^+ A0'
!      write(*,'(10g11.3)')AA0(1:10,1:10)
!      call cpu_time(finish)
!      print '("Time = ",g13.5," seconds.")',finish-start

!      call cpu_time(start)
!      Ai=A0
!      call pseudoinv(Ai,A0inv,ndim,nrsamp,truncation)
!      call dgemm('N','N',nrsamp,nrsamp,ndim,1.0,A0inv,nrsamp,A0,ndim,0.0,AA0,nrsamp)
!      write(*,'(a)')'Checking pseudo inversion A0^+ A0'
!      write(*,'(10g11.3)')AA0(1:10,1:10)
!      call cpu_time(finish)
!      print '("Time = ",g13.5," seconds.")',finish-start

      do i=1,maxienksit
         print *,' '
         print '(a,i3)','Iteration: ',i
         YY(1,:)=(Yi(1,:) - sum(Yi(1,1:nrsamp))/real(nrsamp)) / n1

         if (ndim < nrsamp-1 .and. beta /= 0.0 .and. lcyyreg) then
            print *,'Activating AAi projection for Y'
            allocate (AAi(nrsamp,nrsamp))
            do k=1,ndim
               Ai(k,:)=( Ei(k,:) - sum(Ei(k,1:nrsamp))/real(nrsamp) )/n1
            enddo
            call aaprojection(Ai,AAi,ndim,nrsamp,1.0)
            YY=matmul(YY,AAi)
            deallocate(AAi)
         endif

         if (i == 1) then
            H(1,:)=Dens(1,:) - Yi(1,:)             
            S=YY

         else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Computing S1 = Yi * Ai^+ * A0 (for testing and diagnostics)
!            do k=1,ndim
!               Ai(k,:)=( Ei(k,:) - sum(Ei(k,1:nrsamp))/real(nrsamp) )/n1
!            enddo
!            call pseudoinv(Ai,Ainv,ndim,nrsamp,truncation)
!            call dgemm('N','N',1,ndim,nrsamp,1.0,YY,1,Ainv,nrsamp,0.0,YAinv,1)
!            call dgemm('N','N',1,nrsamp,ndim,1.0,YAinv,1,A0,ndim,0.0,S1,1)
!            print '(a,10g16.8)','S1 = YY * Ai^+ * A0     =',S1(1,1:10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Computing S2= Yi * (Ai^+ * Ai) Omega_i^+  
            do n=1,nrsamp
               aveW(n)=sum(W(n,1:nrsamp))/real(nrsamp) 
               WW(n,1:nrsamp)=W(n,1:nrsamp)-aveW(n)
            enddo
            WW=WW/n1
            WW=transpose(WW)
            do m=1,nrobs ! Loop over measurements
               bb(:)=YY(m,:)
               xx(:)=bb(:)
               do k=1,1000
                  xxold(:)=xx(:)
                  xx(:) = bb(:) - matmul(WW,xx)
                  diffx= maxval(abs(xx-xxold))
                  if (diffx < 0.000000001) exit
                  if (k==1000) print *,'Linear solver not converged'
               enddo
               S2(m,:)=xx(:)
            enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Standard LU with multiple rhs
!           call cpu_time(start)
!            WW=transpose(W)/n1
!            do n=1,nrsamp
!               WW(n,n)=WW(n,n)+1.0
!            enddo
!            YT=transpose(YY)
!            call dgesv(nrsamp,1,WW,nrsamp,ipiv,YT,nrsamp,info)
!            if (info > 0) stop 'dgesv singular'
!            S2=transpose(YT)
!            print '(a,10g13.5)','S2c=',S2(1,1:10)
!           call cpu_time(finish)
!           print '("Time = ",f6.3," seconds.")',finish-start
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            print '(a,10g16.8)','S2 = Yi * Omega_i^+     =',S2(1,1:10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if (IEnKS_S == 1 ) then
               S=S1
            elseif (IEnKS_S == 2 ) then
               S=S2
            else
               stop 'incorrect IEnKS_S'
            endif               

            open(10,file='ss.dat')
            write(10,*)'TITLE = "S"'
            write(10,*)'VARIABLES = "n" "S1=Yi*Ai^+*A0" "S2=Yi*Omega_i^+" "S3=S2*A0^+*A0" "S2-S1" "S3-S2"' 
            write(10,*)'ZONE T= "Samples" F=POINT, I=',nrsamp
            do n=1,nrsamp
               write(10,'(i8,6g17.8)')n,S1(1,n),S2(1,n),S3(1,n),S2(1,n)-S1(1,n),S3(1,n)-S1(1,n)
            enddo
            close(10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            H(1,:)=Dens(1,:) - Yi(1,:)             
            call dgemm('N','N',1,nrsamp,nrsamp,1.0,S,1,W,nrsamp,1.0,H,1)
         endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Compute C
         C(1,1)=Cdd
         call dgemm('N','T',1,1,nrsamp,1.0,S,1,S,1,1.0,C,1)
         S(1,:)=S(1,:)/C(1,1) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Update W
         WW=W
         call dgemm('T','N',nrsamp,nrsamp,1,gamma_IEnKS,S,1,H,1,1.0-gamma_IEnKS,W,nrsamp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Compute solution Ei=E0+matmul(A0,W) 
         Ei=E0
         call dgemm('N','N',ndim,nrsamp,nrsamp,1.0,A0,ndim,W,nrsamp,1.0,Ei,ndim)

!    Compute new model prediction
         do n=1,nrsamp
            Yi(1,n)=func(Ei(1,n),beta,funcmode)
            if (ndim>1) Yi(1,n)=Yi(1,n)+Ei(2,n)
         enddo

         samples(:,1,2)=Ei(1,:)
         if (ndim > 1) samples(:,2,2)=Ei(2,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Convergence test
         diffW=maxval(abs(W-WW))
         print '(a,g13.5,a,g13.5)','Diff 1=',diffW

! dumping iterations
         call getcaseid(caseid,'IEnKS',alphageo,nmda,esamp,gradient,beta,sigw,i)
         call tecmargpdf('x',Ei(1,:),nrsamp,caseid,xa,xb,nx)

         if (diffW < 0.0000000001) then
            print '(a,i3)','Exiting at iteration: ',i
            exit
         endif

      enddo

!    End of iteration loop
      samples(:,1,2)=Ei(1,:)
      samples(:,2,2)=Ei(2,:)


!    Recomputing ysamp with some noise for nicer plotting
      if (sigw < sigq) then
         do n=1,nrsamp
            ysamp(n)=ysamp(n)+sigq*normal()
         enddo
      endif
      call getcaseid(caseid,'IEnKS',alphageo,nmda,esamp,gradient,beta,sigw,0)
      call tecmargpdf('x',Ei(1,1:nrsamp),nrsamp,caseid,xa,xb,nx)
      if (ndim > 1) call tecmargpdf('q',Ei(2,1:nrsamp),nrsamp,caseid,qa,qb,nx)
      call tecmargpdf('y',Yi(1,1:nrsamp),nrsamp,caseid,ya,yb,ny)
      write(*,'(a)')'IEnKS analysis completed'



   endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! E n S T E I N !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (lenstein) then
   write(*,'(a)')'++++++++++++++++++++++++++++++++++++++++++++++'
!    EnStein update
      write(*,'(a)',advance='yes')'EnStein analysis...'
      allocate(sf(2),si(2,nrsamp),grad(2,nrsamp),newgrad(2,nrsamp))
!      allocate(kern(nrsamp,nrsamp))

         xlength = 2.0/3.14159265 ; print '(a,f13.5)','xlength=',xlength
!        xlength = 1.0            ; print '(a,f13.5)','xlength=',xlength

      do n=1,nrsamp
!        Initial guess
         iconv(n)=0
         xsamp(n)=xsampini(n)
         qsamp(n)=qsampini(n)
         ysamp(n)=func(xsamp(n),beta,funcmode)+qsamp(n)
!        Prior in cost function 
         sf(1)=x0   
         sf(2)=0.0
      enddo
      call cov(Cxx,Cyy,Cqq,Cyx,Cqy,Cqx,xsamp,ysamp,qsamp,nrsamp)
      Czz(1,1)=Cxx; Czz(2,2)=Cqq; Czz(1,2)=Cqx; Czz(2,1)=Cqx
      CIzz(1,1)=Cqq; CIzz(2,2)=Cxx; CIzz(1,2)=-Cqx; CIzz(2,1)=-Cqx
      CIzz=CIzz/(Cxx*Cqq-Cqx*Cqx)

      call getcaseid(caseid,'EnSTEIN',alphageo,nmda,esamp,gradient,beta,sigw,0)
      sumconv=0
      do i=1,50 !maxiesit
         if (mod(i,10) == 0) then
            write(*,'(a,i4,a)')'i=',i,'...'
         endif
         call cov(Pxx,Pyy,Pqq,Pyx,Pqy,Pqx,xsamp,ysamp,qsamp,nrsamp)
         Pzz(1,1)=Pxx; Pzz(2,2)=Pqq; Pzz(1,2)=Pqx; Pzz(2,1)=Pqx
         PIzz(1,1)=Pqq; PIzz(2,2)=Pxx; PIzz(1,2)=-Pqx; PIzz(2,1)=-Pqx
         PIzz=PIzz/(Pxx*Pqq-Pqx*Pqx)
         Pzy(1,1)=Pyx; Pzy(2,1)=Pqy
         Pyz(1,1)=Pyx; Pyz(1,2)=Pqy

         do n=1,nrsamp
            si(1,n)=xsamp(n)
            si(2,n)=qsamp(n)
         enddo

         do n=1,nrsamp
! Prior term
            call dgemm('N','N',ndim,1,ndim,1.0,CIzz,ndim,si(:,n)-sf(:),ndim,0.0,grad(:,n),ndim)
! C_d^{-1} (y_j - d) with unperturbed data
            fac=(ysamp(n)-d)/sigo**2   
! Data term in EnStein  G=Pxx^{-1} Pxy
            call dgemm('N','N',ndim,1,ndim,fac,PIzz,ndim,Pzy,ndim,1.0,grad(:,n),ndim)
! Data term in Stein with analytic G
!            grad(1,n)=grad(1,n)+dfunc(si(1,n),beta,funcmode)*fac
!            grad(2,n)=grad(2,n)+fac
         enddo


         do n=1,nrsamp
            newgrad(:,n)=0.0
            do k=1,100
               call random_number(tmp)
               j=nint(tmp*real(nrsamp-1)+1.0)
               if (k==1) j=n
               newgrad(:,n)=newgrad(:,n)+&
                       steinkern(si(:,n),si(:,j),xlength,2)*  &
                       (grad(:,j) + (2.0*(si(:,j) - si(:,n))/xlength))
            enddo
            newgrad(:,n)=newgrad(:,n)/real(100.0)
         enddo

         do n=1,nrsamp
            si(:,n) = si(:,n) - 0.25*newgrad(:,n)
            xsamp(n)=si(1,n)
            qsamp(n)=si(2,n)
            ysamp(n)=func(xsamp(n),beta,funcmode) + qsamp(n)
         enddo

         if (mod(i,i)==0) call tecmargpdf('x',xsamp,nrsamp,caseid,xa,xb,nx,i)

      enddo
      samples(:,1,4)=xsamp(:)
      samples(:,2,4)=qsamp(:)

!    Recomputing ysamp with some noise for nicer plotting
      if (sigw < sigq) then
         do n=1,nrsamp
            ysamp(n)=ysamp(n)+sigq*normal()
         enddo
      endif
      call tecmargpdf('x',xsamp,nrsamp,caseid,xa,xb,nx)
      call tecmargpdf('y',ysamp,nrsamp,caseid,ya,yb,ny)
      call tecmargpdf('q',qsamp,nrsamp,caseid,qa,qb,nx)
!      call tecpdf(x,y,nx,ny,xsamp,ysamp,nrsamp,xa,ya,dx,dy,caseid)
      write(*,'(a)')'Stein analysis completed'
   endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! POST PROCESSING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!   open(10,file='samples.dat')
!   write(10,*)'TITLE = "Samples"'
!   write(10,*)'VARIABLES = "n" "ES" "IES" "IEnKS" "IEnKS-IES" "IEnKS-ES" "IES-ES"'
!   write(10,*)'ZONE T= "Samples" F=POINT, I=',nrsamp
!   do n=1,nrsamp
!      write(10,'(i8,6g17.8)')n,samples(n,1,0),samples(n,1,1),samples(n,1,2),&
!                               samples(n,1,2)-samples(n,1,1),&
!                               samples(n,1,2)-samples(n,1,0),&
!                               samples(n,1,1)-samples(n,1,0)
!   enddo
!   close(10)

   do n=1,min(ndim,2)
   do j=0,4
      call avevar(samples(1:nrsamp,n,j),nrsamp,ave,var)
      if (var > 0.0) write(*,'(4a,2g16.8,a,8g16.8)')method(j),' ',variable(n),' :',ave,sqrt(var),'  Samples:',samples(1:8,n,j)
   enddo
   enddo
end program 
