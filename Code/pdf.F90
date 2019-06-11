program iterative_smoothers
   use mod_inistat
   use mod_xyqgrid
   use m_inipdfs
   use m_marginalpdf
   use m_set_random_seed2
   use m_avevar
   use m_func
   use m_cyyreg
   use m_iniens
   use m_es
   use m_ies
   use m_iese
   use m_sies
   use m_esmda
   use m_enstein
   use m_tecpdf
   use m_tecmargpdf
   use m_printcostf
   use m_tecjointpdf
   use m_tecfunc
   use m_teccostens
   use m_normal
   use m_integrals
   implicit none
   integer nrsamp,esamp

   integer, parameter :: nrmethods=6
   character(len=7) :: method(0:nrmethods) =['INI    ','ES     ','IES    ','SIES   ','ESMDA  ','EnSTEIN','IESE   ']
   logical          :: lactive(0:nrmethods)=(/.true.  ,.true.   ,.true.   ,.true.   ,.true.   ,.true.   ,.true.   /)
   character(len=1) :: variable(1:2)=['x','q']

   real, allocatable :: xsampini(:), xsamp(:)
   real, allocatable :: qsampini(:), qsamp(:)
   real, allocatable :: ysampini(:), ysamp(:)
   real, allocatable :: dpert(:)
   real, allocatable :: samples(:,:,:) 

   integer i,j,n
   real ave,var
   character(len=2) ca
   character(len=40) caseid
   logical ladjoints

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call set_random_seed2()
   print *,'Normal number', normal()

   open(10,file='infile.in')
      read(10,*)esamp       ; print '(a,i4,i12)',   'number  of samples 10^x    :',esamp,10**esamp
      read(10,'(a)')ca      
      if (ca /= '#1') then
         print *,'#1: error in infile.in'
         stop
      endif
      read(10,*)xa            ; print '(a,f10.3)',    'xa                         :',xa
      read(10,*)xb            ; print '(a,f10.3)',    'xb                         :',xb
      read(10,*)qa            ; print '(a,f10.3)',    'qa                         :',qa
      read(10,*)qb            ; print '(a,f10.3)',    'qb                         :',qb
      read(10,*)ya            ; print '(a,f10.3)',    'ya                         :',ya
      read(10,*)yb            ; print '(a,f10.3)',    'yb                         :',yb
      read(10,'(a)')ca
      if (ca /= '#2') then
         print *,'#2: error in infile.in'
         stop
      endif
      read(10,*)x0            ; print '(a,f10.3)',    'prior estimate of x        :',x0
      read(10,*)siga          ; print '(a,f10.3)',    'err std dev of prior       :',siga
      read(10,*)d             ; print '(a,f10.3)',    'measurement of y           :',d
      read(10,*)sigo          ; print '(a,f10.3)',    'measurement std dev        :',sigo ; cdd=sigo**2
      read(10,*)sigw          ; print '(a,f10.3)',    'model std dev              :',sigw
      read(10,'(a)')ca
      if (ca /= '#3') then
         print *,'#3: error in infile.in'
         stop
      endif
      read(10,*)beta          ; print '(a,f10.3)',    'model parameter beta       :',beta
      read(10,*)funcmode      ; print '(a,tr7,i3)',   'function to use            :',funcmode
      read(10,'(a)')ca
      if (ca /= '#4') then
         print *,'#4: error in infile.in'
         stop
      endif
      read(10,*)lcyyreg       ; print '(a,tr10,l1)',  'Regression for Cyy         :',lcyyreg
      read(10,*)ladjoints     ; print '(a,tr10,l1)',  'Adjoint model sens         :',ladjoints
      read(10,*)lactivepdf    ; print '(a,tr10,l1)',  'Activate pdf printing      :',lactivepdf
      read(10,'(a)')ca
      if (ca /= '#5') then
         print *,'#5: error in infile.in'
         stop
      endif
      read(10,*)lmda          ; print '(a,tr10,l1)',  'Run MDA                    :',lmda; lactive(4)=lmda
      read(10,*)nmda          ; print '(a,tr7,i3)',   'number of mda iterations   :',nmda
      read(10,*)alphageo      ; print '(a,f10.3)',    'geometrical alpha value    :',alphageo
      read(10,'(a)')ca
      if (ca /= '#6') then
         print *,'#6: error in infile.in'
         stop
      endif
      read(10,*)lies          ; print '(a,tr10,l1)',  'Run IES                    :',lies ; lactive(2)=lies
      read(10,*)maxiesit      ; print '(a,i5)',       'maxiesit                   :',maxiesit
      read(10,*)gamma_ies     ; print '(a,f10.3)',    'gamma_ies                  :',gamma_ies
      read(10,*)IESv          ; print '(a,i1)',       'IESv                       :',IESv
      read(10,'(a)')ca
      if (ca /= '#7') then
         print *,'#7: error in infile.in'
         stop
      endif
      read(10,*)lsies         ; print '(a,tr10,l1)',  'Run SIES                   :',lsies ; lactive(3)=lsies
      read(10,*)maxsiesit     ; print '(a,i5)',       'maxiesit                   :',maxsiesit
      read(10,*)gamma_sies    ; print '(a,f10.3)',    'gamma_sies                 :',gamma_sies
      read(10,'(a)')ca
      if (ca /= '#8') then
         print *,'#8: error in infile.in'
         stop
      endif
      read(10,*)lenstein      ; print '(a,tr10,l1)',  'Run EnSTEIN                :',lenstein; lactive(5)=lenstein
      read(10,*)maxensteinit  ; print '(a,i5)',       'maxiesit                   :',maxensteinit
      read(10,*)gamma_enstein ; print '(a,f10.3)',    'gamma_sies                 :',gamma_enstein
   close(10)

   if (ladjoints) then
      if (sigw > 0.0 ) stop 'Must have sigw=0.0 with ladjoints=true'
      if (lies .and. IESv /= 1 )  stop 'Must have IESv=1 with ladjoints=true'
      lesadjoint=.true.
      lesmdaadjoint=.true.
      liesadjoint=.true.
      if (lsies) lsies=.false.
      if (lenstein) lenstein=.false.
   endif

   nrsamp=10**esamp
   nrsamp=1*nrsamp

   allocate(qsampini(nrsamp), xsampini(nrsamp) , ysampini(nrsamp), dpert(nrsamp))
   allocate(samples(nrsamp,2,0:nrmethods)); samples=0.0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Theroretical values for statistical moements (trustat.dat)
!   call integrals(maxiesit)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Grid domain for plotting pdfs in x, q, and y
   call xyqgrid()

   call compute_prior()

   call compute_likelihood()

   call compute_costfunction()

   call compute_uncond_jointpdf(esamp)

   call compute_cond_jointpdf(esamp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Ensemble initialization
   call iniens(samples(1:nrsamp,1:2,0),xsampini,qsampini,ysampini,dpert,nrsamp,esamp)

! ES 
   call es(samples(1:nrsamp,1:2,1),xsampini,qsampini,dpert,nrsamp,esamp)

! IES
   if (lies) then
      call ies(samples(1:nrsamp,1:2,2),xsampini,qsampini,nrsamp,esamp,dpert)
      call iese(samples(1:nrsamp,1:2,6),xsampini,qsampini,nrsamp,esamp,dpert)
   endif

! Subspace IES
   if (lsies) then
      call sies(samples(1:nrsamp,1:2,3),xsampini,qsampini,nrsamp,esamp,dpert)
   endif

! ESMDA
   if (lmda) then
      call esmda(samples(1:nrsamp,1:2,4),xsampini,qsampini,nrsamp,esamp)
   endif

! EnSTEIN
   if (lenstein) then
      call enstein(samples(1:nrsamp,1:2,5),xsampini,qsampini,nrsamp,esamp)
   endif

! POST PROCESSING
!   call avevar(dpert,nrsamp,ave,var)
!   if (var > 0.0) write(*,'(4a,2g16.8,a,8g16.8)')'MEAS   ',' ',variable(1),' :',ave,sqrt(var),'  Samples:',dpert(1:8)
   do n=1,2
   do j=1,nrmethods
      call avevar(samples(1:nrsamp,n,j),nrsamp,ave,var)
      if (var > 0.0) write(*,'(4a,2g16.8,a,8g16.8)')method(j),' ',variable(n),' :',ave,sqrt(var),'  Samples:',samples(1:8,n,j)
   enddo
   enddo



   allocate(xsamp(nrsamp), qsamp(nrsamp), ysamp(nrsamp))
   do j=0,nrmethods
      if (lactive(j)) then
         print '(a,a)','printing :',trim(method(j))
         if (ladjoints) then
            if (trim(method(j)) == 'ES')    method(j)='ESAD   '
            if (trim(method(j)) == 'IES')   method(j)='IESAD  '
            if (trim(method(j)) == 'ESMDA') method(j)='ESMDAD '
         endif
         do i=1,nrsamp
            xsamp(i)=samples(i,1,j)
            qsamp(i)=samples(i,2,j)
            ysamp(i)=func(xsamp(i),ysamp(i))
            if (sigw < sigq) ysamp(i)=ysamp(i)+sigq*normal()
         enddo
         if (trim(method(j)) == 'ESMDA') then
            call getcaseid(caseid,method(j),alphageo,nmda,esamp,sigw,0)
         else
            call getcaseid(caseid,method(j),-1.0,-1,esamp,sigw,0)
         endif
         call tecpdf(x,y,nx,ny,xsamp,ysamp,nrsamp,xa,ya,dx,dy,caseid)
         call tecmargpdf('x',xsamp,nrsamp,caseid,xa,xb,nx)
         call tecmargpdf('y',ysamp,nrsamp,caseid,ya,yb,ny)
         if (sigw > 0.0) call tecmargpdf('q',qsamp,nrsamp,caseid,qa,qb,nx)
      endif
   enddo
   deallocate(xsamp, qsamp, ysamp)

! Printing some cost functions and the ES solution before and after (see also dead code in ies
!   call printcostf(samples(1:nrsamp,1,0),samples(1:nrsamp,1,1),dpert,nrsamp)

end program 
