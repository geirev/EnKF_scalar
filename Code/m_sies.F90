module m_sies
use mod_inistat
use mod_xyqgrid
use m_cyyreg
use m_cov
use m_aaprojection
use m_getcaseid
use m_func
use m_normal
use m_tecmargpdf
use m_tecsampini
implicit none
logical, save :: lsies         ! Run SIES or not
integer, save :: maxsiesit     ! Maximum number of iterations
real,    save :: gamma_sies    ! steplength
integer, parameter :: nrobs=1
integer, parameter :: ndim=2
real, save :: sies_eps=0.0000001
contains 
   subroutine sies(samples,xf,qf,dpert,nrsamp,esamp)
   integer, intent(in)  :: nrsamp            ! Number of samples
   integer, intent(in)  :: esamp             ! Number of samples nrsamp=10^esamp for plotting
   real,    intent(in)  :: xf(nrsamp)        ! Prior samples of x
   real,    intent(in)  :: qf(nrsamp)        ! Prior samples of q
   real,    intent(in)  :: dpert(nrsamp)     ! The perturbed measurement ensemble
   real,    intent(out) :: samples(nrsamp,2) ! Returns posterior samples

   real, allocatable :: xsamp(:)
   real, allocatable :: qsamp(:)
   real, allocatable :: ysamp(:)
   integer, allocatable :: iconv(:)

   real, allocatable :: YAinv(:,:)
   real, allocatable :: W(:,:)   ! update coefficients
   real, allocatable :: WW(:,:)   ! update coefficients from previous iterate
   real, allocatable :: Yi(:,:)   ! predicted measurements 
   real, allocatable :: YY(:,:)  ! predicted measurement anomalies 
   real, allocatable :: S(:,:)   ! predicted measurement anomalies Y times W^+
   real, allocatable :: YT(:,:),ST(:,:),STO(:,:)
   real, allocatable :: Dens(:,:)   ! Ensemble of perturbed measurements
   real, allocatable :: H(:,:)     ! "Innovation "
   real, allocatable :: E0(:,:)   ! Initial model ensemble 
   real, allocatable :: Ei(:,:)  ! model ensemble iteration i
   real, allocatable :: A0(:,:)   ! Initial model ensemble perturbation
   real, allocatable :: Ai(:,:)   ! model ensemble perturbation at iteration i
   real, allocatable :: Ainv(:,:)   !  Model ensemble perturbation at iteration i
   real, allocatable :: A0inv(:,:)   !  Model ensemble perturbation at iteration i
   real, allocatable :: AAi(:,:)   !  Ai^+Ai
   real, allocatable :: aveW(:)   
   real, allocatable :: xx(:),xxold(:),bb(:)
   integer, allocatable :: ipiv(:)

   integer n,i,k,m
   real n1,diffW,diffx,C(1,1)
   character(len=40) caseid

   allocate(xsamp(nrsamp))
   allocate(qsamp(nrsamp))
   allocate(ysamp(nrsamp))
   allocate(iconv(nrsamp)) 


   allocate (W(nrsamp,nrsamp), WW(nrsamp,nrsamp))
   allocate (Yi(1,nrsamp))
   allocate (YY(1,nrsamp))
   allocate (YAinv(1,ndim))
   allocate (S(1,nrsamp))
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
   allocate(xx(nrsamp),bb(nrsamp), xxold(nrsamp))
   allocate (ipiv(nrsamp))

   write(*,'(a)')'++++++++++++++++++++++++++++++++++++++++++++++'
   write(*,'(a)')'SIES analysis...'

   n1=sqrt(real(nrsamp-1))
   do n=1,nrsamp
      Dens(1,n)=dpert(n)
      E0(1,n)=xf(n)
      E0(2,n)=qf(n)
      Yi(1,n)=func(E0(1,n),E0(2,n))
   enddo

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

   do i=1,maxsiesit
      if (mod(i,1) == 0) then
         write(*,'(a,i4)',advance='yes')'iteration=',i
      endif

      YY(1,:)=(Yi(1,:) - sum(Yi(1,1:nrsamp))/real(nrsamp)) / n1


      if (ndim < nrsamp-1 .and. beta /= 0.0 .and. lcyyreg) then
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
!    Computing S = Yi * Ai^+ * A0 (for testing and diagnostics)
!            do k=1,ndim
!               Ai(k,:)=( Ei(k,:) - sum(Ei(k,1:nrsamp))/real(nrsamp) )/n1
!            enddo
!
!            call pseudoinv(Ai,Ainv,ndim,nrsamp,truncation)
!            call dgemm('N','N',1,ndim,nrsamp,1.0,YY,1,Ainv,nrsamp,0.0,YAinv,1)
!            call dgemm('N','N',1,nrsamp,ndim,1.0,YAinv,1,A0,ndim,0.0,S,1)
!            print '(a,10g16.8)','S = YY * Ai^+ * A0     =',S(1,1:10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Computing S= Yi * (Ai^+ * Ai) Omega_i^+  
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
            S(m,:)=xx(:)
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
!            S=transpose(YT)
!            print '(a,10g13.5)','Sc=',S(1,1:10)
!           call cpu_time(finish)
!           print '("Time = ",f6.3," seconds.")',finish-start
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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
      call dgemm('T','N',nrsamp,nrsamp,1,gamma_sies,S,1,H,1,1.0-gamma_sies,W,nrsamp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Compute solution Ei=E0+matmul(A0,W) 
      Ei=E0
      call dgemm('N','N',ndim,nrsamp,nrsamp,1.0,A0,ndim,W,nrsamp,1.0,Ei,ndim)

!    Compute new model prediction
      do n=1,nrsamp
         Yi(1,n)=func(Ei(1,n),Ei(2,n))
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Convergence test
      diffW=maxval(abs(W-WW))

! dumping iterations
      call getcaseid(caseid,'SIES   ',-1.0,-1,esamp,sigw,i)
      call tecmargpdf('x',Ei(1,:),nrsamp,caseid,xa,xb,nx)

      if (diffW < sies_eps) then
         print '(a,i3)','Exiting at iteration: ',i
         exit
      endif

   enddo

   samples(:,1)=Ei(1,:)
   samples(:,2)=Ei(2,:)

   deallocate(xsamp,ysamp,qsamp,iconv)
   deallocate(W, WW )
   deallocate(Yi, YY, YAinv, S, YT,ST,STO, Dens, H, E0, Ei, A0)
   deallocate(Ai, Ainv, A0inv, aveW, xx,bb, xxold, ipiv) 

   write(*,'(a)')'SIES analysis completed'
   write(*,'(a)')'++++++++++++++++++++++++++++++++++++++++++++++'
   write(*,'(a)')

end subroutine
end module
