module m_ies
use mod_inistat
use mod_xyqgrid
use m_cyyreg
use m_cov
use m_getcaseid
use m_func
use m_normal
use m_tecmargpdf
use m_tecsampini
implicit none
logical, save :: lies                ! Run IES or not
logical, save :: liesadjoint=.false. ! Run IES with adjoint sensitivites
integer, save :: maxiesit            ! maximumn number of iterations
real,    save :: gamma_ies           ! step length used in IES
integer, save :: IESv                ! version of IES implementation
real, save :: ies_eps=0.000001

contains 
   subroutine ies(samples,xf,qf,dpert,nrsamp)
   implicit none
   integer, intent(in)  :: nrsamp            ! Number of samples
   real,    intent(out) :: samples(nrsamp,2) ! Returns posterior samples
   real,    intent(in)  :: xf(nrsamp)        ! Prior samples of x
   real,    intent(in)  :: qf(nrsamp)        ! Prior samples of q
   real,    intent(in)  :: dpert(nrsamp)     ! The perturbed measurement ensemble

   real, allocatable :: xsamp(:)
   real, allocatable :: qsamp(:)
   real, allocatable :: ysamp(:)
   integer, allocatable :: iconv(:)
!   integer :: nrits=100
!   real, allocatable :: costens(:,:)
!   real, allocatable :: costite(:,:)
!   real, allocatable :: xsampit(:,:)

   integer n,i,sumconv
   real grad1,grad2(2)
   real Cxx,Cyy,Cqq,Cxy,Cqy,Cqx
   real Czz(2,2),CIzz(2,2)
   real Pzz(2,2),PIzz(2,2)
   real zi(2),zf(2)
   real Pzy(2,1),Pyz(1,2)
   real Pxx,Pyy,Pqq,Pxy,Pqy,Pqx
   real dj(2),ddj,dg(1,2)
!   character(len=40) caseid

   allocate(xsamp(nrsamp))
   allocate(qsamp(nrsamp))
   allocate(ysamp(nrsamp))
   allocate(iconv(nrsamp)) 
!   allocate(costens(nx,10))
!   allocate(costite(10,nrits))
!   allocate(xsampit(10,nrits))

   write(*,'(a)')'++++++++++++++++++++++++++++++++++++++++++++++'
   write(*,'(a)',advance='yes')'IES analysis'

   if (liesadjoint .and. IESv /= 1) then
      print *,'liesadjoint only implemented for IESv=1'
      stop
   endif

   do n=1,nrsamp
      iconv(n)=0
      xsamp(n)=xf(n)
      qsamp(n)=qf(n)
      ysamp(n)=func(xsamp(n),qsamp(n))
   enddo
   call cov(Cxx,Cyy,Cqq,Cxy,Cqy,Cqx,xsamp,ysamp,qsamp,nrsamp)
   Czz(1,1)=Cxx; Czz(2,2)=Cqq; Czz(1,2)=Cqx; Czz(2,1)=Cqx
   CIzz(1,1)=Cqq; CIzz(2,2)=Cxx; CIzz(1,2)=-Cqx; CIzz(2,1)=-Cqx
   if (Cqq > 0.0) CIzz=CIzz/(Cxx*Cqq-Cqx*Cqx)

   sumconv=0
   do i=1,maxiesit
      if (mod(i,1) == 0) then
         write(*,'(a,i4,a)',advance='no')'i=',i,'...'
         print '(a,i10,i10,f10.2)','converged realizations=', sumconv,nrsamp,real(100)*(real(sumconv)/real(nrsamp))
      endif
      call cov(Pxx,Pyy,Pqq,Pxy,Pqy,Pqx,xsamp,ysamp,qsamp,nrsamp)

      Pzz(1,1)=Pxx; Pzz(2,2)=Pqq; Pzz(1,2)=Pqx; Pzz(2,1)=Pqx
      Pzy(1,1)=Pxy; Pzy(2,1)=Pqy
      Pyz(1,1)=Pxy; Pyz(1,2)=Pqy
      PIzz(1,1)=Pqq; PIzz(2,2)=Pxx; PIzz(1,2)=-Pqx; PIzz(2,1)=-Pqx
      if (Pqq > 0.0) PIzz=PIzz/(Pxx*Pqq-Pqx*Pqx)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (IESv == 1) then
!        Old method (Evensen 2018b paper) but using correct Cyy 

         do n=1,nrsamp
            if (iconv(n) > 0) cycle 
            ddj= Pyy + Cdd 
            if (liesadjoint) then
               dg(1,:)=dfunc(xsamp(n),qsamp(n))
               dj(1)=(xsamp(n)-xf(n))/siga**2 + dg(1,1)*(ysamp(n)-dpert(n))/sigo**2
               dj(2)=(qsamp(n)-qf(n))/sigw**2 + dg(1,2)*(ysamp(n)-dpert(n))/sigo**2  
            else
               dg=matmul(Pyz,PIzz)
               dj(1)=(xsamp(n)-xf(n))/Cxx  + dg(1,1)*(ysamp(n)-dpert(n))/Cdd
               dj(2)=(qsamp(n)-qf(n))/Cqq  + dg(1,2)*(ysamp(n)-dpert(n))/Cdd
            endif


            xsamp(n)=xsamp(n) - gamma_ies*dj(1)/ddj
            qsamp(n)=qsamp(n) - gamma_ies*dj(2)/ddj
            ysamp(n)=func(xsamp(n),qsamp(n)) 

            if (abs(dj(1)) < ies_eps) iconv(n)=1
         enddo

      elseif (IESv == 3 .and. Cqq > 0.0) then
! Matrix version (takes into account initial correlations between x and q

         do n=1,nrsamp
            if (iconv(n) > 0) cycle ! do nothing for converged realizations

            zf(1)=xf(n)
            zf(2)=qf(n)

            zi(1)=xsamp(n)
            zi(2)=qsamp(n)

            grad2 = matmul(Pzz,matmul(CIzz,zi-zf)) &
                  - matmul(Pzy , ( matmul(Pyz, matmul(CIzz,zi-zf)) - ysamp(n) + dpert(n) )) / (Pyy + Cdd)

            zi(:) = zi(:) - gamma_ies*grad2(:)

            xsamp(n)=zi(1)
            qsamp(n)=zi(2)
            ysamp(n)=func(xsamp(n),qsamp(n))
   
            if ((abs(grad2(1)) < ies_eps) .and. (abs(grad2(2)) < ies_eps)) iconv(n)=1
         enddo

      elseif (IESv == 3 .and. Cqq == 0.0) then
!        Matrix version without model erros
         do n=1,nrsamp
            if (iconv(n) > 0) cycle 

            grad1 = (Pxx/Cxx)*(xsamp(n) - xf(n)) &
                      -  (Pxy/(Pyy + Cdd))*( (Pxy/Cxx)*(xsamp(n) - xf(n)) - ysamp(n) + dpert(n)) 

            xsamp(n) = xsamp(n) - gamma_ies*grad1
            qsamp(n)=0.0
            ysamp(n)=func(xsamp(n),qsamp(n)) 
   
            if (abs(grad1) < ies_eps) iconv(n)=1

         enddo
      else
         print *,'No matching options for IESv (1,3):',IESv
         stop
      endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      sumconv=sum(iconv(1:nrsamp))
      if (sumconv == nrsamp) then
         write(*,'(a,i4,a)',advance='no')'i=',i,'...'
         print '(a,i10,i10,f10.2)','converged realizations=', sumconv,nrsamp,real(100)*(real(sumconv)/real(nrsamp))
         print '(a)','Exiting IES iterations'
         exit
      endif

! Dumping costfunction stuff for plotting the updates for the first nrits iterates of IES
!      if (i < nrits) then
!         do j=1,10
!            xsampit(j,i+1)=xsamp(j)
!            costite(j,i+1)= (xsamp(j)-xf(j))**2/siga**2 + (func(xsamp(j),qsamp(j))-dpert(j))**2/Cdd    
!         enddo
!      endif

   enddo
   samples(:,1)=xsamp(:)
   samples(:,2)=qsamp(:)

! Dumping costfunction stuff
!   call tecsampini('sampiniIES',costite,xsampit,min(nrits,i),10)

   deallocate(xsamp,ysamp,qsamp,iconv)
   write(*,'(a)')'IES analysis completed'
   write(*,'(a)')'++++++++++++++++++++++++++++++++++++++++++++++'
   write(*,'(a)')

end subroutine
end module
