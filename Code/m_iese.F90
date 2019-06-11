module m_iese
use mod_inistat
use mod_xyqgrid
use m_cyyreg
use m_cov
use m_getcaseid
use m_func
use m_iescostf
use m_normal
use m_tecmargpdf
use m_tecsampini
use m_tecpdf
implicit none
logical, save :: liese=.true.         ! Run IESE or not
integer, save :: maxieseit=400        ! maximumn number of iterations
real,    save :: gamma_iese=0.4       ! step length used in IESE
real,    save :: iese_eps=0.0000001

contains 
   subroutine iese(samples,xf,qf,nrsamp,esamp,dpert)
   implicit none
   integer, intent(in)  :: nrsamp            ! Number of samples
   integer, intent(in)  :: esamp             ! Number of samples nrsamp=10^esamp for plotting
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
   real Cxx,Cyy,Cqq,Cyx,Cqy,Cqx
   real Czz(2,2),Pzz(2,2),CIzz(2,2),zi(2),zf(2),Pzy(2,1),Pyz(1,2)
   real pxx,Pyy,pqq,pyx,pqy,pqx
   real djx,ddjx,dgx
!   character(len=40) caseid

   allocate(xsamp(nrsamp))
   allocate(qsamp(nrsamp))
   allocate(ysamp(nrsamp))
   allocate(iconv(nrsamp)) 
!   allocate(costens(nx,10))
!   allocate(costite(10,nrits))
!   allocate(xsampit(10,nrits))

   write(*,'(a)')'++++++++++++++++++++++++++++++++++++++++++++++'
   write(*,'(a)',advance='yes')'IESE analysis'

   do n=1,nrsamp
      iconv(n)=0
      xsamp(n)=xf(n)
      qsamp(n)=qf(n)
      ysamp(n)=func(xsamp(n),qsamp(n))
   enddo
   call cov(Cxx,Cyy,Cqq,Cyx,Cqy,Cqx,xsamp,ysamp,qsamp,nrsamp)
   Czz(1,1)=Cxx; Czz(2,2)=Cqq; Czz(1,2)=Cqx; Czz(2,1)=Cqx
   CIzz(1,1)=Cqq; CIzz(2,2)=Cxx; CIzz(1,2)=-Cqx; CIzz(2,1)=-Cqx
   if (Cqq > 0.0) CIzz=CIzz/(Cxx*Cqq-Cqx*Cqx)

   sumconv=0
   do i=1,maxieseit
      if (mod(i,1) == 0) then
         write(*,'(a,i4,a)',advance='no')'i=',i,'...'
         print *,'converged realizations=', sumconv,nrsamp,real(100*sumconv)/real(nrsamp)
      endif
      call cov(Pxx,Pyy,Pqq,Pyx,Pqy,Pqx,xsamp,ysamp,qsamp,nrsamp)
      if (lcyyreg) Pyy=cyyreg(Pxx,Pqq,Pyx,Pqy,Pqx)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Matrix version (takes into account initial correlations between x and q
      Pzz(1,1)=Pxx; Pzz(2,2)=Pqq; Pzz(1,2)=Pqx; Pzz(2,1)=Pqx
      Pzy(1,1)=Pyx; Pzy(2,1)=Pqy
      Pyz(1,1)=Pyx; Pyz(1,2)=Pqy

      do n=1,nrsamp
         if (iconv(n) > 0) cycle ! do nothing for converged realizations

         zf(1)=xf(n)
         zf(2)=qf(n)

         zi(1)=xsamp(n)
         zi(2)=qsamp(n)
         dpert(n)=zi(2)

         hess=CIzz + (matmul(Pyz,PIzz)  I  ) Cdd
         grad2 = matmul(Pzz,matmul(CIzz,zi-zf)) &
               - matmul(Pzy , ( matmul(Pyz, matmul(CIzz,zi-zf)) - ysamp(n) + dpert(n) )) / (Pyy + Cdd)

         zi(:) = zi(:) - gamma_iese*grad2(:)

         xsamp(n)=zi(1)
         qsamp(n)=zi(2)
         ysamp(n)=func(xsamp(n),qsamp(n))

         if ((abs(grad2(1)) < iese_eps) .and. (abs(grad2(2)) < iese_eps)) iconv(n)=1
      enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      sumconv=sum(iconv(1:nrsamp))
      if (sumconv == nrsamp) then
         write(*,'(a,i4,a)',advance='no')'i=',i,'...'
         print *,'converged realizations=', sumconv,nrsamp,real(100*sumconv)/real(nrsamp)
         print *,'Exiting IESE iterations'
         exit
      endif

   enddo
   samples(:,1)=xsamp(:)
   samples(:,2)=qsamp(:)


   deallocate(xsamp,ysamp,qsamp,iconv)
   write(*,'(a)')'IESE analysis completed'
   write(*,'(a)')'++++++++++++++++++++++++++++++++++++++++++++++'
   write(*,'(a)')

end subroutine
end module
