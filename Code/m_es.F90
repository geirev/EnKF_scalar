module m_es
use mod_inistat
use mod_xyqgrid
use m_cyyreg
use m_cov
use m_getcaseid
use m_func
use m_normal
use m_tecmargpdf
implicit none
logical, save :: les=.true.               ! Run ES or not
logical, save :: lesadjoint=.false.       ! Run ES with adjoint sensitivites
contains
subroutine es(samples,xf,qf,dpert,nrsamp)
   integer, intent(in)  :: nrsamp
   real,    intent(out) :: samples(nrsamp,2)
   real,    intent(in)  :: xf(nrsamp) 
   real,    intent(in)  :: qf(nrsamp) 
   real,    intent(in)  :: dpert(nrsamp) 

   real, allocatable :: xsamp(:)
   real, allocatable :: qsamp(:)
   real, allocatable :: ysamp(:)

   integer i
   real Cxx,Cyy,Cqq,Cxy,Cqy,Cqx
!   real Czz(2,2)
   real dg(2) ! model sensitivity wrt x and q

   allocate(xsamp(nrsamp))
   allocate(qsamp(nrsamp))
   allocate(ysamp(nrsamp))

   write(*,'(a)')'++++++++++++++++++++++++++++++++++++++++++++++'
   write(*,'(a)')'ES analysis...'
   do i=1,nrsamp
      xsamp(i)=xf(i)
      qsamp(i)=qf(i)
      ysamp(i)=func(xsamp(i),qsamp(i))
   enddo

   call cov(Cxx,Cyy,Cqq,Cxy,Cqy,Cqx,xsamp,ysamp,qsamp,nrsamp)

   if (lcyyreg) cyy=cyyreg(Cxx,Cqq,Cxy,Cqy,Cqx)

   do i=1,nrsamp
      if (lesadjoint) then
         dg=dfunc(xsamp(i),qsamp(i))
         ! Czz * G^T
         Cxy=siga**2 * dg(1) +     0.0 * dg(2)
         Cqy=    0.0 * dg(1) + sigw**2 * dg(2) 

         ! G * Czz * G^T
         Cyy=dg(1)*Cxy + dg(2)*Cqy
      endif

      xsamp(i)=xsamp(i) + (Cxy/(Cyy+cdd))*(dpert(i)-ysamp(i))
      qsamp(i)=qsamp(i) + (Cqy/(Cyy+cdd))*(dpert(i)-ysamp(i))
      ysamp(i)=func(xsamp(i),qsamp(i))

!      if (updatemode==0) then
!         ysamp(i)=ysamp(i)+ (cyy/(cyy+cdd))*(dpert(i)-ysamp(i))
!      else
!      endif

   enddo

   samples(:,1)=xsamp(:)
   samples(:,2)=qsamp(:)
   deallocate(xsamp,ysamp,qsamp)
   write(*,'(a)')'ES analysis completed'
   write(*,'(a)')'++++++++++++++++++++++++++++++++++++++++++++++'
   write(*,'(a)')
end subroutine
end module
