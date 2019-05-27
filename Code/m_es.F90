module m_es
use mod_inistat
use mod_xyqgrid
use m_cyyreg
use m_cov
use m_getcaseid
use m_func
use m_normal
use m_tecmargpdf
use m_tecpdf
implicit none
logical, save :: les=.true.               ! Run IES or not
logical, save :: lesadjoint=.false. ! Run IES with adjoint sensitivites
contains
subroutine es(samples,xf,qf,dpert,nrsamp,esamp)
   integer, intent(in)  :: nrsamp
   integer, intent(in)  :: esamp
   real,    intent(out) :: samples(nrsamp,2)
   real,    intent(in)  :: xf(nrsamp) 
   real,    intent(in)  :: qf(nrsamp) 
   real,    intent(in)  :: dpert(nrsamp) 

   real, allocatable :: xsamp(:)
   real, allocatable :: qsamp(:)
   real, allocatable :: ysamp(:)

   integer i
   real Cxx,Cyy,Cqq,Cyx,Cqy,Cqx
   real dgx
!   character(len=40) caseid

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

   call cov(Cxx,Cyy,Cqq,Cyx,Cqy,Cqx,xsamp,ysamp,qsamp,nrsamp)

   if (lcyyreg) cyy=cyyreg(Cxx,Cqq,Cyx,Cqy,Cqx)

   do i=1,nrsamp
      if (lesadjoint) then
         dgx=dfunc(xsamp(i),qsamp(i))
         cyx=dgx*siga**2
         cyy=dgx*siga**2*dgx
      endif

      xsamp(i)=xsamp(i) + (cyx/(cyy+cdd))*(dpert(i)-ysamp(i))
      qsamp(i)=qsamp(i) + (cqy/(cyy+cdd))*(dpert(i)-ysamp(i))
      ysamp(i)=func(xsamp(i),qsamp(i))

!      if (updatemode==0) then
!         ysamp(i)=ysamp(i)+ (cyy/(cyy+cdd))*(dpert(i)-ysamp(i))
!      else
!      endif

   enddo

!   if (sigw < sigq) then
!      do i=1,nrsamp
!         ysamp(i)=ysamp(i)+sigq*normal()
!      enddo
!   endif
!   call getcaseid(caseid,'ES',-1.0,-1,esamp,sigw,0)
!   call tecpdf(x,y,nx,ny,xsamp,ysamp,nrsamp,xa,ya,dx,dy,caseid)
!   call tecmargpdf('x',xsamp,nrsamp,caseid,xa,xb,nx)
!   call tecmargpdf('y',ysamp,nrsamp,caseid,ya,yb,ny)
!   call tecmargpdf('q',qsamp,nrsamp,caseid,qa,qb,nx)
!
   samples(:,1)=xsamp(:)
   samples(:,2)=qsamp(:)
   deallocate(xsamp,ysamp,qsamp)
   write(*,'(a)')'ES analysis completed'
   write(*,'(a)')'++++++++++++++++++++++++++++++++++++++++++++++'
   write(*,'(a)')
end subroutine
end module
