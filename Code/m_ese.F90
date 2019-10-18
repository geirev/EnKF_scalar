module m_ese
use mod_inistat
use mod_xyqgrid
use m_cyyreg
use m_cov
use m_getcaseid
use m_func
use m_normal
use m_tecmargpdf
implicit none
logical, save :: lese=.false.         ! Run ESe or not
logical, save :: leseadjoint=.false. ! Run ESe with adjoint sensitivites
contains
subroutine ese(samples,xf,qf,dpert,nrsamp)
   integer, intent(in)  :: nrsamp
   real,    intent(out) :: samples(nrsamp,2)
   real,    intent(in)  :: xf(nrsamp) 
   real,    intent(in)  :: qf(nrsamp) 
   real,    intent(in)  :: dpert(nrsamp) 

   real, allocatable :: xsamp(:)
   real, allocatable :: qsamp(:)
   real, allocatable :: ysamp(:)

   integer i
   real Cxx,Cyy,Cqq,Cyx,Cqy,Cqx,Cxy
   real Czz(2,2),CIzz(2,2)
   real Czy(2,1),Cyz(1,2)
   real G(1,2),Pyy(1,1)

   allocate(xsamp(nrsamp))
   allocate(qsamp(nrsamp))
   allocate(ysamp(nrsamp))

   write(*,'(a)')'++++++++++++++++++++++++++++++++++++++++++++++'
   write(*,'(a)')'ESE analysis...'
   do i=1,nrsamp
      xsamp(i)=xf(i)
      qsamp(i)=qf(i)
      ysamp(i)=func(xsamp(i),qsamp(i))
   enddo

   call cov(Cxx,Cyy,Cqq,Cyx,Cqy,Cqx,xsamp,ysamp,qsamp,nrsamp)
   if (Cqq == 0.0) then
      print *,'ese: Cqq must be larger than zero'
      return
   endif
   Czz(1,1)=Cxx;  Czz(2,2)=Cqq;  Czz(1,2)=Cqx;   Czz(2,1)=Cqx
   CIzz(1,1)=Cqq; CIzz(2,2)=Cxx; CIzz(1,2)=-Cqx; CIzz(2,1)=-Cqx
   CIzz=CIzz/(Cxx*Cqq-Cqx*Cqx)

   Cyz(1,1)=Cyx; Cyz(1,2)=Cqy
   G=matmul(Cyz,CIzz)
  ! G(1,2)=G(1,2)-1.0    ! erroneous alternative G-I_d 
   Pyy=matmul(matmul(G,Czz),transpose(G))
   Cyy=Pyy(1,1)
   
   Czy=matmul(Czz,transpose(G))
   cxy=Czy(1,1)
   cqy=Czy(2,1)

   do i=1,nrsamp
      xsamp(i)=xsamp(i) + (cyx/(Cyy+cdd))*(dpert(i)-ysamp(i))
      qsamp(i)=qsamp(i) + (cqy/(Cyy+cdd))*(dpert(i)-ysamp(i))
      ysamp(i)=func(xsamp(i),qsamp(i))
   enddo

   samples(:,1)=xsamp(:)
   samples(:,2)=qsamp(:)
   deallocate(xsamp,ysamp,qsamp)
   write(*,'(a)')'ESE analysis completed'
   write(*,'(a)')'++++++++++++++++++++++++++++++++++++++++++++++'
   write(*,'(a)')
end subroutine
end module
