module m_printcostf
use mod_inistat
use mod_xyqgrid
use m_teccostens
use m_tecsampini
use m_func
implicit none
integer, parameter :: nrits=99

contains
subroutine printcostf(xsampini,xsamp,dpert,nrsamp)
! Printing the first member's cost functions 
   integer, intent(in) :: nrsamp
   real, intent(in) :: xsampini(nrsamp)
   real, intent(in) :: xsamp(nrsamp)
   real, intent(in) :: dpert(nrsamp)

   real, allocatable :: costens(:,:)
   real, allocatable :: costite(:,:)
   real, allocatable :: xsampit(:,:)
   integer i,n
   if (sigw > 0.0) return

   allocate(costens(nx,10))
   do n=1,10
      do i=1,nx
         costens(i,n)=(x(i)-xsampini(n))**2/siga**2 + (func(x(i),0.0)-dpert(n))**2/cdd    
      enddo
   enddo
   call teccostens('costens',costens,x,nx,10,'x')

   allocate(costite(10,nrits))
   allocate(xsampit(10,nrits))
   do n=1,10
      xsampit(n,1)=xsampini(n)
      costite(n,1)= (func(xsampit(n,1),0.0)-dpert(n))**2/cdd    
      xsampit(n,2)=xsamp(n)
      costite(n,2)= (xsampit(n,2)-xsampini(n))**2/siga**2 + (func(xsampit(n,2),0.0)-dpert(n))**2/cdd    
   enddo
   call tecsampini('sampiniES',costite,xsamp,2,10)
   deallocate(costens)
   deallocate(costite)
   deallocate(xsampit)
end subroutine
end module
