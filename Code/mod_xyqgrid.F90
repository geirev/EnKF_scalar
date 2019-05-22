module mod_xyqgrid
   implicit none
   integer, parameter :: nx=50
   integer, parameter :: nq=50
   integer, parameter :: ny=50

   real, save :: x(nx), y(ny), q(nq)
   real, save :: margx(nx), margy(ny), margq(nq)
   real, save :: datum(ny), prior(nx), pdf(nx,ny), cost(nx)
   
   real, save :: xa,xb,ya,yb,qa,qb 
   real, save :: dx,dy,dq
   contains
   subroutine xyqgrid(x,q,y)
      real,    intent(out) :: x(nx),q(nq),y(ny)
      integer i,j


      dx=(xb-xa)/real(nx-1); print '(a,f12.5)','dx=',dx
      dy=(yb-ya)/real(ny-1); print '(a,f12.5)','dy=',dy
      dq=(qb-qa)/real(nx-1); print '(a,f12.5)','dq=',dq

      do i=1,nx
         x(i)=xa + real(i-1)*dx
      enddo

      do j=1,ny
         y(j)=ya + real(j-1)*dy
      enddo

      do i=1,nq
         q(i)=qa + real(i-1)*dq
      enddo
   end subroutine
end module
