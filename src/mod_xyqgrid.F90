module mod_xyqgrid
use mod_inistat
use m_func
implicit none
integer, parameter :: nx=100
integer, parameter :: nq=50
integer, parameter :: ny=50

real, save :: x(nx), y(ny), q(nq)
real, save :: margx(nx), margy(ny), margq(nq)
real, save :: datum(ny), priorx(nx), priorq(nq), pdf(nx,ny), cost(nx)

integer, parameter :: nxx=500
integer, parameter :: nyy=500
integer, parameter :: nqq=500
real, save :: xx(nxx), yy(nyy), qq(nxx)
real, save :: dxx
real, save :: dyy
real, save :: dqq

real, allocatable :: array(:,:)

real, save :: xa,xb,ya,yb,qa,qb 
real, save :: dx,dy,dq
contains
subroutine xyqgrid()
   integer i,j
   real yfield(nx,nq)

   if (xa==xb) then
      xa=x0 - 5.0*siga
      xb=x0 + 5.0*siga
      print '(a,2f10.4)','xa and xb set to :',xa,xb
   endif
    
   dx=(xb-xa)/real(nx-1); print '(a,f12.5)','dx=',dx
    
   do i=1,nx
      x(i)=xa + real(i-1)*dx
   enddo

   if (qa==qb) then
      qa=-5.0*max(sigq,sigw)
      qb= 5.0*max(sigq,sigw)
      print '(a,2f10.4)','qa and qb set to :',qa,qb
   endif
    
   dq=(qb-qa)/real(nq-1); print '(a,f12.5)','dq=',dq
    
   do i=1,nq
      q(i)=qa + real(i-1)*dq
   enddo

   if (ya==yb) then
      do j=1,nq
      do i=1,nx
         yfield(i,j)=func(x(i),q(j))
      enddo
      enddo
      ya=minval(yfield)
      yb=maxval(yfield)
      print '(a,2f10.4)','ya and yb set to :',ya,yb
   endif
    
   dy=(yb-ya)/real(ny-1); print '(a,f12.5)','dy=',dy
    
   do j=1,ny
      y(j)=ya + real(j-1)*dy
   enddo


! Fine grid for analytic pdfs
   dxx=(xb-xa)/real(nxx-1); print '(a,f12.5)','dxx=',dxx
    
   do i=1,nxx
      xx(i)=xa + real(i-1)*dxx
   enddo
    
   dqq=(qb-qa)/real(nqq-1); print '(a,f12.5)','dqq=',dqq
    
   do i=1,nqq
      qq(i)=qa + real(i-1)*dqq
   enddo
    
   dyy=(yb-ya)/real(nyy-1); print '(a,f12.5)','dyy=',dyy
    
   do j=1,nyy
      yy(j)=ya + real(j-1)*dyy
   enddo

end subroutine
end module
