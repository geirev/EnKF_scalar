module m_marginalpdf
contains
subroutine marginalpdf(pdf,margx,margy,nx,ny,dx,dy)
implicit none
integer, intent(in) :: nx
integer, intent(in) :: ny
real, intent(in)    :: pdf(nx,ny)
real, intent(out)   :: margx(nx)
real, intent(out)   :: margy(ny)
real, intent(in)    :: dx
real, intent(in)    :: dy

real sump
integer i,j


! Marginal pdf for X
   do i=1,nx
      margx(i)=sum(pdf(i,1:ny))
   enddo
   sump=sum(margx(:))*dx
   margx=margx/sump

! Marginal pdf for Y
   do j=1,ny
      margy(j)=sum(pdf(1:nx,j))
   enddo
   sump=sum(margy(:))*dy
   margy=margy/sump

end subroutine
end module
