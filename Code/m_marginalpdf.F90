module m_marginalpdf
contains
subroutine marginalpdf(pdf,margx,margy,nx,ny,x,y,dx,dy)
implicit none
integer, intent(in) :: nx
integer, intent(in) :: ny
real, intent(in)    :: pdf(nx,ny)
real, intent(in)    :: x(nx)
real, intent(in)    :: y(ny)
real, intent(out)   :: margx(nx)
real, intent(out)   :: margy(ny)
real, intent(in)    :: dx
real, intent(in)    :: dy

real sump
integer loc1(1)
integer loc2(2)
integer i,j

! Mode of joint pdf
!   loc2=maxloc(pdf)
!   mode(1)=x(loc2(1))+0.5*dx
!   mode(2)=y(loc2(2))+0.5*dy

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


! Mean of pdf
!   mean(1:2)=0.0
!   do i=1,nx
!      mean(1)=mean(1)+x(i)*margx(i)*dx
!   enddo
!   do j=1,ny
!      mean(2)=mean(2)+y(j)*margy(j)*dy
!   enddo

end subroutine
end module
