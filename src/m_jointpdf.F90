module m_jointpdf
logical, save :: lactivepdf
contains
   subroutine jointpdf(pdf,x,y,nx,ny,xsamp,ysamp,nrsamp,xa,ya,dx,dy,caseid)
   implicit none
   integer,          intent(in) :: nrsamp
   integer,          intent(in) :: nx
   integer,          intent(in) :: ny
   real,             intent(in) :: x(nx)
   real,             intent(in) :: y(ny)
   real,             intent(in) :: xsamp(nx)
   real,             intent(in) :: ysamp(ny)
   real,             intent(in) :: xa
   real,             intent(in) :: ya
   real,             intent(in) :: dx
   real,             intent(in) :: dy
   character(len=40),intent(in) :: caseid

   real,             intent(out) :: pdf(nx,ny)

   integer i,ival,jval
   real sump
   if (.not.lactivepdf) return
! Joint pdf (only used for plotting)
   pdf=0.0
   do i=1,nrsamp
      ival= nint((xsamp(i)-xa)/dx) + 1
      jval= nint((ysamp(i)-ya)/dy) + 1
      if (0 < ival .and. ival < nx+1 .and. 0 < jval .and. jval < ny+1) then
         pdf(ival,jval)=pdf(ival,jval)+1.0
      else
!         print *,'jointpdf: not in interval',ival,xsamp(i),jval,ysamp(i)
      endif
   enddo
   sump=sum(pdf(:,:))*dx*dy
   pdf=pdf/sump
   
!   call tecjointpdf(pdf,x,y,nx,ny,caseid)
end subroutine
end module
