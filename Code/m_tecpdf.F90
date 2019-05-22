module m_tecpdf
logical, save :: lactivepdf
contains
   subroutine tecpdf(x,y,nx,ny,xsamp,ysamp,nrsamp,xa,ya,dx,dy,caseid)
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

   real pdf(nx,ny)

   integer m,i,j,ival,jval
   logical lopen
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
   
   open(10,file='pdf'//trim(caseid)//'.dat',status='unknown')
      write(10,*)'TITLE = "PDF for '//trim(caseid)//'"'
      write(10,*)'VARIABLES = "i-index" "j-index" "x" "y" "pdf"'
      write(10,'(a,i5,a,i5,a)')' ZONE T="'//trim(caseid)//'"  F=BLOCK, I=',nx,', J=',ny,', K=1'
      write(10,'(30I4)')((i,i=1,nx),j=1,ny)
      write(10,'(30I4)')((j,i=1,nx),j=1,ny)
      write(10,900) ((x(i),i=1,nx),j=1,ny)
      write(10,900) ((y(j),i=1,nx),j=1,ny)
      write(10,900) ((pdf(i,j),i=1,nx),j=1,ny)
   close(10)

   900 format(10(1x,e16.9))
end subroutine
end module
