module m_tecjointpdf
contains
   subroutine tecjointpdf(pdf,x,y,nx,ny,caseid)
   implicit none
   integer,          intent(in) :: nx
   integer,          intent(in) :: ny
   real,             intent(in) :: pdf(nx,ny)
   real,             intent(in) :: x(nx)
   real,             intent(in) :: y(ny)
   character(len=*),intent(in) :: caseid

   integer i,j

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
