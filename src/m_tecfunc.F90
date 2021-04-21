module m_tecfunc
contains
subroutine tecfunc(fname,y,x,nx,cxvar,czone)
   implicit none
   character(len=*), intent(in) :: fname
   character(len=*), intent(in) :: cxvar
   character(len=*), intent(in) :: czone
   integer,          intent(in) :: nx
   real,             intent(in) :: x(nx)
   real,             intent(in) :: y(nx)
   integer i

   open(10,file=trim(fname)//'.dat',status='unknown')
      write(10,*)'TITLE = "'//trim(fname)//'"'
      write(10,*)'VARIABLES = "'//trim(cxvar)//'" "Marg pdf" "Marg pdfsh"'
      write(10,*)'ZONE T= "'//trim(czone)//'" F=POINT, I=',nx
      do i =1,nx
         write(10,'(3f13.5)')x(i),y(i),y(i)
      enddo
   close(10)
end subroutine
end module
