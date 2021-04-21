module m_teccostens
contains
subroutine teccostens(fname,y,x,nx,nr,cxvar)
   implicit none
   character(len=*), intent(in) :: fname
   character(len=*), intent(in) :: cxvar
   integer,          intent(in) :: nx
   integer,          intent(in) :: nr
   real,             intent(in) :: x(nx)
   real,             intent(in) :: y(nx,nr)
   integer i
   character(len=2) tag2

   open(10,file=trim(fname)//'.dat',status='unknown')
      write(10,*)'TITLE = "'//trim(fname)//'"'
      write(10,'(a,a,a)',advance='no')'VARIABLES = "'//trim(cxvar)//'" '
      do i=1,min(99,nr)
         write(tag2,'(i2.2)')i
         write(10,'(a,a)',advance='no')'"f'//tag2//'" '
      enddo
      write(10,*)' '
      write(10,*)'ZONE T= "Ensemble cost functions" F=POINT, I=',nx
      do i =1,nx
         write(10,'(100f13.5)')x(i),y(i,1:min(99,nr))
      enddo
   close(10)
end subroutine
end module
