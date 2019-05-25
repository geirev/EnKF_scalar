module m_tecsampini
contains
subroutine tecsampini(fname,cost,x,nit,nr)
   implicit none
   character(len=*), intent(in) :: fname
   integer,          intent(in) :: nit
   integer,          intent(in) :: nr
   real,             intent(in) :: x(nr,nit)
   real,             intent(in) :: cost(nr,nit)
   integer i,j

   open(10,file=trim(fname)//'.dat',status='unknown')
      write(10,*)'TITLE = "'//trim(fname)//'"'
      write(10,'(a)')'VARIABLES = "x" "cost"'

      do i =1,nr
         write(10,*)'ZONE T= "Cost for iteration x" F=POINT, I=',nit
         do j=1,nit
            write(10,'(2f13.5)')x(i,j),cost(i,j)
         enddo
      enddo
   close(10)
end subroutine
end module
