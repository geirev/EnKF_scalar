   module m_tecstatIES
   use m_getcaseid
   contains
   subroutine tecstatIES(fname,esmethod,ave,var,skew,curt,nrsamp,iterations,beta,sigw)
   implicit none
   character(len=*), intent(in) :: fname
   character(len=*), intent(in) :: esmethod
   integer,          intent(in) :: nrsamp
   integer,          intent(in) :: iterations
   real,             intent(in) :: ave(0:iterations)
   real,             intent(in) :: var(0:iterations)
   real,             intent(in) :: skew(0:iterations)
   real,             intent(in) :: curt(0:iterations)
   real,             intent(in) :: beta
   real,             intent(in) :: sigw
   integer i
   character(len=80) totfname
   character(len=40) caseid
   i=0
   call getcaseid(caseid,esmethod,1.0,1,nrsamp,sigw,i)

   totfname(:)=' '
   totfname=trim(fname)//trim(caseid)//'.dat'
   print *,'totfnam:',trim(totfname)

   open(10,file=trim(totfname),status='unknown')
      write(10,*)'TITLE = "IES statisitcs"'
      write(10,*)'VARIABLES = "Iteration" "Mean" "Variance" "Skewness" "Kurtosis"'
      write(10,'(a,i5)')' ZONE T="'//trim(caseid)//'" F=POINT, I= ',iterations+1
      do i=0,iterations
         write(10,'(i5,6g13.4)')i,ave(i),var(i),skew(i),curt(i)
      enddo
   close(10)

   end subroutine
   end module
