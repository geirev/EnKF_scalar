   module m_tecstat
   use m_getcaseid
   implicit none
   contains
   subroutine tecstat(fname,esmethod,nmda,alpha,ave,var,skew,curt,nrsamp,alphageo,beta,sigw)
   implicit none
   character(len=*), intent(in) :: fname
   character(len=*), intent(in) :: esmethod
   integer,          intent(in) :: nrsamp
   integer,          intent(in) :: nmda
   real,             intent(in) :: alphageo
   real,             intent(in) :: beta
   real,             intent(in) :: sigw
   real,             intent(in) :: alpha(nmda)
   real,             intent(in) :: ave(0:nmda)
   real,             intent(in) :: var(0:nmda)
   real,             intent(in) :: skew(0:nmda)
   real,             intent(in) :: curt(0:nmda)
   integer i
   real salpha
   character(len=80) totfname
   character(len=40) caseid


   i=0
   call getcaseid(caseid,esmethod,alphageo,nmda,nrsamp,sigw,i)
   totfname(:)=' '
   totfname=trim(fname)//trim(caseid)//'.dat'
   print *,'totfnam:',trim(totfname)

   open(10,file=trim(totfname),status='unknown')
      write(10,*)'TITLE = "MDA statistics"'
      write(10,*)'VARIABLES = "Iteration" "sumAlpha" "Alpha" "Mean" "Variance" "Skewness" "Kurtosis"'
      write(10,'(a,i5)')' ZONE T="'//trim(caseid)//'" F=POINT, I= ',nmda+1
      salpha=0.0
      do i=0,nmda
         if (i>0) salpha=salpha+1.0/alpha(i)
         write(10,'(i5,6g13.4)')i,salpha,alpha(max(i,1)),ave(i),var(i),skew(i),curt(i)
      enddo
   close(10)

   end subroutine
   end module
