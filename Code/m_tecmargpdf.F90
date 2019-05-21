module m_tecmargpdf
contains
subroutine tecmargpdf(id,qsamp,nrsamp,caseid,xmin,xmax,nx,it)
   use m_getcaseid
   implicit none
   integer, intent(in)  :: nrsamp
   real,    intent(in)  :: qsamp(nrsamp)
   integer, intent(in)  :: nx
   real,    intent(in)  :: xmin
   real,    intent(in)  :: xmax
   integer, optional    :: it
   character(len=*),   intent(in)  :: id
   character(len=*), intent(in) :: caseid


   character(len=80) :: fname
   real                 :: pdf(nx)
   real qmin,qmax,dx,ave,var,sumpdf
   integer i,ival
   character(len=3) tag3

   qmin=minval(qsamp); qmin=max(qmin,xmin)
   qmax=maxval(qsamp); qmax=min(qmax,xmax)
   qmin=xmin
   qmax=xmax
!   print *,'margq:',qmin,qmax,trim(id),trim(caseid)
   if (qmax == qmin) then
!      print *,'returns'
      return
   endif
   
   dx=(qmax-qmin)/real(nx-1)
   ave=0.0
   var=0.0
   do i =1,nrsamp
      ave=ave+qsamp(i)
      var=var+qsamp(i)**2
   enddo
   ave=ave/real(nrsamp)
   var=var/real(nrsamp) - ave**2
!   print '(a,3f10.3)','Statistics of marginal pdf for q (ave, var, std) :',ave,var,sqrt(var)

   pdf(:)=0.0
   do i =1,nrsamp
      ival= nint((qsamp(i)-qmin)/dx) + 1
      if (0 < ival .and. ival < nx+1) then
         pdf(ival)=pdf(ival)+1.0
      endif 
   enddo
   sumpdf=sum(pdf(:))*dx
   pdf=pdf/sumpdf



   fname(:)=' '
   if (present(it)) then
      write(tag3,'(i3.3)')it
      fname='marg'//trim(id)//trim(caseid)//'_'//tag3//'.dat'
   else
      fname='marg'//trim(id)//trim(caseid)//'.dat'
   endif
   open(10,file=trim(fname),status='unknown')
      write(10,*)'TITLE = "marginal PDFs"'
      write(10,*)'VARIABLES = "'//trim(id)//'" "Marg pdf"'
      if (present(it)) then
      write(10,*)'ZONE T= "',trim(caseid)//'_'//tag3,'" F=POINT, I=',nx
      else
      write(10,*)'ZONE T= "',trim(caseid)//'" F=POINT, I=',nx
      endif
      do i =1,nx
         write(10,'(2f15.7)')qmin+real(i-1)*dx,pdf(i)
      enddo
   close(10)
    
end subroutine tecmargpdf
end module m_tecmargpdf

