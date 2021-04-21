module m_tecmargpdf
contains
subroutine tecmargpdf(id,xsamp,nrsamp,caseid,xmin,xmax,nx,it)
   use m_getcaseid
   use mod_shapiro
   implicit none
   integer, intent(in)  :: nrsamp
   real,    intent(in)  :: xsamp(nrsamp)
   integer, intent(in)  :: nx
   real,    intent(in)  :: xmin
   real,    intent(in)  :: xmax
   integer, optional    :: it
   character(len=*),   intent(in)  :: id
   character(len=*), intent(in) :: caseid


   character(len=80) :: fname
   real                 :: pdf(nx)
   real                 :: pdfsh(nx)
   real xxmin,xxmax,dx,sumpdf
   integer i,ival
   character(len=3) tag3
   integer, parameter :: nshapiro=4
   real sh(0:nshapiro)

   if (xmax == xmin) then
      xxmin=minval(xsamp)
      xxmax=maxval(xsamp)
   else
      xxmin=xmin
      xxmax=xmax
   endif
   if (xxmax == xxmin) return
   
   dx=(xxmax-xxmin)/real(nx-1)

   pdf(:)=0.0
   do i =1,nrsamp
      ival= nint((xsamp(i)-xxmin)/dx) + 1
      if (0 < ival .and. ival < nx+1) then
         pdf(ival)=pdf(ival)+1.0
      endif 
   enddo
   call shfact(nshapiro,sh)
   call shfilt(nshapiro,sh,nx,pdf,1,pdfsh,1,nshapiro)
   sumpdf=sum(pdf(:))*dx
   pdf=pdf/sumpdf
   sumpdf=sum(pdfsh(:))*dx
   pdfsh=pdfsh/sumpdf



   fname(:)=' '
   if (present(it)) then
      write(tag3,'(i3.3)')it
      fname='marg'//trim(id)//trim(caseid)//'_'//tag3//'.dat'
   else
      fname='marg'//trim(id)//trim(caseid)//'.dat'
   endif
   open(10,file=trim(fname),status='unknown')
      write(10,*)'TITLE = "marginal PDFs"'
      write(10,*)'VARIABLES = "'//trim(id)//'" "Marg pdf" "Marg pdfsh"'
      if (present(it)) then
      write(10,*)'ZONE T= "',trim(caseid)//'_'//tag3,'" F=POINT, I=',nx
      else
      write(10,*)'ZONE T= "',trim(caseid)//'" F=POINT, I=',nx
      endif
      do i =1,nx
         write(10,'(3f15.7)')xxmin+real(i-1)*dx,pdf(i),pdfsh(i)
      enddo
   close(10)
    
end subroutine tecmargpdf
end module m_tecmargpdf

