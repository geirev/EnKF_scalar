module m_getcaseid
contains
subroutine getcaseid(caseid,esmethod,alphageo,nmda,esamp,gradient,beta,sigw,it)
   character(len=*), intent(in) :: esmethod
   real, intent(in)    :: alphageo
   real, intent(in)    :: beta
   real, intent(in)    :: sigw
   integer, intent(in) :: gradient
   integer, intent(in) :: nmda
   integer, intent(in) :: esamp
   integer, intent(in) :: it
   character(len=40), intent(out) :: caseid 

   character(len=3) calphageo
   character(len=3) cnmda
   character(len=3) csamp
   character(len=3) zgrad
   character(len=3) csigw
   character(len=1) cenks
   character(len=2) ln
   character(len=2) cit
 
   write(calphageo,'(f3.1)')alphageo
   write(csigw,'(f3.1)')sigw

   write(cnmda,'(i3.3)')nmda

   write(cit,'(i2.2)')it

   csamp(:)=' '
!   esamp=nint(LOG10(real(nrsamp)))
   if (esamp < 10) then
      write(csamp(1:1),'(i1.1)')esamp
   else
      write(csamp(1:2),'(i2.2)')esamp
   endif

   if (beta == 0.0) then
      ln(1:2)='_L'
      ln(1:2)=' '
   else
      ln(1:2)=' '
   endif

   if (gradient == 0) then
      zgrad(1:3)='ANA'
   else
      zgrad(1:3)='ENS'
   endif

!!!   csigw(1:3)=zgrad(1:3)

   caseid(:)=' '

   select case (trim(esmethod))
   case('MDA')
      if (it > 0) then
         caseid=esmethod(1:5)//'_'//trim(csamp)//'_'//cnmda//'_'//csigw//'_'//cit 
      else
         caseid=esmethod(1:3)//'_'//trim(csamp)//'_'//cnmda//'_'//csigw
      endif
!!      caseid=esmethod(1:3)//'_'//trim(csamp)//'_'//csigw//'_'//cnmda//'_'//calphageo
   case('IES')
         caseid=esmethod(1:3)//'_'//trim(csamp)//'_'//csigw
   case('IEnKS')
      if (it > 0) then
         caseid=esmethod(1:5)//'_'//trim(csamp)//'_'//csigw//'_'//cit 
      else
         caseid=esmethod(1:5)//'_'//trim(csamp)//'_'//csigw 
      endif
   case('STEIN')
      caseid=esmethod(1:5)//'_'//trim(csamp)//'_'//csigw 
   case('EnSTEIN')
      caseid=esmethod(1:7)//'_'//trim(csamp)//'_'//csigw 
   case('ES')
      caseid=esmethod(1:2)//'_'//trim(csamp)//'_'//csigw
   case('INI')
      caseid=esmethod(1:3)//'_'//trim(csamp)//'_'//csigw
   case('PDFJ')
      caseid=esmethod(1:4)//'_'//trim(csamp)//'_'//csigw
   case('PDFC')
      caseid=esmethod(1:4)//'_'//trim(csamp)//'_'//csigw
   case default
      print *,'Invalid esmethod in getcaseid:',trim(esmethod)
      stop
   end select
   print *,'caseid:',trim(caseid),'+++'
end subroutine
end module
