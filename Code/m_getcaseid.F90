module m_getcaseid
contains
subroutine getcaseid(caseid,esmethod,alphageo,nmda,esamp,sigw,it)
   character(len=7), intent(in) :: esmethod
   real, intent(in)    :: alphageo
   real, intent(in)    :: sigw
   integer, intent(in) :: nmda
   integer, intent(in) :: esamp
   integer, intent(in) :: it
   character(len=40), intent(out) :: caseid 

   character(len=3) calphageo
   character(len=3) cnmda
   character(len=3) csamp
   character(len=3) csigw
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

   caseid(:)=' '
   caseid=trim(esmethod)//'_'//trim(csamp)//'_'//csigw
   if (nmda > 0) caseid=trim(caseid)//'_'//cnmda
   if (it > 0) caseid=trim(caseid)//'_'//cit 

end subroutine
end module
