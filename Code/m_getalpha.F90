module m_getalpha
contains
real function getalpha(n,nmda,alphageo)
   implicit none
   integer, intent(in) :: n
   integer, intent(in) :: nmda
   real,    intent(in) :: alphageo
   real, allocatable :: alpha(:)
   real alphasum
   integer i

   allocate(alpha(nmda))

   alpha(1)=1000.0
   alphasum=1.0/alpha(1)
   do i=2,nmda
      alpha(i)=alpha(i-1)/alphageo
      alphasum=alphasum+1.0/alpha(i)
   enddo
   alpha(:)=alphasum*alpha(:) 
!   write(*,*)'Alpha: ',alpha(1:nmda)
   getalpha=alpha(n)

   deallocate(alpha)
end function getalpha
end module m_getalpha
