module m_normal
contains
real function normal()
!  Returns a random value N(0,1)
   implicit none
   real work1,work2
   real, parameter   ::  pi=3.141592653589
   integer i
   do i=1,100
      call random_number(work1)
      call random_number(work2)
      normal= sqrt(-2.0*log(work1))*cos(2.0*pi*work2)
      if (abs(normal) < 4.2) exit
   enddo

end function normal
end module m_normal
