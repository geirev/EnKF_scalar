module m_avevar
contains
!************************************************************
!* Given an array of Data of length n, this routine returns * 
!* its mean ave, average deviation adev, standard deviation *
!* sdev, variance var, skewness skew, and kurtosis curt.    *
!************************************************************ 
subroutine avevar(datavector,nrsamp,ave,var)
   implicit none
   integer, intent(in)  :: nrsamp
   real,    intent(in)  :: datavector(nrsamp)
   real,    intent(out) :: ave,var
   integer j
   real p,s

   if(nrsamp.le.1) then
      print *,' Nrsamp must be at least 2!'
      stop
   endif

   s=sum(datavector(1:nrsamp))
   ave=s/real(nrsamp)
   var=0.0
   do j=1,nrsamp
      s=datavector(j)-ave
      p=s*s
      var=var+p
   enddo
   var=var/(real(nrsamp-1))
end
end module m_avevar
