module m_moments
contains
!************************************************************
!* Given an array of Data of length n, this routine returns * 
!* its mean ave, average deviation adev, standard deviation *
!* sdev, variance var, skewness skew, and kurtosis curt.    *
!************************************************************ 
subroutine moments(datavector,nrsamp,ave,adev,sdev,var,skew,curt)
   implicit none
   integer, intent(in)  :: nrsamp
   real,    intent(in)  :: datavector(nrsamp)
   real,    intent(out) :: ave,adev,sdev,var,skew,curt
   integer j
   real p,s

   if(nrsamp.le.1) then
      print *,' Nrsamp must be at least 2!'
      stop
   endif
   s=sum(datavector(1:nrsamp))
   ave=s/real(nrsamp)
   adev=0.0
   var=0.0
   skew=0.0
   curt=0.0
   do j=1,nrsamp
      s=datavector(j)-ave
      adev=adev+abs(s)
      p=s*s
      var=var+p
      p=p*s
      skew=skew+p
      p=p*s
      curt=curt+p
   enddo
   adev=adev/real(nrsamp)
   var=var/(real(nrsamp-1))
   sdev=sqrt(var)
   if(var.ne.0.0) then
      skew=skew/(real(nrsamp)*sdev**3)
      curt=curt/(real(nrsamp)*var**2)-3.0
   else
      print *,' No skew or kurtosis when zero variance.'
   end if
return
end subroutine
end module
