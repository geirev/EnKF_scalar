module m_integrals 
use m_func
contains
subroutine integrals(x0,d,siga,sigo,beta,funcmode,maxiesit)
   implicit none
   real, intent(in) :: x0
   real, intent(in) :: d
   real, intent(in) :: siga
   real, intent(in) :: sigo
   real, intent(in) :: beta
   integer, intent(in) :: funcmode
   integer, intent(in) :: maxiesit

   integer, parameter :: n=10000
   real :: xmin=-5.0
   real :: xmax=5.0
   real :: intf=0.0
   real :: intfx1=0.0
   real :: intfx2=0.0
   real :: dummy=0.0
   integer i
   real x
   real y
   real dx
   real jxfunc
   real :: mu=0.0
   real sig
   real :: sig2=0.0
   real :: skew=0.0
   real :: kurt=0.0

   dx=(xmax-xmin)/real(n-1)

   do i=1,n-1
      x=xmin+real(i-1)*dx+0.5*dx
      y=func(x,beta,funcmode)
      jxfunc=(x-x0)**2/siga**2 + (y-d)**2/sigo**2
      intf=intf+dx*exp(-0.5*jxfunc)
      mu=mu+dx*x*exp(-0.5*jxfunc)
   enddo
   mu=mu/intf


   do i=1,n-1
      x=xmin+real(i-1)*dx+0.5*dx
      y=func(x,beta,funcmode)
      jxfunc=(x-x0)**2/siga**2 + (y-d)**2/sigo**2
      sig2=sig2+dx*(x-mu)**2*exp(-0.5*jxfunc)
      skew=skew+dx*(x-mu)**3*exp(-0.5*jxfunc)
      kurt=kurt+dx*(x-mu)**4*exp(-0.5*jxfunc)
   enddo
   sig2=sig2/intf
   sig=sqrt(sig2)

   skew=skew/intf
   skew=skew/sig**3

   kurt=kurt/intf
   kurt=kurt/sig**4-3.0

   print '(a,5f12.4)','Bayes stats: ',mu,sig2,skew,kurt

   open(10,file='truestat.dat',status='unknown')
      write(10,*)'TITLE = "True statisitcs"'
      write(10,*)'VARIABLES = "Iteration" "sumAlpha" "Alpha" "Mean" "Variance" "Skewness" "Kurtosis"'
      write(10,'(a,i5)')' ZONE T="TRUESTAT" F=POINT, I= ',2
         write(10,'(i5,6g13.4)')1,0.0,0.0,mu,sig2,skew,kurt
         write(10,'(i5,6g13.4)')2,1.0,1.0,mu,sig2,skew,kurt
   close(10)

   open(10,file="IEStruestat.dat",status='unknown')
      write(10,*)'TITLE = "True statisitcs"'
      write(10,*)'VARIABLES = "Iteration" "Mean" "Variance" "Skewness" "Kurtosis"'
      write(10,'(a,i5)')' ZONE T="TRUESTAT" F=POINT, I= ',2
         write(10,'(i5,6g13.4)')0,mu,sig2,skew,kurt
         write(10,'(i5,6g13.4)')maxiesit,mu,sig2,skew,kurt
   close(10)
end subroutine 
end module 
