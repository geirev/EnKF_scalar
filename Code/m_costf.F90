module m_costf
use m_func
contains

real function jxfunc(x,x0,siga,d,sigo,beta,funcmode)
   implicit none
   real, intent(in) :: x
   real, intent(in) :: x0
   real, intent(in) :: siga
   real, intent(in) :: d
   real, intent(in) :: sigo
   real, intent(in) :: beta
   integer, intent(in) :: funcmode
  
   jxfunc=(x-x0)**2/siga**2 + (func(x,beta,funcmode)-d)**2/sigo**2
   
end function jxfunc

real function djxfunc(x,x0,siga,d,sigo,beta,funcmode)
   implicit none
   real, intent(in) :: x
   real, intent(in) :: x0
   real, intent(in) :: siga
   real, intent(in) :: d
   real, intent(in) :: sigo
   real, intent(in) :: beta
   integer, intent(in) :: funcmode

   djxfunc=sigo**2*(x-x0) + siga**2*(func(x,beta,funcmode)-d)*dfunc(x,beta,funcmode)

end function djxfunc

real function ddjxfunc(x,x0,siga,d,sigo,beta,funcmode)
   implicit none
   real, intent(in) :: x
   real, intent(in) :: x0
   real, intent(in) :: siga
   real, intent(in) :: d
   real, intent(in) :: sigo
   real, intent(in) :: beta
   integer, intent(in) :: funcmode

   ddjxfunc=sigo**2 + siga**2*dfunc(x,beta,funcmode)**2

end function ddjxfunc

end module m_costf
