module m_cyyreg
implicit none
logical, save :: lcyyreg
contains
real function cyyreg(Cxx,Cqq,Cyx,Cqy,Cqx)
   implicit none
   real, intent(in) :: Cxx,Cqq,Cyx,Cqy,Cqx
   real Pzz(2,2),PIzz(2,2),Pzy(2,1),Pyz(1,2),Pyymat(1,1),G

   if (Cqq == 0.0) then
         G=Cyx/Cxx
         Cyyreg=G*Cxx*G
   else
      Pzz(1,1)=Cxx; Pzz(2,2)=Cqq; Pzz(1,2)=Cqx; Pzz(2,1)=Cqx
      PIzz(1,1)=Cqq; PIzz(2,2)=Cxx; PIzz(1,2)=-Cqx; PIzz(2,1)=-Cqx
      PIzz=PIzz/(Cxx*Cqq-Cqx*Cqx)

      Pzy(1,1)=Cyx; Pzy(2,1)=Cqy
      Pyz(1,1)=Cyx; Pyz(1,2)=Cqy

      Pyymat=matmul(Pyz,matmul(PIzz,Pzy))
      Cyyreg=Pyymat(1,1)
   endif
end function
end module
