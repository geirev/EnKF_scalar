module m_omegafact
contains
subroutine omegafact(nrsamp,Omega)
integer, intent(in) :: nrsamp
real,    intent(in) :: Omega(nrsamp,nrsamp)
real VL
real VR(nrsamp,nrsamp)
real WR(nrsamp),WI(nrsamp)
real work(8*nrsamp)
integer info,lwork

   lwork=8*nrsamp
   
   print '(a,i3)','Eigenvalue decomposition of Omega'
   print *,'Omega'
   print '(10f12.6)',Omega(1:10,1:10)

   call dgeev('N','V',nrsamp,Omega,nrsamp,WR,WI,VL,1,VR,nrsamp,work,lwork,info)
   print '(a,i3)','Info=',info

   print *,'eigen values'
   do i=1,nrsamp
      print '(i5,2g13.6)',i,WR(i),WI(i)
   enddo

   print *,'Eigen vectors'
   print '(10f12.6)',VR(1:10,1:10)

end subroutine
end module
