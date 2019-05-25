module m_aaprojection
contains
subroutine aaprojection(A,AA,ndim,nrsamp,truncation)
   implicit none
   integer, intent(in)     :: ndim
   integer, intent(in)     :: nrsamp
   real,    intent(in)     :: A(ndim,nrsamp)
   real,    intent(out)    :: AA(nrsamp,nrsamp)
   real,    intent(in)     :: truncation
   
   real, allocatable :: work(:)
   real, allocatable :: sig0(:)
   real, allocatable :: VT0(:,:)  
   real, allocatable :: V0(:,:)  
   real  U0(1,1)  
   integer ierr
   integer lwork
   integer nrmin
   integer nrsigma,i

   real sigsum,sigsum1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   nrmin=min(ndim,nrsamp) 
!  print '(a,3I10)','   dgesvd: nrmin,ndim,nrsamp:',nrmin,ndim,nrsamp
   lwork=2*MAX(3*MIN(ndim,nrsamp) + MAX(ndim,nrsamp),5*MIN(ndim,nrsamp))
   allocate(work(lwork))
   allocate(sig0(nrmin))
!   allocate(U0(ndim,nrmin))
   allocate(VT0(nrmin,nrsamp))
   allocate(V0(nrsamp,nrmin))


   sig0=0.0
!   print *,'Calling dgesvd'
   call dgesvd('N', 'S', ndim, nrsamp, A, ndim, sig0, U0, 1, VT0, nrmin, work, lwork, ierr)
!   print *,'dgesvd completed'
   if (ierr /= 0) then
      print *,'svdW: ierr from call dgesvd 0= ',ierr; stop
   endif

!   write(*,'(a,10g11.3)')'   dgesvd: sigma= ',sig0(1:min(nrmin,10))

   open(10,file='singularvaluesaa.dat')
      do i=1,nrmin
         write(10,'(i6,g13.5)')i,sig0(i)
      enddo
   close(10)

! Computing significant number of eigenvalues.
   sigsum=0.0
   do i=1,nrmin
      sigsum=sigsum+sig0(i)**2
   enddo
   sigsum1=0.0
   nrsigma=0
   do i=1,nrmin
      if (sigsum1/sigsum < truncation) then
         nrsigma=nrsigma+1
         sigsum1=sigsum1+sig0(i)**2
      endif
   enddo
!   print '(a,i10)','nrsigma:',nrsigma
!   print '(a,i10)','nrmin  :',nrmin
!   print '(a,i10)','ndim   :',ndim
!   print '(a,i10)','nrsamp :',nrsamp

!   write(*,'(a,g12.4)')'Checking VT: 1/sqrt(nrsamp)=',1.0/sqrt(real(nrsamp))
!   write(*,'(10g12.4)')VT0(1:10,1:10)

   if (nrsigma < nrmin) VT0(nrsigma+1:nrmin,:)=0.0
   V0=transpose(VT0)
   call dgemm('N','N',nrsamp,nrsamp,nrmin,1.0,V0,nrsamp,VT0,nrmin,0.0,AA,nrsamp)
  deallocate(work)
  deallocate(sig0)
  deallocate(VT0)
  deallocate(V0)

end subroutine
end module
