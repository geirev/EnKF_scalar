module m_pseudoinv
contains
subroutine pseudoinv(W,Wi,n,nrsamp,truncation)
   implicit none
   integer, intent(in)     :: n
   integer, intent(in)     :: nrsamp
   real,    intent(in)     :: W(n,nrsamp)
   real,    intent(out)    :: Wi(nrsamp,n)
   real,    intent(in)     :: truncation
   
   real, allocatable :: work(:)
   real, allocatable :: sig0(:)
   real, allocatable :: U0(:,:)  
   real, allocatable :: VT0(:,:)  
   real, allocatable :: VS(:,:)  
   integer ierr
   integer lwork
   integer nrmin
   integer nrsigma,i,j
   real :: V(nrsamp,nrsamp)

   real sigsum,sigsum1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   nrmin=min(n,nrsamp) 
!  print '(a,3I10)','   dgesvd: nrmin,n,nrsamp:',nrmin,n,nrsamp
   lwork=2*MAX(3*MIN(n,nrsamp) + MAX(n,nrsamp),5*MIN(n,nrsamp))
   allocate(work(lwork))
   allocate(sig0(nrmin))
   allocate(U0(n,nrmin))
   allocate(VT0(nrmin,nrsamp))

   sig0=0.0
!   print *,'Calling dgesvd'
   call dgesvd('S', 'S', n, nrsamp, W, n, sig0, U0, n, VT0, nrmin, work, lwork, ierr)
!   print *,'dgesvd completed'
   if (ierr /= 0) then
      print *,'svdW: ierr from call dgesvd 0= ',ierr; stop
   endif

!   write(*,'(a,10g11.3)')'   dgesvd: sigma= ',sig0(1:min(nrmin,10))

   open(10,file='singularvalues.dat')
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
      else
         sig0(i:nrmin)=0.0
         exit
      endif
   enddo
!   print '(a,i10,g13.5)','   dgesvd: significant number of eigenvalues= ',nrsigma,sigsum1/sigsum

! pseudo inverse
   do i=1,nrsigma
       sig0(i) = 1.0/sig0(i)
   enddo

!   write(*,'(a,10g11.3)')'   dgesvd: sigma inv= ',sig0(1:min(nrmin,10))

   do i=1,nrsamp
   do j=1,nrmin
      VT0(j,i)=VT0(j,i)*sig0(j)
   enddo
   enddo

  call dgemm('T','T',nrsamp,n,nrmin,1.0,VT0,nrmin,U0,n,0.0,Wi,nrsamp)

  deallocate(work)
  deallocate(sig0)
  deallocate(U0)
  deallocate(VT0)

end subroutine
end module
