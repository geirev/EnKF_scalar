module m_pseudoinv2
contains
subroutine pseudoinv2(W,n,nrsamp,truncation)
   implicit none
   integer, intent(in)     :: n
   integer, intent(in)     :: nrsamp
   real,    intent(inout)  :: W(n,nrsamp)
   real,    intent(in)     :: truncation
   
   real, allocatable :: sig0(:)
   real, allocatable :: U0(:,:)  
   real, allocatable :: VT0(:,:)  
   integer ierr
   integer lwork
   integer nrmin
   real, allocatable :: work(:)
   integer nrsigma,i,j

   real sigsum,sigsum1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   nrmin=min(n,nrsamp)
   allocate(work(20*nrsamp))
   allocate(sig0(nrsamp))
   allocate(U0(nrsamp,nrsamp))
   allocate(VT0(nrsamp,nrsamp))

   lwork=20*nrsamp
   sig0=0.0
!   print *,'Calling dgesvd'
   call dgesvd('A', 'A', nrsamp, nrsamp, W, nrsamp, sig0, U0, nrsamp, VT0, nrsamp, work, lwork, ierr)
!   print *,'dgesvd completed'
   if (ierr /= 0) then
      print *,'svdW: ierr from call dgesvd 0= ',ierr; stop
   endif

   print *,'Checking Sig '
   write(*,'(10g11.3)')sig0(1:min(nrsamp,10))

!   print *,'dumping singular values'
   open(10,file='singularvalues.dat')
      do i=1,nrsamp
         write(10,'(i6,g13.5)')i,sig0(i)
      enddo
   close(10)

! Computing significant number of eigenvalues.
   sigsum=0.0
   do i=1,nrsamp
      sigsum=sigsum+sig0(i)**2
   enddo
   sigsum1=0.0
   nrsigma=0
   do i=1,nrsamp
      if (sigsum1/sigsum < truncation) then
         nrsigma=nrsigma+1
         sigsum1=sigsum1+sig0(i)**2
      else
         sig0(i:nrsamp)=0.0
         exit
      endif
   enddo
   print *,'significant number of eigenvalues:',nrsigma,sigsum1/sigsum

! pseudo inverse
   do i=1,nrsigma
       sig0(i) = 1.0/sig0(i)
   enddo


  print *,'Checking significant Sig '
  write(*,'(10g11.3)')sig0(1:min(nrsamp,10))

  do j=1,nrsamp
  do i=1,nrsamp
     VT0(j,i)=VT0(j,i)*sig0(j)
  enddo
  enddo
  call dgemm('T','T',nrsamp,nrsamp,nrsamp,1.0,VT0,nrsamp,U0,nrsamp,0.0,W,nrsamp)
   
  deallocate(work,sig0,U0,VT0)

end subroutine
end module
