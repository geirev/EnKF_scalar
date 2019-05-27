module m_enstein
use mod_inistat
use mod_xyqgrid
use m_cyyreg
use m_cov
use m_aaprojection
use m_getcaseid
use m_steinkern
use m_func
use m_normal
use m_tecmargpdf
use m_tecsampini
use m_tecpdf
implicit none
logical, save :: lenstein      ! Run SIES or not
integer, save :: maxensteinit  ! Maximum number of iterations
real,    save :: gamma_enstein ! steplength
integer, parameter :: ndim=2

contains 
   subroutine enstein(samples,xf,qf,nrsamp,esamp)
   integer, intent(in)  :: nrsamp            ! Number of samples
   real,    intent(out) :: samples(nrsamp,2) ! Returns posterior samples
   integer, intent(in)  :: esamp             ! Number of samples nrsamp=10^esamp for plotting
   real,    intent(in)  :: xf(nrsamp)        ! Prior samples of x
   real,    intent(in)  :: qf(nrsamp)        ! Prior samples of q

   real, allocatable :: xsamp(:)
   real, allocatable :: qsamp(:)
   real, allocatable :: ysamp(:)
   integer, allocatable :: iconv(:)
   integer n,i,j,k,sumconv
   character(len=40) caseid

   real cxx,cyy,cqq,cyx,cqy,cqx
   real pxx,pyy,pqq,pyx,pqy,pqx
   real, allocatable :: sf(:),si(:,:),grad(:,:),newgrad(:,:)
   real :: xlength=1.0
   real Czz(2,2),Pzz(2,2),PIzz(2,2),CIzz(2,2),Pzy(2,1),Pyz(1,2)
   real fac,tmp

   allocate(xsamp(nrsamp))
   allocate(qsamp(nrsamp))
   allocate(ysamp(nrsamp))
   allocate(iconv(nrsamp)) 

   write(*,'(a)')'++++++++++++++++++++++++++++++++++++++++++++++'
   write(*,'(a)',advance='yes')'EnStein analysis...'
   allocate(sf(2),si(2,nrsamp),grad(2,nrsamp),newgrad(2,nrsamp))

      xlength = 2.0/3.14159265 ; print '(a,f13.5)','xlength=',xlength
!        xlength = 1.0            ; print '(a,f13.5)','xlength=',xlength

   do n=1,nrsamp
!     Initial guess
      iconv(n)=0
      xsamp(n)=xf(n)
      qsamp(n)=qf(n)
      ysamp(n)=func(xsamp(n),qsamp(n))
!     Prior in cost function 
      sf(1)=x0   
      sf(2)=0.0
   enddo
   call cov(Cxx,Cyy,Cqq,Cyx,Cqy,Cqx,xsamp,ysamp,qsamp,nrsamp)
   if (Cqq==0.0) then
      write(*,'(a)',advance='yes')'Only coded for the case with model errors'
      deallocate(xsamp,ysamp,qsamp,iconv)
      return
   endif


   Czz(1,1)=Cxx; Czz(2,2)=Cqq; Czz(1,2)=Cqx; Czz(2,1)=Cqx
   CIzz(1,1)=Cqq; CIzz(2,2)=Cxx; CIzz(1,2)=-Cqx; CIzz(2,1)=-Cqx
   CIzz=CIzz/(Cxx*Cqq-Cqx*Cqx)

   call getcaseid(caseid,'EnSTEIN',-1.0,-1,esamp,sigw,0)
   sumconv=0
   do i=1,maxensteinit
      if (mod(i,10) == 0) then
         write(*,'(a,i4,a)')'i=',i,'...'
      endif
      call cov(Pxx,Pyy,Pqq,Pyx,Pqy,Pqx,xsamp,ysamp,qsamp,nrsamp)
      Pzz(1,1)=Pxx; Pzz(2,2)=Pqq; Pzz(1,2)=Pqx; Pzz(2,1)=Pqx
      PIzz(1,1)=Pqq; PIzz(2,2)=Pxx; PIzz(1,2)=-Pqx; PIzz(2,1)=-Pqx
      PIzz=PIzz/(Pxx*Pqq-Pqx*Pqx)
      Pzy(1,1)=Pyx; Pzy(2,1)=Pqy
      Pyz(1,1)=Pyx; Pyz(1,2)=Pqy

      do n=1,nrsamp
         si(1,n)=xsamp(n)
         si(2,n)=qsamp(n)
      enddo

      do n=1,nrsamp
! Prior term
         call dgemm('N','N',ndim,1,ndim,1.0,CIzz,ndim,si(:,n)-sf(:),ndim,0.0,grad(:,n),ndim)
! C_d^{-1} (y_j - d) with unperturbed data
         fac=(ysamp(n)-d)/sigo**2   
! Data term in EnStein  G=Pxx^{-1} Pxy
         call dgemm('N','N',ndim,1,ndim,fac,PIzz,ndim,Pzy,ndim,1.0,grad(:,n),ndim)
! Data term in Stein with analytic G
!            grad(1,n)=grad(1,n)+dfunc(si(1,n),si(2,n))*fac
!            grad(2,n)=grad(2,n)+fac
      enddo


      do n=1,nrsamp
         newgrad(:,n)=0.0
         do k=1,100
            call random_number(tmp)
            j=nint(tmp*real(nrsamp-1)+1.0)
            if (k==1) j=n
            newgrad(:,n)=newgrad(:,n)+&
                    steinkern(si(:,n),si(:,j),xlength,2)*  &
                    (grad(:,j) + (2.0*(si(:,j) - si(:,n))/xlength))
         enddo
         newgrad(:,n)=newgrad(:,n)/real(100.0)
      enddo

      do n=1,nrsamp
         si(:,n) = si(:,n) - 0.25*newgrad(:,n)
         xsamp(n)=si(1,n)
         qsamp(n)=si(2,n)
         ysamp(n)=func(xsamp(n),qsamp(n))
!         print *,'XXX',xsamp(n),qsamp(n),ysamp(n)
      enddo

!      if (mod(i,i)==0) call tecmargpdf('x',xsamp,nrsamp,caseid,xa,xb,nx,i)

   enddo
   samples(:,1)=xsamp(:)
   samples(:,2)=qsamp(:)

!  Recomputing ysamp with some noise for nicer plotting
!   if (sigw < sigq) then
!      do n=1,nrsamp
!         ysamp(n)=ysamp(n)+sigq*normal()
!      enddo
!   endif
!   call tecmargpdf('x',xsamp,nrsamp,caseid,xa,xb,nx)
!   call tecmargpdf('y',ysamp,nrsamp,caseid,ya,yb,ny)
!   call tecmargpdf('q',qsamp,nrsamp,caseid,qa,qb,nx)
!      call tecpdf(x,y,nx,ny,xsamp,ysamp,nrsamp,xa,ya,dx,dy,caseid)
   deallocate(xsamp,ysamp,qsamp,iconv)
   write(*,'(a)')'Stein analysis completed'
   write(*,'(a)')'++++++++++++++++++++++++++++++++++++++++++++++'
   write(*,'(a)')

end subroutine
end module
