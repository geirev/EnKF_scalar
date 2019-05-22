module m_ies
use mod_inistat
use mod_xyqgrid
use m_cyyreg
use m_cov
use m_getcaseid
use m_func
use m_normal
use m_tecmargpdf
use m_tecsampini
use m_tecpdf
implicit none
logical, save :: lies        ! Run IES or not
integer, save :: maxiesit    ! maximumn number of iterations
real,    save :: gamma_ies   ! step length used in IES
integer, save :: IESv        ! step length used in IES
real, save :: ies_eps=0.0000001

contains 
   subroutine ies(samples,xf,qf,nrsamp,esamp,dpert)
   integer, intent(in)  :: nrsamp            ! Number of samples
   integer, intent(in)  :: esamp             ! Number of samples nrsamp=10^esamp for plotting
   real,    intent(out) :: samples(nrsamp,2) ! Returns posterior samples
   real,    intent(in)  :: xf(nrsamp)        ! Prior samples of x
   real,    intent(in)  :: qf(nrsamp)        ! Prior samples of q
   real,    intent(in)  :: dpert(nrsamp)     ! The perturbed measurement ensemble

   real, allocatable :: xsamp(:)
   real, allocatable :: qsamp(:)
   real, allocatable :: ysamp(:)
   integer, allocatable :: iconv(:)
!   integer :: nrits=100
!   real, allocatable :: costens(:,:)
!   real, allocatable :: costite(:,:)
!   real, allocatable :: xsampit(:,:)

   integer n,i,j,sumconv
   real grad1,grad2(2)
   real Cxx,Cyy,Cqq,Cyx,Cqy,Cqx,alphasum
   real Czz(2,2),Pzz(2,2),PIzz(2,2),CIzz(2,2),zi(2),zf(2),Pzy(2,1),Pyz(1,2)
   real pxx,Pyy,pqq,pyx,pqy,pqx,Pyymat(1,1)
   real pert,jx,djx,djq,ddjx,dgx,dgq
   character(len=40) caseid

   allocate(xsamp(nrsamp))
   allocate(qsamp(nrsamp))
   allocate(ysamp(nrsamp))
   allocate(iconv(nrsamp)) 

!   allocate(costens(nx,10))
!   allocate(costite(10,nrits))
!   allocate(xsampit(10,nrits))

   write(*,'(a)')'++++++++++++++++++++++++++++++++++++++++++++++'
!    IES update
   write(*,'(a,2f13.5)',advance='yes')'IES analysis...d,sigo=',d,sigo
   do n=1,nrsamp
      iconv(n)=0
      xsamp(n)=xf(n)
      qsamp(n)=qf(n)
      ysamp(n)=func(xsamp(n))+qsamp(n)
   enddo
   call cov(Cxx,Cyy,Cqq,Cyx,Cqy,Cqx,xsamp,ysamp,qsamp,nrsamp)
   Czz(1,1)=Cxx; Czz(2,2)=Cqq; Czz(1,2)=Cqx; Czz(2,1)=Cqx
   CIzz(1,1)=Cqq; CIzz(2,2)=Cxx; CIzz(1,2)=-Cqx; CIzz(2,1)=-Cqx
   CIzz=CIzz/(Cxx*Cqq-Cqx*Cqx)

   sumconv=0
   do i=1,maxiesit
      if (i > 10) then
         if (mod(i,10) == 0) then
            write(*,'(a,i4,a)',advance='no')'i=',i,'...'
            print *,'converged realizations=', sumconv,nrsamp,real(100*sumconv)/real(nrsamp)
         endif
      endif
      call cov(Pxx,Pyy,Pqq,Pyx,Pqy,Pqx,xsamp,ysamp,qsamp,nrsamp)
      if (lcyyreg) Pyy=cyyreg(Pxx,Pqq,Pyx,Pqy,Pqx)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (IESv == 1) then
!        Old method (Evensen 2018b paper)
         if (sigw > 0.0) then
            print *,'IESv=1 do not support model errors'
            exit
         endif

         dgx=(Pyx)/(Pxx)
         djx=0.0

         do n=1,nrsamp
            if (iconv(n) > 0) cycle ! do nothing for converged realizations

            djx=cdd*(xsamp(n)-xf(n))  + dgx*Cxx*(ysamp(n)-dpert(n)) 
            ddjx= Pyy + cdd 

            xsamp(n)=xsamp(n) - gamma_ies*djx/ddjx
            ysamp(n)=func(xsamp(n)) 

            if ((abs(djx) < ies_eps) .and. (abs(djq) < ies_eps)) iconv(n)=1
         enddo

      elseif (IESv == 3 .and. Cqq > 0.0) then
! Matrix version (takes into account initial correlations between x and q
         Pzz(1,1)=Pxx; Pzz(2,2)=Pqq; Pzz(1,2)=Pqx; Pzz(2,1)=Pqx
         Pzy(1,1)=Pyx; Pzy(2,1)=Pqy
         Pyz(1,1)=Pyx; Pyz(1,2)=Pqy

         do n=1,nrsamp
            if (iconv(n) > 0) cycle ! do nothing for converged realizations

            zf(1)=xf(n)
            zf(2)=qf(n)

            zi(1)=xsamp(n)
            zi(2)=qsamp(n)

            grad2 = matmul(Pzz,matmul(CIzz,zi-zf)) &
                  - matmul(Pzy , ( matmul(Pyz, matmul(CIzz,zi-zf)) - ysamp(n) + dpert(n) )) / (Pyy + Cdd)

            zi(:) = zi(:) - gamma_ies*grad2(:)

            xsamp(n)=zi(1)
            qsamp(n)=zi(2)
            ysamp(n)=func(xsamp(n)) + qsamp(n)
   
            if ((abs(grad2(1)) < ies_eps) .and. (abs(grad2(2)) < ies_eps)) iconv(n)=1
         enddo

      elseif (IESv == 3 .and. Cqq == 0.0) then
!        Matrix version without model erros
         do n=1,nrsamp
            if (iconv(n) > 0) cycle ! do nothing for converged realizations

            grad1 = (Pxx/Cxx)*(xsamp(n) - xf(n)) &
                      -  (Pyx/(Pyy + Cdd))*( (Pyx/Cxx)*(xsamp(n) - xf(n)) - ysamp(n) + dpert(n)) 

            xsamp(n) = xsamp(n) - gamma_ies*grad1
            ysamp(n)=func(xsamp(n)) 
   
            if (abs(grad1) < ies_eps) iconv(n)=1

         enddo
      else
         print *,'No matching options for IESv (1,3):',IESv
         stop
      endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      sumconv=sum(iconv(1:nrsamp))
      if (sumconv == nrsamp) then
         write(*,'(a,i4,a)',advance='no')'i=',i,'...'
         print *,'converged realizations=', sumconv,nrsamp,real(100*sumconv)/real(nrsamp)
         print *,'Exiting IES iterations'
         exit
      endif

! Dumping costfunction stuff
!      if (i < nrits) then
!         do j=1,10
!            xsampit(j,i+1)=xsamp(j)
!            costite(j,i+1)= (xsamp(j)-xf(j))**2/siga**2 + (func(xsamp(j))-dpert(j))**2/cdd    
!         enddo
!      endif

   enddo
   samples(:,1)=xsamp(:)
   samples(:,2)=qsamp(:)

! Dumping costfunction stuff
!   call tecsampini('sampiniIES',costite,xsampit,min(nrits,i),10)

! Recomputing ysamp with some noise for nicer plotting
   if (sigw < sigq) then
      do n=1,nrsamp
         ysamp(n)=ysamp(n)+sigq*normal()
      enddo
   endif
   call getcaseid(caseid,'IES',-1.0,-1,esamp,sigw,0)
   call tecmargpdf('x',xsamp,nrsamp,caseid,xa,xb,nx)
   call tecmargpdf('y',ysamp,nrsamp,caseid,ya,yb,ny)
   call tecmargpdf('q',qsamp,nrsamp,caseid,qa,qb,nx)
   call tecpdf(x,y,nx,ny,xsamp,ysamp,nrsamp,xa,ya,dx,dy,caseid)

   deallocate(xsamp,ysamp,qsamp,iconv)
   write(*,'(a)')'IES analysis completed'
   write(*,'(a)')'++++++++++++++++++++++++++++++++++++++++++++++'
   write(*,'(a)')

end subroutine
end module
