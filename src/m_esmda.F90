module m_esmda
use mod_inistat
use m_cyyreg
use m_cov
use m_getcaseid
use m_func
use m_getalpha
use m_normal
use m_tecmargpdf
use mod_xyqgrid
implicit none
logical, save :: lmda           ! Run esmda or not
integer, save :: nmda           ! number of mda iterations
real,    save :: alphageo       ! geometrical factor for alpha sequence
logical, save :: lesmdaadjoint=.false. ! Run IES with adjoint sensitivites
contains
subroutine esmda(samples,xsampini,qsampini,nrsamp,esamp)
   integer, intent(in)  :: nrsamp
   integer, intent(in)  :: esamp
   real,    intent(out) :: samples(nrsamp,2)
   real,    intent(in)  :: xsampini(nrsamp) 
   real,    intent(in)  :: qsampini(nrsamp) 


   real, allocatable :: xsamp(:)
   real, allocatable :: qsamp(:)
   real, allocatable :: ysamp(:)
   real, allocatable :: alpha(:)

   integer n,i
   real Cxx,Cyy,Cqq,Cxy,Cqy,Cxq,alphasum
   real pert
   real dg(2) ! model sensitivity wrt x and q

   allocate(xsamp(nrsamp))
   allocate(qsamp(nrsamp))
   allocate(ysamp(nrsamp))
   allocate(alpha(nmda))

   write(*,'(a)')'++++++++++++++++++++++++++++++++++++++++++++++'
   write(*,'(a)')'MDA analysis...'
   do i=1,nrsamp
      xsamp(i)=xsampini(i)
      qsamp(i)=qsampini(i)
      ysamp(i)=func(xsamp(i),qsamp(i))
   enddo

   alphasum=0.0
   do n=1,nmda
      alpha(n)=getalpha(n,nmda,alphageo)
      write(*,'(a,i3,a,i3,a,f13.5,a,i3)',advance='no')'step:',n,' alpha(',n,')=',alpha(n),' nmda=',nmda
      alphasum=alphasum+1.0/alpha(n)

      call cov(Cxx,Cyy,Cqq,Cxy,Cqy,Cxq,xsamp,ysamp,qsamp,nrsamp)
      write(*,'(a,f8.4)',advance='no')'  cyy0=',Cyy

      if (lcyyreg) Cyy=cyyreg(Cxx,Cqq,Cxy,Cqy,Cxq)
!      if (lcyyreg) write(*,'(a,f8.4)',advance='no')'  cyy1=',Cyy

      do i=1,nrsamp
         if (lesmdaadjoint) then
            dg=dfunc(xsamp(i),qsamp(i))
            Cxy=Cxx * dg(1) + Cxq * dg(2)
            Cqy=Cxq * dg(1) + Cqq * dg(2) 
            Cyy=dg(1)*Cxy + dg(2)*Cqy
         endif
         if (i<5) write(*,'(a,i3,3f8.4,a)',advance='no')'  cyy=(',i, Cyy, dg(1), dg(2),')'
         if (i==5) write(*,'(a)')'-'
         pert=sqrt(alpha(n))*sigo*normal()
         xsamp(i)=xsamp(i) + (Cxy/(Cyy+alpha(n)*Cdd))*(d+pert-ysamp(i))
         qsamp(i)=qsamp(i) + (Cqy/(Cyy+alpha(n)*Cdd))*(d+pert-ysamp(i))
         ysamp(i)=func(xsamp(i),qsamp(i))
      enddo
!         call getcaseid(caseid,'MDA',alphageo,nmda,esamp,sigw,n)
!         call tecmargpdf('x',xsamp,nrsamp,caseid,xa,xb,nx)
   enddo

   samples(:,1)=xsamp(:)
   samples(:,2)=qsamp(:)

   write(*,'(a,f8.2)')'ES-MDA analysis completed.  alphasum=',alphasum
   write(*,'(a)')'++++++++++++++++++++++++++++++++++++++++++++++'
   write(*,'(a)')

   deallocate(xsamp,ysamp,qsamp,alpha)

end subroutine
end module
