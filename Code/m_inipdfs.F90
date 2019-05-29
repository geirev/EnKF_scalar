module m_inipdfs
use mod_xyqgrid
use mod_inistat
use m_func
use m_tecfunc
use m_tecjointpdf
use m_getcaseid
use m_marginalpdf
implicit none

contains
subroutine compute_prior()
! Definition of the prior
   integer i
   real sump
   real priorxx(nxx),priorqq(nqq)
   do i=1,nxx
      priorxx(i)=exp(-0.5*(xx(i)-x0)**2/siga**2)
   enddo
   sump=sum(priorxx(:))*dxx
   priorxx=priorxx/sump
   call tecfunc('margx_prior',priorxx,xx,nxx,'x','prior')
   print '(a)','prior pdf printed to margx_prior.dat'

   if (sigw > 0) then
      do i=1,nqq
         priorqq(i)=exp(-0.5*(qq(i))**2/sigw**2)
      enddo
      sump=sum(priorqq(:))*dqq
      priorqq=priorqq/sump
      call tecfunc('margq_prior',priorqq,qq,nqq,'q','prior')
      print '(a)','prior pdf printed to margq_prior.dat'
   endif
end subroutine


subroutine compute_likelihood()
! Likelihood
   integer j
   real sump
   real datumd(nyy)
   do j=1,nyy
      datumd(j)=exp(-0.5*(yy(j)-d)**2/sigo**2)
   enddo
   sump=sum(datumd(:))*dyy
   datumd=datumd/sump
   call tecfunc('margy_datum',datumd,yy,nyy,'y','datum')
   print '(a)','Likelihood printed to margy_datum.dat'
end subroutine


subroutine compute_costfunction()
! Definition of the strong constraint cost function
   integer i,j,k
   real costxq(nxx,nqq),costx(nxx)
   if (sigw == 0.0) then
      k=1
      do i=1,nxx
         costx(i)=(xx(i)-x0)**2/siga**2 + (func(xx(i),qq(k))-d)**2/sigo**2    
      enddo
      call tecfunc('costfx',costx,xx,nxx,'x','')
      print '(a)','Cost function printed to costfx.dat'
   elseif (sigw > 0.0) then
      do j=1,nqq
      do i=1,nxx
         costxq(i,j)=(xx(i)-x0)**2/siga**2 + (qq(j))**2/sigw**2 + (func(xx(i),qq(j))-d)**2/sigo**2    
      enddo
      enddo
      call tecjointpdf(costxq,xx,yy,nxx,nyy,'costxq')
      print '(a)','Cost function 2D printed to costxq.dat'
   endif
end subroutine


subroutine compute_uncond_jointpdf(esamp)
! Definition of the analytical joint unconditional pdf
   integer i,j,k
   integer esamp
   real sump
   character(len=40) caseid
   real pdff(nxx,nyy),pdfff(nqq,nxx,nyy), margxx(nxx), margyy(nyy)
   if (sigw==0.0) then
      k=1
      do j=1,nyy
      do i=1,nxx
         pdff(i,j)=exp( -0.5*(xx(i)-x0)**2/siga**2                  &
                        -0.5*(yy(j)-func(xx(i),qq(k)))**2/sigq**2 )
      enddo
      enddo

   elseif (sigw > 0.0) then
      do j=1,nyy
      do i=1,nxx
      do k=1,nqq
         pdfff(k,i,j)=exp( -0.5*(qq(k))**2/sigw**2                  & 
                           -0.5*(xx(i)-x0)**2/siga**2               &
                           -0.5*(yy(j)-func(xx(i),qq(k)))**2/sigq**2 )
      enddo
      enddo
      enddo

      pdff(:,:)=0.0
      do j=1,nyy
      do i=1,nxx
         pdff(i,j)=sum(pdfff(1:nqq,i,j))
      enddo
      enddo
   endif 

   sump=sum(pdff(:,:))*dxx*dyy
   pdff=pdff/sump
   call getcaseid(caseid,'PRIOR  ',-1.0,-1,esamp,sigw,0)
   call tecjointpdf(pdff,xx,yy,nxx,nyy,caseid)
   call marginalpdf(pdff,margxx,margyy,nxx,nyy,xx,yy,dxx,dyy)
   call tecfunc('margx'//caseid,margxx,xx,nxx,'x','PRIOR')
   call tecfunc('margy'//caseid,margyy,yy,nyy,'y','PRIOR')

end subroutine

subroutine compute_cond_jointpdf(esamp)
! Definition of the analytical joint conditional pdf
   integer i,j,k
   integer esamp
   real sump
   character(len=40) caseid
   real pdff(nxx,nyy),pdfff(nqq,nxx,nyy), margxx(nxx), margyy(nyy)

   if (sigw==0.0) then
      k=1
      do j=1,nyy
      do i=1,nxx
         pdff(i,j)=exp(-0.5*(xx(i)-x0)**2/siga**2                       &
                       -0.5*(yy(j)-func(xx(i),qq(k)))**2/sigq**2        &
                       -0.5*(yy(j)-d)**2/sigo**2)    
      enddo
      enddo

   elseif (sigw > 0.0) then
      do j=1,nyy
      do i=1,nxx
      do k=1,nqq
         pdfff(k,i,j)=exp( -0.5*(qq(k))**2/sigw**2                      & 
                           -0.5*(xx(i)-x0)**2/siga**2                   &
                           -0.5*(yy(j)-func(xx(i),qq(k)))**2/sigq**2    &
                           -0.5*(yy(j)-d)**2/sigo**2)    
      enddo
      enddo
      enddo

      do j=1,nyy
      do i=1,nxx
         pdff(i,j)=sum(pdfff(1:nqq,i,j))
      enddo
      enddo
   endif

   sump=sum(pdff(:,:))*dxx*dyy
   pdff=pdff/sump
   call getcaseid(caseid,'BAYES  ',-1.0,-1,esamp,sigw,0)
   call tecjointpdf(pdff,xx,yy,nxx,nyy,caseid)
   call marginalpdf(pdff,margxx,margyy,nxx,nyy,xx,yy,dxx,dyy)
   call tecfunc('margx'//caseid,margxx,xx,nxx,'x','BAYES')
   call tecfunc('margy'//caseid,margyy,yy,nyy,'y','BAYES')
end subroutine

end module
