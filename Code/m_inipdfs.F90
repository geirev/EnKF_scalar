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
   do i=1,nx
      prior(i)=exp(-0.5*(x(i)-x0)**2/siga**2)
   enddo
   sump=sum(prior(:))*dx
   prior=prior/sump
   call tecfunc('prior',prior,x,nx,'x','')
   print '(a)','prior pdf printed to prior.dat'
end subroutine


subroutine compute_likelihood()
! Likelihood
   integer j
   real sump
   do j=1,ny
      datum(j)=exp(-0.5*(y(j)-d)**2/sigo**2)
   enddo
   sump=sum(datum(:))*dy
   datum=datum/sump
   call tecfunc('datum',datum,y,ny,'y','')
   print '(a)','Likelihood printed to datum.dat'
end subroutine

subroutine compute_costfunction()
! Definition of the strong constraint cost function
   integer i
   do i=1,nx
      cost(i)=(x(i)-x0)**2/siga**2 + (func(x(i))-d)**2/sigo**2    
   enddo
   call tecfunc('costf',cost,x,nx,'x','')
end subroutine


subroutine compute_uncond_jointpdf(esamp)
! Definition of the analytical joint unconditional pdf
   integer i,j
   integer esamp
   real sump
   character(len=40) caseid
   do j=1,ny
   do i=1,nx
      pdf(i,j)=exp( -0.5*(x(i)-x0)**2/siga**2            &
                    -0.5*(y(j)-func(x(i)))**2/max(sigw,sigq)**2 )
   enddo
   enddo
   sump=sum(pdf(:,:))*dx*dy
   pdf=pdf/sump
   call getcaseid(caseid,'PDFJ',-1.0,-1,esamp,sigw,0)
   call tecjointpdf(pdf,x,y,nx,ny,caseid)
   call marginalpdf(pdf,margx,margy,nx,ny,x,y,dx,dy)
   call tecfunc('margx',margx,x,nx,'x',caseid)
   call tecfunc('margy',margy,y,ny,'y',caseid)

end subroutine

subroutine compute_cond_jointpdf(esamp)
! Definition of the analytical joint conditional pdf
   integer i,j
   integer esamp
   real sump
   character(len=40) caseid
   do j=1,ny
   do i=1,nx
      pdf(i,j)=exp(-0.5*(x(i)-x0)**2/siga**2            &
                   -0.5*(y(j)-func(x(i)))**2/max(sigw,sigq)**2  &
                   -0.5*(y(j)-d)**2/sigo**2)    
   enddo
   enddo
   sump=sum(pdf(:,:))*dx*dy
   pdf=pdf/sump
   call getcaseid(caseid,'PDFC',-1.0,1,esamp,sigw,0)
   call tecjointpdf(pdf,x,y,nx,ny,caseid)
   call marginalpdf(pdf,margx,margy,nx,ny,x,y,dx,dy)
end subroutine


subroutine compute_marginals(esamp)
   integer i,j
   integer esamp
   character(len=40) caseid
   margx(:)=0.0
   do i=1,nx
   do j=1,ny
      if (sigw==0.0) then
         margx(i)=margx(i)+exp(-0.5*(x(i)-x0 )**2/siga**2                   &
                               -0.5*(d-func(x(i)))**2/sigo**2)
      else 
         margx(i)=margx(i)+exp(-0.5*(x(i)-x0 )**2/siga**2                   &
                               -0.5*(q(j)-0.0)**2/max(sigw,0.0001)**2  & 
                               -0.5*(d-func(x(i))-q(j))**2/sigo**2)
      endif
   enddo
   enddo
   margx=margx/(sum(margx(:))*dx)
   call getcaseid(caseid,'PDFC',-1.0,1,esamp,sigw,0)
   call tecfunc('margx',margx,x,nx,'x',caseid)
   call tecfunc('margy',margy,y,ny,'y',caseid)
end subroutine

end module
