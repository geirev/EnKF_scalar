module m_iniens
use mod_inistat
use mod_xyqgrid
use m_getcaseid
use m_normal
use m_tecpdf
use m_func
use m_tecmargpdf
use m_moments
implicit none
contains
subroutine iniens(samples,xsampini,qsampini,ysampini,dpert,nrsamp,esamp)
   integer, intent(in)   :: nrsamp,esamp
   real,    intent(out)  :: xsampini(nrsamp),qsampini(nrsamp),ysampini(nrsamp),dpert(nrsamp)
   real,    intent(out)  :: samples(nrsamp,2)

   real avex,adevx,sdevx,varx,skewx,kurtx
   character(len=40) caseid
   integer i

   do i=1,nrsamp
      xsampini(i)=siga*normal()
   enddo
   call moments(xsampini,nrsamp,avex,adevx,sdevx,varx,skewx,kurtx)
   do i=1,nrsamp
      xsampini(i)=xsampini(i)-avex
   enddo

   do i=1,nrsamp
      qsampini(i)=sigw*normal()
   enddo
   call moments(qsampini,nrsamp,avex,adevx,sdevx,varx,skewx,kurtx)
   do i=1,nrsamp
      qsampini(i)=qsampini(i)-avex
   enddo

   do i=1,nrsamp
      xsampini(i)=x0+xsampini(i)
      ysampini(i)=func(xsampini(i)) + qsampini(i)
   enddo

   do i=1,nrsamp
      dpert(i)=d+sigo*normal()
   enddo
   print *,'Sampling done'
   samples(:,1)=xsampini(:)
   samples(:,2)=qsampini(:)

   call getcaseid(caseid,'INI',-1.0,-1,esamp,sigw,0)
   call tecpdf(x,y,nx,ny,xsampini,ysampini,nrsamp,xa,ya,dx,dy,caseid)
   call tecmargpdf('x',xsampini,nrsamp,caseid,xa,xb,nx)
   call tecmargpdf('y',ysampini,nrsamp,caseid,ya,yb,ny)
   call tecmargpdf('q',qsampini,nrsamp,caseid,qa,qb,nx)

end subroutine
end module
