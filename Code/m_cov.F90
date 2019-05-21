module m_cov
contains
subroutine cov(Cxx,Cyy,Cqq,Cxy,Cqy,Cqx,xsamp,ysamp,qsamp,nrsamp)
   integer, intent(in)  :: nrsamp
   real, intent(in)  :: xsamp(nrsamp)
   real, intent(in)  :: ysamp(nrsamp)
   real, intent(in)  :: qsamp(nrsamp)
   real, intent(out) :: Cxx
   real, intent(out) :: Cxy
   real, intent(out) :: Cyy
   real, intent(out) :: Cqy
   real, intent(out) :: Cqq
   real, intent(out) :: Cqx      

   real avex
   real avey
   real aveq
   integer i

   avex=sum(xsamp(1:nrsamp))/real(nrsamp)
   avey=sum(ysamp(1:nrsamp))/real(nrsamp)
   aveq=sum(qsamp(1:nrsamp))/real(nrsamp)
   Cxx=0.0
   Cyy=0.0
   Cxy=0.0
   Cqy=0.0
   Cqq=0.0
   Cqx=0.0
   do i=1,nrsamp
      Cxx=Cxx+(xsamp(i)-avex)*(xsamp(i)-avex)
      Cyy=Cyy+(ysamp(i)-avey)*(ysamp(i)-avey)
      Cxy=Cxy+(xsamp(i)-avex)*(ysamp(i)-avey)
      Cqy=Cqy+(qsamp(i)-aveq)*(ysamp(i)-avey)
      Cqq=Cqq+(qsamp(i)-aveq)*(qsamp(i)-aveq)
      Cqx=Cqx+(qsamp(i)-aveq)*(xsamp(i)-avex)
   enddo
   Cxx=Cxx/real(nrsamp-1)
   Cyy=Cyy/real(nrsamp-1)
   Cxy=Cxy/real(nrsamp-1)
   Cqy=Cqy/real(nrsamp-1)
   Cqq=Cqq/real(nrsamp-1)
   Cqx=Cqx/real(nrsamp-1)

end subroutine
end module
