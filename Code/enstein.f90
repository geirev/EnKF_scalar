!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! E n S T E I N !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (lenstein) then
   write(*,'(a)')'++++++++++++++++++++++++++++++++++++++++++++++'
!    EnStein update
      write(*,'(a)',advance='yes')'EnStein analysis...'
      allocate(sf(2),si(2,nrsamp),grad(2,nrsamp),newgrad(2,nrsamp))
!      allocate(kern(nrsamp,nrsamp))

         xlength = 2.0/3.14159265 ; print '(a,f13.5)','xlength=',xlength
!        xlength = 1.0            ; print '(a,f13.5)','xlength=',xlength

      do n=1,nrsamp
!        Initial guess
         iconv(n)=0
         xsamp(n)=xsampini(n)
         qsamp(n)=qsampini(n)
         ysamp(n)=func(xsamp(n),beta,funcmode)+qsamp(n)
!        Prior in cost function 
         sf(1)=x0   
         sf(2)=0.0
      enddo
      call cov(Cxx,Cyy,Cqq,Cyx,Cqy,Cqx,xsamp,ysamp,qsamp,nrsamp)
      Czz(1,1)=Cxx; Czz(2,2)=Cqq; Czz(1,2)=Cqx; Czz(2,1)=Cqx
      CIzz(1,1)=Cqq; CIzz(2,2)=Cxx; CIzz(1,2)=-Cqx; CIzz(2,1)=-Cqx
      CIzz=CIzz/(Cxx*Cqq-Cqx*Cqx)

      call getcaseid(caseid,'EnSTEIN',alphageo,nmda,esamp,gradient,beta,sigw,0)
      sumconv=0
      do i=1,50 !maxiesit
         if (mod(i,10) == 0) then
            write(*,'(a,i4,a)'),'i=',i,'...'
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
!            grad(1,n)=grad(1,n)+dfunc(si(1,n),beta,funcmode)*fac
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
            ysamp(n)=func(xsamp(n),beta,funcmode) + qsamp(n)
         enddo

         if (mod(i,i)==0) call tecmargpdf('x',xsamp,nrsamp,caseid,xa,xb,nx,i)

      enddo
      samples(:,1,4)=xsamp(:)
      samples(:,2,4)=qsamp(:)

!    Recomputing ysamp with some noise for nicer plotting
      if (sigw < sigq) then
         do n=1,nrsamp
            ysamp(n)=ysamp(n)+sigq*normal()
         enddo
      endif
      call tecmargpdf('x',xsamp,nrsamp,caseid,xa,xb,nx)
      call tecmargpdf('y',ysamp,nrsamp,caseid,ya,yb,ny)
      call tecmargpdf('q',qsamp,nrsamp,caseid,qa,qb,nx)
!      call tecpdf(x,y,nx,ny,xsamp,ysamp,nrsamp,xa,ya,dx,dy,caseid)
      write(*,'(a)')'Stein analysis completed'
   endif
