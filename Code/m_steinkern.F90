module m_steinkern
contains
   real function steinkern(s1,s2,h,ndim)
   integer ndim
   real s1(ndim),s2(ndim),h

   steinkern=exp(-dot_product(s1-s2,s1-s2)/h)

   end function
end module
