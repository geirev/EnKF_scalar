module mod_inistat
   real, save :: x0      ! Prior estimate
   real, save :: d       ! Measurement
   real, save :: siga    ! std dev in prior
   real, save :: sigw    ! std dev in model
   real, save :: sigo    ! std dev in observation
   real, parameter :: sigq= 0.05
   real, save :: cdd
end module

