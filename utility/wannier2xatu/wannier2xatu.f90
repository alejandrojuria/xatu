program wannier2xatu
    use utils
    implicit none

    ! Load FileName from Arguments List
    call LoadSystem()
    call Export2Xatu()
end program wannier2xatu
