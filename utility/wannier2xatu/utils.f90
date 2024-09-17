 module utils
    public :: LoadSystem

    private
        character(len=:), allocatable :: FileName
        !integer, allocatable          :: iRn(:,:)
        !real(kind = 8), allocatable   :: H(:,:,:)
        !VARIAVEL PARA O MOTIF

        !todo
        !w90 hamiltoniana
        !.model hamiltoniana

    contains
        subroutine LoadArguments()
            implicit none
            integer N

            call get_command_argument(1, length = N)
            allocate(character(N) :: FileName)
            call get_command_argument(1, FileName)  
            
            print*, FileName
        end subroutine LoadArguments

        subroutine LoadSystem
            call LoadArguments()

            !!  Implementar a leitura do arquivo wannier90 aqui
        end subroutine

        subroutine Export2Xatu()
        end subroutine
end module utils
