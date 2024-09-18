 module utils
    public :: LoadSystem
    public :: Export2Xatu

    private
        character(len=:), allocatable  :: FileName
        integer                        :: nFock, mSize
        integer, allocatable           :: Degen(:)
        real(kind = 8)                 :: Rn(3, 3)
        integer, allocatable           :: iRn(:,:)
        complex(kind = 8), allocatable :: H(:,:,:)
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
            integer        ::  fp, ii, jj, i, j
            real(kind = 8) ::  R, Im
            
            !read filename by terminal arguments
            call LoadArguments()

            ! open as newunit
            open(action = 'read', file=FileName, newunit = fp)
                read(fp, *)
                read(fp, *) Rn(1, :)
                read(fp, *) Rn(2, :)
                read(fp, *) Rn(3, :)
                read(fp, *) mSize
                read(fp, *) nFock 
           
                ! time to allocate
                allocate(H(nFock, mSize, mSize))
                allocate(degen(nFock))
                allocate(iRn(nFock, 3))

                ! degen read, 15 elements by line
                do i = 1, (nFock / 15)
                    read(fp, *) Degen((i - 1)*15 + 1:(i - 1)*15 + 15)
                enddo
       
                ! Last line of degenerecences
                read(fp, *) Degen((i - 1)*15 + 1:(i - 1)*15 + MOD(nFock, 15))  
                read(fp, *)
                
                ! begin hamiltonian read
                do i = 1, nFock
                    read(fp, *) iRn(i, :)
                    do j = 1, mSize*mSize
                        read(fp, *) ii, jj, R, Im
                        H(i, ii, jj) = complex(R, Im)
                    enddo
                    if (i < nFock) read(fp, *)
                enddo
                close(fp)
        
        end subroutine LoadSystem

        !!  Implementar a escrita do xatu.model aqui
        subroutine Export2Xatu()
        end subroutine
end module utils
