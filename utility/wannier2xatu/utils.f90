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

subroutine Export2Xatu
            character(len=11) :: filename
            integer :: inuit, i, ios
            integer :: norb
            real, dimension(3, 3) :: bravais

            bravais = reshape((/2.49824, 0.00000, 0.00000, -1.24912, 2.16353, 0.00000, 0.0,0.0,0.0/), shape(bravais))
            norb = 8

            filename='stuff.model'
            ! Assign a unit number for file I/O

            ! Open the file for writing (create it if it doesn't exist, overwrite if it does)
            open(NEWUNIT=iunit, file=filename, status='replace', action='write',iostat=ios)
                if (ios /= 0) then
                    print *, "Error opening file:", trim(filename)
                    print *, "I/O status code:", ios
                    stop
                end if

                ! Write "Hello, World!" to the file
                write(iunit, '(A)') '# dimesion'
                ! this should print a number based on... bravais lattice size..? or should always be 2?

                write(iunit, '(A)') '# norbitals'
                write(iunit, '(*(I1,1X))') (1, i=1,norb)
                ! write(inuit, *) norb

                write(inuit, '(A)') ''! blank line
                write(iunit, '(A)') '# bravaislattice'
                ! do i = 1,3
                !     write(inuit, '(10F6.10)') bravais(i, :)
                ! end do
                ! write(inuit, *) bravais(1)


                ! write(inuit, *) ! blank line
                write(iunit, '(A)') '# bravaisvectors'


                ! write(inuit, *) ! blank line
                write(iunit, '(A)') '# motif'


                ! write(inuit, *) ! blank line
                write(iunit, '(A)') '# hamiltonian'


                ! write(inuit, *) ! blank line
                write(iunit, '(A)') '# filling'

                ! Close the file
                close(iunit)
        end subroutine