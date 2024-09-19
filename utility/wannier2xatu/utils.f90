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
            implicit none
            character(len=11) :: outfile
            real*8 :: a1, a2, a3
            
            integer :: iunit, stat, i, j, k

            ! bravais = reshape((/2.49824, 0.00000, 0.00000, -1.24912, 2.16353, 0.00000, 0.0,0.0,0.0/), shape(bravais))
            ! mSize = 8
            ! nFock = 21

            outfile='stuff.model'
            ! Assign a unit number for file I/O

            ! Open the file for writing (create it if it doesn't exist, overwrite if it does)
            open(NEWUNIT=iunit, file=trim(outfile), action='write', status='replace', iostat=stat)
                write(iunit, '(A)') '# dimension'
                write(iunit, '(A)') '2'
                write(iunit, *) ''
                ! write(inuit, *)
                ! this should print a number based on... bravais lattice size..? or should always be 2?

                write(iunit, '(A)') '# norbitals'
                write(iunit, '(*(I1,1X))') (1, i=1,mSize)
                write(iunit, *) ''

                write(iunit, '(A)') '# bravaislattice'
                do i = 1,3
                    write(iunit, *) Rn(i, :)
                end do
                write(iunit, *) ''


                ! ! write(inuit, *) ! blank line
                write(iunit, '(A)') '# bravaisvectors'
                do i=1, nFock
                    a1 = iRn(i,1)*Rn(1,1)+iRn(i,2)*Rn(2,1)+iRn(i,3)*Rn(3,1)
                    a2 = iRn(i,1)*Rn(1,2)+iRn(i,2)*Rn(2,2)+iRn(i,3)*Rn(3,2)
                    a3 = iRn(i,1)*Rn(1,3)+iRn(i,2)*Rn(2,3)+iRn(i,3)*Rn(3,3)
                    write(iunit, *) a1,'    ',a2,'  ',a3
                end do
                write(iunit, *) ''


                ! ! write(inuit, *) ! blank line
                write(iunit, '(A)') '# motif'
                write(iunit, *) ''


                ! ! write(inuit, *) ! blank line
                write(iunit, '(A)') '# hamiltonian'
                do i=1, nFock
                    do j=1, mSize
                        do k=1,mSize
                        write(iunit, '(F20.15, A, F20.15, A)', advance = 'no') real(H(i, j, k)),' ',aimag(H(i, j, k)), 'j    '
                        end do
                        write(iunit, *) ''
                    end do
                    write(iunit, '(A)') '&'
                end do
                write(iunit, *) ''

                ! ! write(inuit, *) ! blank line
                write(iunit, '(A)') '# filling'
                write(iunit, *) ''

            ! Close the file
            ! flush(iunit)
            close(iunit)
        end subroutine
end module