module utils
    public :: LoadSystem
    public :: Export2Xatu

    private
        character(len=:), allocatable  :: FileName
        integer                        :: nFock, mSize
        integer, allocatable           :: Degen(:)
        real(8)                 :: Rn(3, 3)
        integer, allocatable           :: iRn(:,:)
        complex(8), allocatable :: H(:,:,:)
        complex(8), allocatable :: Rhop(:,:,:,:) ! motif

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
            real(8) ::  R, Im
            real(8) ::  a1, a1j, a2, a2j, a3, a3j

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
                allocate(Rhop(3, nFock, mSize, mSize))
                allocate(degen(nFock))
                allocate(iRn(nFock, 3))

                ! degen read, 15 elements by line
                if ((nFock / 15) .gt. 1) then
                    do i = 1, (nFock / 15)
                        read(fp, *) Degen((i - 1)*15 + 1:(i - 1)*15 + 15) 
                    enddo
                end if
                
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
                read(fp,*) ! line skip

                ! begin reading orbital localization -- for motif
                do i = 1, nFock
                    read(fp, *) ! skip iRn
                    do j = 1, mSize*mSize
                        read(fp, *) ii, jj, a1, a1j, a2, a2j, a3, a3j
                        Rhop(1, i, ii, jj) = complex(a1, a1j)
                        Rhop(2, i, ii, jj) = complex(a2, a2j)
                        Rhop(3, i, ii, jj) = complex(a3, a3j)
                    enddo
                    if (i < nFock) read(fp, *)
                enddo
                close(fp)
        
        end subroutine LoadSystem

subroutine Export2Xatu
            implicit none
            character(len=len(FileName)+2) :: outfile ! +2 because (.model=.dat+2)
            integer :: iunit, stat, i, j, k
            integer :: filePos, diag, dimensions
            real(8) :: a1, a2, a3
            logical :: is2D=.True.


            ! Find the position of '.dat' in the filename
            filePos = index(FileName, '.dat')
            if (filePos > 0) then
                ! Replace '.dat' with '.model'
                outfile = FileName(1:filePos-1) // '.model'
            else
                ! If no '.dat' found, just append '.model'
                outfile = FileName // '.model'
            end if

            ! Open the file for writing (create it if it doesn't exist, overwrite if it does)
            open(NEWUNIT=iunit, file=trim(outfile), action='write', status='replace', iostat=stat)
                write(iunit, '(A)') '# dimension'
                do i=1,nFock
                    if (iRn(i,3).ne.0) then
                        is2D = .False.
                    else
                        is2D = .True.  ! 2D system only if all iRn(i, 3) = 0.0
                    end if
                end do

                if (is2D) then
                    write(iunit, *) 2
                    dimensions = 2
                else
                    write(iunit, *) 3
                    dimensions = 3
                end if
            ! ------------------------------------------------------------------------------------ !
                write(iunit, '(A)') '# norbitals'
                write(iunit, '(*(I1,2X))') (1, i=1,mSize)
            ! ------------------------------------------------------------------------------------ !
                write(iunit, '(A)') '# filling'
                write(iunit, *)
            ! ------------------------------------------------------------------------------------ !
                write(iunit, '(A)') '# bravaislattice'
                do i = 1, dimensions
                    write(iunit, *) Rn(i, :)
                end do
            ! ------------------------------------------------------------------------------------ !
                write(iunit, '(A)') '# motif'
                do i=1, nFock
                    if (iRn(i,1).eq.0 .and. iRn(i,2).eq.0 .and. iRn(i,3).eq.0) then
                        diag = i
                    end if
                end do
                do i=1, mSize
                    write(iunit, *) (real(Rhop(k, diag, i,i)),k=1,3), i-1
                end do
            ! ------------------------------------------------------------------------------------ !
                write(iunit, '(A)') '# bravaisvectors'
                do i=1, nFock
                    a1 = iRn(i,1)*Rn(1,1)+iRn(i,2)*Rn(2,1)+iRn(i,3)*Rn(3,1)
                    a2 = iRn(i,1)*Rn(1,2)+iRn(i,2)*Rn(2,2)+iRn(i,3)*Rn(3,2)
                    a3 = iRn(i,1)*Rn(1,3)+iRn(i,2)*Rn(2,3)+iRn(i,3)*Rn(3,3)
                    write(iunit, *) a1,'    ',a2,'  ',a3
                end do
            ! ------------------------------------------------------------------------------------ !
                write(iunit, '(A)') '# hamiltonian'
                do i=1, nFock
                    H(i,:,:) = H(i,:,:)*Degen(i)
                    do j=1, mSize
                        do k=1,mSize
                            write(iunit, '(F20.15, A, F20.15, A)', advance = 'no') real(H(i, j, k)),' ',aimag(H(i, j, k)), 'j    '
                        end do
                        write(iunit, *) ''
                    end do
                    write(iunit, '(A)') '&'
                end do
                write(iunit, *) '#'
            ! ------------------------------------------------------------------------------------ !

            ! flush(iunit)
            close(iunit)
        end subroutine
end module
