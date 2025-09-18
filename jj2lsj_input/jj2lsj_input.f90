program main    
    
    implicit none
    !important - the following needs to be the integer precision of your
    !grasp run... or i thought i coded it another way?
    integer*8,  allocatable :: nevblk(:)
    integer :: numRequest = 100 

    real*8,   allocatable :: energies( : )
    real*8,   allocatable :: energiesForSorting( : )
    real*8,   allocatable :: ave(:)
    !real*8 ,  allocatable :: eigvec(:,:)
    integer,  allocatable :: sortedPointer(:)
    integer,  allocatable :: posIndex(:)
    integer,  allocatable :: blockNumbers(:)
    integer,  allocatable :: indexWriteOut(:,:)
    integer,  allocatable :: blockOffsets (:)

    integer :: thisBlock,index
    !I assume that the energies within each block are in order.. or not? .. i dont have to! 

    !global ground state energy:
    real*8                :: ground

    !global integer sizes, header:
    integer :: NELEC, NCFTOTI, NWI, NVECTOTI, NVECSIZI, NBLOCKI

    !J pi blocks:
    integer :: IB, NCFII, NVECII, IATJP, IASPA
    integer,allocatable :: ibarray(:)
    integer,allocatable :: ncfiiarray(:)
    integer,allocatable :: nveciiarray(:)
    integer,allocatable :: IATJParray(:)
    integer,allocatable :: IASPAarray(:)
    real*8, allocatable :: blockCopy(:) !for energies

    character*6 :: g92
    integer :: offset = 0 ,i,j,jj,Jblock,ii

    open (1,file='grasp0.cm',form='UNFORMATTED')

    read(1)  g92

    READ (1) NELEC, NCFTOTI, NWI, NVECTOTI, NVECSIZI, NBLOCKI
    
    numRequest = min(numRequest,NVECTOTI)

    allocate(ave(NBLOCKI))
    allocate(energies(NVECTOTI))

    allocate(ibarray(NBLOCKI))
    allocate(ncfiiarray(NBLOCKI))
    allocate(nveciiarray(NBLOCKI))
    allocate(IATJParray(NBLOCKI))
    allocate(IASPAarray(NBLOCKI))
    allocate(nevblk(NVECTOTI))
    allocate(blockNumbers(NVECTOTI))
    OFFSET = 0 
    JBLOCK = 1
    DO JBLOCK = 1, NBLOCKI
        READ (1) IB, NCFII, NVECII, IATJP, IASPA
        READ (1)    (nevblk(I+OFFSET),I=1,NVECII)
        READ (1) ave(jblock), (energies(I+OFFSET),I=1,NVECII)
        energies(1+OFFSET : OFFSET + NVECII) =                         &
        energies(1+OFFSET : OFFSET + NVECII) + ave(jblock)
        READ (1) !skip vectors
        do ii = 1,nvecii 
            blockNumbers(ii + offset) = JBLOCK
        end do 
        offset = offset + NVECII
    END DO 

    allocate(sortedPointer(NVECTOTI))
    do i = 1,NVECTOTI 
        sortedPointer(i) = i 
    end do 
    write(0,*) 'found ground value ',ground

    call qsort(energies,NVECTOTI,sortedPointer)
    ground = energies(1)

    !do ii = 1,numRequest
    !    write (*,'(a)') energies(ii), blockNumbers(sortedPointer(ii)), nevblk(sortedPointer(ii))
    !end do 

    allocate(indexWriteOut(NBLOCKI,numRequest))

    allocate(blockOffsets(nblocki))

    blockOffsets = 0
    indexWriteOut = 0 
    do ii = 1, numRequest 
        thisBlock = blockNumbers(sortedPointer(ii))

        blockOffsets(thisBlock) = blockOffsets(thisBlock) + 1
        
        index = blockOffsets(thisBlock)
        !write(0,*) index
        indexWriteOut(thisBlock,index) = nevblk(sortedPointer(ii))
        write(100,*) nevblk(sortedPointer(ii))
        !write (*,'(a)')thisBlock,index,blockOffsets(thisBlock)
        !if (indexWriteOut(thisBlock,index) .eq. 0) print*'yup'
    end do 


    write (*,'(a)') 'grasp0'
    write (*,'(a)') 'y'
    write (*,'(a)') 'n'
    write (*,'(a)') 'n'
    write (*,'(a)') 'n'
    ii=1
    write (*,'(I2)') ii
    write (*,'(100(I2,1X))') indexWriteOut(ii,1:blockOffsets(ii))
    do ii = 2,nblocki
        if (blockOffsets(ii) .gt. 0) then 
           write (*,'(a)')'y'
           write (*,'(I2)') ii
           write (*,'(100(I2,1X))') indexWriteOut(ii,1:blockOffsets(ii))
           if (indexWriteOut(ii,blockOffsets(ii).gt. 99)) then 
                stop 'bigger than 99'
           end if
        end if 
    end do 

    write (*,'(a)')'n'
    write (*,'(a)') '1.0'
    write (*,'(a)') '0.005'
    write (*,'(a)') '0.001'
    write (*,'(a)') 'n'
    write (*,'(a)') 'n'
    write (*,'(a)') '0'

    close(1)

    contains 

        SUBROUTINE qsort(a, n, t)
    !    NON-RECURSIVE STACK VERSION OF QUICKSORT FROM N.WIRTH'S PASCAL
    !    BOOK, 'ALGORITHMS + DATA STRUCTURES = PROGRAMS'.
    !    REAL*8 PRECISION, ALSO CHANGES THE ORDER OF THE ASSOCIATED ARRAY T.
            IMPLICIT NONE
            INTEGER, INTENT(IN)    :: n
            REAL*8, INTENT(INOUT)    :: a(n)
            INTEGER, INTENT(INOUT) :: t(n)
    !    Local Variables
            INTEGER    :: i, j, k, l, r, s, stackl(15), stackr(15), ww
            REAL*8    :: w, x
            s = 1
            stackl(1) = 1
            stackr(1) = n
    !    KEEP TAKING THE TOP REQUEST FROM THE STACK UNTIL S = 0.
     10     CONTINUE
            l = stackl(s)
            r = stackr(s)
            s = s - 1
    !    KEEP SPLITTING A(L), ... , A(R) UNTIL L >= R.
     20     CONTINUE
            i = l
            j = r
            k = (l+r) / 2
            x = a(k)
    !    REPEAT UNTIL I > J.
          DO
          DO
           IF (a(i).LT.x) THEN    ! Search from lower end
          i = i + 1
          CYCLE
           ELSE
           EXIT
           END IF
          END DO
          DO
          IF (x.LT.a(j)) THEN    ! Search from upper end
          j = j - 1
          CYCLE
          ELSE
          EXIT
          END IF
          END DO
          IF (i.LE.j) THEN    ! Swap positions i & j
          w = a(i)
          ww = t(i)
          a(i) = a(j)
          t(i) = t(j)
          a(j) = w
          t(j) = ww
          i = i + 1
          j = j - 1
          IF (i.GT.j) EXIT
          ELSE
          EXIT
          END IF
          END DO
          IF (j-l.GE.r-i) THEN
          IF (l.LT.j) THEN
          s = s + 1
          stackl(s) = l
          stackr(s) = j
          END IF
          l = i
          ELSE
          IF (i.LT.r) THEN
          s = s + 1
          stackl(s) = i
          stackr(s) = r
          END IF
          r = j
          END IF
          IF (l.LT.r) GO TO 20
          IF (s.NE.0) GO TO 10
          RETURN
    END SUBROUTINE qsort


end program 