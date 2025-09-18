program grasp2018_shift

    !shifting energy levels for transition calculations. 
    !this program takes a <state name>.m file,
    !and an input shift.dat 
    !and produces a new <state name>.m file.

    !the shift.dat takes a header of the form 
    !numshift, unit (see below)
    !followed by numshift pairs of 
    !index, newEnergy.
    !where index is the position in the global (unshifted) energy order.

    !the outputted file will have the average block energies, and relative energies
    !adjusted so that each of the shifted levels has energy newEnergy.

    !the blocks are also reordered so that the eigenvectors (and energies) are in their
    !shifted order. 

    implicit none 
    !shifting variables:
    integer :: numshift 
    real*8  :: unit !unit relative to a.u.
    !i.e - what number do I multiply the shifts in shift.dat
    !to get atomic units. I.e - 0.5 for Ryd.
    !i.e I multiply an energy in Ryd by 0.5 to get a.u.
    real*8 :: shift

    !global lists:
    real*8,   allocatable :: energies( : )
    real*8,   allocatable :: energiesForSorting( : )
    real*8,   allocatable :: ave(:)
    real*8 ,  allocatable :: eigvec(:,:)
    integer,  allocatable :: sortedPointer(:)
    integer,  allocatable :: nevblk(:)
    integer,  allocatable :: posIndex(:)

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

    !misc 
    character*6 :: g92
    integer :: offset = 0 ,i,j,jj,Jblock



    open(1,file='final.m',form='UNFORMATTED')

    open(2,file='test.m',form='UNformatted')

    read(1)  g92
    write(2) g92 

    READ (1) NELEC, NCFTOTI, NWI, NVECTOTI, NVECSIZI, NBLOCKI
    write(2)NELEC, NCFTOTI, NWI, NVECTOTI, NVECSIZI, NBLOCKI
    
    !allocate statements
    allocate(eigvec(NCFTOTI,NVECTOTI))
    allocate(nevblk(NVECTOTI))

    allocate(ave(NBLOCKI))
    allocate(energies(NVECTOTI))

    allocate(ibarray(NBLOCKI))
    allocate(ncfiiarray(NBLOCKI))
    allocate(nveciiarray(NBLOCKI))
    allocate(IATJParray(NBLOCKI))
    allocate(IASPAarray(NBLOCKI))

    
    !jblock

    !INPUT 
    jblock = 1
    offset= 0 
    DO JBLOCK = 1, NBLOCKI
        READ (1) IB, NCFII, NVECII, IATJP, IASPA

         ibarray    (JBLOCK)    = IB     
         NCFIIarray(JBLOCK) = NCFII   
         NVECIIarray(JBLOCK)= NVECII  
         IATJParray(JBLOCK) = IATJP  
         IASPAarray(JBLOCK) = IASPA  

        READ (1)    (nevblk(I+OFFSET),I=1,NVECII)
        !PRINT*,(nevblk(I+OFFSET),I=1,NVECII)
        READ (1) ave(jblock), (energies(I+OFFSET),I=1,NVECII)
        print*,NVECII,NVECTOTI
        !add back the average energy for now
        energies(1+OFFSET : OFFSET + NVECII) = energies(1+OFFSET : OFFSET + NVECII) + ave(jblock)
        
        READ (1) ((eigvec(I,J+offset),i=1,ncfii),j=1,NVECII)
        offset = offset + NVECII
    END DO 

    !do the shifting in here... 
    !? ? ? 
    !profit 

    
    energiesForSorting = energies 
    allocate(sortedPointer(NVECTOTI))
    do i = 1,NVECTOTI 
        sortedPointer(i) = i 
    end do 
    ground = minval(energies)
    write(0,*) 'found ground value ',ground




    call qsort(energiesForSorting,NVECTOTI,sortedPointer)

    print*,'hello', sortedPointer(2)

    open(4,file='shift.dat')
    read(4,*) numShift, unit 
    do i= 1,numshift 
        read(4,*) jj , shift
        print*, energies(sortedPointer(jj)),ground + shift * unit 
        energies(sortedPointer(jj)) = ground + shift * unit 
    end do 
    close(4)


    !energies(sortedPointer(2)) = ground + 0.041448 / 2.0 

    !OUTPUT 
    offset= 0 
    print*,NCFIIarray
    DO JBLOCK = 1, NBLOCKI


        IB     = ibarray(JBLOCK)
        NCFII  = NCFIIarray(JBLOCK)
        NVECII = NVECIIarray(JBLOCK)
        IATJP  = IATJParray(JBLOCK)
        IASPA  = IASPAarray(JBLOCK)

        allocate(posIndex(nvecii))
        do i = 1,nvecii
            posIndex(i) = i 
        end do 

        !subtract the average energy back again 
        !possibly recalculate averages here too ... 
        ave(jblock) = sum(energies(1+OFFSET:NVECII+OFFSET)) / NVECII
        do i = 1,nvecii
            energies(I+OFFSET) = energies(I+OFFSET)-ave(jblock)
        end do 
        
        !this code ensures that the level indexing within
        !the blocks are correct. 
        !i.e - that when they come out of the rtransitions
        !their indices are the correct energy order.

        !this will still require a global energy order at some point
        !... 
        !easy way to post process this? 


        allocate(blockCopy(nvecii))
        blockCopy(1:) = energies(1+offset:nvecii + offset)
        call qsort(blockCopy,nvecii,posIndex)


        !write to the new mixing file.
        write (2) IB, NCFII, NVECII, IATJP, IASPA
        write (2) (i,I=1,NVECII) 
        write (2) ave(jblock), (energies(posIndex(I)+OFFSET),I=1,NVECII)
        write (2) ((eigvec(I,posIndex(J)+offset),i=1,ncfii),j=1,NVECII)

        write (3,*) IB, NCFII, NVECII, IATJP, IASPA
        write (3,*) (i,I=1,NVECII) 
        write (3,*) ave(jblock), (energies(posIndex(I)+OFFSET),I=1,NVECII)
        write (3,*) ((eigvec(I,posIndex(J)+offset),i=1,ncfii),j=1,NVECII)

        offset = offset + NVECII

        deallocate(posIndex)
        deallocate(blockCopy)
    END DO 

    close(2)
    close(1)

    deallocate(eigvec)

    stop 

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

end program grasp2018_shift
