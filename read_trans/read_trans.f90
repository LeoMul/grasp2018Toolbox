program read 

      IMPLICIT NONE 
          !global lists:
      real*8,   allocatable :: energies( : )
      real*8,   allocatable :: ave(:)
      integer,  allocatable :: posArray(:)
      integer,  allocatable :: symArray(:)
      integer,  allocatable :: pointer(:)
      real*8,   allocatable :: avalue_array(:,:)
  
      real*8 :: a,f,g,e,s  
      integer :: i,j , iter 
      integer :: iostat , NVECTOTI,counter
  
      integer :: s1,s2,ii1,ii2
      integer :: ii,jj
      logical :: Found 
      character*60 :: aa 
      
      character*2 :: guage 
      character*2 :: f1,f2 
      character*4 :: p1,p2 
      character*4 :: j1,j2 
  
      integer :: sym,pos
  
      call getGlobalOrder
      avalue_array = 0.0d0 
  
      open(1,file='final.final.t')
      Found = .false. 
      
      do iter = 1,4
          read(1,*,iostat=iostat)
      end do 
      iostat = 0 

      !electric transitions. 
      counter = 0
      do while(iostat.ge.0)
          read(1,300,iostat=iostat) f1,i,j1,p1,f2,j,j2,p2 ,e ,& 
          guage,a,f,s
          !print*,f1,i,iostat
          if (f1 .eq. 'Ma' ) then 
              write(0,*) 'Magnetic entered'
              cycle
          end if
          if (f1 .eq. 'El' ) then 
              write(0,*) 'Electric entered'
              cycle            
          end if 
          if (f1 .eq. '==' ) cycle
          if (f1 .eq. 'Up' ) cycle
          if (f1 .eq. 'Le' ) cycle
  
          if (f1 .eq. '  ' ) cycle
  
          if(guage .eq. ' C') then
              read(1,301,iostat=iostat)                      guage,a,f,s
          end if 
  
          !write(0,*) f1
  
  
          s1 = getSymmetry(j1,p1)
          s2 = getSymmetry(j2,p2)
          ii1 = findGlobalIndex(s1,i)
          ii2 = findGlobalIndex(s2,j)
          !print*,ii1,ii2
          avalue_array(ii1,ii2) = avalue_array(ii1,ii2)+ a 
          avalue_array(ii2,ii1) = avalue_array(ii2,ii1)+ a 
          counter = counter + 1
      end do 
    
      NVECTOTI = 67
      do ii = 1,NVECTOTI
          do jj =ii+1,NVECTOTI
              write(*,'(i4,i4,es10.2)' )jj,ii,max(avalue_array(ii,jj),1e-30) 
          end do  
      end do 

    
      close(1)

  300 FORMAT(1X,A2,I3,1X,2A4,A2,I3,1X,2A4,0P,F13.2,A2,1P,  &
           3D13.5)
  301 FORMAT(42X,A2,1P,3D13.5)

        contains 

        function getSymmetry(jstring,parstring) result(sym)
            character*4:: jstring
            character*4:: parstring
            integer     :: sym

            !print*,jstring,parstring
            read(jstring,'(I4)') sym 
            sym = 2*sym+1
            if      (parstring .eq. ' -') then 
                sym = sym * (-1) 
            else if (parstring .ne. ' +') then 
                write(0,*) 'parity failure', parstring
                stop 
            end if 
        end function

        function findGlobalIndex(sym,pos) result(res)
            integer :: res,sym,pos 
            logical :: check1,check2
            do res=1,NVECTOTI 
                check1 = posArray(pointer(res)).eq.pos 
                check2 = symArray(pointer(res)).eq.sym 
                if (check1.and.check2) then 
                    return 
                end if 
            end do 
            
            stop 'symmetry not found'

        end function 

        SUBROUTINE getGlobalOrder 

            IMPLICIT NONE 


            integer :: NELEC, NCFTOTI, NWI, NVECSIZI, NBLOCKI
            integer :: IB, NCFII, NVECII, IATJP, IASPA
            integer :: iter,jj
            integer :: counter

            open(2, file='final.m',form='unformatted')
            read(2)
            read(2) NELEC, NCFTOTI, NWI, NVECTOTI, NVECSIZI, NBLOCKI


            allocate(ave(nblocki))
            allocate(energies(NVECTOTI))
            allocate(posArray(NVECTOTI))
            allocate(symArray(NVECTOTI))
            allocate(pointer(NVECTOTI))
            allocate(avalue_array(NVECTOTI,NVECTOTI))
            print*, NELEC, NCFTOTI, NWI, NVECTOTI, NVECSIZI, NBLOCKI
            counter = 0 
            posArray = 0
            do jj = 1,NVECTOTI 
                pointer(jj) = jj
            end do 
            do iter=1,nblocki
                read(2) IB, NCFII, NVECII, IATJP, IASPA
                read(2) (posArray(jj+counter),jj=1,nvecii)
                read(2) ave(iter), (energies(jj+counter),jj=1,nvecii)
                read(2)
                do jj=1,nvecii 
                    symArray(jj+counter) = IASPA*IATJP
                    energies(jj+counter) = energies(jj+counter)        &
                        + ave(iter)
                end do 
                counter = counter + nvecii 
            end do
            close(2)
            
            call qsort(energies,NVECTOTI,pointer)

            !do jj = 1,NVECTOTI
            !    print*,symArray(pointer(jj)),posArray(pointer(jj))
            !end do 


        END SUBROUTINE

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

            subroutine StripSpaces(string)
    character(len=*) :: string
    integer :: stringLen 
    integer :: last, actual

    stringLen = len (string)
    last = 1
    actual = 1

    do while (actual < stringLen)
        if (string(last:last) == ' ') then
            actual = actual + 1
            string(last:last) = string(actual:actual)
            string(actual:actual) = ' '
        else
            last = last + 1
            if (actual < last) &
                actual = last
        endif
    end do

    end subroutine

end program read 