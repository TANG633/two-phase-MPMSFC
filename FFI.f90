!-------------------------------------------------------------------
!-   free reading procedure                                        -
!-------------------------------------------------------------------
module FFI
  integer:: iord = 11                   ! read input file
  integer:: iomsg = 13                  ! message file
  integer:: iow1 = 14, iow2 = 15        ! write output file
  integer:: iow03 = 101,iow04 = 102     ! plot energy and momentum
  integer:: iow05 = 103                 ! plot the contact force

  character*80:: jobname, FileInp, FileOut
  character*80:: fName, FileCurv, FileAnim

!------------------------------------------------------------------
!-  line_nb - current line number                                 -
!-  nb_word - word left in one line, a pointer of reading position-
!-  nb_read - word read in one line                               -
!------------------------------------------------------------------
  integer:: line_nb, nb_word, nb_read   ! 已经读入行 当前行所包含的数据总数 当前行已经读入的数据
  character(256):: sss                  ! 存储输入数据的当前行
  character(20):: cmd_line(15)          ! 存储当前行读入的有效数据

contains

  subroutine FFIOpen()
!-------------------------------------------------------------------
!-      purpose: Open input/output file and initialize             -
!-------------------------------------------------------------------
    implicit none

    logical:: ext
    character(len=60):: msg

    inquire (file=FileInp, exist=ext)
    
    if (.not. ext) then
        write(msg, '("Input file ", A10, " does not exist .")') FileInp
        write(*,*) msg
        stop 1000
    end if
    !read in
    open(iord, file=FileInp, status='old')
    open(iomsg, file=FileOut, status='unknown')
    !initial
    line_nb = 0
    nb_word = 0
    nb_read = 0

  end subroutine FFIOpen

  subroutine GetString(mystring)
!------------------------------------------------------------------
!-   purpose: Get a string                                        -
!------------------------------------------------------------------
    implicit none
    
    character(256):: mystring
    
    mystring = sss
    nb_word = 0

  end subroutine GetString

  integer function GetInt()
!------------------------------------------------------------------
!-   purpose: Get an integer number                               -
!-   inputs:                                                      -
!-   outputs:                                                     -
!-      GetInt - the integer number                               -
!------------------------------------------------------------------
    implicit none

    do while (nb_word.eq.0)
        nb_word = ReadLine()
        nb_read = 0
    end do

    if (isNumber(cmd_line(nb_read+1))) then
        ! change string to digit
        read(cmd_line(nb_read+1),*) GetInt 
        nb_word = nb_word - 1
        nb_read = nb_read + 1
    else
        print *, '*** Warning *** Non Numeric Field ', cmd_line(nb_read+1), 'on line ', line_nb
        write(iomsg,*), '*** Warning *** Non Numeric Field ', cmd_line(nb_read+1), 'on line ', line_nb
        nb_word = nb_word - 1
        nb_read = nb_read + 1
    end if

  end function GetInt

  real function GetReal()
!------------------------------------------------------------------
!-   purpose: Get a real number                                   -
!-   inputs:                                                      -
!-   outputs:                                                     -
!-      GetReal  - the real number                                -
!------------------------------------------------------------------
    implicit none

    do while (nb_word.eq.0)
        nb_word = ReadLine()
        nb_read = 0
    end do

    if (isNumber(cmd_line(nb_read+1))) then
        read(cmd_line(nb_read+1),*) GetReal
        nb_word = nb_word - 1
        nb_read = nb_read + 1
    else
        print *, '*** Warning *** Non Numeric Field ', cmd_line(nb_read+1), 'on line ', line_nb
        write(iomsg,*) '*** Warning *** Non Numeric Field ', cmd_line(nb_read+1), 'on line ', line_nb
        nb_word = nb_word - 1
        nb_read = nb_read + 1
    end if

  end function GetReal
!从输入文件中读入一个关键字，判断其是否存在于关键词组kw（nbkw）中 存在返回序号 不存在返回-1
  integer function KeyWord(kw,nbkw)
!--------------------------------------------------------------------
!-   purpose: read one line from iord file and identify the keyword -
!-   inputs:                                                        -
!-     kw   - keyword set                                           -
!-     nbkw - number of keyword                                     -
!-  outputs:                                                        -
!-     Keyword - position of inputed keyword                        -
!--------------------------------------------------------------------
    implicit none
    
    character(4):: kw(*)
    integer:: i, nbkw

    KeyWord = -1

    do while (nb_word.eq.0)
        nb_word = ReadLine()
        nb_read = 0
    end do
    
    do i = 1, nbkw
        if (pcomp(cmd_line(nb_read+1),kw(i),4)) then
            KeyWord = i
            nb_word = nb_word - 1
            nb_read = nb_read + 1
            exit
        end if
    end do
    
    if (KeyWord .eq. -1) then
        print *, '*** Warning *** invalid keyword ', cmd_line(nb_read+1), 'on line ', line_nb
        write(iomsg,*) '*** Warning *** invalid keyword ', cmd_line(nb_read+1), 'on line ', line_nb
    end if

  end function KeyWord
!读取输入数据文件中的当前行 存储到sss 提取有效数据到cmd_line 返回该行所包含的数据总数
  integer function ReadLine()
!------------------------------------------------------------------
!-  purpose: read one line from iord file and parse it            -
!-  inputs:                                                       -
!-  outputs:                                                      -
!-     ReadLine - number of word in this line                     -
!-     cmd_line() is filled with words in this line               -
!------------------------------------------------------------------
    implicit none
    
    integer:: i, j, lensss

    ! initialize cmd_line
    cmd_line = ' '
    ! read one line from iord
    read(iord,'(a)',err=901,end=902) sss
    line_nb = line_nb + 1

    ! Strip horizontal tab character (Ctrl-I = ASCII Character 9)
    do i = 1,256
        if(ichar(sss(i:i)).eq.9) sss(i:i) = ' '
    end do
    ! Strip comments
    do i = 1, 256
        if(sss(i:i).eq.'!') then
            sss(i:256) = ' '
            exit
        end if
    end do
    ! Strip leading blanks
    sss = adjustl(sss)
    lensss = len_trim(sss)
    if (lensss.eq.0) then       ! blank line or comment line
        ReadLine = 0
        return
    else
        ReadLine = 1
    end if
    ! Remove extra blanks before separator
    i = 1
    do while (i.lt.lensss)
        do while (sss(i:i).eq.' ' .and. (sss(i+1:i+1).eq.' ' .or. sss(i+1:i+1).eq.','))
            sss(i:lensss-1) = sss(i+1:lensss)
            sss(lensss:lensss) = ' '
            lensss = lensss - 1
        end do
        i = i + 1
    end do
    ! Remove extra blanks after ',' and parse sss to cmd_line
    j = 1
    i = 1
    do while (i.lt.lensss)
       if (sss(i:i).eq.',' .and. sss(i+1:i+1).eq.' ') then
          sss(i+1:lensss-1) = sss(i+2:lensss)
          sss(lensss:lensss) = ' '
          lensss = lensss - 1
       end if
       if (sss(i+1:i+1).eq.',' .or. sss(i+1:i+1).eq.' ') then
          cmd_line(ReadLine) = sss(j:i)
          j = i + 2
          ReadLine = ReadLine + 1
       end if
       i = i + 1
    end do
    cmd_line(ReadLine) = sss(j:lensss)

    return

    ! read error encountered
901 write(iomsg,"('error on reading, line ',i5)") line_nb
    stop 'error on reading'

    ! EOF encountered    
902 write(iomsg,*) 'End of File!'
    stop 'encounter EOF when reading without seeing END'

  end function ReadLine
!比较字符串a b 的前n个字符是否相同 Ignores upper/lower case differences
  logical function pcomp(a,b,n)
!------------------------------------------------------------------
!-  Purpose: Compare character strings for match                  -
!-           Ignores upper/lower case differences.                -
!-  Inputs:                                                       -
!-     a(*)   - Character string 1                                -
!-     b(*)   - Character string 2                                -
!-     n      - Number of characters to compare                   -
!-  Outputs:                                                      -
!-     pcomp  - Flag, true if a = b                               -
!------------------------------------------------------------------
    implicit  none
    
    integer:: n, inc, i, ia, ib
    character:: a*(*), b*(*)

    pcomp = .false.

    ! Compute increment between an upper and lower case letter
    inc = ichar('A') - ichar('a')

    ! Compare for match
    do i = 1,n
        ia = ichar(a(i:i))
        ib = ichar(b(i:i))
        ! Test all permutations of characters for match
        if (ia.ne.ib .and. ia+inc.ne.ib .and. ia.ne.ib+inc ) return
    end do

    pcomp = .true.

  end function pcomp

  logical function isNumber(num)
!------------------------------------------------------------------
!-   purpose: verify a numeric field                              -
!------------------------------------------------------------------
    implicit none
    
    character(20):: num
    character:: a

    a = num(1:1)
    isNumber = .false.
    if(a.eq.'-' .or. (ichar(a).gt.47 .and. ichar(a).lt.58)) isNumber = .true.

  end function isNumber

  subroutine ErrorMsg()
!--------------------------------------------------------------------
!-   Purpose                                                        -
!-      Display error message for input error                       -
!--------------------------------------------------------------------
    implicit none
    
    write(*,*) '*** Error *** : Input file error, line', line_nb
    write(iomsg,*) '*** Error *** : Input file error, line', line_nb
    
  end subroutine ErrorMsg

end module FFI