
c     ============ sub. ropn ===============================

      subroutine ropn (ix,vri,vro)

c     routine for variable output

c     arguments :
c                ix  -  identifier integer
c                vri -  input string
c                vro -  output string

c     notes  :
c             < 'vro' > = <'vri'+ix >

c ----------------------------------------------------------------------

      character*(*) vri
      character*(*) vro
      character*5 tev1
      character*4 tev2
      character*3 tev3
      character*3 aaa
      character*2 aa
      character*1 a

      aaa = '000'
      aa  = '00'
      a   = '0'
      if(ix.lt.10) then
        write(tev1,'(a,a)') vri,aaa
        write(vro,'(a,i1)') tev1,ix
        return
      else if(ix.lt.100) then
        write(tev2,'(a,a)') vri,aa
        write(vro,'(a,i2)') tev2,ix
        return
      else if(ix.lt.1000) then
        write(tev3,'(a,a)') vri,a
        write(vro,'(a,i3)') tev3,ix
        return
      else if(ix.lt.10000) then
        write(vro,'(a,i4)') vri,ix
      else
        write(*,*) ' superato numero massimo di file '
      end if

      return
      end
