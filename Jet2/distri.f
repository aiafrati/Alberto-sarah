
      subroutine distri(a1,a2,rt1,rt2,et1,et2,n1,n2,rft,nn,kk)

      implicit real*8(a-h,o-z)
      dimension rft(*)
* - nn numero pannelli 
      if (kk.eq.1) then ! numero pannelli variabile
        e1  = et1
        e2  = et2
        if (a1.lt.1.d0) then
          rt1 = 0.5d0
          rt2 = 0.5d0
          n1  = int(log( 1.d0-(1.d0-e1)*rt1/a1)/log(e1) )
          n2  = int(log( 1.d0-(1.d0-e2)*rt2/a2)/log(e2) )
        else
          rt2 = 1.d0
          rt1 = 0.d0
          n2  = log( 1.d0-(1.d0-e2)*rt2/a2 )/log(e2)
          n1  = 0
        end if

  100   continue

        if (n2.lt.1) then 
          rt2 = 0.d0
          rt1 = 1.d0
          n2  = 0
          a2  = 0.d0
        else if (n1.lt.1) then
          rt1 = 0.d0
          rt2 = 1.d0
          n1  = 0
          a1  = 0.d0
        end if

cccc        if (n1.gt.0) a1 = rt1*(1.d0-e1)/(1.d0-e1**n1 )
        if (n1.gt.0) call newrap(a1,rt1,e1,n1)
        if (n2.gt.0) a2 = rt2*(1.d0-e2)/(1.d0-e2**n2 )
        nn = n1+n2

        rft(1) = 0.d0
        do i = 1,n1
          rft(i+1) = rft(i) + a1*e1**(i-1)
        enddo

        rft(nn+1) = 1.d0
        do i = 1,n2-1
          rft(nn+1-i) = rft(nn+2-i) - a2*e2**(i-1)
        enddo

        if (n2.le.1.or.n1.le.1) return

        ae1 = a1*e1**(n1-1)
        ae2 = a2*e2**(n2-1)

        if (ae1.gt.1.1*e1*e2*ae2) then
          n1  = n1 - 1
          rt1 = a1*(1.d0-e1**n1)/(1.d0-e1)
          rt2 = 1.d0 - rt1
          go to 100
        else if (ae2.gt.1.1*e1*e2*ae1) then
          n1  = n1 + 1
          rt1 = a1*(1.d0-e1**n1)/(1.d0-e1)
          rt2 = 1.d0 - rt1
          go to 100
        end if

      else 

        if (n1.gt.0.and.n2.gt.0) then
          rt2 = a2*(1.d0-e2**n2)/(1.d0-e2)
          rt1 = 1.d0 - rt2
cccc          a1  = rt1*(1.d0-e1)/(1.d0-e1**n1)     
          call newrap(a1,rt1,e1,n1)
        else if (n1.lt.1) then
          n2 = nn
          a2 = (1.d0-e2)/(1.d0-e2**n2 )
          n1 = 0 
        else if (n2.lt.1) then
          n1  = nn
          n2  = 0 
          rt1 = 1.d0
ccc          a1 = (1.d0-e1)/(1.d0-e1**n1 )
          call newrap(a1,rt1,e1,n1)
        end if
        nn = n1 + n2

        rft(1) = 0.d0
        do i = 1,n1
          rft(i+1) = rft(i) + a1*e1**(i-1)
        enddo

        rft(nn+1) = 1.d0
        do i = 1,n2
          rft(nn+1-i) = rft(nn+2-i) - a2*e2**(i-1)
        enddo

      end if

      return
      end 

c     ++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine newrap(a1,rt1,e1,n1)

      implicit real*8 (a-h,o-z)

      ev1  = e1      
      de1  = ev1
      kitn = 0
      do while (abs(de1).gt.1.d-9)
        q1 = (1.d0-ev1)
        q2 = (1.d0-ev1**n1)
        ff = q2/q1 - rt1/a1
        dd = ( -float(n1)*q1*ev1**(n1-1) + q2 )/q1**2
        de1 = - ff/dd
        ev1 = ev1 - ff/dd
        kitn = kitn + 1
        if (kitn.gt.1000) stop ' NON CONVERGE N-R '
      end do
      e1 = ev1

      return
      end
