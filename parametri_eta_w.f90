!================================================================================

!           Funkcije koje racunaju nepromjenjive parametre eta i w.
!           Ovi parametri su definirani u: Nada Magistarski rad, str. 33

!================================================================================

!Parametri eta definirani su na sljedeci nacin:

!               eta  = 1.0
!                  1

!                                14 ______
!                                --|               k-8.5
!               eta  =  eta    *   | 2000  * (1.05)        ;  k = 2,3,...,15
!                  k       k-1     |             



!       Napomena:
!       ---------
!       Dakle, postoji ukupno 15 parametara eta. Eta je stoga realni vektor
!       duljine 15.
!================================================================================


subroutine vektor_eta(n,eta)
    use funkcije
    implicit none

!   Input
    integer n
!   Output
    real(8),dimension(15) :: eta
!   Local variable
    integer k


    !Tijelo funkcije
    eta(1) = 1.d0
    do k = 2,15
        eta(k) = eta(k-1) * ((2000.d0)**(1.d0/14.d0)) * ((1.05d0)**(k-8.5d0))
    end do

    return
end subroutine


!================================================================================
!Parametri w definirani su na sljedeci nacin:

!                 w  = 0.2
!                  1

!                                23 _______
!                                --| 2500    
!                 w  =    w    *   | -----      ;  k = 2,3,...,24
!                  k       k-1     |  0.2            


!                 w   = oo
!                  25 


!       Napomena:
!       ---------
!       Dakle, postoji ukupno 25 nepromjenjivih parametara  w. w je stoga realni
!       vektor duljine 25 (posljednji clan koji je beskonacan - stavljen je nula
!       jer se kasnije korigira).

!================================================================================


subroutine vektor_w(n,w)
    use funkcije
    implicit none

!   Input
    integer n
!   Output
    real(8),dimension(25) :: w
!   Local variable
    integer k


    !Tijelo funkcije
    w(1) = 0.2d0
    do k = 2,24
        w(k) = w(k-1) * ((2500.d0/0.2d0)**(1.d0/23.d0))
    end do
    w(25) = 0.d0

    return
end subroutine