!================================================================================

!   Subroutina koja odreduje vrijednosti vektora q, a i b potrebnih za
!   odredjivanje vektora Gi i fi. Vidjeti Nada Magistarski rad str. 20 za
!   odredjivanje vektora q, te str. 22 za vektore a i b                                         

!================================================================================

!   Vektor q je definiran na sljedeci nacin:

!                        kapa
!                            i
!               q   =  -----------
!                i       kapa
!                            i+1

!   pri cemu je kapa - kompleksna specificna elektricna vodljivost
!   sredstva, [S/m].

!   Vaktor a definiran je na sljedeci nacin:

!                        1 - q
!                             i
!               a   =  -----------
!                i         2

!   pri cemu je q(i) prethodno definiran.

!   Vektor b je definiran na sljedeci nacin:

!                        1 + q
!                             i
!               b   =  -----------
!                i         2
!   Svi vektori imaju ukupno n-1 elemenata!

subroutine vektori_q_a_b(n,kapa,qi,ai,bi)
    use funkcije
    implicit none

!   Input
    integer n
    complex(8),dimension(:) :: kapa
!   Output variables
    complex(8),dimension(:) :: qi,ai,bi
!   Local variable
    complex(8),parameter :: one = dcmplx(1.0,0.0)
    integer i

    do i = 1,n-1
        qi(i) = (kapa(i))/(kapa(i+1))
        ai(i) = (one - qi(i)) / 2.d0
        bi(i) = (one + qi(i)) / 2.d0
    end do


    return
end subroutine