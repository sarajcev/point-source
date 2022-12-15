!================================================================================

!           Subroutina koja racuna velicinu H koja predstavlja z koordinatu
!           donje granicne plohe za svaki sloj n-slojnog modela.

!================================================================================

!Velicina H skicirana je na donjoj slici:


!       kapa(1)       h(1)=0       H(1)=0               zrak
!       ----------------------------------------------------
!                   ^           |       |       |       tlo
!       kapa(2)     | h(2)      |  H(2) |       |
!                   v           v       |       |
!       --------------------------------|-------|-----------
!                   ^                   | H(3)  |       
!                   |                   |       |
!       kapa(3)     | h(3)              |       |
!                   v                   v       |
!       ----------------------------------------|-----------
!                                               |
!                               .               | H(n-1)
!                               .               |
!                               .               |
!                                               |
!       ----------------------------------------|-----------
!                   ^                           |
!                   |                           |
!       kapa(n-1)   | h(n-1)                    |
!                   v                           v
!       ----------------------------------------------------

!       kapa(n)       h(n)=0                      H(n)=0



!pri cemu su:
!      n - ukupni broj slojeva modela (zrak + viseslojno tlo)
!   h(i) - debljine pojedinih slojeva n-slojnog modela, [m]. Vektor h ima
!          ukupno n elemenata. Ove vrijednosti se zadaju s ulaznim podacima.
!   H(i) - z koordinata donje graniène plohe pojedinog i-tog sloja, [m]. Ovaj
!          vektor ima ukupno n elemenata, pri cemu je H(1) = 0. H(2) je z
!          koordinata prvog sloja u tlu, H(3) drugog u tlu itd. H(n) je nula.
!================================================================================


subroutine vektor_H(n,h_sloja,HD)
    use funkcije
    implicit none

!   Input
    integer n
    real(8),dimension(:) :: h_sloja
!   Output
    real(8),dimension(:) :: HD
!   Local variable
    integer i


    !Tijelo funkcije
    HD(1) = 0.d0
    do i = 2,n-1
        HD(i) = HD(i-1) + h_sloja(i)
    end do
    HD(n) = 0.d0

    return
end subroutine