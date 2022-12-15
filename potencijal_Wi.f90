!================================================================================

!   Subroutina koja racuna prvi dio ukupnog potencijala: Wi. Izraz za ovaj
!   dio potencijala dan je u mom draftu na str. 27.                                     

!================================================================================


!   Potencija uzrokovan tockastim izvorom struje ukopanim na dubini d
!   u viseslojnom tlu dan je sljedecim izrazom:


!               Fi  =  V  + W
!                 i     i    i

!   pri cemu ova subroutina racuna dio potencijala oznacen s W. Ovaj dio
!   potencijala racuna se prema sljedecim izrazima:


!                          15
!                         _____               alfa
!             I           \                       k
!   W  = ------------ *   /      [  ---------------------------  +
!    i    4*PI*kapa       -----      _________________________
!                  i       k=1     -| 2                      2
!                                   |r  +(teta*eta + z - H   )
!                                                 k       i-1

!                  beta
!                      k
!       --------------------------  ]
!        ________________________
!      -| 2                     2
!       |r  + (hi*eta  + H  - z)
!                    k    i


!pri cemu su:
!   teta - koeficijent koji je izracunat u programu: parametri_teta_hi.f90
!   hi   - koeficijent koji je izracunat u programu: parametri_teta_hi.f90
!   alfa  , beta  - koeficjenti koji se racunaju metodom najmanjih kvadrata
!       k       k
!   eta  - poznati koeficjenti
!      k

!   Napomena:
!   ---------
!   Ukoliko se izvor nalazi u zraku (s=1) ili u n-tom sloju s=n, a raspodjela 
!   potencijala se trazi u zemlji (n-ti sloj ili bilo koji drugi) ili u zraku,
!   neki clanovi u gornjim izrazima ce nestati!

!================================================================================



subroutine potencijal_Wi(n,I,sloj,kapa_sloj,alfa,beta,r,z,H,teta,hi,eta,Wi)
    use funkcije
    implicit none
!   Input
    integer n
    complex(8) I
    integer sloj
    complex(8) kapa_sloj
    complex(8),dimension(:) :: alfa,beta
    real(8) r,z
    real(8),dimension(:) :: H
    real(8) teta,hi
    real(8),dimension(:) :: eta
!   Output
    complex(8) Wi
!   Local
    real(8),parameter :: pi = 3.14159265
    integer k


!   Opis ulaznih varijabli:
!       n - ukupni broj slojeva modela (zrak + viseslojno tlo)
!       r - r koordinata tocke u kojoj se racuna potencijal, [m]
!       z - z koordinata tocek u kojoj se racuna potencijal, [m]
!       H - vektor koji sadrzi koordinate donjih granicnih ploha svih slojeva
!       I - jakost struje harmonièkog toèkastog izvora
!       sloj - sloj u kojem se racuna raspodjela potencijala (isl)
!       kapa_sloj - vrijednost specificne kompleksne vodljivosti i-tog sloja (sloja
!       u kojem se racuna raspodjela potencijala)
!       alfa,beta - koeficjenti koji se racunaju metodom najmanjih kvadrata
!       teta - koeficijent koji je izracunat u programu: parametri_teta_hi.f90
!       hi   - koeficijent koji je izracunat u programu: parametri_teta_hi.f90
!       eta - poznati koeficijenti definirani u file-u: parametri_eta_w.f90

!   ------------------------------------------
!           Pocetak proracuna
!   ------------------------------------------

!   Slucaj 1
!   --------
!   Trazi se raspodjela potencijala u zraku (sloj=1)
    if(sloj==1) then
        Wi = dcmplx(0.d0,0.d0)
        do k = 1,15
            Wi = Wi + beta(k) / (dsqrt(r**2+(hi*eta(k)+H(sloj)-z)**2))
        end do

!   Slucaj 2
!   --------
!   Trazi se raspodjela potencijala u n-tom sloju (sloj=n)
    else if(sloj==n) then
        Wi = dcmplx(0.0,0.0)
        do k = 1,15
            Wi = Wi + alfa(k) / (dsqrt(r**2+(teta*eta(k)+z-H(sloj-1))**2))
        end do

!   Slucaj 3
!   --------
!   Trazi se raspodjela potencijala u zemlji ali NE u n-tom sloju
    else if((sloj/=1).and.(sloj/=n)) then
        Wi = dcmplx(0.0,0.0)
        do k = 1,15
            Wi = Wi + alfa(k) / (dsqrt(r**2+(teta*eta(k)+z-H(sloj-1))**2)) + &
            beta(k) / (dsqrt(r**2+(hi*eta(k)+H(sloj)-z)**2))
        end do

    end if

    Wi = (I/(4.d0*pi*kapa_sloj)) * Wi

    return
end subroutine