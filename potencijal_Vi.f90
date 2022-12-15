!================================================================================

!   Subroutina koja racuna prvi dio ukupnog potencijala: Vi. Izraz za ovaj
!   dio potencijala dan je u mom draftu na str. 11.                                     

!================================================================================


!   Potencija uzrokovan tockastim izvorom struje ukopanim na dubini d
!   u viseslojnom tlu dan je sljedecim izrazom:


!               Fi  =  V  + W
!                 i     i    i

!   pri cemu ova subroutina racuna dio potencijala oznacen s V. Ovaj dio
!   potencijala racuna se prema sljedecim izrazima:

!   IZVOR U ZEMLJI (ne u n-tom sloju), RASPODJELA U ZEMLJI (ne u n-tom sloju):

!   za  i <= s:

!                                                       F
!                   I                1                   i-1
!   V  = A  *  ----------- * [ -------------  -  ------------------- +
!    i    i     4*PI*kapa        __________        ________________
!                        i     -| 2       2      -| 2             2
!                               |r + (z-d)        |r  +(z+d-2H   )
!                                                             i-1

!                      F
!                       s
!           + -------------------- ]
!               _________________
!             -| 2              2
!              |r  +(2H - d - z)
!                      s

!   za  i >= s:

!                                                       F
!                   I                1                   s-1
!   V  = A  *  ----------- * [ -------------  -  ------------------- +
!    i    i     4*PI*kapa        __________        ________________
!                        i     -| 2       2      -| 2             2
!                               |r + (z-d)        |r  +(z+d-2H   )
!                                                             s-1

!                      F
!                       i
!           + -------------------- ]
!               _________________
!             -| 2              2
!              |r  +(2H - d - z)
!                      i


!   gdje je:
!       s - redni broj sloja u kojem se nalzi izvor (iso)
!       F  - faktor refleksije na granici z = H  . Ovi su faktori 
!        i                                     i
!            izracunati u file-u: faktori_refleksije.f90.
!       A  - ukupni faktor transmisije od s-tog do i-tog sloja
!        i 
!            Ovi su faktori izracunati u file-u: faktori_transmisije.f90.
!       kapa  - kompleksna specificna elektricna vodljivost, [S/m]. Ove su
!           i
!               vrijednosti izracunate u file-u: proracun_kapa.f90.

!   Napomena:
!   ---------
!   Ukoliko se izvor nalazi u zraku (s=1) ili u n-tom sloju s=n, a raspodjela 
!   potencijala se trazi u zemlji (n-ti sloj ili bilo koji drugi) ili u zraku,
!   neki clanovi u gornjim izrazima ce nestati!

!================================================================================



subroutine potencijal_Vi(n,r,z,d,H,Fr,I,sloj,kapa_sloj,s,Ft_sloj,Vi)
    use funkcije
    implicit none
!   Input
    integer n
    real(8) r,z
    real(8) d
    real(8),dimension(:) :: H
    complex(8),dimension(:) :: Fr
    complex(8) I
    integer sloj
    complex(8) kapa_sloj
    integer s
    complex(8) Ft_sloj
!   Output
    complex(8) Vi
!   Local variables
    real(8),parameter :: pi = 3.141592653589793
    complex(8),parameter :: one = dcmplx(1.d0,0.d0)
    complex(8) temp1,temp2,temp3
    complex(8) potencijal


!   Opis ulaznih vrijednosti:
!       n - ukupni broj slojeva modela (zrak + viseslojno tlo)
!       r - r koordinata tocke u kojoj se racuna potencijal, [m]
!       z - z koordinata tocke u kojoj se racuna potencijal, [m]
!       d - dubina ukopavanja izvora harmonicke struje, [m]
!       H - vektor koji sadrzi koordinate donjih granicnih ploha svih slojeva
!       Fr - vektor faktora refleksije. Oderdjen je u file-u: faktori_refleksije.f90
!       I - jakost struje harmonickog tockastog izvora
!       sloj - sloj u kojem se racuna raspodjela potencijala (isl)
!       kapa_sloj - vrijednost specificne kompleksne vodljivosti i-tog sloja (sloja
!       u kojem se racuna raspodjela potencijala)
!       s - sloj u kojem se nalazi ukopan tockasti izvor struje (iso)
!       Ft_sloj - faktor transmisije sloja u kojem se racuna raspodjela potencijala


!   ------------------------------------------
!           Pocetak proracuna
!   ------------------------------------------

    if((s==1).and.(sloj==1)) then
        temp1 = one / dsqrt(r**2+(z-d)**2)
        temp3 = Fr(s) / dsqrt(r**2+(2.d0*H(s)-d-z)**2)
        potencijal = temp1 + temp3
    else if((s==1).and.(sloj/=1).and.(sloj/=n)) then
        temp1 = one / dsqrt(r**2+(z-d)**2)
        temp3 = Fr(sloj) / dsqrt(r**2+(2.d0*H(sloj)-d-z)**2)
        potencijal = temp1 + temp3
    else if((s==1).and.(sloj==n)) then
        temp1 = one / dsqrt(r**2+(z-d)**2)
        potencijal = temp1
    else if((s==n).and.(sloj==1)) then
        temp1 = one / dsqrt(r**2+(z-d)**2)
        potencijal = temp1
    else if((s==n).and.(sloj/=1).and.(sloj/=n)) then
        temp1 = one / dsqrt(r**2+(z-d)**2)
        temp2 = Fr(sloj-1) / dsqrt(r**2+(z+d-2.d0*H(sloj-1))**2)
        potencijal = temp1 - temp2
    else if((s==n).and.(sloj==n)) then
        temp1 = one / dsqrt(r**2+(z-d)**2)
        temp2 = Fr(s-1) / dsqrt(r**2+(z+d-2.d0*H(s-1))**2)
        potencijal = temp1 - temp2
    else if((s/=1).and.(s/=n).and.(sloj==1)) then
        temp1 = one / dsqrt(r**2+(z-d)**2)
        temp3 = Fr(s) / dsqrt(r**2+(2.d0*H(s)-d-z)**2)
        potencijal = temp1 + temp3
    else if((s/=1).and.(s/=n).and.(sloj/=1).and.(sloj/=n)) then
        if(sloj<=s) then
            temp1 = one / dsqrt(r**2+(z-d)**2)
            temp2 = Fr(sloj-1) / dsqrt(r**2+(z+d-2.d0*H(sloj-1))**2)
            temp3 = Fr(s) / dsqrt(r**2+(2.d0*H(s)-d-z)**2)
            potencijal = temp1 - temp2 + temp3
        else
            temp1 = one / dsqrt(r**2+(z-d)**2)
            temp2 = Fr(s-1) / dsqrt(r**2+(z+d-2.d0*H(s-1))**2)
            temp3 = Fr(sloj) / dsqrt(r**2+(2.d0*H(sloj)-d-z)**2)
            potencijal = temp1 - temp2 + temp3
        end if            
    else if((s/=1).and.(s/=n).and.(sloj==n)) then
        temp1 = one / dsqrt(r**2+(z-d)**2)
        temp2 = Fr(s-1) / dsqrt(r**2+(z+d-2.d0*H(s-1))**2)
        potencijal = temp1 - temp2
    end if

    Vi = ((Ft_sloj*I)/(4.d0*pi*kapa_sloj)) * potencijal

    return
end subroutine