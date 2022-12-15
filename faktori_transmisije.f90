!================================================================================

!           Funkcija koja racuna faktor transmisije A za onaj sloj
!           u kojem se trazi raspodjela potencijala.

!================================================================================


!Faktori transmisije definiraju su na sljedeci nacin:

!       A = (1 + F(s-1))*(1 + F(s-2))*...*(1 + F(i)) ;  za: i < s 

!       A = 1 :  za: i = s

!       A = (1 - F(s))*(1 - F(s+1))*...*(1-F(i-1)) ;  za: i > s 

!pri cemu su:
!   i - trenutna vrijednost sloja za koju se racuna faktor transmisije (isl)
!   s - sloj u kojem se nalazi tockasti izvor struje (iso)
!   F - faktori refleksije, koji su prethodno odredeni. Vektor F je kompleksan
!       i sadrzi n elemenata, pri cemu je n - unupni broj slojeva modela.
!================================================================================


complex(8) function vektor_A(iso,isl,n,F)
!   use funkcije
    implicit none

!   Input
    integer iso,isl,n
    complex(8),dimension(:) :: F
!   Local variables
!   complex(8) prod
    complex(8),parameter :: one = dcmplx(1.d0,0.d0)
    complex(8) A
    integer i
!   Opis varijabli:
!   iso - sloj u kojem se nalazi tockasti izvor struje (s)
!   isl - sloj u kojem se trazi raspodjela potencijala (sloj)


    A = one

    if (isl<iso) then
        do i = isl,iso-1
            A = A * (one + F(i))
        end do

    else if(isl==iso) then
        A = one

    else if (isl>iso) then
        do i = iso,isl-1
            A = A * (one - F(i))
        end do
    end if
    
    vektor_A = A

    return
end function