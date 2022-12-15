!================================================================================

!   Subroutina koja racuna ukupni potencijal u jednoj tocki visesloja, koji
!   stvara tockasti izvor harmonicke struje ukopan na odredjenoj dubini ili
!   smjesten u zraku. Potencijal se racuna uz pomoc analitickih formula.                                        

!================================================================================


subroutine tiuvan(n_sloj,n,Is,d,isl,iso,r,z,p1,p2,ho,kapa,HD,F,potencijal)
    use funkcije
    implicit none

    integer n_sloj
    complex(8) Is
    real(8) d
    integer isl,iso
    real(8) r,z
    complex(8),dimension(:) :: kapa
    real(8),dimension(:) :: HD
    complex(8),dimension(:) :: F
    complex(8) potencijal
    complex(8),dimension(:),allocatable :: Cn
    integer n
    real(8) ho
    integer p1,p2
    real(8) h1,h2
    complex(8) Fi1,Fi2,Fi3,Fi4
    real(8),parameter :: pi = 3.14159265
    complex(8),parameter :: one = dcmplx(1.d0,0.d0)
    complex(8) faktor
    complex(8) temp,temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8
    complex(8) sum,sum1,sum2
    integer k1,k2,k3
    complex(8) C1,C2,C3
    integer i

!   Opis varijabli:
!       n_sloj - ukupni broj slojeva modela viseslojnog tla (zrak + visesloj)
!       d - dubina na kojoj je ukopan izvor struje, [m]. Ukoliko je ova
!           vrijednost negativna, izvor se nalazi u zraku.
!       Is - jakost struje tockastog izvora (re + j*im), [A]
!       freq - rekvencija struje tockastog izvora, [Hz]
!       isl - sloj u kojem se nalazi tocka za koju se racuna potencijal
!       r,z - r i z koordinata doticne tocke, [m].
!       h - vektor duljine n_sloj koji sadrzi debljine svih slojeva n-slojnog
!           modela (zrak + viseslojno tlo), [m]. Debljina n-tog sloja stavlja
!           se da je nula.
!       HD - koordinate donjih granicnih ploha svih slojeva. HD(1)=0 i HD(n)=0.
!            Ovaj realni vektor ima n_sloj elemenata.
!       kapa - kompleksni vektor specificnih elektricnih vodljivosti svih slojeva
!              duljine n_sloj.
!       F - kompleksni vektor faktora refleksija za sve slojeve (duljine n_sloj).
!       iso - sloj u kojem se nalazi izvor harmonicke struje
!       Cn - kompleksni vektor koeficijenata Cn - ukupno n, pri cemu je n ukupni
!            broj slika koje se razmatraju.
!       n - ukupni braoj slika koje se razmatraju. 



!   -----------------------------------------------------------------------------
!           Analiticki proracun raspodjele potencijala za jednu tocku
!   -----------------------------------------------------------------------------

    SELECT CASE(n_sloj)
        CASE(2)     
!           *************************************
!           Dvoslojni model (zrak + homogeno tlo)
!           *************************************
            select case(iso)
                case(1)
!                   ====================================
!                   Izvor se nalazi u prvom sloju (zrak)
!                   ====================================
                    select case(isl)
                        case(1)
!                           -------------------------------------------
!                           Raspodjela potencijala u prvom sloju (zrak)
!                           -------------------------------------------
                            Fi1 = dcmplx(0.d0,0.d0)
                            faktor = Is / (4.d0*pi*kapa(1))
                            Fi1 = one / dsqrt(r**2+(z-d)**2) + F(1) / dsqrt(r**2+(-z+dabs(d))**2)
                            Fi1 = faktor * Fi1
                            potencijal = Fi1
                        case(2)
!                           -------------------------------------------
!                           Raspodjela potencijala u drugom sloju (tlo)
!                           -------------------------------------------
                            Fi2 = dcmplx(0.d0,0.d0)
                            faktor = Is / (4.d0*pi*kapa(2))
                            Fi2 = (one - F(1)) * (one / dsqrt(r**2+(z+dabs(d))**2))
                            Fi2 = faktor * Fi2
                            potencijal = Fi2
                    end select
                case(2)
!                   ====================================
!                   Izvor se nalazi u drugom sloju (tlo)
!                   ====================================
                    select case(isl)
                        case(1)
!                           -------------------------------------------
!                           Raspodjela potencijala u prvom sloju (zrak)
!                           -------------------------------------------
                            Fi1 = dcmplx(0.d0,0.d0)
                            faktor = Is / (4.d0*pi*kapa(1))
                            Fi1 = (one + F(1)) * (one / dsqrt(r**2+(z-d)**2))
                            Fi1 = faktor * Fi1
                            potencijal = Fi1
                        case(2)
!                           -------------------------------------------
!                           Raspodjela potencijala u drugom sloju (tlo)
!                           -------------------------------------------
                            Fi2 = dcmplx(0.d0,0.d0)
                            faktor = Is / (4.d0*pi*kapa(2))
                            Fi2 = one / dsqrt(r**2+(z-d)**2) - F(1) / dsqrt(r**2+(z+d)**2)
                            Fi2 = faktor * Fi2
                            potencijal = Fi2
                    end select
            end select
        CASE(3)
!           **************************************
!           Troslojni model (zrak + dvoslojno tlo)
!           **************************************
!           Proracun koeficijenata Cn
!           -------------------------
            allocate(Cn(n))
            do i = 1,n
                Cn(i) = (-1.d0)**i * F(1)**i * F(2)**i
            end do
!           -------------------------
            select case(iso)
                case(1)
!                   ====================================
!                   Izvor se nalazi u prvom sloju (zrak)
!                   ====================================
                    select case(isl)
                        case(1)
!                           -------------------------------------------
!                           Raspodjela potencijala u prvom sloju (zrak)
!                           -------------------------------------------
                            Fi1 = dcmplx(0.d0,0.d0)
                            faktor = Is / (4.d0*pi*kapa(1))
                            temp = one / dsqrt(r**2+(z-d)**2) + F(1) / dsqrt(r**2+(-z+dabs(d))**2)
                            sum = dcmplx(0.d0,0.d0)
                            do i = 1,n
                                sum = sum + Cn(i) * (one / dsqrt(r**2+(2.d0*i*HD(2)-z+dabs(d))**2))
                            end do
                            sum = ((F(1)**2-one)/F(1)) * sum
                            Fi1 = faktor * (temp + sum)
                            potencijal = Fi1
                        case(2)
!                           -------------------------------------------------------
!                           Raspodjela potencijala u drugom sloju (prvi sloj u tlu)
!                           -------------------------------------------------------
                            Fi2 = dcmplx(0.d0,0.d0)
                            faktor = (Is / (4.d0*pi*kapa(2))) * (one - F(1))
                            temp = one / dsqrt(r**2+(z+dabs(d))**2)
                            sum = dcmplx(0.d0,0.d0)
                            do i = 1,n
                                temp1 = one / dsqrt(r**2+(2.d0*i*HD(2)+z+dabs(d))**2)
                                temp2 = (one/F(1)) * (one / dsqrt(r**2+(2.d0*i*HD(2)-z+dabs(d))**2))
                                sum = sum + Cn(i) * (temp1 - temp2)
                            end do
                            Fi2 = faktor * (temp + sum)
                            potencijal = Fi2
                        case(3)
!                           --------------------------------------------------------
!                           Raspodjela potencijala u trecem sloju (drugi sloj u tlu)
!                           --------------------------------------------------------
                            Fi3 = dcmplx(0.d0,0.d0)
                            faktor = (Is / (4.d0*pi*kapa(3))) * (one - F(1)) * (one - F(2))
                            sum = dcmplx(0.d0,0.d0)
                            do i = 0,n
                                if (i==0) then
                                    sum = sum + one / dsqrt(r**2+(2.d0*i*HD(2)+z+dabs(d))**2)
                                else
                                    sum = sum + Cn(i) * (one / dsqrt(r**2+(2.d0*i*HD(2)+z+dabs(d))**2))
                                end if
                            end do
                            Fi3 = faktor * sum
                            potencijal = Fi3
                    end select
                case(2)
!                   ================================================
!                   Izvor se nalazi u drugom sloju (prvi sloj u tlu)
!                   ================================================
                    select case(isl)
                        case(1)
!                           -------------------------------------------
!                           Raspodjela potencijala u prvom sloju (zrak)
!                           -------------------------------------------
                            Fi1 = dcmplx(0.d0,0.d0)
                            faktor = (Is / (4.d0*pi*kapa(1))) * (one / (one - F(1)))
                            temp = one / dsqrt(r**2+(z-d)**2) - F(1) / dsqrt(r**2+(z+d)**2)
                            sum1 = dcmplx(0.d0,0.d0)
                            sum2 = dcmplx(0.d0,0.d0)
                            do i = 0,n
                                if (i==0) then
                                    temp1 = F(1) / dsqrt(r**2+(2.d0*i*HD(2)-z-d)**2)
                                    temp2 = F(1)**2 / dsqrt(r**2+(2.d0*i*HD(2)-z+d)**2)
                                    sum1 = sum1 + (temp1 - temp2)
                                else
                                    temp1 = F(1) / dsqrt(r**2+(2.d0*i*HD(2)-z-d)**2)
                                    temp2 = F(1)**2 / dsqrt(r**2+(2.d0*i*HD(2)-z+d)**2)
                                    sum1 = sum1 + Cn(i) * (temp1 - temp2)
                                end if
                                if (i==0) then
                                    sum2 = dcmplx(0.d0,0.d0)
                                else
                                    temp3 = -(one/F(1)) * (one / dsqrt(r**2+(2.d0*i*HD(2)-z-d)**2))
                                    temp4 = one / dsqrt(r**2+(2.d0*i*HD(2)-z+d)**2)
                                    sum2 = sum2 + Cn(i) * (temp3 + temp4)
                                end if
                            end do
                            Fi1 = faktor * (temp + sum1 + sum2)
                            potencijal = Fi1
                        case(2)
!                           -------------------------------------------------------
!                           Raspodjela potencijala u drugom sloju (prvi sloj u tlu)
!                           -------------------------------------------------------
                            Fi2 = dcmplx(0.d0,0.d0)
                            faktor = Is / (4.d0*pi*kapa(2))
                            temp = one / dsqrt(r**2+(z-d)**2) - F(1) / dsqrt(r**2+(z+d)**2)
                            sum = dcmplx(0.d0,0.d0)
                            do i = 1,n
                                temp1 = one / dsqrt(r**2+(2.d0*i*HD(2)+z-d)**2)
                                temp2 = F(1) / dsqrt(r**2+(2.d0*i*HD(2)+z+d)**2)
                                temp3 = (one/F(1)) * (one / dsqrt(r**2+(2.d0*i*HD(2)-z-d)**2))
                                temp4 = one / dsqrt(r**2+(2.d0*i*HD(2)-z+d)**2)
                                sum = sum + Cn(i) * (temp1 - temp2 - temp3 + temp4)
                            end do
                            Fi2 = faktor * (temp + sum)
                            potencijal = Fi2
                        case(3)
!                           --------------------------------------------------------
!                           Raspodjela potencijala u trecem sloju (drugi sloj u tlu)
!                           --------------------------------------------------------
                            Fi3 = dcmplx(0.d0,0.d0)
                            faktor = (Is / (4.d0*pi*kapa(3))) * (one - F(2))
                            sum = dcmplx(0.d0,0.d0)
                            do i = 0,n
                                if (i==0) then
                                    temp1 = one / dsqrt(r**2+(2.d0*i*HD(2)+z-d)**2)
                                    temp2 = F(1) / dsqrt(r**2+(2.d0*i*HD(2)+z+d)**2)
                                    sum = sum + (temp1 - temp2)
                                else
                                    temp1 = one / dsqrt(r**2+(2.d0*i*HD(2)+z-d)**2)
                                    temp2 = F(1) / dsqrt(r**2+(2.d0*i*HD(2)+z+d)**2)
                                    sum = sum + Cn(i) * (temp1 - temp2)
                                end if
                            end do
                            Fi3 = faktor * sum
                            potencijal = Fi3
                    end select
                case(3)
!                   =================================================
!                   Izvor se nalazi u trecem sloju (drugi sloj u tlu)
!                   =================================================
                    select case(isl)
                        case(1)
!                           -------------------------------------------
!                           Raspodjela potencijala u prvom sloju (zrak)
!                           -------------------------------------------
                            Fi1 = dcmplx(0.d0,0.d0)
                            faktor = (Is / (4.d0*pi*kapa(1))) * (one + F(1)) * (one + F(2))
                            sum = dcmplx(0.d0,0.d0)
                            do i = 0,n
                                if (i==0) then
                                    temp = one / dsqrt(r**2+(2.d0*i*HD(2)-z+d)**2)
                                    sum = sum + temp
                                else
                                    temp = one / dsqrt(r**2+(2.d0*i*HD(2)-z+d)**2)
                                    sum = sum + Cn(i) * temp
                                end if
                            end do
                            Fi1 = faktor * sum
                            potencijal = Fi1
                        case(2)
!                           -------------------------------------------------------
!                           Raspodjela potencijala u drugom sloju (prvi sloj u tlu)
!                           -------------------------------------------------------
                            Fi2 = dcmplx(0.d0,0.d0)
                            faktor = (Is / (4.d0*pi*kapa(2))) * (one + F(2))
                            sum = dcmplx(0.d0,0.d0)
                            do i = 0,n
                                if (i==0) then
                                    temp1 = one / dsqrt(r**2+(2.d0*i*HD(2)-z+d)**2)
                                    temp2 = F(1) / dsqrt(r**2+(2.d0*i*HD(2)+z+d)**2)
                                    sum = sum + (temp1 - temp2)
                                else
                                    temp1 = one / dsqrt(r**2+(2.d0*i*HD(2)-z+d)**2)
                                    temp2 = F(1) / dsqrt(r**2+(2.d0*i*HD(2)+z+d)**2)
                                    sum = sum + Cn(i) * (temp1 - temp2)
                                end if
                            end do
                            Fi2 = faktor * sum
                            potencijal = Fi2
                        case(3)
!                           --------------------------------------------------------
!                           Raspodjela potencijala u trecem sloju (drugi sloj u tlu)
!                           --------------------------------------------------------
                            Fi3 = dcmplx(0.d0,0.d0)
                            faktor = Is / (4.d0*pi*kapa(3))
                            temp = one / dsqrt(r**2+(z-d)**2) - F(1) / dsqrt(r**2+(z+d)**2)
                            sum = dcmplx(0.d0,0.d0)
                            do i = 1,n
                                temp1 = - F(1) / dsqrt(r**2+(2.d0*i*HD(2)+z+d)**2)
                                temp2 = (one/F(1)) * (one / dsqrt(r**2+(2.d0*i*HD(2)-4.d0*HD(2)+z+d)**2))
                                sum = sum + Cn(i) * (temp1 + temp2)
                            end do
                            Fi3 = faktor * (temp + sum)
                            potencijal = Fi3
                    end select
                !Oslobadjanje vektora Cn
                deallocate(Cn)
            end select
        CASE(4)
!           ******************************************
!           Cetveroslojni model (zrak + troslojno tlo)
!           ******************************************
!           Proracun koeficijenata Cn
!           -------------------------
            h1 = ho * p1
            h2 = ho * p2
            allocate(Cn(n))
            do i = 1,n
                k1 = i-p1
                if(k1<0) then
                    C1 = dcmplx(0.d0,0.d0)
                else if(k1==0) then
                    C1 = one
                else
                    C1 = Cn(k1)
                end if
                k2 = i-p1-p2
                if(k2<0) then
                    C2 = dcmplx(0.d0,0.d0)
                else if(k2==0) then
                    C2 = one
                else
                    C2 = Cn(k2)
                end if
                k3 = i-p2
                if(k3<0) then
                    C3 = dcmplx(0.d0,0.d0)
                else if(k3==0) then
                    C3 = one
                else
                    C3 = Cn(k3)
                end if
                Cn(i) = -F(2)*F(1)*C1 - F(3)*F(1)*C2 - F(2)*F(3)*C3
            end do
!           -------------------------
            select case(iso)
                case(1)
!                   ====================================
!                   Izvor se nalazi u prvom sloju (zrak)
!                   ====================================
                    select case(isl)
                        case(1)
!                           -------------------------------------------
!                           Raspodjela potencijala u prvom sloju (zrak)
!                           -------------------------------------------
                            Fi1 = dcmplx(0.d0,0.d0)
                            faktor = Is / (4.d0*pi*kapa(1))
                            temp = one / dsqrt(r**2+(z-d)**2)
                            sum = dcmplx(0.d0,0.d0)
                            do i = 0,n
                                if (i==0) then
                                    temp1 = F(1) / dsqrt(r**2+(2.d0*i*ho-z+dabs(d))**2)
                                    temp2 = (F(1)*F(2)*F(3)) / dsqrt(r**2+(2.d0*i*ho+2.d0*h2-z+dabs(d))**2)
                                    temp3 = F(2) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(2)-z+dabs(d))**2)
                                    temp4 = F(3) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(3)-z+dabs(d))**2)
                                    sum = sum + (temp1 + temp2 + temp3 + temp4)
                                else
                                    temp1 = F(1) / dsqrt(r**2+(2.d0*i*ho-z+dabs(d))**2)
                                    temp2 = (F(1)*F(2)*F(3)) / dsqrt(r**2+(2.d0*i*ho+2.d0*h2-z+dabs(d))**2)
                                    temp3 = F(2) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(2)-z+dabs(d))**2)
                                    temp4 = F(3) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(3)-z+dabs(d))**2)
                                    sum = sum + Cn(i) * (temp1 + temp2 + temp3 + temp4)
                                end if
                            end do
                            Fi1 = faktor * (temp + sum)
                            potencijal = Fi1
                        case(2)
!                           -------------------------------------------------------
!                           Raspodjela potencijala u drugom sloju (prvi sloj u tlu)
!                           -------------------------------------------------------
                            Fi2 = dcmplx(0.d0,0.d0)
                            faktor = (Is / (4.d0*pi*kapa(2))) * (one - F(1))
                            sum = dcmplx(0.d0,0.d0)
                            do i = 0,n
                                if (i==0) then
                                    temp1 = one / dsqrt(r**2+(2.d0*i*ho+z+dabs(d))**2)
                                    temp2 = (F(2)*F(3)) / dsqrt(r**2+(2.d0*i*ho+2.d0*h2+z+dabs(d))**2)
                                    temp3 = F(2) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(2)-z+dabs(d))**2)
                                    temp4 = F(3) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(3)-z+dabs(d))**2)
                                    sum = sum + (temp1 + temp2 + temp3 + temp4)
                                else
                                    temp1 = one / dsqrt(r**2+(2.d0*i*ho+z+dabs(d))**2)
                                    temp2 = (F(2)*F(3)) / dsqrt(r**2+(2.d0*i*ho+2.d0*h2+z+dabs(d))**2)
                                    temp3 = F(2) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(2)-z+dabs(d))**2)
                                    temp4 = F(3) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(3)-z+dabs(d))**2)
                                    sum = sum + Cn(i) * (temp1 + temp2 + temp3 + temp4)
                                end if
                            end do
                            Fi2 = faktor * sum
                            potencijal = Fi2
                        case(3)
!                           --------------------------------------------------------
!                           Raspodjela potencijala u trecem sloju (drugi sloj u tlu)
!                           --------------------------------------------------------
                            Fi3 = dcmplx(0.d0,0.d0)
                            faktor = (Is / (4.d0*pi*kapa(3))) * (one - F(1)) * (one - F(2))
                            sum = dcmplx(0.d0,0.d0)
                            do i = 0,n
                                if (i==0) then
                                    temp1 = one / dsqrt(r**2+(2.d0*i*ho+z+dabs(d))**2)
                                    temp2 = F(3) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(3)-z+dabs(d))**2)
                                    sum = sum + (temp1 + temp2)
                                else
                                    temp1 = one / dsqrt(r**2+(2.d0*i*ho+z+dabs(d))**2)
                                    temp2 = F(3) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(3)-z+dabs(d))**2)
                                    sum = sum + Cn(i) * (temp1 + temp2)
                                end if
                            end do
                            Fi3 = faktor * sum
                            potencijal = Fi3
                        case(4)
!                           ----------------------------------------------------------
!                           Raspodjela potencijala u cetvrtom sloju (treci sloj u tlu)
!                           ----------------------------------------------------------
                            Fi4 = dcmplx(0.d0,0.d0)
                            faktor = (Is / (4.d0*pi*kapa(4))) * (one - F(1)) * (one - F(2)) * (one - F(3))
                            sum = dcmplx(0.d0,0.d0)
                            do i = 0,n
                                if (i==0) then
                                    temp = one / dsqrt(r**2+(2.d0*i*ho+z+dabs(d))**2)
                                    sum = sum + temp
                                else
                                    temp = one / dsqrt(r**2+(2.d0*i*ho+z+dabs(d))**2)
                                    sum = sum + Cn(i) * temp
                                end if
                            end do
                            Fi4 = faktor * sum
                            potencijal = Fi4
                    end select
                case(2)
!                   ================================================
!                   Izvor se nalazi u drugom sloju (prvi sloj u tlu)
!                   ================================================
                    select case(isl)
                        case(1)
!                           -------------------------------------------
!                           Raspodjela potencijala u prvom sloju (zrak)
!                           -------------------------------------------
                            Fi1 = dcmplx(0.d0,0.d0)
                            faktor = (Is / (4.d0*pi*kapa(1))) * (one + F(1))
                            temp = one / dsqrt(r**2+(-z+d)**2)
                            sum = dcmplx(0.d0,0.d0)
                            do i = 0,n
                                if (i==0) then
                                    temp1 = F(2) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(2)-z-d)**2)
                                    temp2 = F(3) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(3)-z-d)**2)
                                    temp3 = (F(1)*F(2)) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(2)-z+d)**2)
                                    temp4 = (F(1)*F(3)) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(3)-z+d)**2)
                                    sum = sum + (temp1 + temp2 - temp3 - temp4)
                                else
                                    temp1 = F(2) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(2)-z-d)**2)
                                    temp2 = F(3) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(3)-z-d)**2)
                                    temp3 = (F(1)*F(2)) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(2)-z+d)**2)
                                    temp4 = (F(1)*F(3)) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(3)-z+d)**2)
                                    sum = sum + Cn(i) * (temp1 + temp2 - temp3 - temp4)
                                end if
                            end do
                            Fi1 = faktor * (temp + sum)
                            potencijal = Fi1
                        case(2)
!                           -------------------------------------------------------
!                           Raspodjela potencijala u drugom sloju (prvi sloj u tlu)
!                           -------------------------------------------------------
                            Fi2 = dcmplx(0.d0,0.d0)
                            faktor = Is / (4.d0*pi*kapa(2))
                            temp = one / dsqrt(r**2+(z-d)**2) - F(1) / dsqrt(r**2+(z+d)**2)
                            sum = dcmplx(0.d0,0.d0)
                            do i = 0,n
                                if (i==0) then
                                    temp1 = (F(1)*F(2)) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(2)+z-d)**2)
                                    temp2 = (F(1)**2*F(2)) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(2)+z+d)**2)
                                    temp3 = (F(1)*F(3)) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(3)+z-d)**2)
                                    temp4 = (F(1)**2*F(3)) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(3)+z+d)**2)
                                    temp5 = F(2) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(2)-z-d)**2)
                                    temp6 = (F(1)*F(2)) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(2)-z+d)**2)
                                    temp7 = F(3) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(3)-z-d)**2)
                                    temp8 = (F(1)*F(3)) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(3)-z+d)**2)
                                    sum = sum + (-temp1 + temp2 - temp3 + temp4 + temp5 - temp6 + temp7 - temp8)
                                else
                                    temp1 = (F(1)*F(2)) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(2)+z-d)**2)
                                    temp2 = (F(1)**2*F(2)) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(2)+z+d)**2)
                                    temp3 = (F(1)*F(3)) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(3)+z-d)**2)
                                    temp4 = (F(1)**2*F(3)) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(3)+z+d)**2)
                                    temp5 = F(2) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(2)-z-d)**2)
                                    temp6 = (F(1)*F(2)) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(2)-z+d)**2)
                                    temp7 = F(3) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(3)-z-d)**2)
                                    temp8 = (F(1)*F(3)) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(3)-z+d)**2)
                                    sum = sum + Cn(i) * (-temp1 + temp2 - temp3 + temp4 + temp5 - temp6 + temp7 - temp8)
                                end if
                            end do
                            Fi2 = faktor * (temp + sum)
                            potencijal = Fi2
                        case(3)
!                           --------------------------------------------------------
!                           Raspodjela potencijala u trecem sloju (drugi sloj u tlu)
!                           --------------------------------------------------------
                            Fi3 = dcmplx(0.d0,0.d0)
                            faktor = (Is / (4.d0*pi*kapa(3))) * (one - F(2))
                            sum = dcmplx(0.d0,0.d0)
                            do i = 0,n
                                if (i==0) then
                                    temp1 = one / dsqrt(r**2+(2.d0*i*ho+z-d)**2)
                                    temp2 = F(1) / dsqrt(r**2+(2.d0*i*ho+z+d)**2)
                                    temp3 = F(3) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(3)-z-d)**2)
                                    temp4 = (F(1)*F(3)) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(3)-z+d)**2)
                                    sum = sum + (temp1 - temp2 + temp3 - temp4)
                                else
                                    temp1 = one / dsqrt(r**2+(2.d0*i*ho+z-d)**2)
                                    temp2 = F(1) / dsqrt(r**2+(2.d0*i*ho+z+d)**2)
                                    temp3 = F(3) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(3)-z-d)**2)
                                    temp4 = (F(1)*F(3)) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(3)-z+d)**2)
                                    sum = sum + Cn(i) * (temp1 - temp2 + temp3 - temp4)
                                end if
                            end do
                            Fi3 = faktor * sum
                            potencijal = Fi3
                        case(4)
!                           ----------------------------------------------------------
!                           Raspodjela potencijala u cetvrtom sloju (treci sloj u tlu)
!                           ----------------------------------------------------------
                            Fi4 = dcmplx(0.d0,0.d0)
                            faktor = (Is / (4.d0*pi*kapa(4))) * (one - F(2)) * (one - F(3))
                            sum = dcmplx(0.d0,0.d0)
                            do i = 0,n
                                if (i==0) then
                                    temp1 = one / dsqrt(r**2+(2.d0*i*ho+z-d)**2)
                                    temp2 = F(1) / dsqrt(r**2+(2.d0*i*ho+z+d)**2)
                                    sum = sum + (temp1 - temp2)
                                else
                                    temp1 = one / dsqrt(r**2+(2.d0*i*ho+z-d)**2)
                                    temp2 = F(1) / dsqrt(r**2+(2.d0*i*ho+z+d)**2)
                                    sum = sum + Cn(i) * (temp1 - temp2)
                                end if
                            end do
                            Fi4 = faktor * sum
                            potencijal = Fi4
                    end select
                case(3)
!                   =================================================
!                   Izvor se nalazi u trecem sloju (drugi sloj u tlu)
!                   =================================================
                    select case(isl)
                        case(1)
!                           -------------------------------------------
!                           Raspodjela potencijala u prvom sloju (zrak)
!                           -------------------------------------------
                            Fi1 = dcmplx(0.d0,0.d0)
                            faktor = (Is / (4.d0*pi*kapa(1))) * (one + F(1)) * (one + F(2))
                            sum = dcmplx(0.d0,0.d0)
                            do i = 0,n
                                if (i==0) then
                                    temp1 = one / dsqrt(r**2+(2.d0*i*ho-z+d)**2)
                                    temp2 = F(3) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(3)-z-d)**2)
                                    sum = sum + (temp1 + temp2)
                                else
                                    temp1 = one / dsqrt(r**2+(2.d0*i*ho-z+d)**2)
                                    temp2 = F(3) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(3)-z-d)**2)
                                    sum = sum + Cn(i) * (temp1 + temp2)
                                end if
                            end do
                            Fi1 = faktor * sum
                            potencijal = Fi1
                        case(2)
!                           -------------------------------------------------------
!                           Raspodjela potencijala u drugom sloju (prvi sloj u tlu)
!                           -------------------------------------------------------
                            Fi2 = dcmplx(0.d0,0.d0)
                            faktor = (Is / (4.d0*pi*kapa(2))) * (one + F(2))
                            sum = dcmplx(0.d0,0.d0)
                            do i = 0,n
                                if (i==0) then
                                    temp1 = -F(1) / dsqrt(r**2+(2.d0*i*ho+z+d)**2)
                                    temp2 = (F(1)*F(3)) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(3)+z-d)**2)
                                    temp3 = one / dsqrt(r**2+(2.d0*i*ho-z+d)**2)
                                    temp4 = F(3) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(3)-z-d)**2)
                                    sum = sum + (temp1 - temp2 + temp3 + temp4)
                                else
                                    temp1 = -F(1) / dsqrt(r**2+(2.d0*i*ho+z+d)**2)
                                    temp2 = (F(1)*F(3)) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(3)+z-d)**2)
                                    temp3 = one / dsqrt(r**2+(2.d0*i*ho-z+d)**2)
                                    temp4 = F(3) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(3)-z-d)**2)
                                    sum = sum + Cn(i) * (temp1 - temp2 + temp3 + temp4)
                                end if
                            end do
                            Fi2 = faktor * sum
                            potencijal = Fi2
                        case(3)
!                           --------------------------------------------------------
!                           Raspodjela potencijala u trecem sloju (drugi sloj u tlu)
!                           --------------------------------------------------------
                            Fi3 = dcmplx(0.d0,0.d0)
                            faktor = Is / (4.d0*pi*kapa(3))
                            temp = one / dsqrt(r**2+(z-d)**2) + F(3) / dsqrt(r**2+(2.d0*HD(3)-z-d)**2)
                            sum = dcmplx(0.d0,0.d0)
                            do i = 0,n
                                if (i==0) then
                                    temp1 = F(1) / dsqrt(r**2+(2.d0*i*ho+z+d)**2)
                                    temp2 = F(2) / dsqrt(r**2+(2.d0*i*ho-2.d0*HD(2)+z+d)**2)
                                    temp3 = (F(1)*F(3)) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(3)+z-d)**2)
                                    temp4 = (F(2)*F(3)) / dsqrt(r**2+(2.d0*i*ho+2.d0*h2+z-d)**2)
                                    temp5 = (F(1)*F(3)) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(3)-z+d)**2)
                                    temp6 = (F(2)*F(3)) / dsqrt(r**2+(2.d0*i*ho+2.d0*h2-z+d)**2)
                                    temp7 = (F(1)*F(3)**2) / dsqrt(r**2+(2.d0*i*ho+4.0*HD(3)-z-d)**2)
                                    temp8 = (F(2)*F(3)**2) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(3)+2.d0*h2-z-d)**2)
                                    sum = sum + (-temp1 - temp2 - temp3 - temp4 - temp5 - temp6 - temp7 - temp8)
                                else
                                    temp1 = F(1) / dsqrt(r**2+(2.d0*i*ho+z+d)**2)
                                    temp2 = F(2) / dsqrt(r**2+(2.d0*i*ho-2.d0*HD(2)+z+d)**2)
                                    temp3 = (F(1)*F(3)) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(3)+z-d)**2)
                                    temp4 = (F(2)*F(3)) / dsqrt(r**2+(2.d0*i*ho+2.d0*h2+z-d)**2)
                                    temp5 = (F(1)*F(3)) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(3)-z+d)**2)
                                    temp6 = (F(2)*F(3)) / dsqrt(r**2+(2.d0*i*ho+2.d0*h2-z+d)**2)
                                    temp7 = (F(1)*F(3)**2) / dsqrt(r**2+(2.d0*i*ho+4.0*HD(3)-z-d)**2)
                                    temp8 = (F(2)*F(3)**2) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(3)+2.d0*h2-z-d)**2)
                                    sum = sum + Cn(i) * (-temp1 - temp2 - temp3 - temp4 - temp5 - temp6 - temp7 - temp8)
                                end if
                            end do
                            Fi3 = faktor * (temp + sum)
                            potencijal = Fi3
                        case(4)
!                           ----------------------------------------------------------
!                           Raspodjela potencijala u cetvrtom sloju (treci sloj u tlu)
!                           ----------------------------------------------------------
                            Fi4 = dcmplx(0.d0,0.d0)
                            faktor = (Is / (4.d0*pi*kapa(4))) * (one - F(3))
                            sum = dcmplx(0.d0,0.d0)
                            do i = 0,n
                                if (i==0) then
                                    temp1 = F(1) / dsqrt(r**2+(2.d0*i*ho+z+d)**2)
                                    temp2 = F(2) / dsqrt(r**2+(2.d0*i*ho-2.d0*HD(2)+z+d)**2)
                                    temp3 = one / dsqrt(r**2+(2.d0*i*ho+z-d)**2)
                                    temp4 = (F(1)*F(2)) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(2)+z-d)**2)
                                    sum = sum + (-temp1 - temp2 + temp3 + temp4)
                                else
                                    temp1 = F(1) / dsqrt(r**2+(2.d0*i*ho+z+d)**2)
                                    temp2 = F(2) / dsqrt(r**2+(2.d0*i*ho-2.d0*HD(2)+z+d)**2)
                                    temp3 = one / dsqrt(r**2+(2.d0*i*ho+z-d)**2)
                                    temp4 = (F(1)*F(2)) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(2)+z-d)**2)
                                    sum = sum + Cn(i) * (-temp1 - temp2 + temp3 + temp4)
                                end if
                            end do
                            Fi4 = faktor * sum
                            potencijal = Fi4
                    end select
                case(4)
!                   ===================================================
!                   Izvor se nalazi u cetvrtom sloju (treci sloj u tlu)
!                   ===================================================
                    select case(isl)
                        case(1)
!                           -------------------------------------------
!                           Raspodjela potencijala u prvom sloju (zrak)
!                           -------------------------------------------
                            Fi1 = dcmplx(0.d0,0.d0)
                            faktor = (Is / (4.d0*pi*kapa(1))) * (one + F(1)) * (one + F(2)) * (one + F(3))
                            sum = dcmplx(0.d0,0.d0)
                            do i = 0,n
                                if (i==0) then
                                    temp = one / dsqrt(r**2+(2.d0*i*ho-z+d)**2)
                                    sum = sum + temp
                                else
                                    temp = one / dsqrt(r**2+(2.d0*i*ho-z+d)**2)
                                    sum = sum + Cn(i) * temp
                                end if
                            end do
                            Fi1 = faktor * sum
                            potencijal = Fi1
                        case(2)
!                           -------------------------------------------------------
!                           Raspodjela potencijala u drugom sloju (prvi sloj u tlu)
!                           -------------------------------------------------------
                            Fi2 = dcmplx(0.d0,0.d0)
                            faktor = (Is / (4.d0*pi*kapa(2))) * (one + F(2)) * (one + F(3))
                            sum = dcmplx(0.d0,0.d0)
                            do i = 0,n
                                if (i==0) then
                                    temp1 = one / dsqrt(r**2+(2.d0*i*ho+z+d)**2)
                                    temp2 = F(1) / dsqrt(r**2+(2.d0*i*ho-z+d)**2)
                                    sum = sum + (temp1 - temp2)
                                else
                                    temp1 = one / dsqrt(r**2+(2.d0*i*ho+z+d)**2)
                                    temp2 = F(1) / dsqrt(r**2+(2.d0*i*ho-z+d)**2)
                                    sum = sum + Cn(i) * (temp1 - temp2)
                                end if
                            end do
                            Fi2 = faktor * sum
                            potencijal = Fi2
                        case(3)
!                           --------------------------------------------------------
!                           Raspodjela potencijala u trecem sloju (drugi sloj u tlu)
!                           --------------------------------------------------------
                            Fi3 = dcmplx(0.d0,0.d0)
                            faktor = (Is / (4.d0*pi*kapa(3))) * (one + F(3))
                            sum = dcmplx(0.d0,0.d0)
                            do i = 0,n
                                if (i==0) then
                                    temp1 = F(1) / dsqrt(r**2+(2.d0*i*ho+z+d)**2)
                                    temp2 = F(2) / dsqrt(r**2+(2.d0*i*ho-2.d0*HD(2)+z+d)**2)
                                    temp3 = one / dsqrt(r**2+(2.d0*i*ho-z+d)**2)
                                    temp4 = (F(1)*F(2)) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(2)-z+d)**2)
                                    sum = sum + (-temp1 - temp2 + temp3 + temp4)
                                else
                                    temp1 = F(1) / dsqrt(r**2+(2.d0*i*ho+z+d)**2)
                                    temp2 = F(2) / dsqrt(r**2+(2.d0*i*ho-2.d0*HD(2)+z+d)**2)
                                    temp3 = one / dsqrt(r**2+(2.d0*i*ho-z+d)**2)
                                    temp4 = (F(1)*F(2)) / dsqrt(r**2+(2.d0*i*ho+2.d0*HD(2)-z+d)**2)
                                    sum = sum + Cn(i) * (-temp1 - temp2 + temp3 + temp4)
                                end if
                            end do
                            Fi3 = faktor * sum
                            potencijal = Fi3
                        case(4)
!                           ----------------------------------------------------------
!                           Raspodjela potencijala u cetvrtom sloju (treci sloj u tlu)
!                           ----------------------------------------------------------
                            Fi4 = dcmplx(0.d0,0.d0)
                            faktor = Is / (4.d0*pi*kapa(4))
                            temp = one / dsqrt(r**2+(z-d)**2)
                            sum = dcmplx(0.d0,0.d0)
                            do i = 0,n
                                if (i==0) then
                                    temp1 = F(1) / dsqrt(r**2+(2.d0*i*ho+z+d)**2)
                                    temp2 = (F(1)*F(2)*F(3)) / dsqrt(r**2+(2.d0*i*ho-2.d0*h2+z+d)**2)
                                    temp3 = F(2) / dsqrt(r**2+(2.d0*i*ho-2.d0*HD(2)+z+d)**2)
                                    temp4 = F(3) / dsqrt(r**2+(2.d0*i*ho-2.d0*HD(3)+z+d)**2)
                                    sum = sum + (-temp1 - temp2 - temp3 - temp4)
                                else
                                    temp1 = F(1) / dsqrt(r**2+(2.d0*i*ho+z+d)**2)
                                    temp2 = (F(1)*F(2)*F(3)) / dsqrt(r**2+(2.d0*i*ho-2.d0*h2+z+d)**2)
                                    temp3 = F(2) / dsqrt(r**2+(2.d0*i*ho-2.d0*HD(2)+z+d)**2)
                                    temp4 = F(3) / dsqrt(r**2+(2.d0*i*ho-2.d0*HD(3)+z+d)**2)
                                    sum = sum + Cn(i) * (-temp1 - temp2 - temp3 - temp4)
                                end if
                            end do
                            Fi4 = faktor * (temp + sum)
                            potencijal = Fi4
                    end select
                deallocate(Cn)
            end select
    END SELECT

    return
end subroutine