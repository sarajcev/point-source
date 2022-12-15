!================================================================================

!   Subroutina koja racuna ukupni potencijal u jednoj tocki visesloja, koji
!   stvara tockasti izvor harmonicke struje ukopan na odredjenoj dubini ili
!   smjesten u zraku.                                       

!================================================================================


subroutine tiuv(n,freq,Is,d,isl,r,z,kapa,HD,H,F,potencijal)
    use funkcije
    implicit none
!   Input
    integer n
    real(8) freq
    complex(8) Is
    real(8) d
    integer isl
    real(8) r,z
    complex(8),dimension(:) :: kapa
    real(8),dimension(:) :: HD,H
!   Output
    complex(8) potencijal
!   Local variables
    complex(8),dimension(:) :: F
    complex(8) A
    real(8),dimension(15) :: eta
    real(8),dimension(25) :: w
    integer iso
    real(8) teta,hi
    complex(8),dimension(15) :: alfa,beta
    real(8),dimension(25) :: u,lamda
    integer m_E,n_E
    real(8),dimension(15,25) :: E
!   complex(8) fi,Gi
!   complex(8),dimension(:),allocatable :: qi,ai,bi
    complex(8),dimension(25) :: tetan
    complex(8),dimension(25) :: teta_k,hi_k
    complex(8) kapa_sloj
    complex(8) Vi,Wi
    integer i,k

!   Opis ulaznih vrijednosti:
!       n - ukupni broj slojeva modela (zrak +  viselslojno tlo)
!       Is - struja tockastog izvora npr. 1000 + j0 [A]
!       freq - frekvencija struje, [Hz]
!       d - dubina na kojoj je ukopan tockasti izvor struje, [m]. Ukoliko se
!           za d unese negativna vrijednost, smatra se da je izvor u zraku!
!       isl - sloj u kojem se zeli racunati raspodjela potencijala
!       r,z - r i z koordinata tocke u kojoj se zeli racunati raspodjela
!             potencijala, [m]. Koordinatni sustav je cilindrican, sa z
!             koordinatom okrenutom u tlo, pri cemu je z=0 na povrsini tla.
!       epsr - relativne permitivnosti svih slojeva tla. Varijabla epsr je
!              realni vektor duljine n. Pritom epsr(1) predstavlja relativnu
!              permitivnost zraka a epsr(2) do epsr(n) slojeva tla.
!       ro - specifièna elektrièna otpornost svih slojeva modela, [Ohm*m].
!            ro je realni vektor duljine n. ro(1)=0 - zrak, ostlai su za tlo.
!       h_sloj - vektor duljine n koji sadrzi debljine svih slojeva n-slojnog
!                modela tla. h_sloj(1)=0 i h_sloj(n)=0
!   Opis lokalnih varijabli:
!       iso - sloj u kojem se nalazi tockasti izvor struje.
!       F - faktori refleksije za sve slojeve
!       A - faktori transmisije za i-ti sloj
!       eta, w - nepromjenjivi parametri eta i w, repsektivno
!       teta, hi - parametri teta i hi koji se biraju u file-u: parametri.f90
!                  ovisno o broju slojeva, polozaju izvora struje i sl.
!       alfa, beta - koeficijenti za aproksimaciju jezgra-funkcija. Racunaju se
!                    metodom najmanjih kvadrata, mnozenjem s pseudoinverznom
!                    matricom [E]. Koriste se pri proracunu raspodjele potenc.
!       u - vrijednosti 25 uzorkovanih tocaka
!       E - pseudoinverzna matrica. Racuna se u fileu: Pseudoinvezna_matrica.for
!       fi, Gi - elementi fi i Gi koji sluze za
!                dobivanje raspodjele jezgra-funkcija teta i hi
!       tetan - vektor vrijednosti teta_n za n-ti sloj u 25 uzorkovanih tocaka
!               Na osnovu ovog vektora rekurzivno se racunaju vektori 
!               teta i hi (jezgra-funkcije), koje se poslije aproksimiraju
!               pomocu alfa i beta sa 15 exponenc. funkcija.
!       teta_k - vektor elemenata (jezgra funkcije) teta. 
!       hi_k - vektor elemenata (jezgra funkcije) hi. 
!       kapa_sloj - vrijednost kapa za sloj u kojem se racuna raspodjela
!                   potencijala (isl)


!   =============================================================================
!                               Pocetak proracuna
!   =============================================================================

!   Odredjivanje sloja u kojem se nalazi tockasti izvor harmonicke struje
!   -----------------------------------------------------------------------------
!   Ukoliko je d - dubina na kojoj se nalazi izvor negativna vrijednost, rijec je
!   o tockastom izvoru koji se nalazi u zraku, te se pritom postavlja iso=1.
    iso = 1
    do i = 1,n-1
        if(d> HD(i)) then
            iso = iso + 1
        end if
    end do


!   Odredjivanje nepromjenjivih parametara eta i w
!   -----------------------------------------------------------------------------
!   Postoji ukupno 15 parametara eta. eta je stoga realni vektor duljine 15.
!   Postoji ukupno 25 nepromjenjivih parametara  w. w je stoga realni
!   vektor duljine 25 (posljednji clan koji je beskonacan - stavljen je nula
!   jer se kasnije korigira).
    call vektor_eta(n,eta)
    call vektor_w(n,w)


!   Odredjivanje parametara odslikavanja u visesloju: teta i hi
!   -----------------------------------------------------------------------------
    call parametri_teta_hi(iso,isl,d,HD,H,n,teta,hi)


!   Odredjivanje nepoznatih koeficijenata alfa i beta 
!   ako je sloj=1 (isl=1) - zrak ili ako je sloj=n (isl=n) - n-ti sloj
!   -----------------------------------------------------------------------------
    if (isl==1) then
        alfa(:) = dcmplx(0.0,0.0)
    else if(isl==n) then
        beta(:) = dcmplx(0.0,0.0)
    end if


!   -----------------------------------------------------------------------------
!   Odredjivanje nepoznatih koeficijenata *alfa* ako je sloj/=1
!   -----------------------------------------------------------------------------
    if(isl/=1) then
        !Odredjivanje uzorkovanih tocaka
        do k = 1,25
            u(k) = teta * w(k)
        end do

        !Formiranje matrice [E]
        m_E = 15
        n_E = 25
        call invmat(E,m_E,n_E)

        !Odredjivanje vrijednosti lamda
        do k = 1,24
            lamda(k) = 1.d0/u(k)
        end do
        lamda(25) = 0.d0

        !Odredjivanje vrijednosti teta_n u 25 uzorkovanih tocaka
        call vektor_teta_n(n,d,H,HD,lamda,kapa,F,isl,iso,tetan)

        !Odredjivanje vektora teta (pomocu rekurzivnih formula)
        !u 25 uzorkovanih tocaka 
        call vektor_teta_hi(n,d,tetan,kapa,lamda,isl,iso,F,H,HD,teta_k,hi_k)

        !Odredjivanje nepoznatih koef. alfa mnozenjem s pseudoinv. matricom
        alfa = matmul(E,teta_k)
    end if


!   -----------------------------------------------------------------------------
!   Odredjivanje nepoznatih koeficijenata *beta* ako je sloj/=n
!   -----------------------------------------------------------------------------
    if(isl/=n) then
        !Odredjivanje uzorkovanih tocaka
        do k = 1,25
            u(k) = hi * w(k)
        end do

        !Formiranje matrice [E]
        m_E = 15
        n_E = 25
        call invmat(E,m_E,n_E)

        !Odredjivanje vrijednosti lamda
        do k = 1,24
            lamda(k) = 1.d0/u(k)
        end do
        lamda(25) = 0.d0

        !Odredjivanje vrijednosti teta_n u 25 uzorkovanih tocaka
        call vektor_teta_n(n,d,H,HD,lamda,kapa,F,isl,iso,tetan)

        !Odredjivanje vektora hi (pomocu rekurzivnih formula)
        !u 25 uzorkovanih tocaka za sve slojeve
        call vektor_teta_hi(n,d,tetan,kapa,lamda,isl,iso,F,H,HD,teta_k,hi_k)

        !Odredjivanje nepoznatih koef. beta mnozenjem s pseudoinv. matricom
        beta = matmul(E,hi_k)
    end if


!   Odredjivanje faktora transmisije i vrijednosti kapa za sloj u kojem se racuna
!   raspodjela potencijala Vi i Wi (sloj)
!   -----------------------------------------------------------------------------
    A = vektor_A(iso,isl,n,F)
    kapa_sloj = kapa(isl)


!   Proracun velicine Vi (koji je dio ukupnog potencijala)
!   -----------------------------------------------------------------------------
    call potencijal_Vi(n,r,z,d,HD,F,Is,isl,kapa_sloj,iso,A,Vi)


!   Proracun velicine Wi (koji je dio ukupnog potencijala)
!   -----------------------------------------------------------------------------
    call potencijal_Wi(n,Is,isl,kapa_sloj,alfa,beta,r,z,HD,teta,hi,eta,Wi)


!   Odredjivanje ukupnog potencijala u zadanoj tocki
!   -----------------------------------------------------------------------------
    potencijal = Vi + Wi


    return
end subroutine