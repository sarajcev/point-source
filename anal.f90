!   ===============================================================================
!   PROGRAM ZA ANALITICKI PRORACUN RASPODJELE POTENCIJALA USLJED TOCKASTOG
!   IZVORA IZMJENICNE STRUJE U VISESLOJU. OVAJ PROGRAM MOZE UZETI U OBZIR 
!   NAJVISE CETVEROSLOJNI MODEL (ZRAK + TROSLOJNO TLO).
!   ===============================================================================

program tiuv_analitic
    use funkcije
    implicit none

    integer n_sloj,n
    real(8) d
    complex(8) Is
    real(8) re,im,freq
    integer br_toc
    integer isl,iso
    integer,dimension(:),allocatable :: br_sloj
    real(8),dimension(:),allocatable :: svi_r,svi_z
    real(8) r,z
    real(8),dimension(:),allocatable :: epsr,ro
    real(8),dimension(:),allocatable :: h,HD
    complex(8),dimension(:),allocatable :: kapa,F
    complex(8),dimension(:),allocatable :: raspodjela
    complex(8) potencijal
    real(8) mod,kut,time_start,time_end
    real(8) ho
    integer p1,p2
    integer i

!   Opis varijabli:
!       n_sloj - ukupni broj slojeva modela viseslojnog tla (zrak + visesloj)
!       d - dubina na kojoj je ukopan izvor struje, [m]. Ukoliko je ova
!           vrijednost negativna, izvor se nalazi u zraku.
!       Is - jakost struje tockastog izvora (re + j*im), [A]
!       freq - rekvencija struje tockastog izvora, [Hz]
!       br_toc - broj tocaka u kojima ce se istovremeno racunati raspodjela 
!                potencijala
!       br_sloj - vektor koji sadrzi brojeve slojeva u kojima se nalaze tocke
!                 za koje se trazi raspodjela potencijala. Ukoliko se trazi 
!                 raspodjela potencijala u zraku (tocka je u zraku - negativna
!                 z koordinata) stavlja se vrijednost 1, za tocku u prvom sloju
!                 stavlja se vrijednost 2 i tako dalje, za tocku u i-tom sloju
!                 u vektoru br_sloj se stavlja vrijednost i+1!
!       svi_r - vektor duljine br_toc koji sadrzi r-koordinate tocaka za koje
!               se racuna potencijal
!       svi_z - vektor duljine br_toc koji sadrzi z-koordinate tocaka za koje
!               se racuna potencijal
!       epsr - relativne permitivnosti svih slojeva tla. Varijabla epsr je
!              realni vektor duljine n_sloj. Pritom epsr(1) do epsr(n) 
!              predstavlja relativnu permitivnost slojeva 1 do n modela tla.
!       ro - specifièna elektrièna otpornost svih slojeva modela tla [Ohm*m].
!            ro je realni vektor duljine n_sloj. Otpornost prvog sloja (zrak)
!            stavlja se da je jednaka nuli.
!       h - vektor duljine n_sloj koji sadrzi debljine svih slojeva n-slojnog
!           modela tla, [m]. Debljina n-tog sloja stavlja se da je nula.
!       HD - koordinate donjih granicnih ploha svih slojeva. HD(1)=0 i HD(n)=0.
!            Ovaj realni vektor ima n_sloj elemenata.
!       kapa - kompleksni vektor specificnih elektricnih vodljivosti svih slojeva
!              duljine n_sloj.
!       F - kompleksni vektor faktora refleksija za sve slojeve (duljine n_sloj).

!   =============================================================================
!                               Pocetak proracuna
!   =============================================================================

    !Citanje ulaznih podataka
    !------------------------
    open(1,file='input_an.txt')
    !Broj slojeva i dubina izvora
    read(1,*) n_sloj,d

    !Struja
    read(1,*) re,im,freq
    Is = dcmplx(re,im)

    !Vrijednosti epsr i ro za sve slojeve
    allocate(epsr(n_sloj))
    allocate(ro(n_sloj))
    do i = 1,n_sloj
        read(1,*) epsr(i),ro(i)
    end do

    !Vrijednosti h (debljine) svih slojeva
    allocate(h(n_sloj))
    do i = 1,n_sloj
        read(1,*) h(i)
    end do

    !Vrijednosti ho, p1 i p2
    read(1,*) ho,p1,p2
    !Ukupni broj slika koji se razmatra (za sume)
    read(1,*) n

    !Tocke u kojima se racuna raspodjela potencijala
    read(1,*) br_toc
    allocate(br_sloj(br_toc))
    allocate(svi_r(br_toc))
    allocate(svi_z(br_toc))
    do i = 1,br_toc
        read(1,*) br_sloj(i),svi_r(i),svi_z(i)
    end do
    close(1)


    !Tiskanje ulaznih podataka
    !-------------------------
    open(2,file='output_an.txt')
    write(2,'(/,"ANALITICKI PRORACUN POTENCIJALA U ODABRANIM TOCKAMA &
    USLJED",/,"TOCKASTOG IZVORA IZMJENICNE STRUJE U VISESLOJNOM SREDSTVU",/, &
    "(VISESLOJNI MODEL TLA PLUS ZRAK)")')
    write(2,'(//,57("-"))')
    write(2,'(tr20,"ULAZNI PODACI:")')
    write(2,'(57("-"),/)')
    write(2,'("Broj slojeva modela (zrak + viseslojno tlo):",i2)') n_sloj
    write(2,'(/,"Model zrak + viseslojno tlo:")')
    write(2,'("-----------------------------------")')
    write(2,'("Br.",7x,"Epsr",5x,"Ro",9x,"h")')
    write(2,'("sloja",3x,4x,4x,"[Ohm*m]",6x,"[m]")')
    write(2,'("-----------------------------------")')
    do i = 1,n_sloj
        write(2,'(i5,1x,f7.1,3x,f6.2,2x,f8.2)') i,epsr(i),ro(i),h(i)
    end do
    write(2,'("-----------------------------------")')
    write(2,'("Napomena:",/,&
    "Varijable u gornjoj tablici imaju sljedeca znacenja:",/,&
    "epsr - specificna permeabilnost sloja",/,&
    "ro   - specificna elektricna otpornost slojeva",/,&
    "h    - debljina sloja, u metrima.")')
    write(2,'(/,"Tockasti izvor struje:")')
    write(2,'("-----------------------------------")')
    write(2,'("Jakost struje: I =",f10.3," +j",f10.3," [A]")') real(Is),aimag(Is)
    write(2,'("Frekvencija: f =",f10.1," [Hz]")') freq
    write(2,'("Dubina ukopavanja tockastog izvora: d =",f6.1," [m]")') d
    write(2,'(/,"Broj tocaka za koje se racuna potencijal:",i3)') br_toc
    write(2,'("Ukupni broj slika:",i10)') n
    write(2,'("Vrijednosti: ho, p1 i p2:")')
    write(2,'("ho =",f10.2)') ho
    write(2,'("p1 =",i3,5x,"p2 =",i3)') p1,p2
    write(2,'(/,"Geometrija tocaka:")')
    write(2,'("-----------------------------------")')
    write(2,'("Sloj",10x,"r (m)",8x,"z(m)")')
    write(2,'("-----------------------------------")')
    do i = 1,br_toc
        write(2,'(i4,5x,f10.3,2x,f10.3,5x,i5,5x,i5)') br_sloj(i),svi_r(i),svi_z(i)
    end do
    write(2,'("-----------------------------------")')
    write(2,'("Napomena:",/,&
    "Varijable u gornjoj tablici imaju sljedeca znacenja:",/,&
    "sloj  - sloj u kojem se nalazi tocka promatranja",/,&
    "r(m)  - r koordinata doticne tocke, m",/,&
    "z(m)  - z koordinata doticne tocke, m.")')

!   =============================================================================
!               Pocetak proracuna analiticke raspodjele potencijala
!   =============================================================================

!   Proracun kompleksne specificne elektricne vodljivosti (kapa) svih slojeva
!   -----------------------------------------------------------------------------
!   Varijabla kapa je kompleksni vektor duljine n_sloj. kapa(1) je specificna elektr.
!   vodljivost zraka za koju vrijedi: kapa(1) = 0 + j*omega*eps0! kapa(2...n) 
!   predstavlja vodljivosti slojeva tla, [S/m].
    allocate(kapa(n_sloj))
    call proracun_kapa(n_sloj,freq,ro,epsr,kapa)

    deallocate(epsr)
    deallocate(ro)

    !Tiskanje vektora kapa
    !write(*,'("Vektor kapa:")')
    !do i = 1,n
    !   write(*,'(i3,d14.8,2x,d14.8)') i,dreal(kapa(i)),dimag(kapa(i))
    !end do


!   Proracun velicine HD za sve slojeve (z koord. donje granicne plohe)
!   -----------------------------------------------------------------------------
    allocate(HD(n_sloj))
    call vektor_H(n_sloj,h,HD)


!   Proracun faktora refleksije (F) za sve slojeve. 
!   -----------------------------------------------------------------------------
!   Faktor refleksije za sloj -1 i sloj n NE postoji. 
!   Vektor faktora refleksije je kompleksan, duljine n_sloj. 
    allocate(F(n_sloj))
    call vektor_F(n_sloj,kapa,F)

    !Tiskanje vektora refleksije
    !write(*,'("Faktori refleksije (F):")')
    !do i = 1,n
    !   write(*,'(i3,2x,d14.8,2x,d14.8)') i,dreal(F(i)),dimag(F(i))
    !end do

!   Odredjivanje sloja u kojem se nalazi tockasti izvor harmonicke struje
!   -----------------------------------------------------------------------------
!   Ukoliko je d - dubina na kojoj se nalazi izvor negativna vrijednost, rijec je
!   o tockastom izvoru koji se nalazi u zraku, te se pritom postavlja iso=1.
    iso = 1
    do i = 1,n_sloj-1
        if(d> HD(i)) then
            iso = iso + 1
        end if
    end do


!   -----------------------------------------------------------------------------
!           Analiticki proracun potencijala u odabranim tockama
!   ----------------------------------------------------------------------------
    call cpu_time(time_start)
    allocate(raspodjela(br_toc))
    do i = 1,br_toc
        isl = br_sloj(i)
        r = svi_r(i)
        z = svi_z(i)

        !Proracun potencijala u jednoj tocki
        !-----------------------------------
        call tiuvan(n_sloj,n,Is,d,isl,iso,r,z,p1,p2,ho,kapa,HD,F,potencijal)
        raspodjela(i) = potencijal

        !write(*,'("Press Enter to continue...")')
        !read(*,*)
    end do
    call cpu_time(time_end)

    !Tiskanje rezultata proracuna
    !----------------------------
    write(2,'(//,57("-"))')
    write(2,'(tr10,"REZULTATI PRORACUNA - ANALITICKI:")')
    write(2,'(57("-"),/)')
    write(2,'("Potencijali odabranih tocaka (Re + j*Im):",/)')
    write(2,'("    k",10x," Real",15x,"Imag")')
    write(2,'(50("-"))')
    do i = 1,br_toc
        write(2,'(i5,5x,d16.10,5x,d16.10)') i,dreal(raspodjela(i)),dimag(raspodjela(i))
    end do
    write(2,'(/,"Potencijali odabranih tocaka (Mod / kut°):",/)')
    write(2,'("    k",10x,"r [m]",10x,"z [m]",10x," Mod",15x,"Kut")')
    write(2,'(70("-"))')
    do i = 1,br_toc
        re = dreal(raspodjela(i))
        im = dimag(raspodjela(i))
        mod = dsqrt(re**2 + im**2)
        kut = datan2d(re,im)
        write(2,'(i5,5x,f10.2,5x,f10.2,5x,f16.6,5x,f8.3)') i,svi_r(i),svi_z(i),mod,kut
    end do
    write(2,'(/," Ukupno vrijeme trajanja proracuna:",f10.3," [s]")') time_end-time_start
    close(2)

    !Tiskanje rezultata analitickog proracuna - za usporedbu rezultata
    !-----------------------------------------------------------------
    open(1,file='out_an.txt')
    write(1,'(i5)') br_toc
        do i = 1,br_toc
        re = dreal(raspodjela(i))
        im = dimag(raspodjela(i))
        mod = dsqrt(re**2 + im**2)
        kut = datan2d(re,im)
        write(1,'(f16.6,5x,f8.3)') mod,kut
    end do
    close(1)
    
    !Oslobadjanje memorije
    !---------------------
    deallocate(h)
    deallocate(HD)
    deallocate(kapa)
    deallocate(F)
    deallocate(br_sloj)
    deallocate(svi_r)
    deallocate(svi_z)
    deallocate(raspodjela)

    stop
end program