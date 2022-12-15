!================================================================================

!                               GLAVNI PROGRAM

!   NUMERICKI PRORACUN POTENCIJALA USLJED TOCKASTOG IZVORA HARMONICKE STRUJE
!   U VISESLOJNOM SREDSTVU (VISESLOJNA ZEMLJA + ZRAK).

!================================================================================

program glavni_tiuv
    use funkcije
    implicit none

    integer br_toc
    integer n
    complex(8) Is
    real(8) freq,d
    integer,dimension(:),allocatable :: br_sloj
    real(8),dimension(:),allocatable :: svi_r,svi_z 
    real(8),dimension(:),allocatable :: epsr,ro
    real(8),dimension(:),allocatable :: H
    complex(8),dimension(:),allocatable :: raspodjela

    complex(8),dimension(:),allocatable :: kapa
    real(8),dimension(:),allocatable :: HD
    complex(8),dimension(:),allocatable :: F
    real(8) re,im,mod,kut
    integer isl,iso
    real(8) r,z,time_start,time_end
    complex(8) potencijal
    real(8),parameter :: pi = 3.14159265
    complex(8),parameter :: one = dcmplx(1.d0,0.d0)
    complex(8) faktor,temp
    integer i

!   Opis varijabli:
!       n - broj slojeva modela viseslojnog modela (zrak + viseslojno tlo)
!       d - dubina na kojoj je ukopan izvor struje, [m]. Ukoliko je ova
!           vrijednost negativna, izvor se nalazi u zraku!
!       Is - jakost struje tockastog izvora (re + j*im), [A]
!       f - rekvencija struje tockastog izvora, [Hz]
!       br_toc - broj tocaka u kojima ce se istovremeno racunati raspodjela 
!                potencijala
!       br_sloj - vektor koji sadrzi brojeve slojeva u kojim ase nalaze tocke
!                 za koje se trazi raspodjela potencijala. Ukoliko se trazi 
!                 raspodjela potencijala u zraku (tocka je u zraku - negativna
!                 z koordinata) stavlja se vrijednost 1, za tocku u prvom sloju
!                 stavlja se vrijednost 2 i tako dalje, za tocku u i-tom sloju
!                 u vektoru br_sloj se stavlja vrijednost i+1!
!       svi_r - vektor duljine br_toc koji sadrzi r-koordinate tocaka za koje
!               se racuna potencijal
!       svi_z - vektor duljine br_toc koji sadrzi z-koordinate tocaka za koje
!               se racuna potencijal. Za tocke u zraku ova vrijednost je negativna!
!       epsr - relativne permitivnosti svih slojeva tla. Varijabla epsr je
!              realni vektor duljine n. Pritom epsr(1) predstavlja relativnu 
!              permitivnost zraka, a epsr(2) do epsr(n) slojeva modela tla.
!       ro - specifièna elektrièna otpornost svih slojeva modela tla [Ohm*m].
!            ro je realni vektor duljine n.
!       h - vektor duljine n koji sadrzi debljine svih slojeva n-slojnog
!           modela (zrak + viseslojno tlo), [m].


!   =============================================================================
!                               Pocetak proracuna
!   =============================================================================

    open(1,file='input_num.txt')
    !Broj slojeva i dubina izvora
    read(1,*) n,d

    !Struja
    read(1,*) re,im,freq
    Is = dcmplx(re,im)

    !Vrijednosti epsr i ro za sve slojeve
    allocate(epsr(n))
    allocate(ro(n))
    do i = 1,n
        read(1,*) epsr(i),ro(i)
    end do

    !Vrijednosti h (debljine) svih slojeva
    allocate(H(n))
    do i = 1,n
        read(1,*) H(i)
    end do

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
    open(2,file='output_num.txt')
    write(2,'(/,"NUMERICKI PRORACUN POTENCIJALA U ODABRANIM TOCKAMA &
    USLJED",/,"TOCKASTOG IZVORA IZMJENICNE STRUJE U VISESLOJNOM SREDSTVU",/, &
    "(VISESLOJNI MODEL TLA PLUS ZRAK)")')
    write(2,'(//,57("-"))')
    write(2,'(tr20,"ULAZNI PODACI:")')
    write(2,'(57("-"),/)')
    write(2,'("Broj slojeva modela (zrak + viseslojno tlo):",i2)') n
    write(2,'(/,"Model zrak + viseslojno tlo:")')
    write(2,'("-----------------------------------")')
    write(2,'("Br.",7x,"Epsr",5x,"Ro",9x,"h")')
    write(2,'("sloja",3x,4x,4x,"[Ohm*m]",6x,"[m]")')
    write(2,'("-----------------------------------")')
    do i = 1,n
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
    write(2,'(/,"Geometrija tocaka:")')
    write(2,'("-----------------------------------")')
    write(2,'("Sloj",10x,"r (m)",8x,"z(m)")')
    write(2,'("-----------------------------------")')
    do i = 1,br_toc
        write(2,'(i4,5x,f10.3,2x,f10.3)') br_sloj(i),svi_r(i),svi_z(i)
    end do
    write(2,'("-----------------------------------")')
    write(2,'("Napomena:",/,&
    "Varijable u gornjoj tablici imaju sljedeca znacenja:",/,&
    "sloj - sloj u kojem se nalazi tocka promatranja",/,&
    "r(m) - r koordinata doticne tocke, m",/,&
    "z(m) - z koordinata doticne tocke, m.")')


!   =============================================================================
!                           Pocetak numerickog proracuna
!   =============================================================================

!   Proracun kompleksne specificne elektricne vodljivosto svih slojeva
!   -----------------------------------------------------------------------------
!   Varijabla kapa je kompleksni vektor duljine n. kapa(1) je specificna elektr.
!   vodljivost zraka za koju vrijedi: kapa(1) = 0 + j*omega*eps0! kapa(2...n) 
!   predstavlja vodljivosti slojeva tla, [S/m].
    allocate(kapa(n))
    call proracun_kapa(n,freq,ro,epsr,kapa)

    !Tiskanje vektora kapa
    !write(*,'("Vektor kapa:")')
    !do i = 1,n
    !   write(*,'(i3,d14.8,2x,d14.8)') i,dreal(kapa(i)),dimag(kapa(i))
    !end do


!   Proracun velicine H za sve slojeve (z koord. donje granicne plohe
!   -----------------------------------------------------------------------------
    allocate(HD(n))
    call vektor_H(n,H,HD)


!   Proracun faktora refleksije za sve slojeve. 
!   -----------------------------------------------------------------------------
!   za sloj -1 i sloj n faktor refleksije ne postoji. 
!   Vektor faktora refleksije je kompleksan, duljine n. 
    allocate(F(n))
    call vektor_F(n,kapa,F)

    !Tiskanje vektora refleksije
    !write(*,'("Faktori refleksije (F):")')
    !do i = 1,n
    !   write(*,'(i3,2x,d14.8,2x,d14.8)') i,dreal(F(i)),dimag(F(i))
    !end do


    IF (n==2) THEN
!       *************************************************************************
!                               ZRAK + HOMOGENO TLO
!       *************************************************************************
!       Odredjivanje sloja u kojem se nalazi tockasti izvor harmonicke struje
!       -------------------------------------------------------------------------
!       Ukoliko je d - dubina na kojoj se nalazi izvor negativna vrijednost, rijec
!       je o tockastom izvoru koji se nalazi u zraku, te se pritom postavlja iso=1.
        if (d<0.d0) then
            iso = 1
        else
            iso = 2
        end if

!       Proracun potencijala u svim odabranim tockama
!       -------------------------------------------------------------------------
        allocate(raspodjela(br_toc))
        do i = 1,br_toc
            isl = br_sloj(i)
            r = svi_r(i)
            z = svi_z(i)

!           Razmatranje svih slucajeva
            if ((iso==1).and.(isl==1)) then
                faktor = Is / (4.d0*pi*kapa(1))
                temp = one / dsqrt(r**2+(z-d)**2) + F(1) / dsqrt(r**2+(z+d)**2)
                raspodjela(i) = faktor * temp

            else if((iso==1).and.(isl==2)) then
                faktor = (Is / (4.d0*pi*kapa(2))) * (one - F(1))
                temp = one / dsqrt(r**2+(z-d)**2)
                raspodjela(i) = faktor * temp

            else if((iso==2).and.(isl==1)) then
                faktor = (Is / (4.d0*pi*kapa(1))) * (one + F(1))
                temp = one / dsqrt(r**2+(z-d)**2)
                raspodjela(i) = faktor * temp

            else if((iso==2).and.(isl==2)) then
                faktor = Is / (4.d0*pi*kapa(2))
                temp = one / dsqrt(r**2+(z-d)**2) - F(1) / dsqrt(r**2+(z+d)**2)
                raspodjela(i) = faktor * temp
            end if
        end do

        !Tiskanje rezultata proracuna
        !----------------------------
        write(2,'(//,57("-"))')
        write(2,'(tr10,"REZULTATI PRORACUNA - NUMERICKI:")')
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
        close(2)

        !Tiskanje rezultata numerickog proracuna - za usporedbu rezultata
        !-----------------------------------------------------------------
        open(1,file='out_num.txt')
        write(1,'(i5)') br_toc
            do i = 1,br_toc
            re = dreal(raspodjela(i))
            im = dimag(raspodjela(i))
            mod = dsqrt(re**2 + im**2)
            kut = datan2d(re,im)
            write(1,'(f16.6,5x,f8.3)') mod,kut
        end do
        close(1)

        deallocate(raspodjela)
    
    ELSE
!       *************************************************************************
!                               ZRAK + VISESLOJNO TLO
!       *************************************************************************
        call cpu_time(time_start)
        !Proracun potencijala u odabranim tockama
        !----------------------------------------
        allocate(raspodjela(br_toc))
        do i = 1,br_toc
            isl = br_sloj(i)
            r = svi_r(i)
            z = svi_z(i)

            !Proracun potencijala u jednoj tocki - za viseslojno sredstvo
            !============================================================
            call tiuv(n,freq,Is,d,isl,r,z,kapa,HD,H,F,potencijal)
            !============================================================ 
            raspodjela(i) = potencijal

            !write(*,'("Press Enter to continue...")')
            !read(*,*)
        end do
        call cpu_time(time_end)

        !Tiskanje rezultata proracuna
        !----------------------------
        write(2,'(//,57("-"))')
        write(2,'(tr10,"REZULTATI PRORACUNA - NUMERICKI:")')
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
        write(2,'(" Ukupno vrijeme trajanja proracuna:",f10.3," [s]")') time_end-time_start
        close(2)

        !Tiskanje rezultata numerickog proracuna - za usporedbu rezultata
        !-----------------------------------------------------------------
        open(1,file='out_num.txt')
        write(1,'(i5)') br_toc
            do i = 1,br_toc
            re = dreal(raspodjela(i))
            im = dimag(raspodjela(i))
            mod = dsqrt(re**2 + im**2)
            kut = datan2d(re,im)
            write(1,'(f16.6,5x,f8.3)') mod,kut
        end do
        close(1)

        deallocate(raspodjela)
    END IF


    !Oslobadjanje memorije
    deallocate(kapa)
    deallocate(F)
    deallocate(HD)
    deallocate(br_sloj)
    deallocate(svi_r,svi_z)

    stop
end program