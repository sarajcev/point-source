!================================================================================

!           Funkcija koja racuna vrijednosti kapa(i) za sve slojeve
!           n-slojnog modela. Kapa(i) je kompleksna specificna 
!           elektricna vodljivost i-tog sloja, [S/m].

!================================================================================

!Kompleksna specificna elektricna vodljivost i-tog sloja n-slojnog modela tla
!definirana je na sljedeci nacin:

!               kapa(i) = sigma(i) + j*omega*eps(i)

!   pri èemu su:
!       sigma(i) - specifièna elektrièna vodljivost i-tog sloja
!                  modela tla [S/m]. Vrijedi pritom sljedeca 
!                  jednakost:

!                                                1
!                                   sigma(i) = ------
!                                               ro(i)

!                  gdje ro(i) - predstavlja specificnu elektricnu
!                               otpornost i-tog sloja modela tla [Ohm/m].
!       omega - kruzna frekvencija. Jednaka je sljedecem izrazu:

!                                       omega = 2*pi*f

!               pri cemu f - predstavlja frekvenciju, [Hz].
!       eps(i) - relativna permitivnost i-tog sloja tla. Racuna se kao:

!                                       eps(i) = eps0 * epsr(i)

!                gdje su:
!                   eps0 - relativna permitivnost slobodnog prostora.
!                          Ona iznosi:

!                                       eps0 = 8.854e-12 [As/Vm]

!                   epsr(i) - relativna permitivnost i-tog sloja
!                             viseslojnog modela tla, [As/Vm].
!       j - imaginarna jedinica.

!   Napomena:
!   ---------
!   Varijabla kapa je kompleksni vektor duljine n. kapa(1) je specificna elektr.
!   vodljivost zraka. 
!   Varijabla epsr je realni vektor duljine n. Pritom epsr(1) do epsr(n) 
!   predstavlja relativnu permitivnost slojeva 1 do n modela.
!   Varijabla ro je realni vektor duljine n i predstavlja relativne elektricne 
!   otprornosti slojeva (1 do n) viseslojnog modela, [Ohm*m].
!   Varijable n; f; ro(i),i=1,n; epsr(i),i=1,n zadaju se kao ulazni podaci.
!================================================================================

subroutine proracun_kapa(n,f,ro,epsr,kapa)
    use funkcije
    implicit none
!   Input
    integer n
    real(8) f
    real(8),dimension(:) :: ro,epsr
!   Output
    complex(8),dimension(:) :: kapa
!   Local variables
    real(8),parameter :: pi = 3.14159265
    real(8),parameter :: eps0 = 8.854185e-12
    real(8) omega
    real(8) temp,sigma
    integer i


    !Tijelo subroutine
    omega = 2.d0*pi*f
    do i = 1,n
        if (i==1) then
            temp = omega*eps0*epsr(i)
            kapa(i) = dcmplx(0.d0,temp)
        else
            sigma = 1.d0/ro(i)
            temp = omega*eps0*epsr(i)
            kapa(i) = dcmplx(sigma,temp)
        end if
    end do

    return
end subroutine