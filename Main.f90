program Main
    use Basic
    use Wigner
    use SMD
    implicit none
    !Main only accessed input variable
        character*32::JobType
        real*8::q0,p0,sigmaq
!---------- Initialize ----------
    call ShowTime(); call ReadInput(); call Initialize()
!------------- End --------------

!----------- Run job ------------
    select case(JobType)
        case('SMD')
            call Dynamics()
        case('Wigner')
        case default!Throw a warning
            write(*,*)'Program abort: unsupported job type '//trim(adjustl(JobType))
            stop
    end select
!------------- End --------------

!---------- Clean up ------------
    call ShowTime(); write(*,*)'Mission complete'
!------------- End --------------

contains
subroutine ReadInput()
    logical::advance
    open(unit=99,file='SMD.in')
        read(99,*); read(99,*); read(99,*); read(99,*)JobType
        read(99,*); read(99,*)mass
        read(99,*); read(99,*)q0
        read(99,*); read(99,*)p0
        read(99,*); read(99,*)sigmaq
        read(99,*); read(99,*)TotalTime
        read(99,*); read(99,*)dt
        read(99,*); read(99,*)OutputInterval
        read(99,*); read(99,*)SMDEvolutionOrder
        read(99,*); read(99,*)BasisOrder
        read(99,*); read(99,*)PotentialOrder; allocate(PolyPotentialCoeff(PotentialOrder))
        read(99,*); read(99,*)PolyPotentialCoeff
        read(99,*); read(99,*)advance
    close(99)
    if(advance) then
        write(*,*)'Advanced input requested, parameters are set to user specification'
        open(unit=99,file='advance.in',status='old')
            namelist /AdvancedInput/ &
                MaxNCentre,MaxPopDev,MaxImpurity,MaxTRIteration,&
                SMD_FollowStep,SMD_OutputOrder
            read(99,nml=AdvancedInput)
        close(99)
    end if
end subroutine ReadInput

subroutine Initialize()
    integer::i,j
    real*8::sigmap
    call BetterRandomSeed()
    call InitializeBasic(); call InitializeWigner(); call InitializeSMD()
    if(JobType=='SMD') then!Allocate dynamics work space, provide initial value
        !Initial condition: gaussian wave packet
        !Initial time-dependently evolved SMD quantity: accordingly
        sigmap=0.5d0/sigmaq
        SMDquantity(1).Array(0)=q0; SMDquantity(1).Array(1)=p0
        SMDquantity(2).Array(0)=sigmaq; SMDquantity(2).Array(2)=sigmap
        do i=4,SMDEvolutionOrder
            do j=0,i,2
                SMDquantity(i).Array(j)=fct2(i-j-1)*fct2(j-1)/fct(i-j)/fct(j)
            end do
        end do
        !Initial Wigner distribution: accordingly
        NCentre=1
        wigcoeff(1).miuq=q0; wigcoeff(1).miup=p0
        wigcoeff(1).sigmaq=sigmaq; wigcoeff(1).sigmap=sigmap; wigcoeff(1).cor=0d0
        wigcoeff(1).b(1)=1d0; wigcoeff(1).b(2:NLinearCoeffSC)=0d0
        InitialPurity=1d0
    end if
end subroutine Initialize

end program Main