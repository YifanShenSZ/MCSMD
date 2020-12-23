program Main
    use Basic
    use Wigner
    use SMD
    implicit none
    !Main only accessed input variable
    character*32::JobType
    real*8::q0,p0,sigmaq

    call ShowTime()
    call ReadInput()
    call Initialize()

    select case(JobType)
    case("SMD")
        call Dynamics()
    case("Wigner")
        call WignerDistributionPlot()
    case default
        stop "Unknown job type"
    end select

    call ShowTime()
    write(*,*)"Mission success"

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
        read(99,*); read(99,*)qleft
        read(99,*); read(99,*)qright
        read(99,*); read(99,*)dq
        read(99,*); read(99,*)pleft
        read(99,*); read(99,*)pright
        read(99,*); read(99,*)dp
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
    integer::istep,i,j
    real*8::sigmap
    call BetterRandomSeed()
    call InitializeBasic(); call InitializeWigner(); call InitializeSMD()
    select case(JobType)
        case('SMD')!Allocate dynamics work space, provide initial value
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
        case('Wigner')!Read old Wigner coefficient
            open(unit=99,file='t.out',status='old')
                j=0
                do
                    read(99,*,iostat=i); if(i/=0) exit
                    j=j+1
                end do
            close(99)
            allocate(trajwigcoeff(j))
            open(unit=99,file='WignerCoefficient.out',status='old')
                do istep=1,size(trajwigcoeff)
                    read(99,*)trajwigcoeff(istep).NCentre
                    allocate(trajwigcoeff(istep).centre(trajwigcoeff(istep).NCentre))
                    do i=1,trajwigcoeff(istep).NCentre
                        read(99,*)trajwigcoeff(istep).centre(i).miuq,trajwigcoeff(istep).centre(i).miup
                        read(99,*)trajwigcoeff(istep).centre(i).sigmaq,trajwigcoeff(istep).centre(i).cor,trajwigcoeff(istep).centre(i).sigmap
                        allocate(trajwigcoeff(istep).centre(i).b(NLinearCoeffSC))
                        do j=1,NLinearCoeffSC
                            read(99,*)trajwigcoeff(istep).centre(i).b(j)
                        end do
                    end do
                end do
            close(99)
        case default
    end select
end subroutine Initialize

end program Main