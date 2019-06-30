! Q = (q - <q>) / sigma_q
! P = (p - <p>) / sigma_p
module SMD
    use Basic
    use Wigner
    implicit none

!Derived type
    !Example: trajSMD(i).Order(j).Array(m) is < Q^(j-m) * P^m / (j-m)! / m! > at snap shot i
    type SMDTrajectory
        type(d2PArray),allocatable,dimension(:)::Order
    end type SMDTrajectory

!Parameter
    integer::SMD_FollowStep=1000,&!Every how many steps print propagation
             SMD_OutputOrder=2!Output only 1 to SMD_OutputOrder SMD terms

!Global variable
    !Dynamics variable
    type(d2PArray),allocatable,dimension(:)::SMDquantity

!SMD module only variable
    !RK4 for SMD
    real*8,allocatable,dimension(:)::SMD_dtforall,SMD_dtd2,SMD_dtd6
    type(d2PArray),allocatable,dimension(:)::SMD_k1,SMD_k2,SMD_k3,SMD_k4,SMD_temp

contains
subroutine InitializeSMD()
    integer::i
    allocate(SMDquantity(SMDEvolutionOrder))
    do i=1,SMDEvolutionOrder
        allocate(SMDquantity(i).Array(0:i)); SMDquantity(i).Array=0d0
    end do
    !RK4 for SMD
    allocate(SMD_dtforall(SMDEvolutionOrder))
    allocate(SMD_dtd2(SMDEvolutionOrder)); allocate(SMD_dtd6(SMDEvolutionOrder))
    SMD_dtforall=dt; SMD_dtd2=dt/2d0; SMD_dtd6=dt/6d0
    allocate(SMD_k1(SMDEvolutionOrder)); allocate(SMD_k2(SMDEvolutionOrder))
    allocate(SMD_k3(SMDEvolutionOrder)); allocate(SMD_k4(SMDEvolutionOrder))
    do i=1,SMDEvolutionOrder
        allocate(SMD_k1(i).Array(0:i)); allocate(SMD_k2(i).Array(0:i))
        allocate(SMD_k3(i).Array(0:i)); allocate(SMD_k4(i).Array(0:i))
    end do
    allocate(SMD_temp(0:SMDOrder))
    do i=0,SMDOrder
        allocate(SMD_temp(i).Array(0:i))
    end do
end subroutine InitializeSMD

subroutine Dynamics()
    integer::OutputStep,TotalSteps,istep,i,j
    type(d2PArray),allocatable,dimension(:)::SMDquantityold
    !Data storage
        integer::TrajectoryLength
        type(SMDTrajectory),allocatable,dimension(:)::trajSMD
        real*8,allocatable,dimension(:)::trajpurity
    !Job control
        OutputStep=ceiling(OutputInterval/dt)
        dt=OutputInterval/dble(OutputStep)
        TrajectoryLength=ceiling(TotalTime/OutputInterval)
        TotalSteps=TrajectoryLength*OutputStep
    !Allocate data storage
        allocate(trajSMD(0:TrajectoryLength))!SMD
        do j=0,TrajectoryLength
            allocate(trajSMD(j).Order(SMD_OutputOrder))
            do i=1,SMD_OutputOrder
                allocate(trajSMD(j).Order(i).Array(0:i))
            end do
        end do
        forall(i=1:SMD_OutputOrder)
            trajSMD(0).Order(i).Array=SMDquantity(i).Array
        end forall
        allocate(trajwigcoeff(0:TrajectoryLength))!Wigner coefficient
        trajwigcoeff(0).NCentre=NCentre
        allocate(trajwigcoeff(0).centre(NCentre))
        do i=1,NCentre
            allocate(trajwigcoeff(0).centre(i).b(NLinearCoeffSC))
        end do
        trajwigcoeff(0).centre=wigcoeff(1:NCentre)
        allocate(trajpurity(0:TrajectoryLength))!Purity
        trajpurity(0)=purity()
    allocate(SMDquantityold(0:SMDOrder))!Allocate local work space
    do i=0,SMDOrder
        allocate(SMDquantityold(i).Array(0:i))
    end do
    SMDquantityold(0).Array(0)=1d0!Preloop
    forall(i=1:SMDEvolutionOrder)
        SMDquantityold(i).Array=SMDquantity(i).Array
    end forall
    write(*,*)'Total snap shots =',TotalSteps
    write(*,*)'Evolving...'
    do istep=1,TotalSteps
        call SMDRK4(SMDquantityold,SMDquantity,MoyalEOM)
        if(mod(istep,SMD_FollowStep)==0) then!Follow evolving progress
            call ShowTime()
            write(*,*)'evolve time =',istep*dt
            write(*,*)'<x> =',SMDquantityold(1).Array(0)
            write(*,*)'<p> =',SMDquantityold(1).Array(1)
        end if
        if(mod(istep,OutputStep)==0) then!Write trajectory
            j=istep/OutputStep
            forall(i=1:SMD_OutputOrder)!SMD
                trajSMD(j).Order(i).Array=SMDquantity(i).Array
            end forall
            trajwigcoeff(j).NCentre=NCentre!Wigner coefficient
            allocate(trajwigcoeff(j).centre(NCentre))
            do i=1,NCentre
                allocate(trajwigcoeff(j).centre(i).b(NLinearCoeffSC))
            end do
            trajwigcoeff(j).centre=wigcoeff(1:NCentre)
            trajpurity(j)=purity()!Purity
        end if
        do i=1,SMDEvolutionOrder!Check whether overflow
            do j=0,i
                if(isnan(SMDquantity(i).Array(j)).or.dAbs(SMDquantity(i).Array(j))>1d37) exit
            end do
            if(j<=i) exit
        end do
        if(i<=SMDEvolutionOrder) then!Overflow, stop evolution
            TrajectoryLength=istep/OutputStep
            write(*,*)'Overflow at',istep*dt
            exit
        end if
        forall(i=1:SMDEvolutionOrder)!Get ready for next loop
            SMDquantityold(i).Array=SMDquantity(i).Array
        end forall
    end do
    open(unit=99,file='t.out',status='replace')!output
        do i=0,TrajectoryLength
            write(99,*)dble(i)*OutputInterval
        end do
    close(99)
    open(unit=99,file='SMD.out',status='replace')
        do istep=0,TrajectoryLength
            do i=1,SMD_OutputOrder
                do j=0,i
                    write(99,*)trajSMD(istep).Order(i).Array(j)
                end do
            end do
        end do
    close(99)
    open(unit=99,file='WignerCoefficient.out',status='replace')
        do istep=0,TrajectoryLength
            write(99,*)trajwigcoeff(istep).NCentre
            do i=1,trajwigcoeff(istep).NCentre
                write(99,*)trajwigcoeff(istep).centre(i).miuq,trajwigcoeff(istep).centre(i).miup
                write(99,*)trajwigcoeff(istep).centre(i).sigmaq,trajwigcoeff(istep).centre(i).cor,trajwigcoeff(istep).centre(i).sigmap
                do j=1,NLinearCoeffSC
                    write(99,*)trajwigcoeff(istep).centre(i).b(j)
                end do
            end do
        end do
    close(99)
    open(unit=99,file='purity.out',status='replace')
        do i=0,TrajectoryLength
            write(99,*)trajpurity(i)
        end do
    close(99)
end subroutine Dynamics

subroutine MoyalEOM(d,u)!Moyal equation of motion
    type(d2PArray),dimension(SMDEvolutionOrder),intent(inout)::d
    type(d2PArray),dimension(0:SMDOrder),intent(inout)::u
    integer::i,order,m,n,j,jm2
    real*8::q,p,sigmaq,covQP,sigmap,varq,varp,spdsq,dlnsigmaq,dlnsigmap
    real*8::EdV,s,quantum,qtemp,ptemp
    real*8,dimension(0:ForceOrder)::a
    call CutOffScheme(u)
    call Polynomial(EdV,a,u)!Get V'(Q) - <V'>
    q=u(1).Array(0); p=u(1).Array(1); sigmaq=u(2).Array(0); covQP=u(2).Array(1); sigmap=u(2).Array(2)
    varq=sigmaq*sigmaq; varp=sigmap*sigmap; spdsq=sigmap/sigmaq
    !EOM for first 2 order terms
    d(1).Array(0)=u(1).Array(1)/mass
    d(1).Array(1)=-EdV
    d(2).Array(0)=covQP*sigmap/mass
    d(2).Array(1)=-a(1)
    d(2).Array(2)=-a(1)*covQP
    do i=2,ForceOrder!Dependence on higher order terms
        d(2).Array(1)=d(2).Array(1)-a(i)*(i+1)*u(i+1).Array(0)
        d(2).Array(2)=d(2).Array(2)-a(i)*u(i+1).Array(1)
    end do
    dlnsigmaq=d(2).Array(0)/sigmaq; dlnsigmap=d(2).Array(2)/sigmap
    d(2).Array(1)=spdsq/mass+d(2).Array(1)/sigmap-(dlnsigmaq+dlnsigmap)*covQP
    !EOM for 3 and higher order terms
    u(1).Array(0)=0; u(1).Array(1)=0!Transform the first 2 order terms into the uniform form of xi
    u(2).Array(0)=0.5d0; u(2).Array(2)=0.5d0
    do order=3,SMDEvolutionOrder
        !<Q^order> does not have 2nd term
        d(order).Array(0)=1d0/mass*spdsq*u(order).Array(1)-order*dlnsigmaq*u(order).Array(0)
        !Other quantities
        do n=1,order
            m=order-n
            if(n==order) then!<P^order> does not have 1st term
                d(order).Array(n)=-n*dlnsigmap*u(order).Array(n)
            else
                d(order).Array(n)=(n+1)/mass*spdsq*u(order).Array(n+1)-(m*dlnsigmaq+n*dlnsigmap)*u(order).Array(n)
            end if
            !2nd term: force dependence
            s=a(0)*u(order-1).Array(n-1)+a(1)*(m+1)*u(order).Array(n-1)
            do i=2,ForceOrder
                s=s+a(i)*cbn(m+i).Array(m)*u(order+i-1).Array(n-1)
            end do
            d(order).Array(n)=d(order).Array(n)-s/sigmap
            !4th term: quantum effect
            if(n>2) then
                quantum=0d0
                qtemp=varq; ptemp=varp
                do j=1,(n-1)/2
                    s=0d0
                    jm2=j*2
                    do i=jm2,ForceOrder
                        s=s+a(i)*cbn(m+i-jm2).Array(m)*u(m+i-jm2+n-(jm2+1)).Array(n-(jm2+1))
                    end do
                    quantum=quantum+(-1)**(j+1)*(hbar/2d0)**jm2/fct(jm2+1)/qtemp/ptemp*s
                    qtemp=qtemp*varq; ptemp=ptemp*varp
                end do
                d(order).Array(n)=d(order).Array(n)+quantum/sigmap
            end if
        end do
    end do
    u(1).Array(0)=q; u(1).Array(1)=p!Restore the old form of the first 2 order terms
    u(2).Array(0)=sigmaq; u(2).Array(2)=sigmap
end subroutine MoyalEOM

subroutine SMDRK4(old,new,f)!Runge Kutta 4 order, modified for SMD data type
    type(d2PArray),dimension(0:SMDOrder),intent(inout)::old
    type(d2PArray),dimension(SMDEvolutionOrder),intent(inout)::new
    external::f
    integer::i
    call f(SMD_k1,old)
    forall(i=1:SMDEvolutionOrder)
        SMD_temp(i).Array=old(i).Array+SMD_k1(i).Array*SMD_dtd2(i)
    end forall
    call f(SMD_k2,SMD_temp)
    forall(i=1:SMDEvolutionOrder)
        SMD_temp(i).Array=old(i).Array+SMD_k2(i).Array*SMD_dtd2(i)
    end forall
    call f(SMD_k3,SMD_temp)
    forall(i=1:SMDEvolutionOrder)
        SMD_temp(i).Array=old(i).Array+SMD_k3(i).Array*SMD_dtforall(i)
    end forall
    call f(SMD_k4,SMD_temp)
    forall(i=1:SMDEvolutionOrder)
        new(i).Array=old(i).Array+SMD_dtd6(i)*(SMD_k1(i).Array+2d0*SMD_k2(i).Array+2d0*SMD_k3(i).Array+SMD_k4(i).Array)
    end forall
end subroutine SMDRK4

end module SMD