! alpha = (Q + P) / Sqrt(2)
! beta  = (Q - P) / Sqrt(2)
module Wigner
    use Basic
    implicit none

!Derived type
    type WignerCoefficient
        real*8::miuq,miup,sigmaq,sigmap,cor
        real*8,allocatable,dimension(:)::b
    end type WignerCoefficient

!Parameter:
    integer::MaxNCentre=1!Max number of centres allowed
    real*8::MaxPopDev=1d-6,&!Maximum population deviation from 1
            MaxImpurity=1d-6!Maximum purity deviation from initial purity
    integer::MaxTRIteration=1000!Max number of trust region iteration

!Global variable
    integer::NCentre
    type(WignerCoefficient),allocatable,dimension(:)::wigcoeff
    real*8::InitialPurity

!Wigner module only variable
    integer::NLinearCoeffSC!Number of linear coefficients for a single centre
    !Work space for evaluating SMD quantity and fitting Wigner distribution
    real*8,allocatable,dimension(:)::SIGMA,SIGMAk,SIGMAkAlphaBeta,dSIGMAk
    real*8,allocatable,dimension(:,:)::U_Coeff,U_SMDEvolution,U_SMD,M_SMDEvolution,M_SMD,Wk_SMDEvolution,Wk_SMD,dWk
    type(d2PMatrix),allocatable,dimension(:,:)::BlockTemp
    !Purity work space
    type(d2PMatrix),allocatable,dimension(:)::Ak
    type(d2PArray),allocatable,dimension(:)::Bk
    real*8,allocatable,dimension(:)::Ck

contains
subroutine InitializeWigner()
    integer::i,j,m,n,index,indexmrow
    real*8::s
    type(d2PMatrix),allocatable,dimension(:)::Utemp
    type(d2PMatrix),allocatable,dimension(:,:)::Mtemp
    NLinearCoeffSC=(1+BasisOrder+1)*(BasisOrder+1)/2
    !Global variable
    allocate(wigcoeff(MaxNCentre))
    do i=1,MaxNCentre
        allocate(wigcoeff(i).b(NLinearCoeffSC))
    end do
    !Work space for evaluating SMD quantity and fitting Wigner distribution
    allocate(SIGMA(NSMDQuantity))
    allocate(SIGMAk(NSMDQuantity))
    allocate(SIGMAkAlphaBeta(NSMDQuantity))
    SIGMA(1)=1d0
    SIGMAk(1)=1d0
    SIGMAkAlphaBeta(1)=1d0
    allocate(dSIGMAk(NSMDQuantityEvolution))
    dSIGMAk(1)=0d0
    !Canonical transformation matrix from (q-miu_q)/sigma_q,(p-miu_p)/sigma_p to alpha,beta
        allocate(Utemp(2:max(SMDOrder,BasisOrder)))!Prepare
        do j=2,max(SMDOrder,BasisOrder)
            allocate(Utemp(j).Matrix(0:j,0:j))
            do n=0,j
                do m=0,j
                    s=0d0
                    do i=max(0,m+n-j),min(m,n)
                        s=s+cbn(j-m).Array(n-i)*cbn(m).Array(i)*(-1)**i
                    end do
                    Utemp(j).Matrix(m,n)=fct(j-n)*fct(n)/fct(j-m)/fct(m)/(1.4142135623730951d0**j)*s
                end do
            end do
        end do
        allocate(U_Coeff(NLinearCoeffSC,NLinearCoeffSC))!Transform single centre linear coefficient vector
            U_Coeff=0d0
            U_Coeff(1,1)=1d0!0th order
            if(BasisOrder>0) then
                U_Coeff(2,2)= 1d0/1.4142135623730951d0!1st order
                U_Coeff(3,2)= 1d0/1.4142135623730951d0
                U_Coeff(2,3)= 1d0/1.4142135623730951d0
                U_Coeff(3,3)=-1d0/1.4142135623730951d0
                index=4
                do i=2,BasisOrder
                    U_Coeff(index:index+i,index:index+i)=Utemp(i).Matrix
                    index=index+i+1
                end do
            end if
        allocate(U_SMDEvolution(NSMDQuantityEvolution,NSMDQuantityEvolution))!Transform SMD quantities up to SMDEvolutionOrder
            U_SMDEvolution=0d0
            U_SMDEvolution(1,1)=1d0!0th order
            U_SMDEvolution(2,2)= 1d0/1.4142135623730951d0!1st order
            U_SMDEvolution(3,2)= 1d0/1.4142135623730951d0
            U_SMDEvolution(2,3)= 1d0/1.4142135623730951d0
            U_SMDEvolution(3,3)=-1d0/1.4142135623730951d0
            index=4
            do i=2,SMDEvolutionOrder
                U_SMDEvolution(index:index+i,index:index+i)=Utemp(i).Matrix
                index=index+i+1
            end do
        allocate(U_SMD(NSMDQuantity,NSMDQuantity))!Transform SMD quantities
            U_SMD=0d0
            U_SMD(1,1)=1d0!0th order
            U_SMD(2,2)= 1d0/1.4142135623730951d0!1st order
            U_SMD(3,2)= 1d0/1.4142135623730951d0
            U_SMD(2,3)= 1d0/1.4142135623730951d0
            U_SMD(3,3)=-1d0/1.4142135623730951d0
            index=4
            do i=2,SMDOrder
                U_SMD(index:index+i,index:index+i)=Utemp(i).Matrix
                index=index+i+1
            end do
        do i=3,SMDOrder!Clean up
            deallocate(Utemp(i).Matrix)
        end do
        deallocate(Utemp)
    !Linear mapping matrix from single centre linear coefficient vector to single centre SMD quantity
        allocate(Mtemp(0:SMDOrder,0:BasisOrder))!Prepare
        do i=0,BasisOrder
            do j=0,SMDOrder
                allocate(Mtemp(j,i).Matrix(0:j,0:i))
                do n=0,i
                    do m=0,j
                        if(mod(j-m+i-n,2)==0.and.mod(m+n,2)==0) then
                            Mtemp(j,i).Matrix(m,n)=fct2(j-m+i-n-1)*fct2(m+n-1)/fct(j-m)/fct(m)/fct(i-n)/fct(n)
                        else
                            Mtemp(j,i).Matrix(m,n)=0d0
                        end if
                    end do
                end do
            end do
        end do
        allocate(M_SMDEvolution(NSMDQuantityEvolution,NLinearCoeffSC))!Map SMD quantities up to SMDEvolutionOrder
            index=1
            do i=0,BasisOrder
                indexmrow=1
                do j=0,SMDEvolutionOrder
                    M_SMDEvolution(indexmrow:indexmrow+j,index:index+i)=Mtemp(j,i).Matrix
                    indexmrow=indexmrow+j+1
                end do
                index=index+i+1
            end do
        allocate(M_SMD(NSMDQuantity,NLinearCoeffSC))!Map SMD quantities
            index=1
            do i=0,BasisOrder
                indexmrow=1
                do j=0,SMDOrder
                    M_SMD(indexmrow:indexmrow+j,index:index+i)=Mtemp(j,i).Matrix
                    indexmrow=indexmrow+j+1
                end do
                index=index+i+1
            end do
        do i=0,BasisOrder!Clean up
            do j=0,SMDOrder
                deallocate(Mtemp(j,i).Matrix)
            end do
        end do
        deallocate(Mtemp)
    allocate(Wk_SMDEvolution(NSMDQuantityEvolution,NSMDQuantityEvolution))
    allocate(Wk_SMD(NSMDQuantity,NSMDQuantity))
    allocate(dWk(NSMDQuantityEvolution,NSMDQuantityEvolution))
    allocate(BlockTemp(0:SMDOrder,0:SMDOrder))
    do i=0,SMDOrder
        do j=0,SMDOrder
            allocate(BlockTemp(j,i).Matrix(0:j,0:i))
        end do
    end do
    !Purity work space
    allocate(Ak(MaxNCentre))
    allocate(Bk(MaxNCentre))
    allocate(Ck(MaxNCentre))
    do i=1,MaxNCentre
        allocate(Ak(i).Matrix(2,2))
        allocate(Bk(i).Array(2))
    end do
end subroutine InitializeWigner

!Take in the SMD quantities, fill in high order terms
subroutine CutOffScheme(u)
    type(d2PArray),dimension(0:SMDOrder),intent(inout)::u
    integer::i,j,m,n,index,indexwrow
    real*8::q,p,sigmaq,sigmap
    real*8,dimension(NSMDQuantity)::xi
    !First, use time-dependently evolved SMD quantities to fit a Wigner distribution
    call FitWignerDistribution(u(1:SMDEvolutionOrder))
    !Then compute all terms based on the fitted Wigner distribution
    q=u(1).Array(0)!Save the old location and width
    p=u(1).Array(1)
    sigmaq=u(2).Array(0)
    sigmap=u(2).Array(2)
    call EvaluateSMDQuantity(xi,q,p,sigmaq,sigmap)
    if(dAbs(xi(1)-1d0)>MaxPopDev) then
        write(*,*)'Too large population fluctuation:',xi(1)-1d0
    end if
    u(1).Array(0)=xi(2)*sigmaq+q!Replace with fitted location and width
    u(1).Array(1)=xi(3)*sigmap+p
    u(2).Array(0)=dSqrt(2d0*xi(4)*sigmaq*sigmaq+2d0*q*u(1).Array(0)-q*q-u(1).Array(0)*u(1).Array(0))
    u(2).Array(2)=dSqrt(2d0*xi(6)*sigmap*sigmap+2d0*p*u(1).Array(1)-p*p-u(1).Array(1)*u(1).Array(1))
    !Transform other terms to new location and width
    index=2!Construct SIGMA^fit (in SIGMAk), SIGMA is done during calling EvaluateSMDQuantity
    do i=1,SMDOrder
        do j=0,i
            SIGMAk(index)=u(2).Array(0)**(i-j)*u(2).Array(2)**j
            index=index+1
        end do
    end do
    q=q-u(1).Array(0)!Construct W (in Wk_SMD)
    p=p-u(1).Array(1)
    index=1
    do i=0,SMDOrder
        indexwrow=1
        do j=0,SMDOrder
            do n=0,i
                do m=0,j
                    if(max(0,m+i-j)<=n.and.n<=min(m,i)) then
                        BlockTemp(j,i).Matrix(m,n)=q**(j-m-i+n)*p**(m-n)/fct(j-m-i+n)/fct(m-n)
                    else
                        BlockTemp(j,i).Matrix(m,n)=0d0
                    end if
                end do
            end do
            Wk_SMD(indexwrow:indexwrow+j,index:index+i)=BlockTemp(j,i).Matrix
            indexwrow=indexwrow+j+1
        end do
        index=index+i+1
    end do
    xi=matmul(Wk_SMD,SIGMA*xi)/SIGMAk
    u(2).Array(1)=xi(5)
    index=7
    do i=3,SMDOrder
        do j=0,i
            u(i).Array(j)=xi(index)
            index=index+1
        end do
    end do
end subroutine CutOffScheme

!Input:  the time dependent evolved SMD quantities
!Output: the expansion coefficient of Wigner distribution in global variable wigcoeff
subroutine FitWignerDistribution(u)
    type(d2PArray),dimension(SMDEvolutionOrder),intent(inout)::u
    integer::NCoefficient,index,i,j
    real*8::q,p,sigmaq,sigmap!These will be used in contained subroutines
    real*8,dimension(NSMDQuantityEvolution)::xi
    real*8,allocatable,dimension(:)::c!Fit coefficient vector
    !Prepare
        NCoefficient=NCentre*(5+NLinearCoeffSC)
        q=u(1).Array(0)
        p=u(1).Array(1)
        sigmaq=u(2).Array(0)
        sigmap=u(2).Array(2)
        xi(1)=1d0
        xi(2)=0d0
        xi(3)=0d0
        xi(4)=0.5d0
        xi(5)=u(2).Array(1)
        xi(6)=0.5d0
        index=7
        do i=3,SMDEvolutionOrder
            xi(index:index+i)=u(i).Array
            index=index+i+1
        end do
        allocate(c(NCoefficient))
        call wigcoeff2c(c,NCoefficient)
        index=2!Construct SIGMA
        do i=1,SMDEvolutionOrder
            do j=0,i
                SIGMA(index)=sigmaq**(i-j)*sigmap**j
                index=index+1
            end do
        end do
    call TrustRegion(residue,c,NSMDQuantityEvolution,NCoefficient,Jacobian=jacobian,MaxIteration=MaxTRIteration)
    call c2wigcoeff(c,NCoefficient)!Output
    contains
    subroutine residue(r,c,dimdat,dimvar)
        integer,intent(in)::dimdat,dimvar
        real*8,dimension(dimdat),intent(out)::r
        real*8,dimension(dimvar),intent(in)::c
        integer::k,i,j,m,n,index,indexwrow
        real*8::miuqk,miupk,sigmaqk,sigmapk,sigmaalphak,sigmabetak
        call c2wigcoeff(c,dimvar)
        r=0d0
        do k=1,NCentre
            miuqk=wigcoeff(k).miuq-q!Prepare
            miupk=wigcoeff(k).miup-p
            sigmaqk=wigcoeff(k).sigmaq
            sigmapk=wigcoeff(k).sigmap
            sigmaalphak=dSqrt(1d0+wigcoeff(k).cor)
            sigmabetak= dSqrt(1d0-wigcoeff(k).cor)
            index=2!Construct SIGMA_k, SIGMA_k^AlphaBeta
            do i=1,SMDEvolutionOrder
                do j=0,i
                    SIGMAk(index)=sigmaqk**(i-j)*sigmapk**j
                    SIGMAkAlphaBeta(index)=sigmaalphak**(i-j)*sigmabetak**j
                    index=index+1
                end do
            end do
            index=1!Construct W_k
            do i=0,SMDEvolutionOrder
                indexwrow=1
                do j=0,SMDEvolutionOrder
                    do n=0,i
                        do m=0,j
                            if(max(0,m+i-j)<=n.and.n<=min(m,i)) then
                                BlockTemp(j,i).Matrix(m,n)=miuqk**(j-m-i+n)*miupk**(m-n)/fct(j-m-i+n)/fct(m-n)
                            else
                                BlockTemp(j,i).Matrix(m,n)=0d0
                            end if
                        end do
                    end do
                    Wk_SMDEvolution(indexwrow:indexwrow+j,index:index+i)=BlockTemp(j,i).Matrix
                    indexwrow=indexwrow+j+1
                end do
                index=index+i+1
            end do
            r=r+matmul(Wk_SMDEvolution,SIGMAk(1:NSMDQuantityEvolution)*matmul(U_SMDEvolution,SIGMAkAlphaBeta(1:NSMDQuantityEvolution)*matmul(M_SMDEvolution,matmul(U_Coeff,wigcoeff(k).b))))
        end do
        r=r/SIGMA(1:NSMDQuantityEvolution)-xi
    end subroutine residue
    integer function jacobian(Jacob,c,dimdat,dimvar)
        integer,intent(in)::dimdat,dimvar
        real*8,dimension(dimdat,dimvar),intent(out)::Jacob
        real*8,dimension(dimvar),intent(in)::c
        integer::indexj,k,i,j,m,n,index,indexwrow
        real*8::miuqk,miupk,sigmaqk,sigmapk,sigmaalphak,sigmabetak
        real*8,dimension(NLinearCoeffSC)::vectortemp
        real*8,dimension(NLinearCoeffSC,NLinearCoeffSC)::matrixtemp
        call c2wigcoeff(c,dimvar)
        indexj=0
        do k=1,NCentre
            miuqk=wigcoeff(k).miuq-q!Prepare
            miupk=wigcoeff(k).miup-p
            sigmaqk=wigcoeff(k).sigmaq
            sigmapk=wigcoeff(k).sigmap
            sigmaalphak=dSqrt(1d0+wigcoeff(k).cor)
            sigmabetak= dSqrt(1d0-wigcoeff(k).cor)
            index=2!Construct SIGMA_k, SIGMA_k^AlphaBeta
            do i=1,SMDEvolutionOrder
                do j=0,i
                    SIGMAk(index)=sigmaqk**(i-j)*sigmapk**j
                    SIGMAkAlphaBeta(index)=sigmaalphak**(i-j)*sigmabetak**j
                    index=index+1
                end do
            end do
            index=1!Construct W_k
            do i=0,SMDEvolutionOrder
                indexwrow=1
                do j=0,SMDEvolutionOrder
                    do n=0,i
                        do m=0,j
                            if(max(0,m+i-j)<=n.and.n<=min(m,i)) then
                                BlockTemp(j,i).Matrix(m,n)=miuqk**(j-m-i+n)*miupk**(m-n)/fct(j-m-i+n)/fct(m-n)
                            else
                                BlockTemp(j,i).Matrix(m,n)=0d0
                            end if
                        end do
                    end do
                    Wk_SMDEvolution(indexwrow:indexwrow+j,index:index+i)=BlockTemp(j,i).Matrix
                    indexwrow=indexwrow+j+1
                end do
                index=index+i+1
            end do
            !Nonlinear part
            vectortemp=matmul(M_SMDEvolution,matmul(U_Coeff,wigcoeff(k).b))
            index=2!Construct dSIGMA_k^AlphaBeta / dtheta_k
            do i=1,SMDEvolutionOrder
                do j=0,i
                    dSIGMAk(index)=(dble(i-j)*sigmaalphak**(i-j-2)*sigmabetak**(j)-dble(j)*sigmaalphak**(i-j)*sigmabetak**(j-2))/2d0
                    index=index+1
                end do
            end do
            Jacob(:,indexj+5)=dCos(c(indexj+5))*matmul(Wk_SMDEvolution,SIGMAk(1:NSMDQuantityEvolution)*matmul(U_SMDEvolution,dSIGMAk*vectortemp))
            vectortemp=matmul(U_SMDEvolution,SIGMAkAlphaBeta(1:NSMDQuantityEvolution)*vectortemp)
            index=2!Construct dSIGMA_k / dsigma_qk
            do i=1,SMDEvolutionOrder
                do j=0,i
                    dSIGMAk(index)=dble(i-j)*sigmaqk**(i-j-1)*sigmapk**j
                    index=index+1
                end do
            end do
            Jacob(:,indexj+3)=matmul(Wk_SMDEvolution,dSIGMAk*vectortemp)
            index=2!Construct dSIGMA_k / dsigma_pk
            do i=1,SMDEvolutionOrder
                do j=0,i
                    dSIGMAk(index)=dble(j)*sigmaqk**(i-j)*sigmapk**(j-1)
                    index=index+1
                end do
            end do
            Jacob(:,indexj+4)=matmul(Wk_SMDEvolution,dSIGMAk*vectortemp)
            vectortemp=SIGMAk(1:NSMDQuantityEvolution)*vectortemp
            index=1!Construct dW_k / dsigma_qk
            do i=0,SMDEvolutionOrder
                indexwrow=1
                do j=0,SMDEvolutionOrder
                    do n=0,i
                        do m=0,j
                            if(max(0,m+i-j+1)<=n.and.n<=min(m,i)) then
                                BlockTemp(j,i).Matrix(m,n)=miuqk**(j-m-i+n-1)*miupk**(m-n)/fct(j-m-i+n-1)/fct(m-n)
                            else
                                BlockTemp(j,i).Matrix(m,n)=0d0
                            end if
                        end do
                    end do
                    dWk(indexwrow:indexwrow+j,index:index+i)=BlockTemp(j,i).Matrix
                    indexwrow=indexwrow+j+1
                end do
                index=index+i+1
            end do
            Jacob(:,indexj+1)=matmul(dWk,vectortemp)
            index=1!Construct dW_k / dsigma_pk
            do i=0,SMDEvolutionOrder
                indexwrow=1
                do j=0,SMDEvolutionOrder
                    do n=0,i
                        do m=0,j
                            if(max(0,m+i-j)<=n.and.n<=min(m-1,i)) then
                                BlockTemp(j,i).Matrix(m,n)=miuqk**(j-m-i+n)*miupk**(m-n-1)/fct(j-m-i+n)/fct(m-n-1)
                            else
                                BlockTemp(j,i).Matrix(m,n)=0d0
                            end if
                        end do
                    end do
                    dWk(indexwrow:indexwrow+j,index:index+i)=BlockTemp(j,i).Matrix
                    indexwrow=indexwrow+j+1
                end do
                index=index+i+1
            end do
            Jacob(:,indexj+2)=matmul(dWk,vectortemp)
            !Linear part
            matrixtemp=matmul(M_SMDEvolution,U_Coeff)
            forall(i=1:NSMDQuantityEvolution)
                matrixtemp(i,:)=matrixtemp(i,:)*SIGMAkAlphaBeta(i)
            end forall
            matrixtemp=matmul(U_SMDEvolution,matrixtemp)
            forall(i=1:NSMDQuantityEvolution)
                matrixtemp(i,:)=matrixtemp(i,:)*SIGMAk(i)
            end forall
            indexj=indexj+5
            Jacob(:,indexj+1:indexj+NLinearCoeffSC)=matmul(Wk_SMDEvolution,matrixtemp)
            indexj=indexj+NLinearCoeffSC
        end do
        forall(i=1:NSMDQuantityEvolution)
            Jacob(i,:)=Jacob(i,:)/SIGMA(i)
        end forall
        jacobian=0!Return 0
    end function jacobian
end subroutine FitWignerDistribution

!Transformation of wigner distribution parameters between wigcoeff form and vector form c
subroutine wigcoeff2c(c,N)!In c, cor_k will be replaced with sin(theta_k) to avoid exceeding [-1,1]
    integer,intent(in)::N
    real*8,dimension(N),intent(out)::c
    integer::k,index,index2
    index=0
    do k=1,NCentre
        c(index+1)=wigcoeff(k).miuq
        c(index+2)=wigcoeff(k).miup
        c(index+3)=wigcoeff(k).sigmaq
        c(index+4)=wigcoeff(k).sigmap
        c(index+5)=dAsin(wigcoeff(k).cor)
        index=index+5
        index2=index+NLinearCoeffSC
        c(index+1:index2)=wigcoeff(k).b
        index=index2
    end do
end subroutine wigcoeff2c
subroutine c2wigcoeff(c,N)
    integer,intent(in)::N
    real*8,dimension(N),intent(in)::c
    integer::k,index,index2
    index=0
    do k=1,NCentre
        wigcoeff(k).miuq=c(index+1)
        wigcoeff(k).miup=c(index+2)
        wigcoeff(k).sigmaq=c(index+3)
        wigcoeff(k).sigmap=c(index+4)
        wigcoeff(k).cor=dSin(c(index+5))
        index=index+5
        index2=index+NLinearCoeffSC
        wigcoeff(k).b=c(index+1:index2)
        index=index2
    end do
end subroutine c2wigcoeff

subroutine EvaluateSMDQuantity(xi,q,p,sigmaq,sigmap)!Evaluate SMD quantity from global variable wigcoeff
    real*8,dimension(NSMDQuantity),intent(out)::xi
    real*8,intent(in)::q,p,sigmaq,sigmap
    integer::k,i,j,m,n,index,indexwrow
    real*8::miuqk,miupk,sigmaqk,sigmapk,sigmaalphak,sigmabetak
    index=2!Construct SIGMA
    do i=1,SMDOrder
        do j=0,i
            SIGMA(index)=sigmaq**(i-j)*sigmap**j
            index=index+1
        end do
    end do
    xi=0d0!Summation
    do k=1,NCentre
        miuqk=wigcoeff(k).miuq-q!Prepare
        miupk=wigcoeff(k).miup-p
        sigmaqk=wigcoeff(k).sigmaq
        sigmapk=wigcoeff(k).sigmap
        sigmaalphak=dSqrt(1d0+wigcoeff(k).cor)
        sigmabetak= dSqrt(1d0-wigcoeff(k).cor)
        index=2!Construct SIGMA_k, SIGMA_k^AlphaBeta
        do i=1,SMDOrder
            do j=0,i
                SIGMAk(index)=sigmaqk**(i-j)*sigmapk**j
                SIGMAkAlphaBeta(index)=sigmaalphak**(i-j)*sigmabetak**j
                index=index+1
            end do
        end do
        index=1!Construct W_k
        do i=0,SMDOrder
            indexwrow=1
            do j=0,SMDOrder
                do n=0,i
                    do m=0,j
                        if(max(0,m+i-j)<=n.and.n<=min(m,i)) then
                            BlockTemp(j,i).Matrix(m,n)=miuqk**(j-m-i+n)*miupk**(m-n)/fct(j-m-i+n)/fct(m-n)
                        else
                            BlockTemp(j,i).Matrix(m,n)=0d0
                        end if
                    end do
                end do
                Wk_SMD(indexwrow:indexwrow+j,index:index+i)=BlockTemp(j,i).Matrix
                indexwrow=indexwrow+j+1
            end do
            index=index+i+1
        end do
        xi=xi+matmul(Wk_SMD,SIGMAk*matmul(U_SMD,SIGMAkAlphaBeta*matmul(M_SMD,matmul(U_Coeff,wigcoeff(k).b))))
    end do
    xi=xi/SIGMA
end subroutine EvaluateSMDQuantity

real*8 function purity()
    integer::k1,k2
    real*8::temp,sigmaq_2,sigmap_2,sigmaqp_1
    purity=0d0!Prepare
    do k1=1,NCentre
        temp=1d0-wigcoeff(k1).cor*wigcoeff(k1).cor
        sigmaq_2=1d0/wigcoeff(k1).sigmaq/wigcoeff(k1).sigmaq
        sigmap_2=1d0/wigcoeff(k1).sigmap/wigcoeff(k1).sigmap
        sigmaqp_1=1d0/wigcoeff(k1).sigmaq/wigcoeff(k1).sigmap
        Ak(k1).Matrix(1,1)=sigmaq_2!Construct Ak
        Ak(k1).Matrix(2,1)=-wigcoeff(k1).cor*sigmaqp_1
        Ak(k1).Matrix(1,2)=A(k1).Matrix(2,1)
        Ak(k1).Matrix(2,2)=sigmap_2
        Ak(k1).Matrix=A(k1).Matrix/temp
        Bk(k1).Array(1)=wigcoeff(k1).cor*wigcoeff(k1).miup*sigmaqp_1-wigcoeff(k1).miuq*sigmaq_2
        Bk(k1).Array(2)=wigcoeff(k1).cor*wigcoeff(k1).miuq*sigmaqp_1-wigcoeff(k1).miup*sigmap_2
        Bk(k1).Array=-Bk(k1).Array/temp
        Ck(k1)=(wigcoeff(k1).miuq*wigcoeff(k1).miuq*sigmaq_2-2d0*wigcoeff(k1).cor*wigcoeff(k1).miuq*wigcoeff(k1).miup*sigmaqp_1&
               +wigcoeff(k1).miup*wigcoeff(k1).miup*sigmap_2)/temp
    end do
    do k1=1,NCentre
        do k2=k1+1,NCentre
        end do
    end do
    purity=purity/pim2/hbar
end function purity

real*8 function WignerDistribution(q,p)
    real*8,intent(in)::q,p
    integer::k,i,n,index
    real*8::sc,CapQ,CapP,temp
    WignerDistribution=0d0
    do k=1,NCentre
        sc=0d0
        CapQ=(q-wigcoeff(k).miuq)/wigcoeff(k).sigmaq
        CapP=(p-wigcoeff(k).miup)/wigcoeff(k).sigmap
        index=1
        do i=0,BasisOrder
            do n=0,i
                sc=sc+wigcoeff(k).b(index)*CapQ**(i-n)*CapP**n/fct(i-n)/fct(n)
                index=index+1
            end do
        end do
        temp=1d0-wigcoeff(k).cor*wigcoeff(k).cor
        WignerDistribution=WignerDistribution&
            +sc/wigcoeff(k).sigmaq/wigcoeff(k).sigmap/dSqrt(temp)&
            *dExp(-0.5d0/temp*(CapQ*CapQ-2d0*wigcoeff(k).cor*CapQ*CapP+CapP*CapP))
    end do
    WignerDistribution=WignerDistribution/pim2
end function WignerDistribution

end module Wigner