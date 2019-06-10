!Library interface, data storage, basic routine
!
!SMD terms no higher than SMDEvolutionOrder are evolved time-dependently
!SMDEvolutionOrder+1 : SMDOrder are predicted based on evolved terms
!For polynomial potential, SMDOrder = SMDEvolutionOrder + PotentialOrder - 2
module Basic
    use General
    use Mathematics
    use LinearAlgebra
    use NonlinearOptimization
    implicit none

!Programwide accessed input variable
    real*8::mass,TotalTime,dt,OutputInterval
    integer::SMDEvolutionOrder,BasisOrder,PotentialOrder
    real*8,allocatable,dimension(:)::PolyPotentialCoeff

!Global variable
    integer::ForceOrder,SMDOrder,NSMDQuantityEvolution,NSMDQuantity
    real*8,allocatable,dimension(:)::fct,fct2!Store the necessary factorial, double factorial
    type(d2PArray),allocatable,dimension(:)::pmt,cbn!Store the necessary permutation, combination

!Constant
    real*8,parameter::hbar=1d0

contains
subroutine InitializeBasic()
    integer::i,j
    ForceOrder=PotentialOrder-1
    SMDOrder=SMDEvolutionOrder+ForceOrder-1
    NSMDQuantityEvolution=(1+SMDEvolutionOrder+1)*(SMDEvolutionOrder+1)/2
    NSMDQuantity         =(1+SMDOrder         +1)*(SMDOrder         +1)/2
    allocate(fct(0:max(SMDOrder,BasisOrder)))!Store the necessary factorial
    do i=0,max(SMDOrder,BasisOrder)
        fct(i)=dFactorial(i)
    end do
    allocate(fct2(-1:SMDOrder+BasisOrder))!Store the necessary double factorial
    do i=-1,SMDOrder+BasisOrder
        fct2(i)=dFactorial2(i)
    end do
    allocate(pmt(2:PotentialOrder))!Store the necessary permutation
    do i=2,PotentialOrder
        allocate(pmt(i).Array(0:i))
        do j=0,i
            pmt(i).Array(j)=dPermutation(i,j)
        end do
    end do
    allocate(cbn(0:SMDEvolutionOrder+ForceOrder))!Store the necessary combination
    do i=0,SMDEvolutionOrder+ForceOrder
        allocate(cbn(i).Array(0:i))
        do j=0,i
            cbn(i).Array(j)=dCombination(i,j)
        end do
    end do
end subroutine InitializeBasic

!Compute the derivative of polynomial potential V at coordinate expectation value <q>
!Return <V'> = EdV and V'[(q-<q>)/sigma_q] - <V'> = a(i) * [(q-<q>)/sigma_q]^i / i!
subroutine Polynomial(EdV,a,u)
    real*8,intent(out)::EdV
    real*8,dimension(0:ForceOrder),intent(out)::a
    type(d2PArray),dimension(0:SMDOrder),intent(in)::u
    integer::i,j
    real*8::q,sigma
    q=u(1).Array(0)
    sigma=u(2).Array(0)
    do i=0,ForceOrder
        a(i)=fct(i+1)*PolyPotentialCoeff(i+1)
        do j=i+2,PotentialOrder
            a(i)=a(i)+pmt(j).Array(i+1)*PolyPotentialCoeff(j)*q**(j-(i+1))
        end do
        a(i)=a(i)*sigma**i
    end do
    EdV=a(0)+a(2)*0.5d0
    do i=3,ForceOrder
        EdV=EdV+a(i)*u(i).Array(0)
    end do
    a(0)=a(0)-EdV
end subroutine Polynomial

end module Basic