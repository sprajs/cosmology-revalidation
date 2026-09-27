! Synthetic sequence-association proof; no CAMB physics or spline evaluation.
program indexing
    use, intrinsic :: ieee_arithmetic
    implicit none
    integer, parameter :: n=5
    real(8) :: transfer_times(0:n+1), old_knots(n), fixed_knots(n)
    integer :: i
    transfer_times(0)=ieee_value(0d0,ieee_quiet_nan)
    transfer_times(n+1)=-777d0
    do i=1,n
        transfer_times(i)=10d0*i
    end do
    call explicit_shape(transfer_times,n,old_knots)
    call explicit_shape(transfer_times(1:n),n,fixed_knots)
    if (.not.ieee_is_nan(old_knots(1))) error stop 'old first knot not sentinel'
    if (any(old_knots(2:n)/=transfer_times(1:n-1))) error stop 'old shifted knots'
    if (any(fixed_knots/=transfer_times(1:n))) error stop 'fixed knots mismatch'
    if (any(.not.ieee_is_finite(fixed_knots))) error stop 'fixed nonfinite'
    print *, 'PASS: old first knot is sentinel; old remaining knots shift; fixed knots match.'
    print *, 'FIXED_KNOTS',fixed_knots
contains
    subroutine explicit_shape(x,n,y)
        integer, intent(in) :: n
        real(8), intent(in) :: x(n)
        real(8), intent(out) :: y(n)
        y=x
    end subroutine explicit_shape
end program indexing
