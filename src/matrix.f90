module matrix
    use constants, only: dp
    implicit none
    
    private
    
    public  :: matrix_solve_eigenvalue_equation
    
    integer, private, save :: range_lwork  = 1
    integer, private, save :: range_lrwork = 1
    integer, private, save :: range_liwork = 1
    integer, private, save :: bound_lwork  = 1
    integer, private, save :: bound_lrwork = 1
    integer, private, save :: bound_liwork = 1
    integer, private, save :: none_lwork   = 1
    
    interface matrix_solve_eigenvalue_equation
        module procedure matrix_solve_eigenvalue_equation_range, &
                         matrix_solve_eigenvalue_equation_bound, &
                         matrix_solve_eigenvalue_equation_none
    end interface matrix_solve_eigenvalue_equation
            

contains
    recursive subroutine matrix_solve_eigenvalue_equation_range(matrix, eigenvalue, eigenvector, num_found, minimum, maximum)
        complex(dp), intent(in)  :: matrix(:, :)
        real(dp),    intent(in)  :: minimum, maximum
        complex(dp), intent(out) :: eigenvector(size(matrix(:, 1)), size(matrix(1, :)))
        real(dp),    intent(out) :: eigenvalue(size(matrix(:, 1)))
        integer,     intent(out) :: num_found
        complex(dp)              :: copy_matrix(size(matrix(:, 1)), size(matrix(1, :)))
        integer                  :: info
        integer                  :: isuppz(2 * size(matrix(:, 1)))
        complex(dp)              :: work(range_lwork)
        real(dp)                 :: rwork(range_lrwork)
        integer                  :: iwork(range_liwork)
        copy_matrix = matrix
        if (range_lwork < 2 * size(matrix(:, 1)) .or. &
            range_lrwork < 24 * size(matrix(:, 1)) .or. &
            range_liwork < 10 * size(matrix(:, 1))) then
            range_lwork  = -1
            range_lrwork = -1
            range_liwork = -1
            call zheevr('V', 'V', 'U', size(matrix(:, 1)), copy_matrix, size(matrix(:, 1)), minimum, maximum, 0, 0, 0d0, &
                num_found, eigenvalue, eigenvector, &
                size(matrix(:, 1)), isuppz, work, range_lwork, rwork, range_lrwork, iwork, range_liwork, info)
            range_lwork  = int(work(1))
            range_lrwork = int(rwork(1))
            range_liwork = int(iwork(1))
            call matrix_solve_eigenvalue_equation(matrix, eigenvalue, eigenvector, num_found, minimum, maximum)
        else
            call zheevr('V', 'V', 'U', size(matrix(:, 1)), copy_matrix, size(matrix(:, 1)), minimum, maximum, 0, 0, 0d0, &
                num_found, eigenvalue, eigenvector, &
                size(matrix(:, 1)), isuppz, work, range_lwork, rwork, range_lrwork, iwork, range_liwork, info)
        end if
        
    end subroutine matrix_solve_eigenvalue_equation_range
    
    recursive subroutine matrix_solve_eigenvalue_equation_bound(matrix, eigenvalue, eigenvector, num_found, minimum, maximum)
        complex(dp), intent(in)  :: matrix(:, :)
        integer,     intent(in)  :: minimum, maximum
        complex(dp), intent(out) :: eigenvector(size(matrix(:, 1)), size(matrix(1, :)))
        real(dp),    intent(out) :: eigenvalue(size(matrix(:, 1)))
        integer,     intent(out) :: num_found
        complex(dp)              :: copy_matrix(size(matrix(:, 1)), size(matrix(1, :)))
        integer                  :: info
        integer                  :: isuppz(2 * size(matrix(:, 1)))
        complex(dp)              :: work(bound_lwork)
        real(dp)                 :: rwork(bound_lrwork)
        integer                  :: iwork(bound_liwork)
        copy_matrix = matrix
        if (bound_lwork < 2 * size(matrix(:, 1)) .or. &
            bound_lrwork < 24 * size(matrix(:, 1)) .or. &
            bound_liwork < 10 * size(matrix(:, 1))) then
            bound_lwork  = -1
            bound_lrwork = -1
            bound_liwork = -1
            call zheevr('V', 'I', 'U', size(matrix(:, 1)), copy_matrix, size(matrix(:, 1)), 0.0_dp, 0.0_dp, minimum, maximum, &
                0.0_dp, num_found, eigenvalue, eigenvector, &
                size(matrix(:, 1)), isuppz, work, bound_lwork, rwork, bound_lrwork, iwork, bound_liwork, info)
            bound_lwork  = int(work(1))
            bound_lrwork = int(rwork(1))
            bound_liwork = int(iwork(1))
            call matrix_solve_eigenvalue_equation(matrix, eigenvalue, eigenvector, num_found, minimum, maximum)
        else
            call zheevr('V', 'I', 'U', size(matrix(:, 1)), copy_matrix, size(matrix(:, 1)), 0.0_dp, 0.0_dp, minimum, maximum, &
                0.0_dp, num_found, eigenvalue, eigenvector, &
                size(matrix(:, 1)), isuppz, work, bound_lwork, rwork, bound_lrwork, iwork, bound_liwork, info)
        end if
        
    end subroutine matrix_solve_eigenvalue_equation_bound
    
    recursive subroutine matrix_solve_eigenvalue_equation_none(matrix, eigenvalue, eigenvector)
        complex(dp), intent(in)  :: matrix(:, :)
        complex(dp), intent(out) :: eigenvector(size(matrix(:, 1)), size(matrix(1, :)))
        real(dp),    intent(out) :: eigenvalue(size(matrix(:, 1)))
        complex(dp)              :: copy_matrix(size(matrix(:, 1)), size(matrix(1, :)))
        integer                  :: info
        complex(dp)              :: work(none_lwork)
        real(dp)                 :: rwork(3 * size(matrix(:, 1)) - 2)
        copy_matrix = matrix
        if (none_lwork < 2 * size(matrix(:, 1)) - 1) then
            none_lwork  = -1
            call zheev('V', 'U', size(matrix(:, 1)), copy_matrix, size(matrix(:, 1)), eigenvalue, work, none_lwork, rwork, info)
            none_lwork  = int(work(1))
            call matrix_solve_eigenvalue_equation(matrix, eigenvalue, eigenvector)
        else
            call zheev('V', 'U', size(matrix(:, 1)), copy_matrix, size(matrix(:, 1)), eigenvalue, work, none_lwork, rwork, info)
            eigenvector = copy_matrix
        end if
        
    end subroutine matrix_solve_eigenvalue_equation_none

end module matrix
