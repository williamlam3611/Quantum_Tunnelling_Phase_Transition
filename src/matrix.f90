module hqt_matrix
    use hqt_constants, only: hqt_dp
    implicit none
    
    private
    
    public  :: hqt_solve_eigenvalue_equation
    
    integer, private, save :: range_lwork  = 1
    integer, private, save :: range_lrwork = 1
    integer, private, save :: range_liwork = 1
    integer, private, save :: bound_lwork  = 1
    integer, private, save :: bound_lrwork = 1
    integer, private, save :: bound_liwork = 1
    integer, private, save :: none_lwork   = 1
    
    interface hqt_solve_eigenvalue_equation
        module procedure hqt_solve_eigenvalue_equation_range, &
                         hqt_solve_eigenvalue_equation_bound, &
                         hqt_solve_eigenvalue_equation_none
    end interface hqt_solve_eigenvalue_equation
            

contains
    recursive subroutine hqt_solve_eigenvalue_equation_range(matrix, eigenvalue, eigenvector, num_found, minimum, maximum)
        complex(hqt_dp), intent(in)  :: matrix(:, :)
        real(hqt_dp),    intent(in)  :: minimum, maximum
        complex(hqt_dp), intent(out) :: eigenvector(size(matrix(:, 1)), size(matrix(1, :)))
        real(hqt_dp),    intent(out) :: eigenvalue(size(matrix(:, 1)))
        integer,     intent(out)     :: num_found
        complex(hqt_dp)              :: copy_matrix(size(matrix(:, 1)), size(matrix(1, :)))
        integer                      :: info
        integer                      :: isuppz(2 * size(matrix(:, 1)))
        complex(hqt_dp)              :: work(range_lwork)
        real(hqt_dp)                 :: rwork(range_lrwork)
        integer                      :: iwork(range_liwork)
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
            call hqt_solve_eigenvalue_equation(matrix, eigenvalue, eigenvector, num_found, minimum, maximum)
        else
            call zheevr('V', 'V', 'U', size(matrix(:, 1)), copy_matrix, size(matrix(:, 1)), minimum, maximum, 0, 0, 0d0, &
                num_found, eigenvalue, eigenvector, &
                size(matrix(:, 1)), isuppz, work, range_lwork, rwork, range_lrwork, iwork, range_liwork, info)
        end if
        if (info .ne. 0) print *, "Matrix Error! Info code: ", info, " see: "// &
                                  "http://www.netlib.org/lapack/explore-html/d9/dd2/zheevr_8f_source.html"
        
    end subroutine hqt_solve_eigenvalue_equation_range
    
    recursive subroutine hqt_solve_eigenvalue_equation_bound(matrix, eigenvalue, eigenvector, num_found, minimum, maximum)
        complex(hqt_dp), intent(in)  :: matrix(:, :)
        integer,     intent(in)      :: minimum, maximum
        complex(hqt_dp), intent(out) :: eigenvector(size(matrix(:, 1)), size(matrix(1, :)))
        real(hqt_dp),    intent(out) :: eigenvalue(size(matrix(:, 1)))
        integer,     intent(out)     :: num_found
        complex(hqt_dp)              :: copy_matrix(size(matrix(:, 1)), size(matrix(1, :)))
        integer                      :: info, i
        integer                      :: isuppz(2 * size(matrix(:, 1)))
        complex(hqt_dp)              :: work(bound_lwork)
        real(hqt_dp)                 :: rwork(bound_lrwork)
        integer                      :: iwork(bound_liwork)
        copy_matrix = matrix
        if (bound_lwork < 2 * size(matrix(:, 1)) .or. &
            bound_lrwork < 24 * size(matrix(:, 1)) .or. &
            bound_liwork < 10 * size(matrix(:, 1))) then
            bound_lwork  = -1
            bound_lrwork = -1
            bound_liwork = -1
            call zheevr('V', 'I', 'U', size(matrix(:, 1)), copy_matrix, size(matrix(:, 1)), 0.0_hqt_dp, 0.0_hqt_dp, &
                minimum, maximum, 0.0_hqt_dp, num_found, eigenvalue, eigenvector, &
                size(matrix(:, 1)), isuppz, work, bound_lwork, rwork, bound_lrwork, iwork, bound_liwork, info)
            bound_lwork  = int(work(1))
            bound_lrwork = int(rwork(1))
            bound_liwork = int(iwork(1))
            call hqt_solve_eigenvalue_equation(matrix, eigenvalue, eigenvector, num_found, minimum, maximum)
        else
            call zheevr('V', 'I', 'U', size(matrix(:, 1)), copy_matrix, size(matrix(:, 1)), 0.0_hqt_dp, 0.0_hqt_dp, &
                minimum, maximum, 0.0_hqt_dp, num_found, eigenvalue, eigenvector, &
                size(matrix(:, 1)), isuppz, work, bound_lwork, rwork, bound_lrwork, iwork, bound_liwork, info)
        end if
        if (info .ne. 0) print *, "Matrix Error! Info code: ", info, " see: "// &
                                  "http://www.netlib.org/lapack/explore-html/d9/dd2/zheevr_8f_source.html" 
        
    end subroutine hqt_solve_eigenvalue_equation_bound
    
    recursive subroutine hqt_solve_eigenvalue_equation_none(matrix, eigenvalue, eigenvector)
        complex(hqt_dp), intent(in)  :: matrix(:, :)
        complex(hqt_dp), intent(out) :: eigenvector(size(matrix(:, 1)), size(matrix(1, :)))
        real(hqt_dp),    intent(out) :: eigenvalue(size(matrix(:, 1)))
        integer                      :: info
        complex(hqt_dp)              :: work(none_lwork)
        real(hqt_dp)                 :: rwork(3 * size(matrix(:, 1)) - 2)
        eigenvector = matrix
        if (none_lwork < 2 * size(matrix(:, 1)) - 1) then
            none_lwork  = -1
            call zheev('V', 'U', size(matrix(:, 1)), eigenvector, size(matrix(:, 1)), eigenvalue, work, none_lwork, rwork, info)
            none_lwork  = int(work(1))
            call hqt_solve_eigenvalue_equation(matrix, eigenvalue, eigenvector)
        else
            call zheev('V', 'U', size(matrix(:, 1)), eigenvector, size(matrix(:, 1)), eigenvalue, work, none_lwork, rwork, info)
            eigenvector = eigenvector
        end if
        if (info .ne. 0) print *, "Matrix Error! Info code ", info, " see: "// &
                                  "http://www.netlib.org/lapack/explore-html/d6/dee/zheev_8f_source.html"
        
    end subroutine hqt_solve_eigenvalue_equation_none

end module hqt_matrix
