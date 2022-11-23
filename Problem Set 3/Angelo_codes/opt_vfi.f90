! This code solves a stochastic growth model by value function iteration

program  opt_vfi
    implicit none

    !Define useful variables
    real :: norm, lower, upper
    
    integer :: i_A, i_k, i_k_prime, iter, fake_iter       !Indexes and iteration stuff
    integer :: ind_policy, ind_temp

    !Model parameters
    real, parameter :: A=1.0            !Productivity term
    real, parameter :: beta=0.96        !Discount factor
    real, parameter :: alpha=0.4        !Capital share
    real, parameter :: delta=1.d0       !Depreciation

	integer, parameter :: n_k = 1001                           !Capital gridpoints
    integer, parameter :: n_A = 15                           !Capital gridpoints
    integer, parameter :: m = 500                              !Howard's parameter
    integer, parameter :: fundamental_interations = 2          !Initial regular VFI iterations
    integer, parameter :: max_iter = 1500

    !Define grids
	real, dimension(n_A, n_k) :: k_grid
    real, dimension(n_A, n_k) :: k_policy
    real, dimension(n_A, n_k) :: c_policy
    integer, dimension(n_A, n_k) :: policy_index
    real, dimension(n_A, n_k) :: valuefunction
    real, dimension(n_A, n_k) :: V_howard
	real, dimension(n_A, n_k) :: V0=0.0
    real, dimension(n_A, n_k) :: V_old=0.0
    real, dimension(n_A, n_k) :: V_temp=0.0


    real, parameter :: eps=1e-8                 !Ridiculous precision



    !Code parameters
    real :: k, k_prime, c, thevalue, k_ss

    k_ss = (alpha*A / (1/beta + delta -1) ) ** (1/(1-alpha))

    
    do i_k = 1, n_k
        k_grid(i_k) = 0.75*k_ss + (i_k-1)*k_ss/1000
    end do
    !A, beta, alpha, delta, k_ss !, k_lb, k_ub        
    
    
    
    valuefunction = V_old

    do iter = 1, max_iter
        
        !Initial guess
        V0 = valuefunction

        do i_k = 1, n_k !i_k
            k = k_grid(i_k)

            do i_k_prime = 1, n_k !i_k_prime
                
                k_prime = k_grid(i_k_prime)
                !Flow utility
                c = A*k**(alpha) + (1-delta)*k - k_prime
                if ( c<0 ) then
                    c=eps
                end if

                V_temp(i_k_prime) = log(c) + beta*V0(i_k_prime)

                
            end do !i_k_prime

            !Pick optimal policy
            ind_temp = maxloc(V_temp)
            ind_policy = ind_temp
            thevalue = V_temp(ind_policy)

            !Update policy & VF
            valuefunction(i_k) = thevalue
            k_policy(i_k) = k_grid(ind_policy)
            c_policy(i_k) = A*k**alpha + (1-delta)*k -  k_policy(i_k)
            policy_index(i_k) = ind_policy

            
        end do !i_k

        
		do fake_iter = 1, m
            V_howard = log(c_policy) + beta*valuefunction(policy_index)
            valuefunction = V_howard
        end do

        upper = (beta / (1-beta)) * maxval(valuefunction(:)-V0(:))
        lower = (beta / (1-beta)) * minval(valuefunction(:)-V0(:))




        !Check convergence
        norm = upper - lower

        if ( norm<eps ) then
            print *, 'Value function converged'
            print *, 'VF iterations = ', iter
            exit
        else
            print *, 'Iteration =', iter
            print *, 'Error =', norm

        end if

        
    end do
    
    !Restart initial guess 
    V_old = valuefunction

    !Write files
    open(11,file='k_ss.txt',status='replace')
	write(11,*) k_ss
    close(11)

    open(11,file='k_grid.txt',status='replace')
    do i_k=1,n_k
	    write(11,*) k_grid(i_k)
	enddo
    close(11)

    open(11,file='k_policy.txt',status='replace')
    do i_k=1,n_k
	    write(11,*) k_policy(i_k)
	enddo
    close(11)

    open(11,file='c_policy.txt',status='replace')
    do i_k=1,n_k
	    write(11,*) c_policy(i_k)
	enddo
    close(11)

    open(11,file='policy_index.txt',status='replace')
    do i_k=1,n_k
	    write(11,*) policy_index(i_k)
	enddo
    close(11)

    open(11,file='valuefunction.txt',status='replace')
    do i_k=1,n_k
	    write(11,*) valuefunction(i_k)
	enddo
    close(11)


end program  opt_vfi