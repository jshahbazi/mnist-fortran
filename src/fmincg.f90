module fmincg
    use ieee_arithmetic
    use costfunction
    
    implicit none
    
    contains
    !cost               scalar
    !gradient           vector
    !limit              scalar
    !point              scalar
    !search_direction   vector
    !slope              scalar
    

    
real function fmincg2(length, nn_params, input_layer_size, hidden_layer_size, num_labels, inputdata, y, lambda)
    implicit none
    ! Copyright (C) 2001 and 2002 by Carl Edward Rasmussen. Date 2002-02-13
    ! (C) Copyright 1999, 2000 & 2001, Carl Edward Rasmussen
    ! 
    ! Permission is granted for anyone to copy, use, or modify these
    ! programs and accompanying documents for purposes of research or
    ! education, provided this copyright notice is retained, and note is
    ! made of any changes that have been made.
    ! 
    ! These programs and documents are distributed without any warranty,
    ! express or implied.  As the programs were written for research
    ! purposes only, they have not been tested to the degree that would be
    ! advisable in any important application.  All use of these programs is
    ! entirely at the user's own risk.
    !
    ! [ml-class] Changes Made:
    ! 1) Function name and argument specifications
    ! 2) Output display
    !
    ! [John Shahbazian http://fortrandev.wordpress.com] Changes Made:
    ! 1) Ported to Fortran
    ! 2) Change the cost function call to internal.  Replace the
    !    'costandgradient' function to whatever you would like.  It returns
    !    the cost as the result, and the gradient is returned through the first 
    !    argument as an intent(inout) (e.g. 'gradient#').  
    ! 3) Changed the variable names to be readable.
    !    f1 = cost1
    !    df1 = gradient1
    !    s = search_direction
    !    d1 = slope1
    !    z1 = point1
    !    X0 = backup_params
    !    f0 = cost_backup
    !    df0 = gradient_backup    
    ! 4) Optimized for Intel MKL
    
    integer, intent(in) :: length,input_layer_size,hidden_layer_size,num_labels
    real, allocatable, intent(in) :: inputdata(:,:), y(:,:)
    real, intent(in) :: lambda
    real, allocatable, intent(inout) :: nn_params(:,:)


    real :: RHO                            ! a bunch of constants for line searches
    real :: SIG       ! RHO and SIG are the constants in the Wolfe-Powell conditions
    real :: INT    ! don't reevaluate within 0.1 of the limit of the current bracket
    real :: EXT                   ! extrapolate maximum 3 times the current bracket
    integer :: MAXEVALS                         ! max 20 function evaluations per line search
    real :: RATIO                                      ! maximum allowed slope ratio

    real :: mintemp, minstuff, M, A, B
    real :: fX 

    integer :: success
    integer :: i                                          ! zero the run length counter
    integer :: ls_failed                                    ! no previous line search has failed

    real, allocatable :: backup_params(:,:)
    real :: cost1, cost2, cost3, cost_backup
    real, allocatable :: gradient1(:,:), gradient2(:,:), gradient3(:,:), gradient_backup(:,:), tmp(:,:)
    real :: limit
    real :: point1, point2, point3, ptemp
    real, allocatable :: search_direction(:,:), stemp(:,:)
    real :: slope1, slope2, slope3, slope_vector(1,1)
    real :: sqrtnumber
    real :: sd_calc_1(1,1), sd_calc_2(1,1), sd_calc_3(1,1), sd_calc_4
    real, allocatable :: sd_calc_5(:,:)

    integer :: a1,a2,b1,b2
    
    !double precision :: timer1, timer2, average, counter = 0.0

    allocate(gradient1(1:size(nn_params),1:1))
    allocate(gradient2(1:size(nn_params),1:1))
    allocate(gradient3(1:size(nn_params),1:1))
    allocate(gradient_backup(1:size(nn_params),1:1))
    allocate(tmp(1:size(nn_params),1:1))
    allocate(search_direction(1:size(nn_params),1:1))
    allocate(stemp(1:size(nn_params),1:1))
    allocate(backup_params(1:size(nn_params),1:1))
    allocate(sd_calc_5(size(search_direction,1),size(search_direction,2)))

    
    RHO = 0.01                            ! a bunch of constants for line searches
    SIG = 0.5       ! RHO and SIG are the constants in the Wolfe-Powell conditions
    INT = 0.1    ! don't reevaluate within 0.1 of the limit of the current bracket
    EXT = 3.0                    ! extrapolate maximum 3 times the current bracket
    MAXEVALS = 20                         ! max 20 function evaluations per line search
    RATIO = 100                                      ! maximum allowed slope ratio
    
    i=0
    fX = 0.0
    ls_failed = 0

    cost1 = costandgradient(gradient1,nn_params,input_layer_size,hidden_layer_size,num_labels, inputdata, y, lambda)

    i = i + (length<0)                                          ! count epochs?!
    search_direction = -gradient1                               ! search direction is steepest
    !slope_vector = matmul(-transpose(search_direction),search_direction)    ! this is the slope
        a1 = size(search_direction,2)
        a2 = size(search_direction,1)
        b1 = size(search_direction,1)
        b2 = size(search_direction,2)
        call sGEMM('T','N',a1,b2,b1,1.0,-search_direction,a2,search_direction,b1,0.0,slope_vector,a1)
        
    slope1 = slope_vector(1,1)                                  !convert to scalar
    point1  = 1.0/(1.0-slope1)                                  ! initial step is red/(|s|+1)

    
    do while(i < abs(length))                                      
        i = i + (length>0)                                      ! count iterations?!
        backup_params = nn_params
        cost_backup = cost1
        gradient_backup = gradient1                             ! make a copy of current values
        stemp = point1  * search_direction
        nn_params = nn_params + stemp                           ! begin line search
        gradient2 = gradient1

        cost2 = costandgradient(gradient2,nn_params,input_layer_size,hidden_layer_size,num_labels, inputdata, y, lambda)

        i = i + (length<0)                                      ! count epochs?!
        !slope_vector = matmul(transpose(gradient2),search_direction)              ! this is the slope
            a1 = size(gradient2,2)
            a2 = size(gradient2,1)
            b1 = size(search_direction,1)
            b2 = size(search_direction,2)
            call sGEMM('T','N',a1,b2,b1,1.0,gradient2,a2,search_direction,b1,0.0,slope_vector,a1)
                
        
        slope2 = slope_vector(1,1)                              !convert to scalar

        cost3 = cost1
        slope3 = slope1
        point3 = -point1                                        ! initialize point 3 equal to point 1
        if(length > 0)then
            M = MAXEVALS
        else
            M = min(MAXEVALS, -length-i)
        end if
        success = 0
        limit = -1                                              ! initialize quantities
        do
            do while (( (cost2 > (cost1 + (point1 * RHO * slope1))) .or. (slope2 > (-SIG * slope1)) ) .and. (M > 0) )
                limit = point1                                  ! tighten the bracket
                if(cost2 > cost1)then
                    point2 = point3 - (0.5 * slope3 * point3 * point3)/(slope3 * point3 + cost2 - cost3)                 ! quadratic fit
                else
                    A = 6*(cost2 - cost3)/point3 + 3*(slope2 + slope3)                                 ! cubic fit
                    B = 3*(cost3 - cost2) - point3 * (slope3 + 2*slope2)
                    point2 = (sqrt(B*B - A*slope2*point3*point3) - B)/A
                end if
                if(ieee_is_nan(point2) .or. (.not. ieee_is_finite(point2)))then
                    point2 = point3 / 2                         ! if we had a numerical problem then bisect
                end if
                point2 = MAX( MIN(point2, (INT * point3)), ((1.0 - INT) * point3))  ! don't accept too close to limits
                point1  = point1 + point2                       ! update the step
                stemp = point2 * search_direction
                nn_params = nn_params + stemp

                cost2 = costandgradient(gradient2,nn_params,input_layer_size,hidden_layer_size,num_labels, inputdata, y, lambda)

                M = M - 1.0
                i = i + (length<0)                              ! count epochs?!
                !slope_vector = matmul(transpose(gradient2),search_direction)
                    a1 = size(gradient2,2)
                    a2 = size(gradient2,1)
                    b1 = size(search_direction,1)
                    b2 = size(search_direction,2)
                    call sGEMM('T','N',a1,b2,b1,1.0,gradient2,a2,search_direction,b1,0.0,slope_vector,a1)
                                
                slope2 = slope_vector(1,1)                      !convert to scalar
                point3 = point3 - point2                        ! point3 is now relative to the location of point2
            end do


            if ((cost2 > (cost1 + (point1*RHO*slope1)) ) .or. (slope2 > (-SIG * slope1) ))then
                exit                                            ! this is a failure
            elseif (slope2 > (SIG * slope1))then
                success = 1
                exit                                            ! success
            elseif (M == 0)then
                exit                                            ! failure
            end if
            A = 6*(cost2 - cost3)/point3 + 3*(slope2 + slope3)  ! make cubic extrapolation
            B = 3*(cost3 - cost2) - point3*(slope3 + 2*slope2)
            sqrtnumber = (B*B) - A*slope2*point3*point3
            if((.not. ieee_is_normal(sqrtnumber)) .or. sqrtnumber < 0 )then
                if (limit < -0.5) then                          ! if we have no upper limit
                    point2 = point1  * (EXT - 1)                ! the extrapolate the maximum amount
                else
                    point2 = (limit - point1)/2                 ! otherwise bisect
                end if
            else
                point2 = (-slope2 * point3 * point3)/(B + sqrt(sqrtnumber))
                if ((limit > -0.5) .and. ((point2 + point1) > limit))then          ! extraplation beyond max?
                    point2 = (limit - point1)/2                 ! bisect
                elseif ((limit < -0.5) .and. ((point2 + point1) > (point1 * EXT)))then       ! extrapolation beyond limit
                    point2 = point1 * (EXT - 1.0)               ! set to extrapolation limit
                elseif (point2 < (-point3 * INT))then
                    point2 = -point3 * INT
                elseif ((limit > -0.5) .and. (point2 < (limit - point1)*(1.0 - INT)))then   ! too close to limit?
                    point2 = (limit - point1 ) * (1.0 - INT)
                end if
            end if

            cost3 = cost2
            slope3 = slope2
            point3 = -point2               
            point1  = point1  + point2

            stemp = point2 * search_direction
            nn_params = nn_params + stemp                       ! update current estimates

            cost2 = costandgradient(gradient2,nn_params,input_layer_size,hidden_layer_size,num_labels, inputdata, y, lambda)

            M = M - 1.0
            i = i + (length<0)                                  ! count epochs?!
            !slope_vector = matmul(transpose(gradient2),search_direction)
                a1 = size(gradient2,2)
                a2 = size(gradient2,1)
                b1 = size(search_direction,1)
                b2 = size(search_direction,2)
                call sGEMM('T','N',a1,b2,b1,1.0,gradient2,a2,search_direction,b1,0.0,slope_vector,a1)
                            
            slope2 = slope_vector(1,1)                          !convert to scalar
        end do                                                  ! end of line search

        if (success == 1)then                                   ! if line search succeeded
            cost1 = cost2
            fX = cost1
            print '(a, I6, a, f8.5)', "Iteration: ", i, " | Cost: ", cost1

!timer1=secnds(0.0)            
            ! this calculation was split up to make it easier to work with
            !sd_calc_1 = matmul(transpose(gradient2),gradient2)
                a1 = size(gradient2,2)
                a2 = size(gradient2,1)
                b1 = size(gradient2,1)
                b2 = size(gradient2,2)
                call sGEMM('T','N',a1,b2,b1,1.0,gradient2,a2,gradient2,b1,0.0,sd_calc_1,a1)
                                   
            !sd_calc_2 = matmul(transpose(gradient1),gradient2)
                a1 = size(gradient1,2)
                a2 = size(gradient1,1)
                b1 = size(gradient2,1)
                b2 = size(gradient2,2)
                call sGEMM('T','N',a1,b2,b1,1.0,gradient1,a2,gradient2,b1,0.0,sd_calc_2,a1)
                                   
            !sd_calc_3 = matmul(transpose(gradient1),gradient1)
                a1 = size(gradient1,2)
                a2 = size(gradient1,1)
                b1 = size(gradient1,1)
                b2 = size(gradient1,2)
                call sGEMM('T','N',a1,b2,b1,1.0,gradient1,a2,gradient1,b1,0.0,sd_calc_3,a1)
                                   
            sd_calc_4 = (sd_calc_1(1,1) - sd_calc_2(1,1)) / sd_calc_3(1,1)
            sd_calc_5 = sd_calc_4 * search_direction
            search_direction = sd_calc_5 - gradient2
!timer2=secnds(timer1)
!counter = counter + 1.0
!print '(a, f8.6)','timer: ', timer2/counter
            tmp = gradient1
            gradient1 = gradient2
            gradient2 = tmp                                     ! swap derivatives
            !slope_vector = matmul(transpose(gradient1),search_direction)
                a1 = size(gradient1,2)
                a2 = size(gradient1,1)
                b1 = size(search_direction,1)
                b2 = size(search_direction,2)
                call sGEMM('T','N',a1,b2,b1,1.0,gradient1,a2,search_direction,b1,0.0,slope_vector,a1)
                
            
            slope2 = slope_vector(1,1)                          !convert to scalar
            if(slope2 > 0)then                                  ! new slope must be negative
                search_direction = -gradient1                   ! otherwise use steepest direction
                !slope_vector = matmul(-transpose(search_direction),search_direction)
                    a1 = size(search_direction,2)
                    a2 = size(search_direction,1)
                    b1 = size(search_direction,1)
                    b2 = size(search_direction,2)
                    call sGEMM('T','N',a1,b2,b1,1.0,-search_direction,a2,search_direction,b1,0.0,slope_vector,a1)                        
                
                slope2 = slope_vector(1,1)                      !convert to scalar
            end if
            mintemp = slope1/(slope2 - 2.2251D-308)  !2.2551D-308 is min value double precision float
            minstuff = min(RATIO, mintemp)
            point1  = point1  * minstuff                         ! slope ratio but max RATIO
            slope1 = slope2
            ls_failed = 0                                        ! this line search did not fail
        else
            nn_params = backup_params
            cost1 = cost_backup
            gradient1 = gradient_backup                         ! restore point from before failed line search
            if (ls_failed == 1 .or. i > abs(length)) then       ! line search failed twice in a row
                exit                                            ! or we ran out of time, so we give up
            end if
            tmp = gradient1
            gradient1 = gradient2
            gradient2 = tmp                                    ! swap derivatives
            search_direction = -gradient1                      ! try steepest
            !slope_vector = matmul(-transpose(search_direction),search_direction)
                a1 = size(search_direction,2)
                a2 = size(search_direction,1)
                b1 = size(search_direction,1)
                b2 = size(search_direction,2)
                call sGEMM('T','N',a1,b2,b1,1.0,-search_direction,a2,search_direction,b1,0.0,slope_vector,a1)
                    
            slope1 = slope_vector(1,1)                         ! convert to scalar
            point1  = 1.0 / (1.0 - slope1)
            ls_failed = 1                                      ! this line search failed
        end if

    end do

    fmincg2 = fX
end function fmincg2    

    
!real function fmincgfun(length, nn_params, input_layer_size, hidden_layer_size, num_labels, inputdata, y, lambda)
!    implicit none
!    ! Copyright (C) 2001 and 2002 by Carl Edward Rasmussen. Date 2002-02-13
!    ! (C) Copyright 1999, 2000 & 2001, Carl Edward Rasmussen
!    ! 
!    ! Permission is granted for anyone to copy, use, or modify these
!    ! programs and accompanying documents for purposes of research or
!    ! education, provided this copyright notice is retained, and note is
!    ! made of any changes that have been made.
!    ! 
!    ! These programs and documents are distributed without any warranty,
!    ! express or implied.  As the programs were written for research
!    ! purposes only, they have not been tested to the degree that would be
!    ! advisable in any important application.  All use of these programs is
!    ! entirely at the user's own risk.
!    !
!    ! [ml-class] Changes Made:
!    ! 1) Function name and argument specifications
!    ! 2) Output display
!    !
!    ! [John Shahbazian http://fortrandev.wordpress.com] Changes Made:
!    ! 1) Ported to Fortran
!    ! 2) Change the cost function call to internal.  Replace the
!    !    'costandgradient' function to whatever you would like.  It returns
!    !    the cost as the result, and the gradient is returned through the first 
!    !    argument as an intent(inout) (e.g. 'gradient#').  
!    ! 3) Changed the variable names to be readable.
!    !    f1 = cost1
!    !    df1 = gradient1
!    !    s = search_direction
!    !    d1 = slope1
!    !    z1 = point1
!    !    X0 = backup_params
!    !    f0 = cost_backup
!    !    df0 = gradient_backup    
!    
!    integer, intent(in) :: length,input_layer_size,hidden_layer_size,num_labels
!    real, allocatable, intent(in) :: inputdata(:,:), y(:,:)
!    real, intent(in) :: lambda
!    real, allocatable, intent(inout) :: nn_params(:,:)
!
!
!    real :: RHO = 0.01                            ! a bunch of constants for line searches
!    real :: SIG = 0.5       ! RHO and SIG are the constants in the Wolfe-Powell conditions
!    real :: INT = 0.1    ! don't reevaluate within 0.1 of the limit of the current bracket
!    real :: EXT = 3.0                    ! extrapolate maximum 3 times the current bracket
!    integer :: MAXEVALS = 20                         ! max 20 function evaluations per line search
!    real :: RATIO = 100                                      ! maximum allowed slope ratio
!
!    real :: mintemp, minstuff, M, A, B
!    real :: fX = 0.0
!
!    integer :: success
!    integer :: i = 0                                            ! zero the run length counter
!    integer :: ls_failed = 0                                    ! no previous line search has failed
!
!    real, allocatable :: backup_params(:,:)
!    real :: cost1, cost2, cost3, cost_backup
!    real, allocatable :: gradient1(:,:), gradient2(:,:), gradient3(:,:), gradient_backup(:,:), tmp(:,:)
!    real :: limit
!    real :: point1, point2, point3, ptemp
!    real, allocatable :: search_direction(:,:), stemp(:,:)
!    real :: slope1, slope2, slope3, slope_vector(1,1)
!    real :: sqrtnumber
!    real :: sd_calc_1(1,1), sd_calc_2(1,1), sd_calc_3(1,1), sd_calc_4
!    real, allocatable :: sd_calc_5(:,:)
!
!
!    allocate(gradient1(1:size(nn_params),1:1))
!    allocate(gradient2(1:size(nn_params),1:1))
!    allocate(gradient3(1:size(nn_params),1:1))
!    allocate(gradient_backup(1:size(nn_params),1:1))
!    allocate(tmp(1:size(nn_params),1:1))
!    allocate(search_direction(1:size(nn_params),1:1))
!    allocate(stemp(1:size(nn_params),1:1))
!    allocate(backup_params(1:size(nn_params),1:1))
!    allocate(sd_calc_5(size(search_direction,1),size(search_direction,2)))
!
!
!    cost1 = costandgradient(gradient1,nn_params,input_layer_size,hidden_layer_size,num_labels, inputdata, y, lambda)
!
!    i = i + (length<0)                                          ! count epochs?!
!    search_direction = -gradient1                               ! search direction is steepest
!    slope_vector = matmul(-transpose(search_direction),search_direction)    ! this is the slope
!    slope1 = slope_vector(1,1)                                  !convert to scalar
!    point1  = 1.0/(1.0-slope1)                                  ! initial step is red/(|s|+1)
!
!    
!    do while(i < abs(length))                                      
!        i = i + (length>0)                                      ! count iterations?!
!        backup_params = nn_params
!        cost_backup = cost1
!        gradient_backup = gradient1                             ! make a copy of current values
!        stemp = point1  * search_direction
!        nn_params = nn_params + stemp                           ! begin line search
!        gradient2 = gradient1
!
!        cost2 = costandgradient(gradient2,nn_params,input_layer_size,hidden_layer_size,num_labels, inputdata, y, lambda)
!
!        i = i + (length<0)                                      ! count epochs?!
!        slope_vector = matmul(transpose(gradient2),search_direction)              ! this is the slope
!        slope2 = slope_vector(1,1)                              !convert to scalar
!
!        cost3 = cost1
!        slope3 = slope1
!        point3 = -point1                                        ! initialize point 3 equal to point 1
!        if(length > 0)then
!            M = MAXEVALS
!        else
!            M = min(MAXEVALS, -length-i)
!        end if
!        success = 0
!        limit = -1                                              ! initialize quantities
!        do
!            do while (( (cost2 > (cost1 + (point1 * RHO * slope1))) .or. (slope2 > (-SIG * slope1)) ) .and. (M > 0) )
!                limit = point1                                  ! tighten the bracket
!                if(cost2 > cost1)then
!                    point2 = point3 - (0.5 * slope3 * point3 * point3)/(slope3 * point3 + cost2 - cost3)                 ! quadratic fit
!                else
!                    A = 6*(cost2 - cost3)/point3 + 3*(slope2 + slope3)                                 ! cubic fit
!                    B = 3*(cost3 - cost2) - point3 * (slope3 + 2*slope2)
!                    point2 = (sqrt(B*B - A*slope2*point3*point3) - B)/A
!                end if
!                if(ieee_is_nan(point2) .or. (.not. ieee_is_finite(point2)))then
!                    point2 = point3 / 2                         ! if we had a numerical problem then bisect
!                end if
!                point2 = MAX( MIN(point2, (INT * point3)), ((1.0 - INT) * point3))  ! don't accept too close to limits
!                point1  = point1 + point2                       ! update the step
!                stemp = point2 * search_direction
!                nn_params = nn_params + stemp
!
!                cost2 = costandgradient(gradient2,nn_params,input_layer_size,hidden_layer_size,num_labels, inputdata, y, lambda)
!
!                M = M - 1.0
!                i = i + (length<0)                              ! count epochs?!
!                slope_vector = matmul(transpose(gradient2),search_direction)
!                slope2 = slope_vector(1,1)                      !convert to scalar
!                point3 = point3 - point2                        ! point3 is now relative to the location of point2
!            end do
!
!
!            if ((cost2 > (cost1 + (point1*RHO*slope1)) ) .or. (slope2 > (-SIG * slope1) ))then
!                exit                                            ! this is a failure
!            elseif (slope2 > (SIG * slope1))then
!                success = 1
!                exit                                            ! success
!            elseif (M == 0)then
!                exit                                            ! failure
!            end if
!            A = 6*(cost2 - cost3)/point3 + 3*(slope2 + slope3)  ! make cubic extrapolation
!            B = 3*(cost3 - cost2) - point3*(slope3 + 2*slope2)
!            sqrtnumber = (B*B) - A*slope2*point3*point3
!            if((.not. ieee_is_normal(sqrtnumber)) .or. sqrtnumber < 0 )then
!                if (limit < -0.5) then                          ! if we have no upper limit
!                    point2 = point1  * (EXT - 1)                ! the extrapolate the maximum amount
!                else
!                    point2 = (limit - point1)/2                 ! otherwise bisect
!                end if
!            else
!                point2 = (-slope2 * point3 * point3)/(B + sqrt(sqrtnumber))
!                if ((limit > -0.5) .and. ((point2 + point1) > limit))then          ! extraplation beyond max?
!                    point2 = (limit - point1)/2                 ! bisect
!                elseif ((limit < -0.5) .and. ((point2 + point1) > (point1 * EXT)))then       ! extrapolation beyond limit
!                    point2 = point1 * (EXT - 1.0)               ! set to extrapolation limit
!                elseif (point2 < (-point3 * INT))then
!                    point2 = -point3 * INT
!                elseif ((limit > -0.5) .and. (point2 < (limit - point1)*(1.0 - INT)))then   ! too close to limit?
!                    point2 = (limit - point1 ) * (1.0 - INT)
!                end if
!            end if
!
!            cost3 = cost2
!            slope3 = slope2
!            point3 = -point2               
!            point1  = point1  + point2
!
!            stemp = point2 * search_direction
!            nn_params = nn_params + stemp                       ! update current estimates
!
!            cost2 = costandgradient(gradient2,nn_params,input_layer_size,hidden_layer_size,num_labels, inputdata, y, lambda)
!
!            M = M - 1.0
!            i = i + (length<0)                                  ! count epochs?!
!            slope_vector = matmul(transpose(gradient2),search_direction)
!            slope2 = slope_vector(1,1)                          !convert to scalar
!        end do                                                  ! end of line search
!
!        if (success == 1)then                                   ! if line search succeeded
!            cost1 = cost2
!            fX = cost1
!            print *, "Iteration: ", i, " | Cost: ", cost1
!
!            ! this calculation was split up to make it easier to work with
!            sd_calc_1 = matmul(transpose(gradient2),gradient2)
!            sd_calc_2 = matmul(transpose(gradient1),gradient2)
!            sd_calc_3 = matmul(transpose(gradient1),gradient1)
!            sd_calc_4 = (sd_calc_1(1,1) - sd_calc_2(1,1)) / sd_calc_3(1,1)
!            sd_calc_5 = sd_calc_4 * search_direction
!            search_direction = sd_calc_5 - gradient2
!
!            tmp = gradient1
!            gradient1 = gradient2
!            gradient2 = tmp                                     ! swap derivatives
!            slope_vector = matmul(transpose(gradient1),search_direction)
!            slope2 = slope_vector(1,1)                          !convert to scalar
!            if(slope2 > 0)then                                  ! new slope must be negative
!                search_direction = -gradient1                   ! otherwise use steepest direction
!                slope_vector = matmul(-transpose(search_direction),search_direction)
!                slope2 = slope_vector(1,1)                      !convert to scalar
!            end if
!            mintemp = slope1/(slope2 - 2.2251D-308)  !2.2551D-308 is min value double precision float
!            minstuff = min(RATIO, mintemp)
!            point1  = point1  * minstuff                         ! slope ratio but max RATIO
!            slope1 = slope2
!            ls_failed = 0                                        ! this line search did not fail
!        else
!            nn_params = backup_params
!            cost1 = cost_backup
!            gradient1 = gradient_backup                         ! restore point from before failed line search
!            if (ls_failed == 1 .or. i > abs(length))then        ! line search failed twice in a row
!                exit                                            ! or we ran out of time, so we give up
!            end if
!            tmp = gradient1
!            gradient1 = gradient2
!            gradient2 = tmp                                    ! swap derivatives
!            search_direction = -gradient1                      ! try steepest
!            slope_vector = matmul(-transpose(search_direction),search_direction)
!            slope1 = slope_vector(1,1)                         ! convert to scalar
!            point1  = 1.0 / (1.0 - slope1)
!            ls_failed = 1                                      ! this line search failed
!        end if
!
!    end do
!
!    fmincgfun = fX
!end function fmincgfun
    
    
end module fmincg