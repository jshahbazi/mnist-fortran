module costfunction
    use support
    
    implicit none
    
    contains
    
    real function costandgradient(gradient, nn_params,input_layer_size,hidden_layer_size,num_labels, inputdata, y, lambda)

        real, allocatable :: inputdata(:,:), X(:,:), y(:,:), y_representative(:,:), gradient(:,:), ones(:,:)
        real, allocatable :: a_1(:,:),z_2(:,:),a_2(:,:),z_3(:,:),a_3(:,:),hofx(:,:)
        real, allocatable :: delta_2(:,:), temp_delta_2(:,:), delta_3(:,:), gradient_1(:,:), gradient_2(:,:)
        real, allocatable :: Theta1_grad(:,:), Theta2_grad(:,:)
        integer :: input_layer_size,hidden_layer_size,num_labels,m,l,K,i
        integer :: a1,a2,b1,b2,c1,c2
        real :: lambda, temp, J
        real, allocatable :: nn_params(:,:), Theta1(:,:), Theta2(:,:), temptheta1(:,:),temptheta2(:,:)
        real, allocatable :: a_2_ones(:,:), summation(:,:), summation2(:,:), summation3(:,:),z_2_ones(:,:), sigmoid_z_2(:,:)
    
        real :: timer1, timer2, holder1, holder2
        double precision :: average1 = 0.0
        double precision :: average2 = 0.0
        double precision :: average3 = 0.0
        double precision :: counter = 0.0
    
        
        
        m = size(inputdata,1)
        l = size(y,1)
        K = num_labels

        allocate(Theta1(1:hidden_layer_size,1:(input_layer_size+1)))
        allocate(Theta2(1:num_labels,1:(hidden_layer_size+1)))
        allocate(ones(m,1))
        allocate(X(size(inputdata,1),(size(inputdata,2)+1)))        
        allocate(y_representative(l,num_labels))        
        allocate(z_2(m,hidden_layer_size))
        allocate(a_1(size(X,1),size(X,2)))
        allocate(a_2(size(z_2,1),size(z_2,2)))              
        allocate(a_2_ones(size(a_2,1),(size(a_2,2)+1)))      
        allocate(z_3(size(a_2_ones,1),size(Theta2,1)))  !a_2 * Theta2'      
        allocate(temptheta2(size(Theta2,2),size(Theta2,1)))        
        allocate(a_3(size(z_3,1),size(z_3,2)))        
        allocate(summation(size(z_3,1),size(z_3,2)))
        allocate(summation2(size(z_3,1),size(z_3,2)))
        allocate(summation3(size(z_3,1),size(z_3,2)))
        allocate(delta_3(size(a_3,1),size(a_3,2)))        
        allocate(z_2_ones(size(z_2,1),(size(z_2,2)+1)))        
        allocate(sigmoid_z_2(size(z_2_ones,1),size(z_2_ones,2)))        
        allocate(temp_delta_2(size(delta_3,1),size(Theta2,2)))   
        allocate(delta_2(size(temp_delta_2,1),size(sigmoid_z_2,2)))
        allocate(gradient_2(size(delta_3,2),size(a_2_ones,2)))        
        allocate(gradient_1(size(delta_2,2)-1,size(a_1,2)))
        allocate(Theta2_grad(1:size(gradient_2),1))
        allocate(Theta1_grad(1:size(gradient_1),1))        
        
        allocate(temptheta1(size(Theta1,2),size(Theta1,1)))  
        

        Theta1 = reshape(nn_params(1:hidden_layer_size*(input_layer_size + 1),1), (/ hidden_layer_size, (input_layer_size + 1) /))
        Theta2 = reshape(nn_params((1+hidden_layer_size*(input_layer_size + 1)):,1), (/ num_labels, (hidden_layer_size + 1) /))
        
        !print *,shape(Theta1)   !500 x 785
        !print *,shape(Theta2)   !10  x 501
        
        
        ones = 1.0



        
        !timer1=secnds(0.0)


        do i=1,m
            temp = y(i,1)
            do j=1,K
                if(temp == j)then
                    y_representative(i,j) = 1
                else
                    y_representative(i,j) = 0
                end if
            end do
        end do

        X(:,1:1) = ones
        X(:,2:) = inputdata        
        a_1 = X

        z_2 = 0.0
        
        call matrix_multiply(a_1,0,Theta1,1,z_2)


        a_2 = 0.0
        call sigmoid(z_2,a_2)

        
        a_2_ones(:,1:1) = ones
        a_2_ones(:,2:) = a_2
        
        z_3 = 0.0
        call matrix_multiply(a_2_ones,0,Theta2,1,z_3)
        
        call sigmoid(z_3,a_3)

        summation = 0.0
        summation = -y_representative * log(a_3) - (1 - y_representative) * log(1 - a_3)
        
        J = sum(summation)
        J = J / m
        
        

        
        !do concurrent (i = 1:size(delta_3,1))
        !    delta_3(i,:) = a_3(i,:) - y(1,:)
        !end do
        
        delta_3 = a_3 - y_representative

        z_2_ones(:,1:1) = ones
        z_2_ones(:,2:) = z_2   
        

        sigmoid_z_2 = 0.0
        
     
        call sigmoidgradient(z_2_ones,sigmoid_z_2)  


        !temp_delta_2 = matmul(delta_3,Theta2)

!timer1=secnds(0.0)   
        !a1 = size(delta_3,1)
        !a2 = size(delta_3,2)
        !b1 = size(Theta2,1)
        !b2 = size(Theta2,2)
        !call sGEMM('N','N',a1,b2,b1,1.0,delta_3,a1,Theta2,b1,0.0,temp_delta_2,a1)
        call matrix_multiply(delta_3,0,Theta2,0,temp_delta_2)
!timer2=secnds(timer1)
!average1 = average1 + timer2
!counter = counter + 1
!print '(a, f8.4)','matrix_multiply timer: ', average1/counter        
        
        delta_2 =0.0
        !delta_2 = temp_delta_2 * sigmoid_z_2   !elemental multiplication
        call vsMul( size(temp_delta_2), temp_delta_2, sigmoid_z_2, delta_2 )   !elemental multiplication    

        
        gradient_2 = 0.0
        gradient_1 = 0.0
        
        call matrix_multiply(delta_3,1,a_2_ones,0,gradient_2)
        
        call matrix_multiply(delta_2(:,2:),1,a_1,0,gradient_1)        
        
        gradient_2 = gradient_2 / m
        gradient_1 = gradient_1 / m
        
        gradient_2(:,2:) = gradient_2(:,2:) + Theta2(:,2:) * (lambda / m)
        gradient_1(:,2:) = gradient_1(:,2:) + Theta1(:,2:) * (lambda / m)


        Theta2_grad = reshape(gradient_2, (/size(gradient_2),1/))
        Theta1_grad = reshape(gradient_1, (/size(gradient_1),1/))

        gradient(1:size(Theta1_grad),1:1) = Theta1_grad
        gradient((size(Theta1_grad)+1):,1:1) = Theta2_grad
        
        costandgradient = J
        
        
        !print '(f8.4, a, f8.4)',average1/counter,' ',average2/counter
        !timer2=secnds(timer1)
        !average = average*0.9 + timer2*0.1
        !print '(a, f8.4)','costandgradient timer: ', timer2
    end function costandgradient
    
    
    !subroutine sigmoid(input, output)
    !    real, allocatable :: input(:,:), output(:,:)
    !
    !    output = 1.0 / (1.0 + exp(-input))
    !end subroutine sigmoid
    
    
    
    !elemental real function sgm(input)
    !    !real, allocatable :: input(:,:), output(:,:)
    !    real, intent(in) :: input
    !    !real, intent(out) :: output    
    !
    !    sgm = 1.0 / (1.0 + exp(-input))
    !end function sgm
    
    
    
    subroutine sigmoid(input, output)
        real, allocatable :: input(:,:), output(:,:)
    
        call vsexp( size(input), -input, output )
        output = 1.0/(1.0 + output)
    end subroutine sigmoid   
    
    
    
    
    
    
    subroutine sigmoidgradient(input, output)
        real, allocatable :: input(:,:), output(:,:)

        output = 1.0 / (1.0 + exp(-input))
        output = output * (1.0 - output)   
    end subroutine sigmoidgradient    
  
    !subroutine sigmoidgradient2(input, output)
    !    real, allocatable :: input(:,:), output(:,:)
    !    
    !    !call vsexp( size(input), -input, output )
    !    !output = 1.0/(1.0 + output)
    !    call sigmoid(output, output)
    !    output = output * (1.0 - output)   
    !end subroutine sigmoidgradient2    
    

    
end module costfunction