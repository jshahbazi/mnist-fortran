program dir
    use support
    use fmincg
    !use iso_c_binding
    
    implicit none

    !interface
    !    
    !    integer(c_int) function nn_socket(domain, protocol) bind(C,name="nn_socket")
    !        use iso_c_binding
    !        integer(c_int), value, intent(in) :: domain 
    !        integer(c_int), value, intent(in) :: protocol
    !    end function nn_socket
    !
    !end interface
    
    
    real, allocatable :: X(:,:), y(:,:), X_cv(:,:), y_cv(:,:), predictions(:,:), results(:,:)
    integer :: i, j
    integer :: m,n,k,num_labels
    integer :: hidden_layer_size, input_layer_size
    integer temp
    
    real, allocatable :: Theta1(:,:), Theta2(:,:)
    real, allocatable :: nn_params(:), nn_temp1(:), nn_temp2(:)
    real, allocatable :: nn_params2(:,:)

    real :: cost, grad
    real :: accuracy_train, accuracy_cv
    real :: timer1, timer2

    
    
    !!c binding experiment
    !integer(c_int) :: domain, protocol, socket_num, sock2
    !
    !
    !enum, bind(c)
    !    enumerator :: AF_SP = 1
    !    enumerator :: NN_PUB = 32
    !end enum
    !
    !socket_num = nn_socket(AF_SP,NN_PUB)
    !sock2 = nn_socket(AF_SP,NN_PUB)
    !
    !print *,'socket: ',socket_num
    !print *,'socket: ',sock2
    !print *,'AF_SP: ',AF_SP
    
    
    
    
    
    
    !call mkl_set_num_threads(4)
    
    
    m=22000  !number of input examples
    !m=100
    k=784  !input layer size
    n=500  !hidden layer size
    num_labels=10  !output layer size
    
    allocate(X(1:m, 1:k))
    allocate(predictions(size(X,1),1))
    allocate(results(size(X,1),1))

    call read_file (66, 'c:\train_X.csv', m, k, X)
    
    X=X/255
    
    allocate(y(1:m,1:1))
    call read_file (66, 'c:\train_y.csv', m, 1, y)
    
    allocate(Theta1(1:n,1:(k+1)))  !+1 for the 1's
    allocate(Theta2(1:num_labels,1:(n+1))) !+1 for the 1's
    
    !call RANDOM_SEED
    !call RANDOM_NUMBER(Theta1)
    !call RANDOM_NUMBER(Theta2)
    call read_file (66, 'c:\theta1.csv', n, k+1, Theta1)
    call read_file (66, 'c:\theta2.csv', num_labels, n+1, Theta2)

    allocate(nn_params(1:(n*(k+1)+num_labels*(n+1))))
    allocate(nn_temp1(1:(n*(k+1))))
    allocate(nn_temp2(1:(num_labels*(k+1))))

    nn_temp1 = reshape(Theta1, (/ size(Theta1) /))
    nn_temp2 = reshape(Theta2, (/ size(Theta2) /))
    
    nn_params(1:size(nn_temp1))=nn_temp1(:)
    nn_params((size(nn_temp1)+1):) = nn_temp2(:)
    
    allocate(nn_params2(size(nn_params),1))
    nn_params2 = reshape(nn_params, (/SIZE(nn_params), 1/))

    deallocate(nn_temp1)
    deallocate(nn_temp2)
    deallocate(nn_params)

timer1=secnds(0.0)
    
    cost = fmincg2(10000, nn_params2, k, n , num_labels, X, y, 0.01)

timer2=secnds(timer1)
print *,'fmincg completed in: ', timer2, ' seconds.'
    
    hidden_layer_size = n
    input_layer_size = k
    Theta1 = reshape(nn_params2(1:hidden_layer_size*(input_layer_size + 1),1), (/ hidden_layer_size, (input_layer_size + 1) /))
    Theta2 = reshape(nn_params2((1+hidden_layer_size*(input_layer_size + 1)):,1), (/ num_labels, (hidden_layer_size + 1) /))

    call predict(Theta1, Theta2, X, predictions)

    results = (predictions == y)
    accuracy_train = (sum(results)/m) * 100.0
    print '(a,f8.2,a)', "Training accuracy: ",accuracy_train, "%"
    
    !------------------------------------------
    
    m=20000  !number of input examples
    allocate(X_cv(1:m, 1:k))
    allocate(y_cv(1:m,1:1))
    call read_file (66, 'c:\cv_X.csv', m, k, X_cv)
    call read_file (66, 'c:\cv_y.csv', m, 1, y_cv)
    
    X_cv=X_cv/255
    
    call predict(Theta1, Theta2, X_cv, predictions)

    results = (predictions == y_cv)
    accuracy_cv = (sum(results)/m) * 100.0
    print '(a,f8.2,a)', "CV accuracy: ",accuracy_cv, "%"    
    
    
    deallocate(X)
    deallocate(y)
    deallocate(X_cv)
    deallocate(y_cv)
    deallocate(nn_params2)
    deallocate(Theta1)
    deallocate(Theta2)
    deallocate(predictions)
    deallocate(results)
    
    
    print *,'End of program.'
    pause
    

    contains
    

    
    
    subroutine predict(Theta1, Theta2, inputdata, predictions)
        real, allocatable, intent(in) :: Theta1(:,:), Theta2(:,:), inputdata(:,:)
        real, allocatable, intent(inout) :: predictions(:,:)
        real, allocatable :: h1(:,:), h1_ones(:,:), h2(:,:), pre_h1(:,:), pre_h2(:,:), ones(:,:), X(:,:)

        integer :: m
        integer :: a1,a2,b1,b2
    
        m = size(inputdata,1)
    
        allocate(h1(size(inputdata,1),size(Theta1,1)))
        allocate(h1_ones(size(inputdata,1),size(Theta1,1)+1))
        allocate(h2(size(inputdata,1),size(Theta2,1)))
        allocate(pre_h1(size(inputdata,1),size(Theta1,1)))
        allocate(pre_h2(size(inputdata,1),size(Theta2,1)))
        allocate(ones(m,1))
        allocate(X(size(inputdata,1),(size(inputdata,2)+1)))      
    
        ones = 1.0

        X(:,1:1) = ones
        X(:,2:) = inputdata    
    
        call matrix_multiply(X,0,Theta1,1,pre_h1)
        call sigmoid(pre_h1,h1)
    
        h1_ones(:,1:1) = ones
        h1_ones(:,2:) = h1
       
        call matrix_multiply(h1_ones,0,Theta2,1,pre_h2)
        call sigmoid(pre_h2,h2)
    
        predictions = reshape(maxloc(h2,2), (/ m,1 /))
    
    end subroutine predict
    
end program dir

