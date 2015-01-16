module support
    
    
    contains
    
    subroutine representative_matrix(input_matrix, output_matrix)
        real, allocatable, intent(in) :: input_matrix(:,:)
        real, allocatable, intent(inout) :: output_matrix(:,:)
        integer :: i,j
    
        do i=1,size(input_matrix,1)
            temp = input_matrix(i,1)
            do j=1,K
                if(temp == j)then
                    output_matrix(i,j) = 1
                else
                    output_matrix(i,j) = 0
                end if
            end do
        end do
    
    end subroutine representative_matrix
    
    subroutine matrix_multiply(first_matrix, first_transposed, second_matrix, second_transposed, result_matrix)
        real, intent(in) :: first_matrix(:,:), second_matrix(:,:)
        real, intent(inout) :: result_matrix(:,:)
        integer, intent(in) :: first_transposed, second_transposed
        integer :: a1,a2,b1,b2
        
        a1 = size(first_matrix,1)
        a2 = size(first_matrix,2)
        b1 = size(second_matrix,1)
        b2 = size(second_matrix,2)
              
        
        if (first_transposed .AND. second_transposed) then
            if (a1 /= b2) then
                print *, "Problem with matrix sizes: ",a2,"x",a1," * ",b2,"x",b1
            end if
            call sGEMM('T','T',a2,b1,b2,1.0,first_matrix,a1,second_matrix,b1,0.0,result_matrix,a2) 
        else if (.NOT. first_transposed .AND. second_transposed) then
            if (a2 /= b2) then
                print *, "Problem with matrix sizes: ",a1,"x",a2," * ",b2,"x",b1
            end if            
            !a1 = size(a_2_ones,1)
            !a2 = size(a_2_ones,2)
            !b1 = size(Theta2,2)
            !b2 = size(Theta2,1)
            !call sGEMM('N','T',a1,b2,b1,1.0,a_2_ones,a1,Theta2,b2,0.0,z_3,a1)
            call sGEMM('N','T',a1,b1,b2,1.0,first_matrix,a1,second_matrix,b1,0.0,result_matrix,a1)             
        else if (first_transposed .AND. .NOT. second_transposed) then
            if (a1 /= b1) then
                print *, "Problem with matrix sizes: ",a2,"x",a1," * ",b1,"x",b2
            end if            
            !a1 = size(delta_3,2)
            !a2 = size(delta_3,1)
            !b1 = size(a_2_ones,1)
            !b2 = size(a_2_ones,2)
            !call sGEMM('T','N',a1,b2,b1,1.0,delta_3,a2,a_2_ones,b1,0.0,gradient_2,a1)
            call sGEMM('T','N',a2,b2,b1,1.0,first_matrix,a1,second_matrix,b1,0.0,result_matrix,a2)             
        else if (.NOT. first_transposed .AND. .NOT. second_transposed) then
            if (a2 /= b1) then
                print *, "Problem with matrix sizes: ",a1,"x",a2," * ",b1,"x",b2
            end if            
            !a1 = size(delta_3,1)
            !a2 = size(delta_3,2)
            !b1 = size(Theta2,1)
            !b2 = size(Theta2,2)
            !call sGEMM('N','N',a1,b2,b1,1.0,delta_3,a1,Theta2,b1,0.0,temp_delta_2,a1)            
            call sGEMM('N','N',a1,b2,b1,1.0,first_matrix,a1,second_matrix,b1,0.0,result_matrix,a1)   
        end if
            
    end subroutine matrix_multiply
    
    
    subroutine read_file (UnitNum, FileName, NumRows, NumCols, Array )
        implicit none
        
        integer, intent (in) :: UnitNum
        character (len=*), intent (in) :: FileName
        integer, intent (in) :: NumRows, NumCols
        !real, dimension (1:NumRows, 1:NumCols), intent (out) :: Array
        real, allocatable :: Array(:,:)

        character (len=300) :: line
        integer :: i, j

        !allocate(Array(1:NumRows, 1:NumCols))  already allocated interestingly
      
        open (unit=UnitNum, file=FileName, status='old', action='read' )

        do i=1, NumRows
            read (UnitNum, *) (Array (i, j), j=1,NumCols)
        end do

        close (UnitNum)

        return

    end subroutine read_file
  
    
end module support