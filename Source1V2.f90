    use dflib
    !Author : globalpolicy
    !Blog   : c0dew0rth.blogspot.com
    !GitHub : globalpolicy
    !Date   : 06:05PM - 12th August, 2019
    !Purpose: To draw flownet (revised)
    !Changes: 
    !-Head at dam heel is water depth as should be(was assumed 0 before)
    !-Dam was fixed to be the size of a single grid box. Now, width of dam can be specified.
    !-Dam's heel and toe are Dirichlet boundaries i.e. known heads but any node between are treated as unknowns.
    !Parameters : contourInterval, numDys, numDxs, DxsTillDam, DamWidthInDxs, reservoirHead, downstreamBoundaryHead
    
    !Built by observing structure of an example coefficient matrix from my Fortran Simulation notes
    !https://photos.app.goo.gl/BS5C2tp4MNMFp9Dh6 -> example grid
    !https://photos.app.goo.gl/GecBS62BS1hTPNri9 -> matrix construction
    !https://photos.app.goo.gl/WgKvzMWf63CRNBwZA -> final matrix
    !Also, during revision of code, the following problems were taken as a reference:
    !https://photos.app.goo.gl/1p7NWiWeYKTRTBEh6
    !https://photos.app.goo.gl/LF51G7jisALXweCL7
    !Coefficient matrices were constructed for them and solved in Excel:
    !https://drive.google.com/open?id=1I513Eq-tPaADIw_MF5oiTIGo4es4k2f8
    !Great lectures I made my notes off of:
    !https://www.youtube.com/watch?v=FgR3KS51mbA - Prof. Wen Shen, Penn State University
    !https://www.youtube.com/watch?v=i0f2DegLBN8 - Sandip Mazumder
    
    implicit none
    
    integer i,j,col,row
    character(10) char
    real tmp,distance,smallestDistanceIndex,smallestDistance
    integer numDys,numDxs,numUnknowns,DxsTillDam,DamWidthInDxs
    integer(2) viewPortTopLeftx,viewPortTopLefty,viewPortBottomRightx,viewPortBottomRighty,viewPortWidth !viewPortWidth=viewPortHeight
    real reservoirHead,downstreamBoundaryHead
    real A(500,500),B(500,500),X(500,500),AInverse(500,500)
    logical status
    TYPE(wxycoord) wxy
    integer ret
    real H(500,500)
    real EquipotentialX(500,500),EquipotentialY(500,500) !one column for one head, same column's row element from each matrix give the equipotential point x,y
    real searchHead,currentHead,rightHead,upHead,contourInterval
    real interpolatedX,interpolatedY
    integer numberOfContours, contourId
    integer contourPointsCounters(500) !stores the number of equipotential points obtained for each contour
    
    contourInterval=10
    
    numDys=3 !total grid boxes along Y axis, enclosed within boundaries top and bottom
    numDxs=4 !total grid boxes along X axis, enclosed within boundaries left and right
    DxsTillDam=2
    DamWidthInDxs=2
    
    reservoirHead=100
    downstreamBoundaryHead=50
    
    viewPortWidth = 400 !width(=height) of the viewport
    viewPortTopLeftx = 200
    viewPortTopLefty = 80
    
    numUnknowns=numDys*numDxs+(DamWidthInDxs-1)
    numberOfContours=reservoirHead/contourInterval + 1
    
    !generate coefficient matrix A
    !the rows, columns and nodes in the comments are referencing the full domain grid including the boundaries(https://photos.app.goo.gl/nToWZ9p3znyGxvwy8), not the array/matrix used here
    do i=1,numUnknowns !row index
        do j=1,numUnknowns !column index
            if(i.EQ.j) then
            
                A(i,j)=-4
                
                !for immediate right element of diagonal's row (right blob of Laplace computational stencil)
                if(i.LT.numDxs*numDys.AND.MODULO(i,numDxs).EQ.1) then !if current node is on u/s boundary
                    A(i,j+1)=2
                elseif(i.LE.numDxs*numDys.AND.MODULO(i,numDxs).EQ.0) then !if current node is on second-last column (before d/s boundary) and below the surface row
                    A(i,j+1)=0
                    B(i,1)=B(i,1)-downstreamBoundaryHead
                elseif(i.GT.numDxs*numDys.AND.i.EQ.numUnknowns) then !if current node is the final unknown i.e. surface node right beneath dam, right before toe of dam(d/s boundary)
                    A(i,j+1)=0
                    B(i,1)=B(i,1)-0
                else
                    A(i,j+1)=1
                endif
                
                !for immediate left element of diagonal's row (left blob of Laplace computational stencil)
                if(j.GT.1) then !if current node is not the bottom-left one
                    if(i.LT.numDxs*numDys.AND.MODULO(i,numDxs).EQ.1) then !if current node is on u/s boundary
                        A(i,j-1)=0
                    elseif(DamWidthInDxs.GE.2.AND.i.EQ.numDxs*numDys+1) then !if current node is right beneath the dam, one node to right of heel
                        A(i,j-1)=0
                        B(i,1)=B(i,1)-reservoirHead
                    else
                        A(i,j-1)=1
                    endif
                endif
                
                !for far right element of diagonal's row (top blob of Laplace computational stencil)
                if(i.LE.numDxs) then !if current node is in the bottom row
                    A(i,j+numDxs)=2
                elseif(j+numDxs.LE.numDxs*numDys) then !if current row is below the second-top row
                    A(i,j+numDxs)=1
                elseif(j.LE.numDxs*numDys) then !if current row is second-top row(row immediately below surface)
                    if(j.LE.numUnknowns-(DamWidthInDxs-1+numDxs-DxsTillDam-1)) then !if this node is u/s of heel of dam, heel inclusive
                        B(i,1)=B(i,1)+(-reservoirHead)
                    elseif(j.GT.numDxs*numDys-(DamWidthInDxs-1).AND.j.LE.numDxs*numDys) then !if node is between(exclusive) heel and toe
                        A(i,j+(DamWidthInDxs-1))=1
                    endif
                else !if current row is topmost row(unknown nodes directly beneath the dam)
                    B(i,1)=B(i,1)+0
                endif
                
                !for far left element of diagonal's row (bottom blob of Laplace computational stencil)
                if(j.GT.numDxs) then
                    if(j.LE.numDxs*numDys) then
                        A(i,j-numDxs)=1
                    else
                        A(i,j-(DamWIdthInDxs-1))=2
                    endif
                endif
            
            endif
        end do   
    end do
    
    !generated coefficient matrix A written to file for inspection
    open(1,file='A.txt')
    do i=1,numUnknowns
        write(1,'(F5.1,$)') (A(i,j),j=1,numUnknowns)
        write(1,*)
    end do
    close(1)
    
    !generated B vector written to file for inspection
    open(1,file='B.txt')
    do i=1,numUnknowns
        write(1,*) B(i,1)
    end do
    close(1)
    
    call Inverse(A,AInverse,numUnknowns)
    X=matmul(AInverse,B)
    
    !calculated X vector
    open(1,file='X.txt')
    do i=1,numUnknowns
        write(1,*) X(i,1)
    end do
    close(1)
    
    !setup viewport
    viewPortBottomRightx = int2(viewPortTopLeftx+real(numDxs)/real(numDys)*viewPortWidth)
    viewPortBottomRighty = int2(viewPortTopLefty+viewPortWidth)
    call SETVIEWPORT(200,80,viewPortBottomRightx,viewPortBottomRighty)
    status = SETCOLOR(4) !red
    status = FLOODFILL(51,81,4) !flood fill viewport
    status = SETWINDOW(.TRUE.,real(0,8),real(0,8),real(numDxs,8),real(numDys,8))
    status = SETCOLOR(1)

    
    !draw horizontal lines from bottom to top
    do i=0,numDys
        call MOVETO_W(0.0_8,real(i,8),wxy)
        ret = LINETO_W(real(numDxs,8),real(i,8))
    end do
    
    !draw vertical lines from left to right
    do i=0,numDxs
        call MOVETO_W(real(i,8),0.0_8,wxy)
        ret = LINETO_W(real(i,8),real(numDys,8))
    end do
    
    !create Head matrix for whole domain from calculated heads at nodes and the boundary condition values
    !populate non-last row(below surface) of Head matrix from the computed X vector
    do i=1,numDys
        do j=1,numDxs+1
            if(j.EQ.numDxs+1) then !if final column
                H(i,j)=downstreamBoundaryHead
            else
                H(i,j)=X((i-1)*numDxs+j,1)
            endif
        end do
    end do
    !populate the last row (corresponding to the surface heads)
    i=1
    do j=1,numDxs+1
        if(j.LE.DxsTillDam+1) then
            H(numDys+1,j)=reservoirHead
        elseif(j.GT.DxsTillDam+1.AND.j.LT.numDxs+1) then
            H(numDys+1,j)=X(numDxs*numDys+i,1)
            i=i+1
        else
            H(numDys+1,j)=0
        endif
    end do
    
    !output H matrix to file for inspection
    !observe that it is of the form (taking the referenced example at the top)
    !           h1           h2             h3          h4    downstreamBoundaryHead
    !           h5           h6             h7          h8    downstreamBoundaryHead
    !           h9          h10             h11         h12   downstreamBoundaryHead
    !   reservoirHead   reservoirHead  reservoirHead    h13     0
    !where h1,h2,... are computed heads from X vector
    !also observe that H matrix has exactly the same number of elements as nodes in our whole domain grid(i.e. numDys+1 by numDxs+1) 
    open(1,file='H.txt')
    do i=1,numDys+1
        write(1,'(F5.1,$,"  ")') (H(i,j),j=1,numDxs+1)
        write(1,*)
    end do
    close(1)
            
    contourId=0
    do searchHead=0,reservoirHead,contourInterval
        contourId=contourId+1
        !iterate through grid nodes one column at a time, left to right, looking for same potential(searchHead) right and top
        !note that correspondence between our grid coordinate and H matrix's element is
        !Head at coordinate (col,row) of grid  = (row+1,col+1) element of H matrix
        !Makes use of the fact that in a row, head decreases from left to right.
        !Depending on column on grid, in a column, head decreases top to bottom for col<damX and increases for col>=damX
        do col=0,numDxs
10            do row=0,numDys
                currentHead=H(row+1,col+1)
                if(currentHead.NE.searchHead) then
                    if(col<numDxs) then !look toward east node
                        rightHead=H(row+1,col+2)
                        if(searchHead.LE.currentHead.AND.searchHead.GE.rightHead) then
                            interpolatedX=col+1/(rightHead-currentHead)*(searchHead-currentHead)
                            contourPointsCounters(contourId)=contourPointsCounters(contourId)+1
                            EquipotentialX(contourPointsCounters(contourId),contourId)=interpolatedX
                            EquipotentialY(contourPointsCounters(contourId),contourId)=row
                        endif

                    endif
                    if(row<numDys) then !look toward north node
                        upHead=H(row+2,col+1)
                        if(col<DxsTillDam) then !if we're left of dam
                            !head increases towards up
                            if(searchHead.GE.currentHead.AND.searchHead.LE.upHead) then
                                interpolatedY=row+1/(upHead-currentHead)*(searchHead-currentHead)
                                contourPointsCounters(contourId)=contourPointsCounters(contourId)+1
                                EquipotentialY(contourPointsCounters(contourId),contourId)=interpolatedY
                                EquipotentialX(contourPointsCounters(contourId),contourId)=col
                            endif
                        else !if we're at or right of dam's heel
                            !head decreases towards up
                            if(searchHead.LE.currentHead.AND.searchHead.GE.upHead) then
                                interpolatedY=row+1/(upHead-currentHead)*(searchHead-currentHead)
                                contourPointsCounters(contourId)=contourPointsCounters(contourId)+1
                                EquipotentialY(contourPointsCounters(contourId),contourId)=interpolatedY
                                EquipotentialX(contourPointsCounters(contourId),contourId)=col
                            endif
                        endif
                    endif
                else
                    !we're standing at our searchHead!!
                    contourPointsCounters(contourId)=contourPointsCounters(contourId)+1
                    EquipotentialY(contourPointsCounters(contourId),contourId)=row
                    EquipotentialX(contourPointsCounters(contourId),contourId)=col
                endif
            end do
        end do
    end do
    

    !!----DEPRECATED EQUIPOTENTIAL POINT SORTING METHOD----
    !!sort equipotential points in ascending (descending works too) order of x coordinate for proper equipotential line
    !!sketching (left to right in ascending, right to left in descending). 
    !!EDIT : Equipotential curves that are improper functions (two y's for a x) don't behave nicely with
    !!this method though. Only a couple of curves, generally towards the center of the flownet are as such
    !!in my experience.
    !do contourId=1,numberOfContours
    !    do i=1,contourPointsCounters(contourId)
    !        do j=i,contourPointsCounters(contourId)
    !            if(EquipotentialX(j,contourId).LT.EquipotentialX(i,contourId)) then
    !                tmp=EquipotentialY(i,contourId)
    !                EquipotentialY(i,contourId)=EquipotentialY(j,contourId)
    !                EquipotentialY(j,contourId)=tmp
    !                
    !                tmp=EquipotentialX(i,contourId)
    !                EquipotentialX(i,contourId)=EquipotentialX(j,contourId)
    !                EquipotentialX(j,contourId)=tmp
    !            endif     
    !        end do
    !    end do
    !end do
    
    !equipotential points sorting for proper equipotential line sketching
    !first sort equipotential points in descending order of Y coordinate(topmost to bottom-most)
    !then sort based on shortest consecutive distance.
    !better than sorting by x coordinate but not totally devoid of some weird artifacts
    do contourId=1,numberOfContours
        !sort descending in Y coordinate
        do i=1,contourPointsCounters(contourId)
            do j=i,contourPointsCounters(contourId)
                if(EquipotentialY(j,contourId).GT.EquipotentialY(i,contourId)) then
                    !swap Y coordinate
                    tmp=EquipotentialY(i,contourId)
                    EquipotentialY(i,contourId)=EquipotentialY(j,contourId)
                    EquipotentialY(j,contourId)=tmp
                    !swap X coordinate
                    tmp=EquipotentialX(i,contourId)
                    EquipotentialX(i,contourId)=EquipotentialX(j,contourId)
                    EquipotentialX(j,contourId)=tmp
                endif     
            end do
        end do
        !sort based on shortest consecutive distance (basically ordering the shortest path downward, downward coz the array is descending inY)
        do i=1,contourPointsCounters(contourId)-1
            smallestDistance=10000 !very large number
            do j=i+1,contourPointsCounters(contourId)
                !distance between i and j point
                distance = (EquipotentialX(i,contourId)-EquipotentialX(j,contourId))**2 + (EquipotentialY(i,contourId)-EquipotentialY(j,contourId))**2
                if(distance.LT.smallestDistance) then
                    smallestDistance=distance
                    smallestDistanceIndex=j
                endif
            end do
            !swap i+1 point and point at smallestDistanceIndex
            !X
            tmp = EquipotentialX(i+1,contourId)
            EquipotentialX(i+1,contourId) = EquipotentialX(smallestDistanceIndex,contourId)
            EquipotentialX(smallestDistanceIndex,contourId) = tmp
            !Y
            tmp = EquipotentialY(i+1,contourId)
            EquipotentialY(i+1,contourId) = EquipotentialY(smallestDistanceIndex,contourId)
            EquipotentialY(smallestDistanceIndex,contourId) = tmp
        end do    
    end do
    
    !write equipotential points to file for observation
    open(1,file='Isobars.txt')
        do contourId=1,numberOfContours
            write(1,*) (contourId-1)*contourInterval
            write(1,'("(",F5.2,",",F5.2,")      ",$)') (EquipotentialX(i,contourId),EquipotentialY(i,contourId),i=1,contourPointsCounters(contourId))
            write(1,*) ''
        end do 
    close(1)
    
    !draw equipotential lines
    do contourId=1,numberOfContours
        call MOVETO_W(real(EquipotentialX(1,contourId),8),real(EquipotentialY(1,contourId),8),wxy)
        
        !draw a circle on the equipotential point(not necessary)
        status = SETCOLOR(6)
        ret = ELLIPSE_W($GFILLINTERIOR,real(EquipotentialX(1,contourId)-0.015,8),real(EquipotentialY(1,contourId)+0.015,8),real(EquipotentialX(1,contourId)+0.015,8),real(EquipotentialY(1,contourId)-0.015,8))
        status = SETCOLOR(2)
        
        do i=2,contourPointsCounters(contourId)
            !draw a circle on the equipotential point(not necessary)
            status = SETCOLOR(6)
            ret = ELLIPSE_W($GFILLINTERIOR,real(EquipotentialX(i,contourId)-0.015,8),real(EquipotentialY(i,contourId)+0.015,8),real(EquipotentialX(i,contourId)+0.015,8),real(EquipotentialY(i,contourId)-0.015,8))
            status = SETCOLOR(2)
            
            !avoids jumping lines
            call GETCURRENTPOSITION_W(wxy)
            if(abs(wxy.wx-EquipotentialX(i,contourId)).GT.4.OR.abs(wxy.wy-EquipotentialY(i,contourId)).GT.4) then
                call MOVETO_W(real(EquipotentialX(i,contourId),8),real(EquipotentialY(i,contourId),8),wxy)
            endif     
            
            ret = LINETO_W(real(EquipotentialX(i,contourId),8),real(EquipotentialY(i,contourId),8))
            
        end do
    end do
    
    !draw one grid length horizontal line on grid top representing base of dam
    status = SETCOLOR(8)
    call MOVETO_W(real(DxsTillDam,8),real(numDys,8),wxy)
    ret = LINETO_W(real(DxsTillDam+DamWidthInDxs,8),real(numDys,8))

    
    !write heads on grid nodes
    !note that correspondence between our grid coordinate and H matrix's element is
    !head at coordinate (col,row) of grid  = (row+1,col+1) element of H matrix
    status = SETCOLOR(8)
    status = INITIALIZEFONTS()
    status = SETFONT('t''Arial''h12w8pvb')
    !loop for except the bottom and right boundaries and bottom-right node (else they're clipped)
    do i=0,numDxs-1
        do j=1,numDys
            call MOVETO_W(real(i,8),real(j,8),wxy)      
            write(char,'(F6.2)') H(j+1,i+1)
            call OUTGTEXT(char)
        end do
    end do
    !for bottom boundary
    do i=0,numDxs-1
        call MOVETO_W(real(i,8),real(0,8)+0.1,wxy)
        write(char,'(F6.2)') H(0+1,i+1)
        call OUTGTEXT(char)
    end do
    !for right boundary
    do i=1,numDys
        call MOVETO_W(real(numDxs,8)-0.3,real(i,8),wxy)
        write(char,'(F6.2)') H(i+1,numDxs+1)
        call OUTGTEXT(char)
    end do
    !for bottom right node
    call MOVETO_W(real(numDxs,8)-0.3,real(0,8)+0.1,wxy)
    write(char,'(F6.2)') H(0+1,numDxs+1)
    call OUTGTEXT(char)
    
    end program
    
    
    
    
    
    
    
    
subroutine Inverse(X,XI,size)
dimension X(500,500),XI(500,500),originalX(500,500)
integer rowindex,colindex,size,tmpi,tmpj,i,j
real pivotelement

originalX=X

do i=1,size
    do j=1,size
        if(i.EQ.j) then
            XI(i,j)=1
        else
            XI(i,j)=0 !this is for good measure coz by default, XI may have been changed by the caller itself though it is 0 by default
        endif
    end do
end do


do rowindex=1,size !loop over rows of X
    
    pivotelement=X(rowindex,rowindex) !this is crucial!
    
    !divide the row by pivot element
    do colindex=1,size
        X(rowindex,colindex)=X(rowindex,colindex)/pivotelement
    end do
    do colindex=1,size
        XI(rowindex,colindex)=XI(rowindex,colindex)/pivotelement
    end do
    
    
    !perform row operations to make non-pivot elements of this column zero
    do i=1,size !i gives the row number being manipulated
        pivotelement=X(i,rowindex) !pivot element for the currently manipulated row. this i crucial
        
        if(i.NE.rowindex) then !row operations are for non-current rows only
            
            do j=1,size !column index for ith row of X matrix
                X(i,j)=X(i,j)-pivotelement*X(rowindex,j)
            end do
            do j=1,size !column index for ith row of XI matrix
                XI(i,j)=XI(i,j)-pivotelement*XI(rowindex,j)
            end do
            
            
        endif
    end do
    

    
end do

!restore X to original X since upon inversion completion, X has been reduced to an Identity matrix. good practice
X=originalX

end subroutine