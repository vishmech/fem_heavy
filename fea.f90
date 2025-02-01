module fea

    !! This module contains procedures that are common to FEM-analysis
    
    implicit none
    save
    private
    public :: displ, initial, buildload, buildstiff, enforce, recover

contains

!
!--------------------------------------------------------------------------------------------------
!
    subroutine initial
        
        !! This subroutine is mainly used to allocate vectors and matrices
        
        use fedata
        use link1
        use plane42
               
!        Hint for continuum elements:
!        integer, parameter :: mdim = 8 
!        integer, dimension(mdim) :: edof

        ! This subroutine computes the number of global equation,
        ! half bandwidth, etc and allocates global arrays.

        ! Calculate number of equations
        
        integer :: e, bwww(ne),max_nn,min_nn,y_max_nn,x_min_nn 
        

        neqn = 2*nn
        if (.not. banded) then
            allocate (kmatr(neqn, neqn))
        else
        ! Calculate Bandwidth
        do e = 1, ne
            max_nn=maxval(element(e)%ix)
            min_nn=minval(element(e)%ix)
            y_max_nn=max_nn*2
            x_min_nn=min_nn*2 - 1
            bwww(e)=y_max_nn - x_min_nn + 1
        end do      
        bw=maxval(bwww)
        print*,'Bandwidth=',bw,'neqn=',neqn    
        allocate (kmatr(bw, neqn))
        end if
        

        nngg=4   ! total no. gauss points = 2*2 = 4
        cyc=13        
        total_incr=2000 ! total no. increments

        ! field variables are stored for each gauss points and time increments

        allocate (p(neqn), del_p(neqn),del_p_corrected1(neqn), del_p_corrected2(neqn), &
                 p_final(neqn), d(neqn), del_d(neqn), R_int(neqn))
        allocate (strain(ne, 3, nngg, total_incr), stress(ne, 3, nngg, total_incr))
        allocate (EQP_strain(ne, nngg, total_incr), vonstress(ne, nngg, total_incr))
        allocate (yield_stress(ne, nngg, total_incr),del_lamda(ne, nngg, total_incr))
        allocate (plastic_check(ne, nngg, total_incr))


        ! Initial stress and strain
        strain = 0
        stress = 0 
        EQP_strain=0
        vonstress = 0  
        p=0
        del_p=0
        d=0
        del_d=0
        p_final=0
        R_int=0


        ! plastic_check=0 ELASTIC
        ! plastic_check=1 PLASTIC

        plastic_check=0
	print*, ne
        

    end subroutine initial
!
!--------------------------------------------------------------------------------------------------
!
    subroutine displ

        !! This subroutine calculates displacements

        use fedata
        use numeth
        use processor
        integer :: e,i,idof,j,eleme
        real(wp), dimension(:), allocatable :: plotval
        !real(wp) :: del_p_corrected(neqn),R_int(neqn)
        real(wp) :: hhh(neqn)
    


      OPEN(unit=1,FILE='res1')
      OPEN(unit=2,FILE='res2')

        ! initialize yield stress for 1st step

        do e=1,ne
        yield_stress(e,:,1)=mprop(element(e)%mat)%sgy
        del_lamda(e,:,1)=0.0
        EQP_strain(e,:,1)=0.0
        end do


        ! Build load-vector in P_final
        call buildload
        

        d=0
        R_int=0
        p=0
        del_p_corrected2=p-R_int
        del_p=p_final/(total_incr)
        

        !!! START INCREMENT !!!
        do n_incr=2,total_incr
            

               

            ! Build stiffness matrix
            call buildstiff

            del_p_corrected1=del_p+del_p_corrected2    
           ! del_p_corrected1=del_p

            ! Remove rigid body modes
            ! Enforcing boundary cond on K and del_p_corrected1

                 if (.not. banded) then
                    do i = 1, nb
                        idof = int(2*(bound(i,1)-1) + bound(i,2))
                        
                     ! p(1:neqn) = p(1:neqn) - kmatr(1:neqn, idof) * bound(i, 3) ! for non zero pres.  displacement
                        del_p_corrected1(idof) = bound(i, 3)
                        kmatr(1:neqn, idof) = 0.0_wp
                        kmatr(idof, 1:neqn) = 0.0_wp
                        kmatr(idof, idof) = 1.0_wp
                    end do
                 else
                        do i = 1, nb
                            idof = int(2*(bound(i,1)-1) + bound(i,2))
        
                            if (abs(bound(i, 3)) .gt. 10.0**(-5)) then

                            print*, 'Doesnt WORK with NON-ZERO DISPLACEMENT BC, ADD CODE' 
                            stop
                            end if
                            
                         !!    p(1:neqn) = p(1:neqn) - K_bound(1:neqn)* bound(i, 3)
                            
                            del_p_corrected1(idof) = bound(i, 3)                   
                            kmatr(1:bw, idof) = 0
                            do j=1,bw
                            if ((idof-j).lt.0) exit
                            kmatr(j,idof-j+1)=0.0_wp
                            end do
                            kmatr(1, idof) = 1.0_wp
                        end do
                end if


            ! Solving del_d

            hhh=del_p_corrected1

            if (.not. banded) then
            ! Factor stiffness matrix
                call factor(kmatr)
            ! Solve for displacement vector
                call solve(kmatr,hhh)
            else    
            ! Factor stiffness matrix
                call bfactor(kmatr)
            ! Solve for displacement vector
                call bsolve(kmatr,hhh)        
            end if  

            ! Transfer results
            del_d(1:neqn) = hhh(1:neqn)
            
            d(1:neqn)=d(1:neqn)+del_d(1:neqn)

            ! Recover stress
            call recover


           ! p=p+del_p

            p=p+del_p


            
            del_p_corrected2=p-R_int

            ! eleme=1

            ! WRITE(1,3)n_incr,EQP_strain(eleme, 1, n_incr),yield_stress(eleme, 1, n_incr),vonstress(eleme, 1, n_incr), & 
            !             strain(eleme,1,1, n_incr), &
            !             strain(eleme,2,1, n_incr),strain(eleme,3,1, n_incr),stress(eleme,1,1, n_incr), &
            !              stress(eleme,2,1, n_incr),stress(eleme,3,1, n_incr)
            ! 3 FORMAT(i8,3x,e17.5,e17.5,e17.5,e17.5,e17.5,e17.5,e17.5,e17.5,e17.5)
               ! 511 node is at 5,5

              WRITE(2,4)n_incr,d(196),p(196),d(195),p(195)
             4 FORMAT(i8,3x,e17.5,e17.5,e17.5,e17.5)
         !  print*,stress(eleme,1,2, n_incr)

           call plot( undeformed + deformed,device=1,wait=.false.,title="Deformed",iter=n_incr )
            call plot( elements , device=2,wait=.false.,eval=vonstress(:, 1, n_incr), title="von Mises stress", &
            iter=n_incr, legend=.true. ) 
            call plot( elements , device=3,wait=.false.,eval=EQP_strain(:, 1, n_incr), title="EQP_strain", &
            iter=n_incr, legend=.true. ) 
        !    call plot( elements , device=4,wait=.false.,eval=stress(:, 2,1, n_incr), title="stress22", &
        !    iter=n_incr, legend=.true. ) 
        !    call plot( elements , device=5,wait=.false.,eval=stress(:, 3,1, n_incr), title="stress12", &
        !    iter=n_incr, legend=.true. ) 



        end do
            allocate(plotval(ne))
            plotval=0

            do e=1,ne

                plotval(e)=(vonstress(e, 1, total_incr)+vonstress(e, 2, total_incr)+vonstress(e, 3, total_incr)+&
                                        vonstress(e, 4, total_incr))/4.0_wp

            end do

            call plot( elements , device=matlab,eval=plotval(:), title="von Mises stress", &
            legend=.true. ) 

            plotval=0

            do e=1,ne

                plotval(e)=(EQP_strain(e, 1, total_incr)+EQP_strain(e, 2, total_incr)+EQP_strain(e, 3, total_incr)+&
                            EQP_strain(e, 4, total_incr))/4.0_wp

            end do
            call plot( elements , device=matlab,eval=plotval(:), title="EQP_strain", &
            legend=.true. ) 


            call plot( elements , device=matlab,eval=vonstress(:, 1, total_incr), title="von Mises stress1", &
            legend=.true. ) 
            call plot( elements , device=matlab,eval=vonstress(:, 2, total_incr), title="Evon Mises stress2", &
            legend=.true. ) 
            call plot( elements , device=matlab,eval=vonstress(:, 3, total_incr), title="von Mises stress3", &
            legend=.true. ) 
            call plot( elements , device=matlab,eval=vonstress(:, 4, total_incr), title="von Mises stress4", &
            legend=.true. )


            call plot( elements , device=matlab,eval=EQP_strain(:, 1, total_incr), title="EQP_strain1", &
            legend=.true. ) 
            call plot( elements , device=matlab,eval=EQP_strain(:, 2, total_incr), title="EQP_strain2", &
            legend=.true. ) 
            call plot( elements , device=matlab,eval=EQP_strain(:, 3, total_incr), title="EQP_strain3", &
            legend=.true. ) 
            call plot( elements , device=matlab,eval=EQP_strain(:, 4, total_incr), title="EQP_strain4", &
            legend=.true. )
             call plot( elements , device=matlab,eval=stress(:, 1,1, total_incr), title="stress11", &
             legend=.true. ) 
            call plot( elements , device=matlab,eval=stress(:, 2,1, total_incr), title="stress22", &
            legend=.true. ) 
            call plot( elements , device=matlab,eval=stress(:, 3,4, total_incr), title="stress12-4GP", &
            legend=.true. ) 
            call plot( elements , device=matlab,eval=strain(:, 1,1, total_incr), title="strain11", &
            legend=.true. ) 
            call plot( elements , device=matlab,eval=strain(:, 2,1, total_incr), title="strain22", &
            legend=.true. ) 
            call plot( elements , device=matlab,eval=strain(:, 3,1, total_incr), title="strain12", &
            legend=.true. ) 

    end subroutine displ
!
!--------------------------------------------------------------------------------------------------
!
    subroutine buildload

        !! This subroutine builds the global load vector

        use fedata
        use plane42

        integer :: i,j,e
        integer :: nen,eface
        integer, dimension(mdim) :: edof
        real(wp), dimension(mdim) :: xe
        real(wp), dimension(mdim) :: re
        real(wp) :: fe,thk
        ! Build load vector
        p_final(1:neqn) = 0
        do i = 1, np         
            select case(int(loads(i, 1)))
            case( 1 )
                ! Build nodal load contribution
                if(int(loads(i,3))==1) then
                   ! p( int(2*loads(i,2))-1)=loads(i,4);
                    p_final(int(2*loads(i,2))-1)=p_final(int(2*loads(i,2))-1)+loads(i,4)
                end if
                if(int(loads(i,3))==2) then
                  !  p( int(2*loads(i,2)))=loads(i,4);
                    p_final( int(2*loads(i,2)))=p_final(int(2*loads(i,2)))+loads(i,4)
                end if                
            case( 2 )
            	! Build uniformly distributed surface (pressure) load contribution
                e=int(loads(i,2))
                eface=int(loads(i,3))
                fe=loads(i,4)   
                thk = mprop(element(e)%mat)%thk             
                ! Find coordinates and degrees of freedom
                nen = element(e)%numnode
                do j = 1, nen
                    xe(2*j-1) = x(element(e)%ix(j),1)
                    xe(2*j  ) = x(element(e)%ix(j),2)
                    edof(2*j-1) = 2 * element(e)%ix(j) - 1  
                    edof(2*j)   = 2 * element(e)%ix(j)
                end do
                !
                call plane42_re(xe, eface, fe, thk, re)
                p_final(edof)=p_final(edof)+re
            case default
                print *, 'ERROR in fea/buildload'
                print *, 'Load type not known'
                stop
            end select
        end do
    end subroutine buildload
!
!--------------------------------------------------------------------------------------------------
!
    subroutine buildstiff

        !! This subroutine builds the global stiffness matrix from
        !! the local element stiffness matrices

        use fedata
        use link1
        use plane42

        integer :: e, i, j
        integer :: nen 
! Hint for system matrix in band form:
        integer :: irow, icol
        integer, dimension(mdim) :: edof
        real(wp), dimension(mdim) :: xe,de
        real(wp), dimension(mdim, mdim) :: ke 
! Hint for modal analysis:
!        real(wp), dimension(mdim, mdim) :: me 
        real(wp) :: young, area, sgy, n_expo 
! Hint for modal analysis and continuum elements:
        real(wp) :: nu, dens, thk

        ! Reset stiffness matrix
        if (.not. banded) then
            kmatr = 0
        else
            do irow=1,bw
               do icol=1,neqn
                  kmatr(irow,icol)=0
               end do
            end do
            
        end if

        do e = 1, ne
          
            ! Find coordinates and degrees of freedom
            nen = element(e)%numnode
            do i = 1, nen
                 xe(2*i-1) = x(element(e)%ix(i),1)
                 xe(2*i  ) = x(element(e)%ix(i),2)
                 edof(2*i-1) = 2 * element(e)%ix(i) - 1  
                 edof(2*i)   = 2 * element(e)%ix(i)
                 de(2*i-1) = d(edof(2*i-1))
                 de(2*i)   = d(edof(2*i))
            end do

            ! Gather material properties and find element stiffness matrix
            select case( element(e)%id )
            case( 1 )
                 young = mprop(element(e)%mat)%young
                 area  = mprop(element(e)%mat)%area
                 call link1_ke(xe, young, area, ke)
            case( 2 )
                 young = mprop(element(e)%mat)%young
                 area  = mprop(element(e)%mat)%area
                 nu    = mprop(element(e)%mat)%nu
                 thk   = mprop(element(e)%mat)%thk
                 sgy   = mprop(element(e)%mat)%sgy
                 n_expo= mprop(element(e)%mat)%n_expo
                 
                 call plane42_ke(xe, de, young, nu, thk, sgy, n_expo, stress(e,:,:,n_incr-1), & 
                                yield_stress(e,:,n_incr-1), plastic_check(e,:,n_incr-1),ke)
            end select

            ! Assemble into global matrix
            if (.not. banded) then
                do i = 1, 2*nen
                    do j = 1, 2*nen
                        kmatr(edof(i), edof(j)) = kmatr(edof(i), edof(j)) + ke(i, j)
                    end do
                end do
            else
                 do i = 1, 2*nen
                    do j = 1, 2*nen 
                     if (edof(i) .ge. edof(j)) then     
                        kmatr(1+edof(i)-edof(j),edof(j))=kmatr(1+edof(i)-edof(j),edof(j))+ ke(i, j)
                     end if
                    end do
                end do 
                

            end if 
        end do
    end subroutine buildstiff
!
!--------------------------------------------------------------------------------------------------
!
    subroutine enforce

        !! This subroutine enforces the support boundary conditions

        use fedata

        integer :: i, idof, j
        real(wp) :: penal,toler2
        real(wp), dimension(neqn) :: K_bound
        toler2=10.0**(-5)
        
        ! Correct for supports
        if (.not. banded) then
            if (.not. penalty) then
                do i = 1, nb
                    idof = int(2*(bound(i,1)-1) + bound(i,2))
                    
                    p(1:neqn) = p(1:neqn) - kmatr(1:neqn, idof) * bound(i, 3)
                    p(idof) = bound(i, 3)
                    kmatr(1:neqn, idof) = 0
                    kmatr(idof, 1:neqn) = 0
                    kmatr(idof, idof) = 1
                end do
            else
                penal = penalty_fac*maxval(kmatr)
                do i = 1, nb
                    idof = int(2*(bound(i,1)-1) + bound(i,2))
                    kmatr(idof, idof) = kmatr(idof, idof) + penal
                    p(idof) = penal * bound(i, 3)  
                end do  
            end if
        else
            if (.not. penalty) then
            !!  K_bound(1:neqn)=0
                do i = 1, nb
                    idof = int(2*(bound(i,1)-1) + bound(i,2))
!
                    if (abs(bound(i, 3)) .gt. toler2) then

                       print*, 'Doesnt WORK with NON-ZERO DISPLACEMENT BC, ADD CODE' 
                       stop
                    end if
                    
                    
                !!    do j=1,bw
                !!    if ((idof-j).lt.0) exit
                !!       K_bound(idof+j-1)=kmatr(j, idof)
                !!       K_bound(idof-j+1)=kmatr(j,idof-j+1)                                                         
                !!    end do
                    
                   
                    
                !!    p(1:neqn) = p(1:neqn) - K_bound(1:neqn)* bound(i, 3)
                    
                    p(idof) = bound(i, 3)
                    
                    
                    kmatr(1:bw, idof) = 0
                    do j=1,bw
                    if ((idof-j).lt.0) exit
                    kmatr(j,idof-j+1)=0
                    end do
                    kmatr(1, idof) = 1
                end do
            else        
        
            print *, 'ERROR in fea/enforce'
            print *, 'Penalty in Band form not implemented -- you need to add your own code here'
            stop
            
            end if
            
            
        end if
    end subroutine enforce
!
!--------------------------------------------------------------------------------------------------
!
    subroutine recover

        !! This subroutine recovers the element stress, element strain, 
        !! and nodal reaction forces

        use fedata
        use link1
        use plane42

        integer :: e, i, nen
        integer :: edof(mdim)
        real(wp), dimension(mdim) :: xe, del_de,de, re_int
        real(wp) :: young, area 
        real(wp):: nu, thk,sgy,n_expo
        real(wp) :: stress_dumm(3,4),strain_dumm(3,4),yield_stress_dumm(4), &
                        EQP_strain_dummy(4),del_lamda_dumm(4)
        integer :: plastic_check_dumm(4)

        ! Reset force vector
        R_int = 0
        do e = 1, ne

          
            ! Find coordinates etc...
            nen = element(e)%numnode
            do i = 1,nen
                xe(2*i-1) = x(element(e)%ix(i), 1)
                xe(2*i)   = x(element(e)%ix(i), 2)
                edof(2*i-1) = 2 * element(e)%ix(i) - 1
                edof(2*i)   = 2 * element(e)%ix(i)
                del_de(2*i-1) = del_d(edof(2*i-1))
                del_de(2*i)   = del_d(edof(2*i))
                de(2*i-1) = d(edof(2*i-1))
                de(2*i)   = d(edof(2*i))
            end do

            ! Find stress and strain
            select case( element(e)%id )
            case( 1 )
                print*,"Currently, TRUSS Suppressed! Check?"
                stop
                ! young = mprop(element(e)%mat)%young
                ! area  = mprop(element(e)%mat)%area
                ! call link1_ke(xe, young, area, ke)
                ! p(edof(1:2*nen)) = p(edof(1:2*nen)) + matmul(ke(1:2*nen,1:2*nen), de(1:2*nen))
                ! call link1_ss(xe, de, young, estress, estrain)
                ! stress(e, 1:3) = estress
                ! strain(e, 1:3) = estrain
                ! principal_stress(e, 1:3)=estress
            case( 2 )
                young = mprop(element(e)%mat)%young
                area  = mprop(element(e)%mat)%area
                nu    = mprop(element(e)%mat)%nu
                thk   = mprop(element(e)%mat)%thk
                sgy   = mprop(element(e)%mat)%sgy
                n_expo   = mprop(element(e)%mat)%n_expo
               

               ! call plane42_ke(xe, de, young, nu, thk, sgy, n_expo, stress(e,:,:,n_incr), & 
               !                    yield_stress(e,:,n_incr), ke)
               ! p(edof(1:2*nen)) = p(edof(1:2*nen)) + matmul(ke(1:2*nen,1:2*nen), de(1:2*nen))

                stress_dumm=stress(e,:,:,n_incr-1)
                yield_stress_dumm=yield_stress(e,:,n_incr-1)
                strain_dumm=strain(e,:,:,n_incr-1)
                EQP_strain_dummy=EQP_strain(e,:,n_incr-1)
                del_lamda_dumm=del_lamda(e,:,n_incr-1)
                plastic_check_dumm=plastic_check(e,:,n_incr-1)
                call plane42_ss(xe, del_de, young, nu, thk, sgy, n_expo, re_int, stress_dumm, &
                    yield_stress_dumm,strain_dumm, EQP_strain_dummy,del_lamda_dumm,plastic_check_dumm)

                        !Transfer values

                stress(e,:,:,n_incr)=stress_dumm
                yield_stress(e,:,n_incr)=yield_stress_dumm
                strain(e,:,:,n_incr)=strain_dumm
                EQP_strain(e,:,n_incr)=EQP_strain_dummy
                del_lamda(e,:,n_incr)=del_lamda_dumm
                plastic_check(e,:,n_incr)=plastic_check_dumm

                R_int(edof(1:2*nen))=R_int(edof(1:2*nen)) + re_int(1:2*nen)

                do i=1,4
                    call cal_von_stress(stress(e,:,i,n_incr),vonstress(e,i,n_incr))
                end do

           
            end select
        end do
    end subroutine recover



end module fea
