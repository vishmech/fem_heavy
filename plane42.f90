module plane42

    !! This module contains subroutines specific to the plane42 element 
    !!        
    !! The plane42 element has 4 nodes. Each node has 2 degrees-of-freedom,
    !! namely, displacement along the \(x\)- and \(y\)-coordinate directions.
    !!
    !! The nodes are numbered counter-clockwise as indicated below. The
    !! figure also shows how the element edges are labelled. For example,
    !! edge 1 is the edge between element node 1 and 2.
    !! It also shows how Gauss Points are Numbered. 
    !! IMPORTANT: Field variables at GPs are stored according to this sequence!!
    !!
    !!        N4         E3          N3
    !!          o-------------------o 
    !!          |                   |
    !!          |    GP4      GP3   |
    !!          |                   |
    !!       E4 |                   | E2
    !!          |                   |
    !!          |                   |
    !!          |    GP1      GP2   |
    !!          |                   |
    !!          o-------------------o
    !!        N1         E1          N2
    !!             
    !!
    !! `N1` = element node 1, `N2` = element node 2, etc  
    !! `E1` = element face 1, `E2` = element face 2, etc  
    
    use types
    implicit none
    save
    
    private
    public :: plane42_ke, plane42_re, plane42_ss, cal_von_stress

contains

    subroutine plane42_ke(xe, de, young, nu, thk, sgy, n_expo, estress, eyield_stress, eplastic_check, ke)

        !! This subroutine constructs the stiffness matrix for
        !! a rectangular 4-noded quad element.

        real(wp), intent(in) :: young
            !! Young's Modulus for this element
        real(wp), intent(in) :: nu
            !! Poisson's Ratio for this element
        real(wp), intent(in) :: thk
            !! Thickness of this element
        real(wp), intent(in) :: sgy
            !! Initial Yield Stress    
        real(wp), intent(in) :: n_expo
            !! Hardening Exponent 
        real(wp), dimension(:), intent(in) :: xe, de
            !! Nodal coordinates of this element in undeformed configuration
            !!
            !! * `xe(1:2)` = \((x,y)\)-coordinates of element node 1
            !! * `xe(3:4)` = \((x,y)\)-coordinates of element node 1
            !! * `xe(5:6)` = \((x,y)\)-coordinates of element node 2
            !! * `xe(7:8)` = \((x,y)\)-coordinates of element node 2
            !!
            !! See also [[plane42rect]]
        real(wp), intent(in) :: estress(3,4),eyield_stress(4)
        integer, intent(in) :: eplastic_check(4)

        real(wp), dimension(:,:), intent(out) :: ke
            !! Stiffness matrix

        real(wp) :: cmat(3,3),C_ep(3,3),bmat(3,8),jac(2,2),n(2,8),fact,weight(2), Et,&
                   gconst,detjac,xi,eta,von_mises_stress_in_GP,yield_toler,h_hard, dF_dsigma(1,3)
        integer :: ng,i,j, GP
    
        yield_toler=1.0E-9

        ! build constitutive matrix (plane stress)
        cmat = 0
        fact = young/(1-nu**2)
        cmat(1,1) = fact
        cmat(1,2) = fact*nu
        cmat(2,1) = fact*nu
        cmat(2,2) = fact
        cmat(3,3) = fact*(1-nu)/2

        ! ng=no. of gauss points in each direction = 2 
        ng=2
        weight(1)=1.0_wp
        weight(2)=1.0_wp

        gconst=1.0_wp/sqrt(3.0_wp)


        ke=0

        
        do i=1,ng
        do j=1,ng

            if(i.eq.1 .and. j.eq.1) then
                xi=-gconst  !! GP1
                eta=-gconst
                GP=1
            end if

            if(i.eq.2 .and. j.eq.1) then
                xi=gconst   !! GP2
                eta=-gconst
                GP=2
            end if

            if(i.eq.2 .and. j.eq.2) then
                xi=gconst   !! GP3
                eta=gconst
                GP=3
            end if

            if(i.eq.1 .and. j.eq.2) then
                xi=-gconst  !! GP4
                eta=gconst
                GP=4
            end if

            call shape(xe,xi,eta,n,bmat,jac,detjac)

            
            
            C_ep=0
           
            if (eplastic_check(GP) .eq. 0 ) then

                ! ELASTIC
 
                ke=ke+ weight(i)*weight(j)*thk*detjac*matmul(matmul(transpose(bmat),cmat),bmat)

            else if (eplastic_check(GP) .eq. 1 ) then

   
                  ! PLASTIC
                call cal_von_stress(estress(:,GP),von_mises_stress_in_GP)
                dF_dsigma(1,1)= (2.0_wp*estress(1,GP) - estress(2,GP))/(2.0_wp*von_mises_stress_in_GP)  
                dF_dsigma(1,2)= (2.0_wp*estress(2,GP) - estress(1,GP))/(2.0_wp*von_mises_stress_in_GP)  
                dF_dsigma(1,3)= (6.0_wp*estress(3,GP))/(2.0_wp*von_mises_stress_in_GP)  

                
                !Et=young/500.0
                Et=young/40.0_wp
                !h_hard=young*n_expo*((eyield_stress(GP)/sgy)**((n_expo-1.0_wp)/n_expo))

                h_hard=Et*young/(young-Et)

                C_ep=cmat-  matmul(matmul(cmat,transpose(dF_dsigma)),matmul(dF_dsigma,cmat))/( &
                dot_product( matmul(dF_dsigma(1,:),cmat),dF_dsigma(1,:) )  + h_hard)

                ! print*,(C_ep)
                ! print*, 'JJJJ GP=', GP
                ! print*,cmat


                           
                ke=ke+ weight(i)*weight(j)*thk*detjac*matmul(matmul(transpose(bmat),C_ep),bmat)

            end if


        end do
        end do


    end subroutine plane42_ke
!
!--------------------------------------------------------------------------------------------------
!
    subroutine plane42_re(xe, eface, fe, thk, re)

        !! This subroutine computes the element load vector due
        !! to surface traction (traction is always perpendicular
        !! to element face).

        integer, intent(in) :: eface
            !! Element face where traction (pressure) is applied

        real(wp), intent(in) :: fe
            !! Value of surface traction (pressure)
        real(wp), intent(in) :: thk
            !! Thickness of this element
        real(wp), dimension(:), intent(in) :: xe
            !! Nodal coordinates of this element in undeformed configuration (see also [[plane42rect_ke]])
        real(wp), intent(out) :: re(8)
        real(wp) :: gconst
            !! Element force vector
            !!
            !! * `re(1:2)` = \((f_x^1, f_y^1)\) force at element node 1 in \(x\)- and y-direction
            !! * `re(3:4)` = \((f_x^2, f_y^2)\) force at element node 1 in \(x\)- and y-direction
            !! * etc...
        real(wp) :: n(2,8),bmat(3,8),jac(2,2),detjac,jac_part(2),pressure
        real(wp),allocatable,dimension (:) :: xi,eta,weight
        integer :: ng,i,j


        ng=2


        if (ng .eq. 1) then
            allocate (xi(1))
            allocate (eta(1))
            allocate (weight(1))

            xi=0
            eta=0
            weight=2

        else if (ng .eq. 2) then
            allocate (xi(2))
            allocate (eta(2))
            allocate (weight(2))

            gconst=1.0_wp/sqrt(3.0_wp)

            xi(1)=-gconst
            xi(2)=gconst
            eta(1)=-gconst
            eta(2)=gconst

            weight(1)=1
            weight(2)=1
        else if (ng .eq. 3) then
            allocate (xi(3))
            allocate (eta(3))
            allocate (weight(3))

            gconst=sqrt(0.6_wp)

            xi(1)=-gconst
            xi(2)=0.0_wp
            xi(3)=gconst
            eta(1)=-gconst
            eta(2)=0.0_wp
            eta(3)=gconst

            weight(1)=5.0_wp/9.0_wp
            weight(2)=5.0_wp/9.0_wp
            weight(3)=8.0_wp/9.0_wp
        end if



        re=0

        ! Sign Convention for Pressure (+ve into the cell)
        if (eface .eq. 1) then 
            pressure=fe
            deallocate(eta)
            allocate(eta(1))          
            eta=-1.0_wp
            do i=1,ng
                
                call shape(xe,xi(i),eta(1),n,bmat,jac,detjac)
                jac_part(1)=-jac(1,2)
                jac_part(2)=jac(1,1)
                re=re+ weight(i)*thk*pressure* matmul(transpose(n),jac_part)

            end do

        else if (eface .eq. 3) then
            pressure=fe
            deallocate(eta)
            allocate(eta(1)) 
            eta=1.0_wp
            do i=1,ng
                
                call shape(xe,xi(i),eta(1),n,bmat,jac,detjac)
                jac_part(1)=jac(1,2)
                jac_part(2)=-jac(1,1)
                re=re+ weight(i)*thk*pressure* matmul(transpose(n),jac_part)

            end do

        else if (eface .eq. 2) then 
            pressure=fe
            deallocate(xi)
            allocate(xi(1)) 
            xi=1.0_wp 
            do i=1,ng
                
                call shape(xe,xi(1),eta(i),n,bmat,jac,detjac)
                jac_part(1)=-jac(2,2)
                jac_part(2)=jac(2,1)
                re=re+ weight(i)*thk*pressure* matmul(transpose(n),jac_part)

            end do


        else if (eface .eq. 4) then
            pressure=fe
            deallocate(xi)
            allocate(xi(1)) 
            xi=-1.0_wp 
            do i=1,ng
                
                call shape(xe,xi(1),eta(i),n,bmat,jac,detjac)
                jac_part(1)=jac(2,2)
                jac_part(2)=-jac(2,1)
                re=re+ weight(i)*thk*pressure* matmul(transpose(n),jac_part)

            end do

        else
          print *, 'ERROR in plane42rect/plane42rect_re'
          print *, 'Element Face incorrect, not 1,2,3 or 4'
          stop
        end if
 
    end subroutine plane42_re
!
!--------------------------------------------------------------------------------------------------

    subroutine plane42_ss(xe, del_de, young, nu, thk, sgy, n_expo, re_int, estress, eyield_stress, & 
                      estrain, eEQP_strain, del_lamda_e, eplastic_check)

        !! This subrotuine computes the element stress and strain (The location inside the element
        !! where stress and and strain is evaluated, is defined inside the subroutine).

        real(wp), intent(in) :: young, sgy, n_expo
            !! Young's Modulus for this element
        real(wp), intent(in) :: nu , thk
            !! Poisson's Ratio for this element
        real(wp), dimension(:), intent(in)  :: xe
            !! Nodal coordinates of this element in undeformed configuration (see also [[plane42rect_ke]])
        real(wp), dimension(:), intent(in)  :: del_de
        real(wp), dimension(:), intent(out)  :: re_int,  del_lamda_e
            !! Displacement field
            !!
            !! * `de(1:2)` = displacement of element node 1 along \(x\)- and \(y\)-axis, respectively
            !! * `de(3:4)` = displacement of element node 2 along \(x\)- and \(y\)-axis, respectively
            !! * etc...
        real(wp), dimension(:,:), intent(inout) :: estress,estrain
        real(wp), dimension(:), intent(inout) ::eyield_stress, eEQP_strain
        integer, dimension(:), intent(out) :: eplastic_check

            !! Stress at a point inside the element
            !!
            !! * `estress(1)` = \(\sigma_{11}\)
            !! * `estress(2)` = \(\sigma_{22}\)
            !! * `estress(3)` = \(\sigma_{12}\)
            !! Strain at a point inside the element
            !!
            !! * `estrain(1)` = \(\epsilon_{11}\)
            !! * `estrain(2)` = \(\epsilon_{22}\)
            !! * `estrain(3)` = \(\epsilon_{12}\)
        real(wp) :: fact,n(2,8), bmat(3, 8), cmat(3, 3),jac(2,2),detjac,xi, eta
        real(wp) ::  pstrain(3),gconst,von_mises_stress_in_GP,yield_toler,h_hard, Et,&
                     dF_dsigma(1,3),del_estrain(3),del_estress(3),del_sgy,weight(2), &
                     von_mises_stress_trial_in_GP, sigma_trial_in_GP(3),del_eplastic_strain(3,1)
        integer :: ng,i,j, GP

        ng=2
        gconst=1.0_wp/sqrt(3.0_wp)
        yield_toler=1.0E-9


            ! Build constitutive matrix (plane stress)
            cmat = 0
            fact = young/(1-nu**2)
            cmat(1,1) = fact
            cmat(1,2) = fact*nu
            cmat(2,1) = fact*nu
            cmat(2,2) = fact
            cmat(3,3) = fact*(1-nu)/2

        weight(1)=1.0_wp
        weight(2)=1.0_wp

        re_int=0

        do i=1,ng
        do j=1,ng

            if(i.eq.1 .and. j.eq.1) then
                xi=-gconst  !! GP1
                eta=-gconst
                GP=1
            end if

            if(i.eq.2 .and. j.eq.1) then
                xi=gconst   !! GP2
                eta=-gconst
                GP=2
            end if

            if(i.eq.2 .and. j.eq.2) then
                xi=gconst   !! GP3
                eta=gconst
                GP=3
            end if

            if(i.eq.1 .and. j.eq.2) then
                xi=-gconst  !! GP4
                eta=gconst
                GP=4
            end if



            ! Build strain-displacement matrix


            call shape(xe,xi,eta,n,bmat,jac,detjac)

            ! Compute element strain
            del_estrain = matmul(bmat, del_de)


            sigma_trial_in_GP=estress(:,GP)+ matmul(cmat,del_estrain)

            call cal_von_stress(sigma_trial_in_GP,von_mises_stress_trial_in_GP)

            if ((von_mises_stress_trial_in_GP - eyield_stress(GP)) .lt. yield_toler) then

                ! ELASTIC
                del_lamda_e(GP)=0.0_wp
                del_estress=matmul(cmat,del_estrain)                
                del_sgy=0.0_wp
                eplastic_check(GP)=0


            else if ((von_mises_stress_trial_in_GP - eyield_stress(GP)) .ge. yield_toler) then

                 ! PLASTIC
                call cal_von_stress(estress(:,GP),von_mises_stress_in_GP)
                dF_dsigma(1,1)= (2.0_wp*estress(1,GP) - estress(2,GP))/(2.0_wp*von_mises_stress_in_GP) 
                dF_dsigma(1,2)= (2.0_wp*estress(2,GP) - estress(1,GP))/(2.0_wp*von_mises_stress_in_GP)  
                dF_dsigma(1,3)= (6.0_wp*estress(3,GP))/(2.0_wp*von_mises_stress_in_GP)  

                !Et=young/500.0
                Et=young/40.0_wp
                !h_hard=young*n_expo*((eyield_stress(GP)/sgy)**((n_expo-1.0_wp)/n_expo))
                h_hard=Et*young/(young-Et)

                del_lamda_e(GP)= dot_product(   matmul(dF_dsigma(1,:),cmat), del_estrain)/( &
                        dot_product( matmul(dF_dsigma(1,:),cmat),dF_dsigma(1,:) )  + h_hard)



                    if (del_lamda_e(GP) .le. yield_toler) then 
                       

                        del_estress=matmul(cmat,del_estrain)
                        del_sgy=0.0_wp
                        del_lamda_e(GP)=0.0_wp
                        eplastic_check(GP)=0
                        

                    else if (del_lamda_e(GP) .gt. yield_toler) then 

                        del_eplastic_strain= del_lamda_e(GP)*transpose(dF_dsigma)

                        del_estress=matmul(cmat,(del_estrain-del_eplastic_strain(:,1)))
                        del_sgy=h_hard*del_lamda_e(GP)
                        eplastic_check(GP)=1
                      !stop

                    end if 

                    
            end if

            ! Compute element stress
            estress(:,GP)=estress(:,GP)+del_estress
            estrain(:,GP)=estrain(:,GP)+del_estrain
            eyield_stress(GP)=eyield_stress(GP)+del_sgy
            eEQP_strain(GP)=eEQP_strain(GP)+del_lamda_e(GP)

            re_int=re_int + weight(i)*weight(j)*thk*detjac*matmul(transpose(bmat),estress(:,GP))

            


        end do
        end do
            ! call cal_von_stress(estress(:,1),von_mises_stress_in_GP)
            ! print*,eyield_stress(1),von_mises_stress_in_GP

!   (plane stress problem)
        !   pstress=0
  !!      pstrain=0
        
  !!      pstrain(1)=(estrain(1)+estrain(2))/2  + sqrt(((estrain(1)-estrain(2))/2)**2 + (estrain(3)/2)**2)
  !!      pstrain(2)=(estrain(1)+estrain(2))/2  - sqrt(((estrain(1)-estrain(2))/2)**2 + (estrain(3)/2)**2)
  !!      pstrain(3)=-nu*(estress(1)+estress(2))/young

        ! pstress(1)=(estress(1)+estress(2))/2  + sqrt(((estress(1)-estress(2))/2)**2 + (estress(3))**2)
        ! pstress(2)=(estress(1)+estress(2))/2  - sqrt(((estress(1)-estress(2))/2)**2 + (estress(3))**2)
        ! sinnn= -2.0_wp *estress(3)/(pstress(1)-pstress(2))
        ! cosss= (estress(1)-estress(2)) /(pstress(1)-pstress(2))
        ! pstress(3)= (atan2(sinnn,cosss)/2.0_wp) * 180.0_wp/pi
      !  print*,pi


        
    end subroutine plane42_ss
!
!
!
    subroutine shape(xe,xi,eta,n,bmat,jac,detjac)

        real(wp), dimension(:), intent(in) :: xe
        real(wp), intent(in) :: xi,eta
        real(wp), dimension(:,:), intent(out) :: n
        real(wp), dimension(:,:), intent(out) :: bmat
        real(wp), dimension(:,:), intent(out) :: jac
        real(wp), intent(out) :: detjac

        real(wp), dimension(2,2) :: gamma
        real(wp), dimension(4,4) :: gamma_tilde
        real(wp), dimension(4,8) :: N_tilde
        real(wp), dimension(3,4) :: L
        real(wp), dimension(4) :: dN_dxi,dN_deta,xe_X,xe_Y,nn
        integer :: i 
        
        L=0

        L(1,1)=1
        L(2,4)=1
        L(3,2)=1
        L(3,3)=1

        nn(1)=(1.0_wp-xi)*(1.0_wp-eta)/4.0_wp
        nn(2)=(1.0_wp+xi)*(1.0_wp-eta)/4.0_wp
        nn(3)=(1.0_wp+xi)*(1.0_wp+eta)/4.0_wp
        nn(4)=(1.0_wp-xi)*(1.0_wp+eta)/4.0_wp

        dN_dxi(1)=-(1.0_wp-eta)/4.0_wp
        dN_dxi(2)=(1.0_wp-eta)/4.0_wp
        dN_dxi(3)=(1.0_wp+eta)/4.0_wp
        dN_dxi(4)=-(1.0_wp+eta)/4.0_wp

        dN_deta(1)=-(1.0_wp-xi)/4.0_wp
        dN_deta(2)=-(1.0_wp+xi)/4.0_wp
        dN_deta(3)=(1.0_wp+xi)/4.0_wp
        dN_deta(4)=(1.0_wp-xi)/4.0_wp

        do i=1,4
            xe_X(i)=xe(2*i-1)
            xe_Y(i)=xe(2*i)
        end do

        jac(1,1)=dot_product(dN_dxi,xe_X)
        jac(1,2)=dot_product(dN_dxi,xe_Y)
        jac(2,1)=dot_product(dN_deta,xe_X)
        jac(2,2)=dot_product(dN_deta,xe_Y)
        
        detjac=jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1)
        
        gamma(1,1)=jac(2,2)/detjac
        gamma(1,2)=-jac(1,2)/detjac
        gamma(2,1)=-jac(2,1)/detjac
        gamma(2,2)=jac(1,1)/detjac

        gamma_tilde=0
        gamma_tilde(1:2,1:2)=gamma
        gamma_tilde(3:4,3:4)=gamma

        N_tilde=0
        n=0
        do i=1,4

            N_tilde(1,2*i-1)=dN_dxi(i)
            N_tilde(2,2*i-1)=dN_deta(i)

            N_tilde(3,2*i)=dN_dxi(i)
            N_tilde(4,2*i)=dN_deta(i)            

            n(1,2*i-1)=nn(i)
            n(2,2*i)=nn(i)

        end do

        bmat=matmul(matmul(L,gamma_tilde),N_tilde)
    end subroutine shape

    ! subroutine hardening_law(sig_y0,epsilon_0,n_expo,eqps,sig_y)

    !     real(wp), intent(in) :: sig_y0,epsilon_0,n_expo,eqps
    !     real(wp), intent(out) :: sig_y

    !     sig_y=sig_y0*(1.0_wp +  eqps/epsilon_0)**n_expo
    ! end subroutine hardening_law

    subroutine cal_von_stress(stress_in_GP,von_mises_stress_in_GP)

        real(wp), intent(in) :: stress_in_GP(3)
        real(wp), intent(out) :: von_mises_stress_in_GP

        von_mises_stress_in_GP=sqrt(stress_in_GP(1)**2 + stress_in_GP(2)**2 + &
                    3.0_wp*stress_in_GP(3)**2 - stress_in_GP(1)*stress_in_GP(2)) 
    end subroutine cal_von_stress

end module plane42
