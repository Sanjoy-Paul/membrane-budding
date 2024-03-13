      Module module_Initialize_system
      Integer,Allocatable,Dimension(:,:) :: Vneighlist
      Contains
!---------------------------------------------------------------------------------------------------------------------------
!             Subroutine to make a patch containing epsins
!---------------------------------------------------------------------------------------------------------------------------
        Subroutine Initialize_PatchOfProteins(num_proteins)                                                                                     ! Makes a patch with scalar field 1 and phase field 1 
        Use module_datastruct
        Implicit None
        Integer :: i,j,ncount,ncount1,verlist(nver),rnumber,num_proteins
        Integer :: jj1
        Real(Kind=8) :: ran2,splo(3,1)
        Real,Dimension(nver) :: sec23list
        Integer:: i1,j_old,jj,sjj
        Integer ::ier,num_rings,nbuffer,buffsize,hexnline,sec23num,sec23
        Real(Kind=8),Dimension(:),Allocatable ::pern,Kparstr,Kperstr,cparstr,cperstr
        Real(Kind=8):: bud_area,kappa_unscaled,hexbl
        Real(Kind=8), Allocatable ::surfaceprop(:,:),energyprop(:,:)
        Character(100) ::Cnname,fname,start_state,surffile,energyfile,syscmd,nem_profile


        ver(:)%splo(1,1)=0.0   ; ver(:)%splo(2,1)=0.0 
        ver(:)%splo(3,1)=0.0   ; ver(:)%phase=2 ;   verlist=0  
        ver(:)%cpar=0.0        ; ver(:)%cper=0.0 
        rnumber=Nint(nver*ran2(seed))+1                                                                                                 ! The patch is centered around a rnumberomly chosen vertex rnumber 




       fname='./parameters_circular.in'
        OPEN(01,File=Trim(Adjustl(fname)))
! Read the parameters from a  file 
        READ(01,*),kappa_unscaled,pr
        Allocate(pern(Nphase))
        READ(01,*),pern(1:Nphase)
! Percentage of nematic field of kind A and B  
        Allocate(Kparstr(Nphase))
        READ(01,*),Kparstr(1:Nphase)
! Kparallel flag for each phase 
        Allocate(Kperstr(Nphase))
        READ(01,*),Kperstr(1:Nphase)
! Kperp flag for each phase 
        Allocate(cparstr(Nphase))
        READ(01,*),cparstr(1:Nphase)
! c0 parallel for each phase
        Allocate(cperstr(Nphase))
        READ(01,*),cperstr(1:Nphase)
! c0 perp for each phase
!        Allocate(nempar(Nphase,Nphase))
!        Allocate(unempar(Nphase,Nphase))
!        Do i=1,Nphase,1
!        READ(01,*),nempar(i:i,:)
!        Enddo
!        Read(01,*),start_state
! starting state 'MEMBRANE' or 'RESTART' or 'HEXAHEDRON'
!        Read(01,*),hexnline,hexbl
! number of discretization and bondlength                                
!                Read(01,*),itime,ftime,dataint,confint
!                Read(01,*),nem_profile,bud_area,num_rings       
        Close(01)



! For s % spacer the following line will be -
!       sec23num=Nint(num_proteins*(s/100))
        sec23num=Nint(num_proteins*0.15) 
        sec23=1
 
        If(rnumber.Gt.nver) rnumber=nver
        ncount=1 ; ncount1=1 ; verlist(ncount1)=rnumber ;
        j_old = rnumber              
        ver(rnumber)%phase=1  
         Do While((ncount .Lt. num_proteins) .And. (ncount1.Lt.nver))
         Do i=1,ver(verlist(ncount1))%nonei
         j=ver(verlist(ncount1))%vneipt(i)
         If(ver(j)%phase .Eq. 2)Then
         ncount=ncount+1 ; ver(j)%phase=1 ; verlist(ncount)=j
         splo = reshape((/1-2*ran2(Seed),1-2*ran2(Seed),zero/),(/3,1/))
         ver(j)%splo = splo/Sqrt(Sum(splo**2))
         ver(j)%cpar=Cparstr(1) ; j_old = j 
         Endif
         Enddo
         ncount1=ncount1+1
         Enddo



        j_old=verlist(ncount); sec23list(1) = j_old
        Do While(sec23 .Lt. sec23num) 
!        ncount=ncount+1 
!        If(ncount .Lt. num_proteins)Then
        
        rnumber=Nint(ran2(Seed)*sec23num)+1
        
        if(ver(verlist(rnumber))%phase .eq. 1) Then        
        ver(verlist(rnumber))%splo = reshape((/0,0,0/),(/3,1/))
        ver(verlist(rnumber))%phase=3
        ver(verlist(rnumber))%cpar=0.0
        Print*,"Here is a Sec23 containing vertex and number of spacer",verlist(ncount), sec23
        sec23list(sec23) = verlist(ncount)
        sec23=sec23+1
        Endif

!        Endif
        Enddo


         Print*,'rnumber,sec23',num_proteins,sec23
        End Subroutine Initialize_PatchOfProteins

!---------------------------------------------------------------------------------------------------------------------------
!                 Subroutine to distribute the epsins at rnumberom locations
!---------------------------------------------------------------------------------------------------------------------------
        Subroutine Initialize_Proteins_Randomly(num_proteins)                                                                                   ! Makes a patch with scalar field 1 and phase field 1 
        Use module_datastruct
        Implicit None
        Integer :: ncount,rnumber,num_proteins
        Real(Kind=8) :: ran2,splo(3,1)
        ver(:)%splo(1,1)=0.0   ; ver(:)%splo(2,1)=0.0 
        ver(:)%splo(3,1)=0.0
        ver(:)%phase=2 ; ncount=0
        Do While((ncount.Lt. num_proteins))
          rnumber=Nint(ran2(Seed)*nver)+1
          If(rnumber.Gt.nver) rnumber=nver
          If(ver(rnumber)%phase .Eq. 2)Then
          ver(rnumber)%phase = 1
          ncount=ncount+1
          splo=reshape((/1.-2.*ran2(Seed),1.-2.*ran2(Seed),zero/),(/3,1/))
          ver(rnumber)%splo=splo/Sqrt(Sum(splo**2))
          Endif
        Enddo
        End Subroutine Initialize_Proteins_Randomly
!---------------------------------------------------------------------------------------------------------------------------
!                              Initialize  vector field in an annulus
!---------------------------------------------------------------------------------------------------------------------------
        Subroutine Initialize_Protein_In_Annulus(Loc_Nei,vicinity)                                                                  ! Makes a patch with scalar field 1 and phase field 1 
        Use module_datastruct ; Use module_dataformat
        Implicit None
        Integer :: i,ncount,rnumber,Loc_nei
        Real(Kind=8) :: ran2,veclo(3,1)
        Integer :: vicinity

        rnumber=Nint(nver*ran2(seed))+1                                                                                                ! The patch is centered around a rnumberomly chosen vertex rnumber 
        If(rnumber.Gt.nver) rnumber=nver
        Call Complete_neighbourhood(rnumber,'RETURN_LIST')
        ver(:)%phase=2 ;  ncount=0
        ver(:)%splo(1,1)=0.0   ; ver(:)%splo(2,1)=0.0 
        ver(:)%splo(3,1)=0.0
        Do i=1,nver,1
        If(Vneighlist(i,2).Ge.Loc_nei .And. Vneighlist(i,2).Le.Loc_nei+vicinity )Then
        ver(vneighlist(i,1))%phase=1 
        veclo=reshape((/1.-2*ran2(Seed),1.-2*ran2(seed),zero/),(/3,1/))
        ver(Vneighlist(i,1))%splo=veclo/Sqrt(Sum(veclo**2))
        Endif
        Enddo
        DeAllocate(Vneighlist)
        End Subroutine Initialize_Protein_In_Annulus

!---------------------------------------------------------------------------------------------------------------------------
!                       Initialize  vector field in an annulus with num_proteins  nematics and area patch_area
!---------------------------------------------------------------------------------------------------------------------------
        Subroutine Initialize_Annulus_Area_and_Number(patch_area,num_proteins)                                                      ! Makes a patch with scalar field 1 and phase field 1 
        Use module_datastruct ; Use module_dataformat
        Implicit None
        Integer :: i,j,ncount,rnumber,num_proteins
        Real(Kind=8) :: ran2,veclo(3,1),patch_area,area_comp
        Integer :: min_ring,max_ring,area_flag,num_count
     
	area_comp=0.0 ; area_flag=0 ; num_count=0

        rnumber=Nint(nver*ran2(seed))+1                                                                                                ! The patch is centered around a rnumberomly chosen vertex rnumber 
        If(rnumber.Gt.nver) rnumber=nver
        Call Complete_neighbourhood(rnumber,'RETURN_LIST')
        ver(:)%phase=2 ;  ncount=0
        ver(:)%splo(1,1)=0.0   ; ver(:)%splo(2,1)=0.0 ;  ver(:)%splo(3,1)=0.0

	min_ring=Minval(Vneighlist(:,2)) ; max_ring=Maxval(Vneighlist(:,2))

	j=min_ring

	Do While(num_count .Lt. num_proteins)
		Do i=1,nver,1
			If(Vneighlist(i,2) .Eq. j) Then
				If(area_flag .Eq.0) Then
				area_comp=area_comp+ver(i)%totarea
				if (area_comp .Ge. patch_area) area_flag=1
				Else
			        ver(vneighlist(i,1))%phase=1 
			        veclo=reshape((/1.-2*ran2(Seed),1.-2*ran2(seed),zero/),(/3,1/))
			        ver(Vneighlist(i,1))%splo=veclo/Sqrt(Sum(veclo**2))
				num_count=num_count+1
				If(num_count .Eq. num_proteins) Return
				Endif
			Endif
		Enddo
		j=j+1
	Enddo

        DeAllocate(Vneighlist)
        End Subroutine Initialize_Annulus_Area_and_Number

!---------------------------------------------------------------------------------------------------------------------------
!                       Initialize  vector field in an annulus with num_proteins  nematics and area patch_area
!---------------------------------------------------------------------------------------------------------------------------
        Subroutine Initialize_Annulus_Area_and_rings(patch_area,num_rings)                                                      ! Makes a patch with scalar field 1 and phase field 1 
        Use module_datastruct ; Use module_dataformat
        Implicit None
        Integer :: i,j,ncount,rnumber,num_rings,num_proteins
        Real(Kind=8) :: ran2,veclo(3,1),patch_area,area_comp
        Integer :: min_ring,max_ring,area_flag,num_count
     
	area_comp=0.0 ; area_flag=0 ; num_count=0 ; num_proteins = 0

        rnumber=Nint(nver*ran2(seed))+1                                                                                                ! The patch is centered around a rnumberomly chosen vertex rnumber 
        If(rnumber.Gt.nver) rnumber=nver
        Call Complete_neighbourhood(rnumber,'RETURN_LIST')
        ver(:)%phase=2 ;  ncount=0
        ver(:)%splo(1,1)=0.0   ; ver(:)%splo(2,1)=0.0 ;  ver(:)%splo(3,1)=0.0

	min_ring=Minval(Vneighlist(:,2)) ; max_ring=Maxval(Vneighlist(:,2))

	j=min_ring

	Do While(num_count .Le. num_rings)
		Do i=1,nver,1
			If(Vneighlist(i,2) .Eq. j) Then
				If(area_flag .Eq.0) Then
				area_comp=area_comp+ver(i)%totarea
				if (area_comp .Ge. patch_area) Then
				       	area_flag=1 ; exit 
				Endif
				Else
			        ver(vneighlist(i,1))%phase=1 
			        veclo=reshape((/1.-2*ran2(Seed),1.-2*ran2(seed),zero/),(/3,1/))
			        ver(Vneighlist(i,1))%splo=veclo/Sqrt(Sum(veclo**2))
				num_proteins = num_proteins+1
				Endif
			Endif
		Enddo
		If(area_flag.Eq.1) num_count=num_count+1
	       	j=j+1
	Enddo
	mp%nadd = num_proteins
	mp%nmem = nver -mp%nadd

        DeAllocate(Vneighlist)
        End Subroutine Initialize_Annulus_Area_and_rings

!---------------------------------------------------------------------------------------------------------------------------
!                       Initialize  vector field in an annulus with num_proteins  nematics and area patch_area
!---------------------------------------------------------------------------------------------------------------------------
        Subroutine Initialize_Patch_Specified_Area(patch_area)                                                                      ! Makes a patch with scalar field 1 and phase field 1 
        Use module_datastruct ; Use module_dataformat
        Implicit None
        Integer :: i,j,ncount,rnumber,num_proteins
        Real(Kind=8) :: ran2,veclo(3,1),patch_area,area_comp
        Integer :: min_ring,max_ring,area_flag,num_count
     
		area_comp=0.0 ; area_flag=0 ; num_count=0 ; num_proteins=0

        rnumber=Nint(nver*ran2(seed))+1                                                                                                ! The patch is centered around a rnumberomly chosen vertex rnumber 
        If(rnumber.Gt.nver) rnumber=nver
        Call Complete_neighbourhood(rnumber,'RETURN_LIST')
        ver(:)%phase=2 ;  ncount=0
        ver(:)%splo(1,1)=0.0   ; ver(:)%splo(2,1)=0.0 ;  ver(:)%splo(3,1)=0.0

	min_ring=Minval(Vneighlist(:,2)) ; max_ring=Maxval(Vneighlist(:,2))

	j=min_ring

	Do While(area_flag .Eq. 0)
		Do i=1,nver,1
			If(Vneighlist(i,2) .Eq. j) Then
				If(area_flag .Eq.0) Then
				area_comp=area_comp+ver(i)%totarea
			        ver(vneighlist(i,1))%phase=1 
			        veclo=reshape((/1.-2*ran2(Seed),1.-2*ran2(seed),zero/),(/3,1/))
			        ver(Vneighlist(i,1))%splo=veclo/Sqrt(Sum(veclo**2))
				num_proteins = num_proteins+1
				If (area_comp .Ge. patch_area) Then
				area_flag=1 ; exit 
				Endif
				Endif
			Endif
		Enddo
		If(area_flag.Eq.1) num_count=num_count+1
	       	j=j+1
	Enddo

	mp%nadd = num_proteins
	mp%nmem = nver -mp%nadd
        DeAllocate(Vneighlist)
        End Subroutine Initialize_Patch_Specified_Area
	
	
	

!---------------------------------------------------------------------------------------------------------------------------
!                             Compute the curvature of an Epsin in its neighbourhood after being moved
!---------------------------------------------------------------------------------------------------------------------------
        Subroutine Complete_neighbourhood(v,retflag)                    
        Use module_datastruct
        Implicit None
        Integer :: v
        Character(Len=9),Allocatable :: accounted(:)
        Character(Len=9), Intent(In), Optional :: retflag
        Integer :: i,j,ncount,ncount1
        Integer,Allocatable :: neilist(:,:)
        ver(:)%neigh=0

        Allocate(neilist(nver,2)) ; Allocate(accounted(nver))
        accounted='NOTACC' ; neilist=0
        neilist(1,1)=v ; neilist(1,2)=0                                 							    ! The first entry is the vertex itself  
        ncount=ver(v)%nonei+1
        neilist(2:ncount,1)=ver(v)%vneipt(1:ver(v)%nonei)                							    ! Make a list with vertex v as a starting point             
        neilist(2:ncount,2)=1                                                                                                       ! Set the number of the neighbourhood we are accounting for 
        Forall (i=1:ncount) accounted(neilist(i,1))='ACCOUNTED'                                                                     ! this maintains a list that avoids overcounting            
        Forall (i=2:ncount) ver(neilist(i,1))%neigh=neilist(i,2)                                                                    ! this maintains a list that avoids overcounting            

        ncount1=2
        Do While(ncount.Lt.nver)
           Do i=1,ver(neilist(ncount1,1))%nonei,1
              j=ver(neilist(ncount1,1))%vneipt(i)
              If(accounted(j).Eq.'NOTACC')Then
              ncount=ncount+1 
              neilist(ncount,1)=j
              neilist(ncount,2)=neilist(ncount1,2)+1
              ver(j)%neigh=neilist(ncount,2)
              accounted(j)='ACCOUNTED'
              Endif
           Enddo
           ncount1=ncount1+1
        Enddo

        If(Present(retflag)) Then                                                                                                   ! Return the list only if asked for 
        Allocate(Vneighlist(nver,2))        
        Vneighlist=neilist
        Endif

        DeAllocate(neilist) ; DeAllocate(accounted)
        End Subroutine Complete_neighbourhood      

      End Module module_Initialize_system

!---------------------------------------------------------------------------------------------------------------------------
!                             Initialize peptide arrangemnet along a
!                             line (series arrangement)
!
!---------------------------------------------------------------------------------------------------------------------------



