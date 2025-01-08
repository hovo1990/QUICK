!---------------------------------------------------------------------!
! Created by Hovakim Grabski on 12/08/2024                            !
!                                                                     !
! Copyright (C) 2020-2021 Merz lab                                    !
! Copyright (C) 2020-2021 Götz lab                                    !
!                                                                     !
! This Source Code Form is subject to the terms of the Mozilla Public !
! License, v. 2.0. If a copy of the MPL was not distributed with this !
! file, You can obtain one at http://mozilla.org/MPL/2.0/.            !
!_____________________________________________________________________!

#include "util.fh"


module quick_qcschema_module

    ! -- TODO use json_module
    use,intrinsic :: iso_fortran_env, only: wp => real64
    use json_module


    implicit none
    private

    public :: quick_qcschema
    public :: initializeExportQC, finalizeExportQC, exportCoordinatesQC, exportBasisQC, exportMOQC, &
              exportSCFQC, exportOPTQC

    type quick_qcschema_type
      integer :: iQCSchemaFile

      ! true if this a geometry optimization
      logical :: opt
      
      ! number of atoms in molecule
      integer :: natom

      ! atom symbols
      character(len=2), allocatable :: atom_symbol(:)

      ! number of scf iterations to converge
      integer, dimension(:), allocatable :: nscf_snapshots

      ! scf energy during each iteration
      double precision,dimension(:,:),allocatable :: e_snapshots

      ! geometry during optimization
      double precision,dimension(:,:,:),allocatable :: xyz_snapshots

      ! counter to keep track of number of snapshots
      integer :: iexport_snapshot
    
      ! temporary vector for reorganizing mo coefficients
      double precision, dimension(:), allocatable :: reord_mo_vec

      ! -- TODO initiliazing json objects
      type(json_core) :: json
      type(json_value), pointer :: p

       ! counter to keep track of number of snapshots
      integer :: testy

  

    end type quick_qcschema_type

    type (quick_qcschema_type),save:: quick_qcschema

    interface initializeExportQC
        module procedure initialize_qcschema
    end interface initializeExportQC

    interface finalizeExportQC
        module procedure finalize_qcschema
    end interface finalizeExportQC

    interface exportCoordinatesQC
        module procedure write_coordinates
    end interface exportCoordinatesQC

    interface exportBasisQC
        module procedure write_basis_info
    end interface exportBasisQC

    interface exportMOQC
        module procedure write_mo
    end interface exportMOQC

    interface exportSCFQC
        module procedure write_scf
    end interface exportSCFQC

    interface exportOPTQC
        module procedure write_opt
    end interface exportOPTQC

contains

subroutine write_coordinates(self, ierr)

    use quick_molspec_module, only: quick_molspec, xyz, natom
    use quick_constants_module, only : symbol, BOHRS_TO_A
    use json_module
    implicit none
    type (quick_qcschema_type), intent(inout) :: self
    type(json_value),pointer :: inp
    type(json_value),pointer :: geometry
    integer, intent(out) :: ierr
    integer :: i, j, k
    integer :: testC
    integer :: nx, ny, nz
    real, allocatable :: array_slice(:,:)
    real, allocatable :: one_dim_reshaped(:)
    integer :: total_elements

    ! -- TODO simple example
    ! integer, dimension(3, 3) :: two_dim_array = reshape((/ 1, 2, 3, 4, 5, 6, 7, 8, 9 /), [3, 3])
    ! integer, dimension(9) :: one_dim_array

    
    ! ! Calculate the total number of elements in the 2D array
    ! total_elements = size(two_dim_array)

    ! ! Convert 2D to 1D using reshape
    ! one_dim_array = reshape(two_dim_array, [total_elements])
    
    ! -- TODO working code, for debug
    ! call self%json%create_object(inp,'writeCoordinates')
    ! call self%json%add(self%p, inp) !add it to the root


    call self%json%create_object(inp,'molecule')
    call self%json%add(self%p, inp) !add it to the root


    ! ! ! -- * add some data to inputs:
    call self%json%add(inp, 'symbols', self%atom_symbol)
    ! call self%json%add(inp, 'tutusf', 1.1_wp)


    ! nullify(inp)  !don't need this anymore

    ! -- TODO write atomic coordinates
    if (self%opt) then
        ! in case of geometry optimization write last stored geometry
        ! we need to do this because optimizers may return the geometry
        ! for the next step which may be stored in xyz
        k = self%iexport_snapshot - 1
        print *, 'Stop Here'
        array_slice = self%xyz_snapshots(:, :, k)
        total_elements = size(array_slice)


        allocate ( one_dim_reshaped (total_elements) )  
        one_dim_reshaped = reshape(array_slice, [total_elements])

        ! Loop through the array and assign each value with 10 digits after the decimal point
        ! do i = 1, size(one_dim_reshaped )
        !     one_dim_reshaped(i) = nint(one_dim_reshaped(i) * 1.0d10) / 1.0d10
        !  end do

        call self%json%add(inp, 'geometry', one_dim_reshaped )

        deallocate(one_dim_reshaped )
      else
        ! if it's a single point calculation we can use xyz
        ! we can't use xyz_snapshots because they have not been populated
        write(self%iQCSchemaFile, '(3(2x,F20.10))') (xyz(j,i),j=1,3)
     endif
    testC = 1

    ! -- * write atomic labels and coordinates
    ! write(self%iQCSchemaFile, '("[Atoms] (AU)")')
    ! do i=1,natom
    !     if (self%opt) then
    !        ! in case of geometry optimization write last stored geometry
    !        ! we need to do this because optimizers may return the geometry
    !        ! for the next step which may be stored in xyz
    !        k = self%iexport_snapshot - 1
    !     !    write(self%iQCSchemaFile, '(3(2x,F20.10))') (self%xyz_snapshots(j,i,k),j=1,3)
    !        write(*, '(3(2x,F20.10))') (self%xyz_snapshots(j, i, k), j = 1, 3)
    !      else
    !        ! if it's a single point calculation we can use xyz
    !        ! we can't use xyz_snapshots because they have not been populated
    !        write(self%iQCSchemaFile, '(3(2x,F20.10))') (xyz(j,i),j=1,3)
    !     endif
    ! enddo

end subroutine write_coordinates

subroutine write_basis_info(self, ierr)

    use quick_basis_module, only: quick_basis, nshell, nbasis, ncontract
    use quick_molspec_module, only: natom
    use json_module

    implicit none
    type (quick_qcschema_type), intent(inout) :: self
    integer, intent(out) :: ierr
    integer :: iatom, ishell, ibas, iprim, nprim, j, ishell_idx
    logical :: print_gto
    double precision :: val_gccoeff, xnorm

    ! write basis function information
    write(self%iQCSchemaFile, '("[GTO] (AU)")')
    do iatom=1, natom
        write(self%iQCSchemaFile, '(2x, I5)') iatom

        ! s basis functions
        do ishell=1, nshell
            if(quick_basis%katom(ishell) .eq. iatom) then
                nprim = quick_basis%kprim(ishell)
                if(quick_basis%ktype(ishell) .eq. 1) then
                    write(self%iQCSchemaFile, '(2x, "s", 4x, I2)') nprim
                    do iprim=1, nprim
                        ishell_idx=quick_basis%ksumtype(ishell)
                        write(self%iQCSchemaFile, '(2E20.10)') &
                        quick_basis%gcexpo(iprim,ishell_idx), quick_basis%unnorm_gccoeff(iprim,ishell_idx) 
                    enddo                    
                endif
            endif
        enddo

        ! s, p basis functions of sp shell
        do ishell=1, nshell
            if(quick_basis%katom(ishell) .eq. iatom) then
                nprim = quick_basis%kprim(ishell)
                if(quick_basis%ktype(ishell) .eq. 4) then
                    write(self%iQCSchemaFile, '(2x, "s", 4x, I2)') nprim
                    do iprim=1, nprim
                        ishell_idx=quick_basis%ksumtype(ishell)
                        write(self%iQCSchemaFile, '(2E20.10)') &
                        quick_basis%gcexpo(iprim,ishell_idx), quick_basis%unnorm_gccoeff(iprim,ishell_idx)
                    enddo
                    write(self%iQCSchemaFile, '(2x, "p", 4x, I2)') nprim
                    do iprim=1, nprim
                        ishell_idx=quick_basis%ksumtype(ishell)
                        write(self%iQCSchemaFile, '(2E20.10)') &
                        quick_basis%gcexpo(iprim,ishell_idx), (quick_basis%unnorm_gccoeff(iprim,ishell_idx+1))
                    enddo
                endif
            endif
        enddo
        
        ! p, d, and f basis functions
        do ishell=1, nshell
            if(quick_basis%katom(ishell) .eq. iatom) then
                nprim = quick_basis%kprim(ishell)
                print_gto=.false.
                if(quick_basis%ktype(ishell) .eq. 3) then
                    print_gto=.true.
                    write(self%iQCSchemaFile, '(2x, "p", 4x, I2)') nprim
                elseif(quick_basis%ktype(ishell) .eq. 6) then
                    print_gto=.true.
                    write(self%iQCSchemaFile, '(2x, "d", 4x, I2)') nprim
                elseif(quick_basis%ktype(ishell) .eq. 10) then
                    print_gto=.true.
                    write(self%iQCSchemaFile, '(2x, "f", 4x, I2)') nprim
                endif
                
                do iprim=1, nprim
                    if(print_gto) then
                        ishell_idx=quick_basis%ksumtype(ishell)
                        write(self%iQCSchemaFile, '(2E20.10)') &
                        quick_basis%gcexpo(iprim,ishell_idx), quick_basis%unnorm_gccoeff(iprim,ishell_idx)
                    endif
                enddo
            endif
        enddo

        write(self%iQCSchemaFile, '("")')
    enddo

end subroutine write_basis_info

subroutine write_mo(self, ierr)

    use quick_basis_module, only: quick_basis, nbasis
    use quick_calculated_module, only: quick_qm_struct
    use quick_scratch_module
    use quick_molspec_module, only: quick_molspec
    use quick_method_module, only: quick_method

    implicit none
    type (quick_qcschema_type), intent(inout) :: self

    ! type(json_core), intent(inout) :: jsontodo
    ! type(json_value), pointer, intent(inout) :: ptodo
    ! type(json_value), pointer, intent(inout) :: inptodo


    integer, intent(out) :: ierr    
    integer :: i, j, k, neleca, nelecb
    character(len=5) :: lbl1
    double precision :: occnum, occval

    if(.not. allocated(self%reord_mo_vec)) allocate(self%reord_mo_vec(nbasis))

    write(self%iQCSchemaFile, '("[MO]")')

    if(.not.quick_method%unrst) then
        neleca = quick_molspec%nElec/2
        occval = 2.0d0
    else
        neleca = quick_molspec%nElec
        nelecb = quick_molspec%nElecb
        occval = 1.0d0
    endif

    do i=1, nbasis
        if(neleca .gt. 0 ) then
            occnum=occval
            neleca=neleca-1
        else
            occnum=0.0d0
        endif

        write(lbl1,'(I5)') i
        write(self%iQCSchemaFile, '(A11)') "  Sym= a"//trim(adjustl(lbl1))
        write(self%iQCSchemaFile, '("  Ene= ", E20.10)') quick_qm_struct%E(i)
        write(self%iQCSchemaFile, '(2x, "Spin= Alpha" )') 

        ! write orbital occupation numbers
        write(self%iQCSchemaFile, '("  Occup= ", F10.8)') occnum 

        ! reorder molecular orbital coefficients        
         call reorder_mo_coeffs(quick_qm_struct%co, quick_basis%KLMN, nbasis, i, self%reord_mo_vec, ierr)

        ! write molecular orbital coefficients  
        do j=1, nbasis
            write(self%iQCSchemaFile, '(I4,F15.10)') j, self%reord_mo_vec(j)
        enddo
    enddo

    if(quick_method%unrst) then
        do i=1, nbasis
            if(nelecb .gt. 0 ) then
                occnum=occval
                nelecb=nelecb-1
            else
                occnum=0.0d0
            endif
    
            write(lbl1,'(I5)') i
            write(self%iQCSchemaFile, '(A11)') "  Sym= b"//trim(adjustl(lbl1))
            write(self%iQCSchemaFile, '("  Ene= ",E20.10)') quick_qm_struct%Eb(i)
            write(self%iQCSchemaFile, '(2x, "Spin= Beta" )')
    
            ! write orbital occupation numbers
            write(self%iQCSchemaFile, '("  Occup= ", F10.8)') occnum

            ! reorder molecular orbital coefficients        
            call reorder_mo_coeffs(quick_qm_struct%cob, quick_basis%KLMN, nbasis, i, self%reord_mo_vec, ierr)
    
            ! write molecular orbital coefficients        
            do j=1, nbasis
                write(self%iQCSchemaFile, '(I4,F15.10)') j, self%reord_mo_vec(j)
            enddo
        enddo
    endif

    if(allocated(self%reord_mo_vec)) deallocate(self%reord_mo_vec)

end subroutine write_mo

subroutine write_scf(self, ierr)

    implicit none
    type (quick_qcschema_type), intent(inout) :: self
    integer, intent(out) :: ierr
    integer :: i, j
    character(len=9) :: label

    write(self%iQCSchemaFile, '("[SCFCONV]")')

    do i=1, self%iexport_snapshot-1
       if (i == 1) then
          label = "scf-first"
       else if (i == (self%iexport_snapshot - 1) ) then
          label = "scf-last"
       else
          label = "scf"
       end if
       write(self%iQCSchemaFile, '(2x, a, 2x, I3, 2x, "THROUGH", 2x, I3)') trim(label), 1, self%nscf_snapshots(i)
        do j=1, self%nscf_snapshots(i)
            write(self%iQCSchemaFile, '(2x, E16.10)') self%e_snapshots(j,i)
        enddo
    enddo

end subroutine write_scf

subroutine write_opt(self, ierr)

    use quick_constants_module, only : BOHRS_TO_A
    implicit none
    type (quick_qcschema_type), intent(inout) :: self
    integer, intent(out) :: ierr
    integer :: i, j, k
    character(len=8) :: lbl1
    character(len=2) :: lbl2

    write(self%iQCSchemaFile, '("[GEOCONV]")')

    write(self%iQCSchemaFile, '("energy")')
    do i=1, self%iexport_snapshot-1
        write(self%iQCSchemaFile, '(2x, E16.10)') self%e_snapshots(self%nscf_snapshots(i),i)
    enddo

    write(self%iQCSchemaFile, '("[GEOMETRIES] (XYZ)")')

    do k=1, self%iexport_snapshot-1
       write(self%iQCSchemaFile, '(2x, I5)') self%natom
       write(self%iQCSchemaFile, '("")')
        do i=1,self%natom
           write(self%iQCSchemaFile,'(A2, 2x, 3F14.6)') &
                self%atom_symbol(i), (self%xyz_snapshots(j,i,k)*BOHRS_TO_A,j=1,3)
        enddo    
    enddo

end subroutine write_opt

subroutine initialize_qcschema(self, ierr)
    
    use quick_files_module, only : iQCSchemaFile, qcSchemaFileName
    use quick_method_module, only: quick_method
    use quick_molspec_module, only: natom, quick_molspec
    use quick_constants_module, only : symbol




    implicit none
    type (quick_qcschema_type), intent(inout) :: self

    integer, intent(out) :: ierr
    type(json_value),pointer :: inp

    integer :: i, dimy

    self%iQCSchemaFile = iQCSchemaFile
    self%iexport_snapshot=1
    self%natom = natom



    ! ! -- TODO initialize the json object
    ! ! -- * initialize the class
    call self%json%initialize()

    ! -- * initialize the structure:
    call self%json%create_object(self%p,'')

!    self%json = json
!    self%p = p


    dimy = 1
    if (quick_method%opt) then
       self%opt = .true.
       dimy = quick_method%iopt
    else
       self%opt = .false.
    end if

    ! -- * allocate memory
    if(.not. allocated(self%atom_symbol)) allocate(self%atom_symbol(natom))
    if(.not. allocated(self%nscf_snapshots)) allocate(self%nscf_snapshots(quick_method%iscf))
    if(.not. allocated(self%e_snapshots)) allocate(self%e_snapshots(quick_method%iscf, dimy))
    if(.not. allocated(self%xyz_snapshots)) allocate(self%xyz_snapshots(3, natom, dimy))

    ! -- * store atom symbols
    do i = 1, natom
       self%atom_symbol(i) = symbol(quick_molspec%iattype(i))
    end do

    self%testy = 1
    ! -- * open file
    ! call quick_open(self%iQCSchemaFile,qcSchemaFileName,'U','F','R',.false.,ierr)

    ! write(self%iQCSchemaFile, '("[QCSchema Format]")')


    ! -- TODO this is working code
    ! ! -- * add an "inputs" object to the structure:
    ! call self%json%create_object(inp,'inputs')
    ! call self%json%add(self%p, inp) !add it to the root

    ! ! ! -- * add some data to inputs:
    ! call self%json%add(inp, 't0', 0.1_wp)
    ! call self%json%add(inp, 'tf', 1.1_wp)

    ! nullify(inp)  !don't need this anymore

    ! ! -- * write the file:
    ! call self%json%print(self%p, qcSchemaFileName)

    ! ! -- * cleanup:
    ! call self%json%destroy(self%p)
    ! if (self%json%failed()) stop 1


end subroutine initialize_qcschema

subroutine finalize_qcschema(self, ierr)

    use quick_files_module, only : iQCSchemaFile, qcSchemaFileName
    implicit none
    type (quick_qcschema_type), intent(inout) :: self
    integer, intent(out) :: ierr
    ! type(json_value),pointer :: inp


    ! -- TODO working code, for debug
    ! call self%json%create_object(inp,'finalizeQCSchema')
    ! call self%json%add(self%p, inp) !add it to the root

    ! ! ! -- * add some data to inputs:
    ! call self%json%add(inp, 'tutus0', 0.1_wp)
    ! call self%json%add(inp, 'tutusf', 1.1_wp)

    ! nullify(inp)  !don't need this anymore


    ! -- * write the file:
    call self%json%print(self%p, qcSchemaFileName)

    ! -- * cleanup:
    call self%json%destroy(self%p)
    if (self%json%failed()) stop 1


    ! deallocate memory
    if(allocated(self%atom_symbol)) deallocate(self%atom_symbol)
    if(allocated(self%nscf_snapshots)) deallocate(self%nscf_snapshots)
    if(allocated(self%e_snapshots)) deallocate(self%e_snapshots)
    if(allocated(self%xyz_snapshots)) deallocate(self%xyz_snapshots)


    self%testy = 10
    ! close file
    ! close(self%iQCSchemaFile)


end subroutine finalize_qcschema

subroutine reorder_mo_coeffs(co, KLMN, nbasis, i, reord_mo_vec, ierr)

    implicit none
    double precision, intent(in) :: co(:,:)
    integer, intent(in) :: KLMN(:,:)
    integer, intent (in) :: nbasis, i 
    double precision, intent(inout) :: reord_mo_vec(:)
    integer, intent(out) :: ierr
    integer :: j
   
    j=1
    do while (j <= nbasis)
        reord_mo_vec(j)=co(j,i)
        ! order d functions. xx is the first, followed by yy, zz, xy, xz, yz
        if(KLMN(1,j) .eq. 2) then
            reord_mo_vec(j+1)=co(j+2,i)
            reord_mo_vec(j+2)=co(j+5,i)
            reord_mo_vec(j+3)=co(j+1,i)
            reord_mo_vec(j+4)=co(j+3,i)
            reord_mo_vec(j+5)=co(j+4,i)
            j=j+5
        endif

        ! order f functions. xxx is the first, followed by yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
        if(KLMN(1,j) .eq. 3) then
            reord_mo_vec(j+1)=co(j+3,i) 
            reord_mo_vec(j+2)=co(j+9,i) 
            reord_mo_vec(j+3)=co(j+2,i) 
            reord_mo_vec(j+4)=co(j+1,i) 
            reord_mo_vec(j+5)=co(j+4,i)
            reord_mo_vec(j+6)=co(j+7,i)
            reord_mo_vec(j+7)=co(j+8,i)
            reord_mo_vec(j+8)=co(j+6,i)
            reord_mo_vec(j+9)=co(j+5,i)
            j=j+9
        endif

        j=j+1
    enddo 
    
end subroutine reorder_mo_coeffs

end module quick_qcschema_module
