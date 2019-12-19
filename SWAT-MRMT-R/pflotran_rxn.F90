!=======================================================================
! PFLOTRAN v2.0 LA-CC-09-047
!=======================================================================

!Copyright 2009. Los Alamos National Security, LLC. This material was produced under U.S. 
!Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is operated 
!by Los Alamos National Security, LLC for the U.S. Department of Energy. The U.S. Government has 
!rights to use, reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS 
!NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE 
!USE OF THIS SOFTWARE.  If software is modified to produce derivative works, such modified software 
!should be clearly marked, so as not to confuse it with the version available from LANL.
!Additionally, this library is free software; you can redistribute it and/or modify it under the 
!terms of the GNU Lesser General Public License as published by the Free Software Foundation; 
!either version 2.1 of the License, or (at your option) any later version. Accordingly, this 
!library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even 
!the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser 
!General Public License for more details.

! Send all bug reports/questions/comments to:
!
! Peter C. Lichtner
! Los Alamos National Laboratory
! Earth and Environmental Sciences
! EES-16, MS: D469
! (505) 667-3420
! lichtner@lanl.gov
! Los Alamos, NM

! or

! Glenn E. Hammond
! Sandia National Laboratories
! Applied Systems Analysis & Research
! 413 Cherry Blossom Lp
! Richland, WA 99352
! (505) 235-0665
! gehammo@sandia.gov

!=======================================================================


module BatchChem

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: BatchChemInitializeReactions, &
            BatchChemProcessConstraints

contains

! ************************************************************************** !

subroutine BatchChemInitializeReactions(option, input, reaction)

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Reaction_module
  use Reaction_Aux_module
  use Reaction_Database_module
  use Option_module
  use Input_Aux_module
  use String_module

  implicit none

  type(option_type), pointer :: option
  type(input_type), pointer :: input
  type(reaction_type), pointer :: reaction
  character(len=MAXSTRINGLENGTH) :: string

  ! check for a chemistry block in the  input file
  string = "CHEMISTRY"
  call InputFindStringInFile(input, option, string)
  if (.not.InputError(input)) then
    ! found a chemistry block, initialize the chemistry.

    ! NOTE(bja): ReactionInit() only does a first pass through the
    ! input file to check for a select items
    call ReactionInit(reaction, input, option)
    ! rewind the input file to prepare for the second pass
    call InputFindStringInFile(input, option, string)
    ! the second pass through the input file to read the remaining blocks
    call ReactionReadPass2(reaction, input, option)
  else
     ! TODO(bja): no chemistry block --> fatal error
  endif
    
  if (associated(reaction)) then
    if (reaction%use_full_geochemistry) then
       call DatabaseRead(reaction, option)
       call BasisInit(reaction, option)    
    else
      ! NOTE(bja): do we need this for the batch chemistry driver?

      ! turn off activity coefficients since the database has not been read
      reaction%act_coef_update_frequency = ACT_COEF_FREQUENCY_OFF
      allocate(reaction%primary_species_print(option%ntrandof))
      reaction%primary_species_print = PETSC_TRUE
    endif
  endif

end subroutine BatchChemInitializeReactions

! ************************************************************************** !

subroutine BatchChemProcessConstraints(option, input, reaction, &
     global_auxvars, rt_auxvars, material_auxvars, transport_constraints, &
     constraint_coupler)

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Reaction_module
  use Reaction_Aux_module
  use Reaction_Database_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_class
  use Transport_Constraint_module
  use Option_module
  use Input_Aux_module
  use String_module

  implicit none


  type(option_type), pointer :: option
  type(input_type), pointer :: input
  type(reaction_type), pointer :: reaction
  character(len=MAXSTRINGLENGTH) :: string
  type(global_auxvar_type), pointer :: global_auxvars
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars
  class(material_auxvar_type), pointer :: material_auxvars

  character(len=MAXWORDLENGTH) :: card
  character(len=MAXWORDLENGTH) :: word
  type(tran_constraint_type), pointer :: tran_constraint
  type(tran_constraint_list_type), pointer :: transport_constraints
  type(tran_constraint_coupler_type), pointer :: constraint_coupler 
  PetscBool :: use_prev_soln_as_guess
  PetscInt :: num_iterations
  

  !
  ! read the constraints...
  !

  ! look through the input file
  rewind(input%fid)        
  do
    call InputReadPflotranString(input, option)
    if (InputError(input)) exit

    call InputReadWord(input, option, word, PETSC_FALSE)
    call StringToUpper(word)
    card = trim(word)

    option%io_buffer = 'pflotran card:: ' // trim(card)
    call printMsg(option)

    select case(trim(card))
      case('CONSTRAINT')
        if (.not.associated(reaction)) then
          option%io_buffer = 'CONSTRAINTs not supported without CHEMISTRY.'
          call printErrMsg(option)
        endif
        tran_constraint => TranConstraintCreate(option)
        call InputReadWord(input, option, tran_constraint%name, PETSC_TRUE)
        call InputErrorMsg(input, option, 'constraint', 'name') 
        call printMsg(option, tran_constraint%name)
        call TranConstraintRead(tran_constraint, reaction, input, option)
        call TranConstraintAddToList(tran_constraint, transport_constraints)
        nullify(tran_constraint)

      case default
         ! do nothing
    end select
  enddo

  !
  ! process constraints
  !
  num_iterations = 0
  use_prev_soln_as_guess = PETSC_FALSE
  tran_constraint => transport_constraints%first
  ! NOTE(bja): we only created one set of global and rt auxvars, so if
  ! there is more than one constratint in the input file, they will be
  ! over written.
  do 
     if (.not. associated(tran_constraint)) exit
     ! initialize constraints
     option%io_buffer = "initializing constraint : " // tran_constraint%name
     call printMsg(option)
     call ReactionProcessConstraint(reaction, &
                                    tran_constraint%name, &
                                    tran_constraint%aqueous_species, &
                                    tran_constraint%free_ion_guess, &
                                    tran_constraint%minerals, &
                                    tran_constraint%surface_complexes, &
                                    tran_constraint%colloids, &
                                    tran_constraint%immobile_species, &
                                    option)

     ! link the constraint to the constraint coupler
     constraint_coupler%constraint_name = tran_constraint%name
     constraint_coupler%aqueous_species => tran_constraint%aqueous_species
     constraint_coupler%free_ion_guess => tran_constraint%free_ion_guess
     constraint_coupler%minerals => tran_constraint%minerals
     constraint_coupler%surface_complexes => tran_constraint%surface_complexes
     constraint_coupler%colloids => tran_constraint%colloids
     constraint_coupler%global_auxvar => global_auxvars
     constraint_coupler%rt_auxvar => rt_auxvars
     
     ! equilibrate
     option%io_buffer = "equilibrate constraint : " // tran_constraint%name
     call printMsg(option)
     call ReactionEquilibrateConstraint(rt_auxvars, global_auxvars, &
                                        material_auxvars, reaction, &
                                        tran_constraint%name, &
                                        tran_constraint%aqueous_species, &
                                        tran_constraint%free_ion_guess, &
                                        tran_constraint%minerals, &
                                        tran_constraint%surface_complexes, &
                                        tran_constraint%colloids, &
                                        tran_constraint%immobile_species, &
                                        num_iterations, &
                                        use_prev_soln_as_guess, &
                                        option)
     call ReactionPrintConstraint(constraint_coupler, reaction, option)
     tran_constraint => tran_constraint%next
  enddo

end subroutine BatchChemProcessConstraints


end module BatchChem

module pflotran_rxn
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Reaction_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_class
  use Reaction_Database_module
  use Option_module
  use Input_Aux_module
  use String_module
  
  use Transport_Constraint_module
  use PFLOTRAN_Constants_module

  use BatchChem
  !use parm
  implicit none
  public :: pflotran_rxn_init, &
            pflotran_rxn_react, &
            pflotran_rxn_cleanup
  type(reaction_type), public, pointer :: reaction
  type(option_type), public, pointer :: option
  type(input_type), pointer :: input

  type(global_auxvar_type), public, pointer :: global_auxvars
  type(reactive_transport_auxvar_type), public, pointer :: rt_auxvars
  class(material_auxvar_type), public, pointer :: material_auxvars

contains
! ************************************************************************** !
subroutine pflotran_rxn_init( )
  
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Reaction_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_class
  use Reaction_Database_module
  use Option_module
  use Input_Aux_module
  use String_module
  
  use Transport_Constraint_module
  use PFLOTRAN_Constants_module

  use BatchChem
  use parm, only: nstep, mch
  use parm, only: hhimmobile, cimmobile
  use parm, only: qs_v,qs_l,ts_v,ts_l,xsarea,ch_l2,hharea !related to nexss calculation

  implicit none


  PetscErrorCode :: ierr
  PetscBool :: option_found  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH) :: filename_out
  !type(reaction_type), pointer :: reaction
  !type(option_type), pointer :: option
  type(input_type), pointer :: input

  !type(global_auxvar_type), pointer :: global_auxvars
  !type(reactive_transport_auxvar_type), pointer :: rt_auxvars
  !class(material_auxvar_type), pointer :: material_auxvars

  character(len=MAXWORDLENGTH) :: card
  character(len=MAXWORDLENGTH) :: word
  type(tran_constraint_list_type), pointer :: transport_constraints
  type(tran_constraint_coupler_type), pointer :: constraint_coupler
! 
  PetscInt :: istart, iend, iendaq
  PetscInt :: immobile_start, immobile_end
  PetscInt :: step, max_steps, num_iterations
  PetscReal, dimension(:), allocatable :: tran_xx_p
  PetscReal, dimension(:), allocatable :: cimmobile_init
  PetscReal :: t_max, dt
  PetscInt :: jrch,i,itm_tmp
  PetscReal :: xv_mu,xl_mu, x_sigma, x_qv, x_ql
!
  option => OptionCreate()
  option%fid_out = OUT_UNIT


  call MPI_Init(ierr)
  option%global_comm = MPI_COMM_WORLD
  call MPI_Comm_rank(MPI_COMM_WORLD, option%global_rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, option%global_commsize, ierr)
  call MPI_Comm_group(MPI_COMM_WORLD, option%global_group, ierr)
  option%mycomm = option%global_comm
  option%myrank = option%global_rank
  option%mycommsize = option%global_commsize
  option%mygroup = option%global_group

  ! check for non-default input filename
  option%input_filename = "pflotran.in"
  string = '-pflotranin'
  call InputGetCommandLineString(string, option%input_filename, option_found, option)

  string = '-output_prefix'
  call InputGetCommandLineString(string, option%global_prefix, option_found, option)

  PETSC_COMM_WORLD = MPI_COMM_WORLD
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr);CHKERRQ(ierr)

  input => InputCreate(IN_UNIT, option%input_filename, option)

  filename_out = trim(option%global_prefix) // trim(option%group_prefix) // &
                 '.out'

  if (option%myrank == option%io_rank .and. option%print_to_file) then
    open(option%fid_out, file=filename_out, action="write", status="unknown")
  endif

  !
  ! manual initialization...
  !
  option%nphase = 1
  option%liquid_phase = 1
  option%reference_water_density = 998.2

  call BatchChemInitializeReactions(option, input, reaction)

  !
  ! create the storage containers
  !
  ! NOTE(bja) : batch chem --> one cell

  ! global_auxvars --> cell by cell temperature, pressure, saturation, density
  allocate(global_auxvars)
  call GlobalAuxVarInit(global_auxvars, option)

  ! rt_auxvars --> cell by cell chemistry data
  allocate(rt_auxvars)
  call RTAuxVarInit(rt_auxvars, reaction, option)

  ! material_auxvars --> cell by cell material property data
  allocate(material_auxvars)
  call MaterialAuxVarInit(material_auxvars, option)
  
  material_auxvars%porosity = 1.d0
  material_auxvars%volume= 1.d0

  ! assign default state values
  global_auxvars%pres = option%reference_pressure
  global_auxvars%temp = option%reference_temperature
  ! global_auxvars%den_kg = option%reference_water_density
  ! NOTE(bja): option%ref_density = 0.0, so we set it manually. This is a Bad Thing(TM)
  global_auxvars%den_kg = 998.2
  global_auxvars%sat = option%reference_saturation  

  ! create the constraint list
  !allocate(transport_constraints)
  !call TranConstraintInitList(transport_constraints)
  !allocate(constraint_coupler)
  !constraint_coupler => TranConstraintCouplerCreate(option)

  !call BatchChemProcessConstraints(option, input, reaction, &
  !                                 global_auxvars, rt_auxvars, &
  !                                 material_auxvars, transport_constraints, &
  !                                 constraint_coupler)
! initialize - should be in a file in the future
! cimmobile - in each channel, hhimmobile - subtimestep xxxx within storage zone
! 
  allocate(cimmobile(reaction%immobile%nimmobile,mch))
  allocate(hhimmobile(reaction%immobile%nimmobile,nstep))
  allocate(cimmobile_init(reaction%immobile%nimmobile))
  open(9105,file='hz_init_conc.dat',status='old')
  read(9105,*) cimmobile_init(1:reaction%immobile%nimmobile)
  do i=1,mch
    cimmobile(:,i) = cimmobile_init(:)
  enddo
!  do i=1,nstep
    hhimmobile(:,1) = cimmobile_init(:)
!  enddo
! nexss related input
  open(5000,file='nxs2swat.txt',status='unknown')
!
! itm_tmp: number of sub-hyporheic zones
! x_sigma: standard deviation of residence time in hr

  read(5000,*) itm_tmp, x_sigma

  allocate(qs_v(mch,itm_tmp))
  allocate(qs_l(mch,itm_tmp))
  allocate(ts_l(mch,itm_tmp))
  allocate(ts_v(mch,itm_tmp))
  allocate(xsarea(mch,itm_tmp))
!mean
  do i=1,1
  do jrch=1,mch-1 !one reach not simulated in swat
! qs_v: vertical exchange rate (m^3/s)
! qs_l: lateral exchange rate (m^3/s)
! ts_v: vertical residence time (hr)
! ts_l: lateral residence time (hr)
! reach cross sectional area (m^2), not used
    read(5000,*) qs_v(jrch,i),qs_l(jrch,i),ts_v(jrch,i),ts_l(jrch,i),xsarea(jrch,i)
  enddo
  enddo
  do jrch = 1,mch-1 
    x_qv = qs_v(jrch,1)
    x_ql = qs_l(jrch,1)
    qs_v(jrch,:) = x_qv/itm_tmp !evenly assign flux
    qs_l(jrch,:) = x_ql/itm_tmp
  enddo
  if(itm_tmp>1) then
  do jrch = 1, mch-1
     xv_mu = max(ts_v(jrch,1),1.d-20)
     xl_mu = max(ts_l(jrch,1),1.d-20)
!     do i=1,itm_tmp
!       if(xv_mu > 1.e4 .or. xv_mu == 0.d0) then
!        ts_v(jrch,i) = xv_mu
!       else
        !call transfer_coe1(itm_tmp,dlog(xv_mu),x_sigma,ts_v(jrch,:)) 
        call transfer_coex(itm_tmp,dlog(xv_mu),x_sigma,ts_v(jrch,:)) 
!       endif
!       if(xl_mu > 1.e4 .or. xl_mu == 0.d0) then
!        ts_l(jrch,i) = xl_mu
!       else
         !call transfer_coe1(itm_tmp,dlog(xl_mu),x_sigma,ts_l(jrch,:)) 
         call transfer_coex(itm_tmp,dlog(xl_mu),x_sigma,ts_l(jrch,:)) 
!       endif
!     enddo  
!test uniform 
!    ts_v(jrch,:) = xv_mu
!    ts_l(jrch,:) = xl_mu
  enddo
  endif

  ts_v(:,:) = ts_v(:,:)*3600.d0 !s
  ts_l(:,:) = ts_l(:,:)*3600.d0 !s
!print *,'ts_v',ts_v
!print *,'ts_l',ts_l
!stop
  close(5000)
  !rt_auxvars%immobile(1:reaction%immobile%nimmobile) = 1.d-5
end subroutine pflotran_rxn_init


subroutine pflotran_rxn_react(c_swat_sp, swat_params, dt, volume,cimmob,num_mrmt)
  
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Reaction_module
  use Reaction_Aux_module
  use Reactive_Transport_aux_module
  use Global_Aux_module
  use Material_Aux_class
  use Reaction_Database_module
  use Option_module
  use Input_Aux_module
  use String_module
  
  use Transport_Constraint_module
  use PFLOTRAN_Constants_module

  use BatchChem
  !use parm

  implicit none


  PetscErrorCode :: ierr
  PetscBool :: option_found  
  character(len=MAXSTRINglength) :: string
  character(len=MAXSTRINglength) :: filename_out
  !type(reaction_type), pointer :: reaction
  !type(option_type), pointer :: option
  !type(input_type), pointer :: input

  !type(global_auxvar_type), pointer :: global_auxvars
  !type(reactive_transport_auxvar_type), pointer :: rt_auxvars
  !class(material_auxvar_type), pointer :: material_auxvars

  character(len=MAXWORDLength) :: card
  character(len=MAXWORDLength) :: word
  type(tran_constraint_list_type), pointer :: transport_constraints
  type(tran_constraint_coupler_type), pointer :: constraint_coupler
! 
  PetscInt :: istart, iend, iendaq
  PetscInt :: immobile_start, immobile_end
  PetscInt :: step, max_steps, num_iterations
  PetscInt :: jrch
  PetscReal, dimension(reaction%ncomp) :: tran_xx_p
  PetscReal, dimension(reaction%ncomp) :: tran_work
  PetscReal :: t_max, dt, timed
  PetscReal :: Cno3, src_no3, cimmob(reaction%immobile%nimmobile)
  PetscReal :: molality_to_molarity, volume, day2sec
  PetscReal, dimension(12) :: c_swat_sp
  PetscReal, dimension(23+100) :: swat_params
  PetscInt :: i
  PetscInt :: icvg
  PetscReal :: dt_work, dt_sum
  PetscReal :: dt_acc = 1.25d0 ! accelerator
  PetscInt :: icut
  PetscInt :: cnt_mrmt, num_mrmt
!
!  timed = 0.d0
! convert hour to sec
  
  dt = dt * 3600.d0
  dt_work = dt
  
!
  istart = 1
  iend = reaction%ncomp
  iendaq = istart + reaction%naqcomp - 1
  if (reaction%immobile%nimmobile > 0) then
      immobile_start = reaction%offset_immobile + 1
      immobile_end = reaction%offset_immobile + reaction%immobile%nimmobile
  endif
! mg/L to mol/L
  rt_auxvars%total(1,1) = max(c_swat_sp(1), 1.d-15)*1.d-3/14.d0 !nh4
  rt_auxvars%total(2,1) = max(c_swat_sp(2), 1.d-15)*1.d-3/32.d0 !o2
  rt_auxvars%total(3,1) = max(c_swat_sp(3), 1.d-15)*1.d-3/14.d0 !no3
  rt_auxvars%total(4,1) = max(c_swat_sp(4), 1.d-15)*1.d-3/14.d0 !no2
  rt_auxvars%total(5,1) = max(c_swat_sp(5), 1.d-15)*1.d-3/14.d0 !orgn
! kept in mg
  rt_auxvars%total(6,1) = max(c_swat_sp(6), 1.d-15) !algae
  rt_auxvars%total(7,1) = max(c_swat_sp(7), 1.d-15)*1.d-3/30.974d0 !orgp
  rt_auxvars%total(8,1) = max(c_swat_sp(8), 1.d-15)*1.d-3/30.974d0 !solp
! mg O2 for CBOD (biochemical oxygen demand)
  rt_auxvars%total(9,1) = max(c_swat_sp(9), 1.d-15)*1.d-3/32.d0 !o2
!  rt_auxvars%total(10,1) = rt_auxvars%total(9,1) ! doc
!Fang - fix concentration
!  rt_auxvars%total(10,1) = 1.1e-4 ! doc
  rt_auxvars%total(10,1) = max(c_swat_sp(10),1.d-15)*1.d-3/12.d0  !doc 
  rt_auxvars%total(11,1) = max(c_swat_sp(11), 1.d-15)*1.d-3/44.d0 !co2
  rt_auxvars%total(12,1) = max(c_swat_sp(12), 1.d-15)*1.d-3/14.d0 !n2
! molarlity - mol/kg.  molarity - mol/volume => convert to mol/L
  molality_to_molarity = global_auxvars%den_kg(1)*1.d-3
! convert to mol/kg
  rt_auxvars%pri_molal(1:5) = rt_auxvars%total(1:5,1)/molality_to_molarity 
  rt_auxvars%pri_molal(6) = rt_auxvars%total(6,1)
  rt_auxvars%pri_molal(7:) = rt_auxvars%total(7:,1)/molality_to_molarity 
  tran_xx_p(istart:iendaq) =  rt_auxvars%total(istart:iendaq, 1)
! assign storage zone conc
  if (reaction%immobile%nimmobile > 0) then
      rt_auxvars%immobile(1:reaction%immobile%nimmobile) = &
        cimmob(1:reaction%immobile%nimmobile)
!Fang - assign HZ concentration at 2011 first day
!      if(option%time == 3600*365*2*24.0) then
!        rt_auxvars%immobile(6) = 3.5e-4
!      endif
      tran_xx_p(immobile_start:immobile_end) = &
         rt_auxvars%immobile(1:reaction%immobile%nimmobile)
  endif
  tran_work(:) = tran_xx_p(:)
!    i = reaction%nauxiliary - 23 - 4 ! 4 = two sets of nexss parameters
    i = reaction%nauxiliary - 23 - 4*num_mrmt ! two sets of nexss parameters
! stoichiometric coefficients
    rt_auxvars%auxiliary_data(i+1:i+7) = swat_params(1:7)
! rate constants
    rt_auxvars%auxiliary_data(i+8) = swat_params(8) /3600.d0 !rho 1/s
    rt_auxvars%auxiliary_data(i+9) = swat_params(9) /3600.d0 !mu
    rt_auxvars%auxiliary_data(i+10) = swat_params(10) /3600.d0 !beta1
    rt_auxvars%auxiliary_data(i+11) = swat_params(11) /3600.d0 !beta2
    rt_auxvars%auxiliary_data(i+12) = swat_params(12) /3600.d0 !beta4
    rt_auxvars%auxiliary_data(i+13) = swat_params(13) /3600.d0 !k1
    rt_auxvars%auxiliary_data(i+14) = swat_params(14) /3600.d0 !k3
    rt_auxvars%auxiliary_data(i+15) = swat_params(15) /3600.d0 !sig1OVd
    rt_auxvars%auxiliary_data(i+16) = swat_params(16) /3600.d0 !sig5
    rt_auxvars%auxiliary_data(i+17) = swat_params(17) /3600.d0 !sig4
    rt_auxvars%auxiliary_data(i+18) = swat_params(18) * 1.d-3/30.974d0/3600.d0 !sig2OVd, mol/L/s
    rt_auxvars%auxiliary_data(i+19) = swat_params(19) /3600.d0 !k2
    rt_auxvars%auxiliary_data(i+20) = swat_params(20) * 1.d-3/32.d0 !eq, mol/L
    rt_auxvars%auxiliary_data(i+21) = swat_params(21) *1.d-3/32.d0/3600.d0 !k4d, mol/l/s
    rt_auxvars%auxiliary_data(i+22) = swat_params(22) /3600.d0 !beta3
    rt_auxvars%auxiliary_data(i+23) = swat_params(23) * 1.d-3/14.d0/3600.d0 !sig3Ovd, mol/L/s
    !num_mrmt = (size(swat_params)-23)/4 
    do cnt_mrmt=1,num_mrmt !mrmt in each storage zone
      rt_auxvars%auxiliary_data(i+23+cnt_mrmt) = swat_params(23+cnt_mrmt) !nexss transfer rate in vertical
      rt_auxvars%auxiliary_data(i+23+num_mrmt+cnt_mrmt) = swat_params(23+num_mrmt+cnt_mrmt) !transfer rate in lateral
      rt_auxvars%auxiliary_data(i+23+2*num_mrmt+cnt_mrmt) = swat_params(23+2*num_mrmt+cnt_mrmt) !a/as in vertical (ratio of cross section area)
      rt_auxvars%auxiliary_data(i+23+3*num_mrmt+cnt_mrmt) = swat_params(23+3*num_mrmt+cnt_mrmt) !a/as in lateral
   enddo
!! convert to m^3
!    material_auxvars%volume = volume/1000.d0
    material_auxvars%volume = volume
  ! stepping
    
!   if (reaction%immobile%nimmobile > 0) then
!      tran_xx_p(immobile_start:immobile_end) = &
!         rt_auxvars%immobile(1:reaction%immobile%nimmobile)
!    endif
    dt_sum = 0.d0
    icut = 0
    do
      option%tran_dt = dt
      if(dt < 0.1) then
          print *, 'Time step is too small in pflotran_rxn_react: stop: ',dt
          stop
      endif

      call RReact(rt_auxvars,global_auxvars, &
                material_auxvars, &
                tran_xx_p(istart:iend), &
                num_iterations,reaction,option,icvg)
      dt_sum = dt_sum + dt
      if(icvg == 0 .and. dt_sum == dt_work) exit
      if(icvg == 0) then 
        dt = min( dt*(1.d0+dt_acc), dt_work-dt_sum)
        tran_xx_p(istart:iendaq) = rt_auxvars%total(istart:iendaq, 1)
        if (reaction%immobile%nimmobile > 0) then
         tran_xx_p(immobile_start:immobile_end) = rt_auxvars%immobile
        endif
        tran_work(:) = tran_xx_p(:)
        icut = 0
        call RUpdateKineticstate(rt_auxvars, &
                             global_auxvars, &
                             material_auxvars, &
                             reaction,option)
      else
        dt_sum = dt_sum - dt
        dt = dt/2.d0
        tran_xx_p(:) = tran_work(:)
        if(icut > 8) then
          print *, 'Too many time step cuts in pflotran_rxn_react: stop: ',icut
          stop
        endif
        icut = icut + 1
       
      endif
    enddo


    ! set primary dependent var back to free-ion molality
    tran_xx_p(istart:iendaq) = rt_auxvars%total(istart:iendaq, 1)

    C_swat_sp(1:12) = tran_xx_p(istart:iendaq)

    if (reaction%immobile%nimmobile > 0) then
      tran_xx_p(immobile_start:immobile_end) = rt_auxvars%immobile
    endif
    if (reaction%immobile%nimmobile > 0) then
      cimmob(1:reaction%immobile%nimmobile) = tran_xx_p(immobile_start:immobile_end)
    end if

end subroutine pflotran_rxn_react

subroutine pflotran_rxn_cleanup
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Reaction_module
  use Reaction_Aux_module
  use Reactive_Transport_aux_module
  use Global_Aux_module
  use Material_Aux_class
  use Reaction_Database_module
  use Option_module
  use Input_Aux_module
  use String_module
  
  use Transport_Constraint_module
  use PFLOTRAN_Constants_module

  use BatchChem

  implicit none

  PetscErrorCode :: ierr
  PetscBool :: option_found  
  character(len=MAXSTRINglength) :: string
  character(len=MAXSTRINglength) :: filename_out
  type(reaction_type), pointer :: reaction
  type(option_type), pointer :: option
  type(input_type), pointer :: input

  type(global_auxvar_type), pointer :: global_auxvars
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars
  class(material_auxvar_type), pointer :: material_auxvars

  character(len=MAXWORDLength) :: card
  character(len=MAXWORDLength) :: word
  type(tran_constraint_list_type), pointer :: transport_constraints
  type(tran_constraint_coupler_type), pointer :: constraint_coupler
! 
  PetscInt :: istart, iend, iendaq
  PetscInt :: immobile_start, immobile_end
  PetscInt :: step, max_steps, num_iterations
  PetscReal, dimension(:), allocatable :: tran_xx_p
  PetscReal :: t_max, dt
!
  call TranConstraintCouplerdestroy(constraint_coupler)
  call TranConstraintDestroylist(transport_constraints)
  ! FIXME(bja) : causes error freeing memory.
  !call RTAuxVarDestroy(rt_auxvars)
  !call GlobalAuxVarDestroy(global_auxvars)
  call ReactionDestroy(reaction,option)
  call GlobalAuxVarStrip(global_auxvars)
  deallocate(global_auxvars)
  nullify(global_auxvars)
  call RTAuxVarStrip(rt_auxvars)
  deallocate(rt_auxvars)
  nullify(rt_auxvars)
  call MaterialAuxVarStrip(material_auxvars)
  deallocate(material_auxvars)
  nullify(material_auxvars)
  call InputDestroy(input)
  call OptionDestroy(option)
  call PetscFinalize(ierr);chkerrQ(ierr)
  call MPI_Finalize(ierr)
!
end subroutine pflotran_rxn_cleanup

        subroutine transfer_coe1(ME,mu,sigma,am)        !lognormal distribution
        implicit none
        integer ME,i,j
        double precision mu,sigma,am(ME)
        double precision aa,x1,x2,f1,f2,dxx,rtbis,xmid,fmid
!        double precision prob
!
          write(*,*) 'Intermediate multi-rate parameters'
        aa=0.5d0/dble(ME)
        do i=1,ME
                x1=1.d-29
                x2=x1+1.d29
                f1=prob(x1,mu,sigma)-aa
                f2=prob(x2,mu,sigma)-aa
                do while (f2.lt.0.d0)
                        x2=x1+0.5*(x2-x1)
                        f2=prob(x2,mu,sigma)
                enddo
                rtbis=x1
                dxx=(x2-x1)
                do j=1,500
                        dxx=0.5d0*dxx
                        xmid=rtbis+dxx
                        fmid=prob(xmid,mu,sigma)-aa
                        if(fmid.le.0.d0) rtbis=xmid
                        if(dabs(dxx).lt.1.d-50) goto 10
                enddo
   10           continue
                write(*,*) i,rtbis,prob(rtbis,mu,sigma)-aa,j
                am(i)=rtbis
                aa=aa+1.d0/dble(ME)
        enddo
       return
       end
        
       function prob(x,mu,sigma)
        implicit none
        integer i,j
        double precision mu,sigma,x
        double precision a1,error,prob
!
        a1=(dlog(x)-mu)/(dsqrt(2.d0)*sigma)
        prob = 0.5d0*(1.d0+erf(a1))
        return
        end

        subroutine transfer_coex(ME,mu,sigma,am)        !exponentialdistribution
        implicit none
        integer ME,i,j
        double precision mu,sigma,am(ME)
        double precision aa,x1,x2,f1,f2,dxx,rtbis,xmid,fmid
!        double precision prob
!
          write(*,*) 'Intermediate multi-rate parameters'
        aa=0.5d0/dble(ME)
        do i=1,ME
                x1=1.d-29
                x2=x1+1.d29
                f1=prob_exp(x1,mu)-aa
                f2=prob_exp(x2,mu)-aa
                do while (f2.lt.0.d0)
                        x2=x1+0.5*(x2-x1)
                        f2=prob_exp(x2,mu)
                enddo
                rtbis=x1
                dxx=(x2-x1)
                do j=1,500
                        dxx=0.5d0*dxx
                        xmid=rtbis+dxx
                        fmid=prob_exp(xmid,mu)-aa
                        if(fmid.le.0.d0) rtbis=xmid
                        if(dabs(dxx).lt.1.d-50) goto 10
                enddo
   10           continue
!                write(*,*) i,rtbis,prob_exp(rtbis,mu,sigma)-aa,j
                am(i)=rtbis
                aa=aa+1.d0/dble(ME)
        enddo
       return
       end


        function prob_exp(x,mu)
        implicit none
        integer i,j
        double precision mu,sigma,x
        double precision a1,error,prob_exp
!
        
        prob_exp = 1.0-exp(-x/mu)
        return
        end
end module pflotran_rxn
