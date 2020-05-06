module Reaction_Sandbox_Cyber_Mrmt_class

  use Reaction_Sandbox_Base_class
  
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  
  use PFLOTRAN_Constants_module
  use Utility_module

  implicit none
  
  private
  
#include "petsc/finclude/petscsys.h"

  PetscInt, parameter :: DOC_MASS_STORAGE_INDEX = 1
  PetscInt, parameter :: NH4_MASS_STORAGE_INDEX = 2
  PetscInt, parameter :: O2_MASS_STORAGE_INDEX = 3
  PetscInt, parameter :: NO3_MASS_STORAGE_INDEX = 4
  PetscInt, parameter :: NO2_MASS_STORAGE_INDEX = 5
  PetscInt, parameter :: CO2_MASS_STORAGE_INDEX = 6
  PetscInt, parameter ::  MRMT_T = 200
  PetscInt :: MRMT = 0
type myptr
    PetscReal, pointer :: ptr
end type myptr

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_cyber_mrmt_type
    PetscInt :: nh4_id
    PetscInt :: o2_id
    PetscInt :: no3_id
    PetscInt :: no2_id
    PetscInt :: n2_id
    PetscInt :: doc_id
    PetscInt :: biomass_id
    PetscInt :: co2_id
    PetscInt :: a_id
    PetscInt :: dop_id
    PetscInt :: cbod_id
    PetscInt :: orgn_id
    PetscInt :: orgp_id
    PetscInt, dimension(MRMT_T) :: nh4_st_id
    PetscInt, dimension(MRMT_T) :: o2_st_id
    PetscInt, dimension(MRMT_T) :: no3_st_id
    PetscInt, dimension(MRMT_T) :: no2_st_id
    PetscInt, dimension(MRMT_T) :: n2_st_id
    PetscInt, dimension(MRMT_T) :: doc_st_id
    PetscInt, dimension(MRMT_T) :: biomass_st_id
    PetscInt, dimension(MRMT_T) :: co2_st_id
    PetscInt :: carbon_consumption_species_id
    character(len=MAXWORDLENGTH) :: carbon_consumption_species
    PetscReal :: f1
    PetscReal :: f2
    PetscReal :: f3
!    PetscReal :: f_act   ! fraction of active biomass
    PetscReal :: k_deg   ! biomass degradation rate
    PetscReal :: k1      ! nitrate rate constant
    PetscReal :: k2      ! nitrite rate constant
    PetscReal :: k3      ! oxygen rate constant
    PetscReal :: Kd1
    PetscReal :: Ka1
    PetscReal :: Kd2
    PetscReal :: Ka2
    PetscReal :: Kd3
    PetscReal :: Ka3
    PetscReal :: Ki2
    PetscReal :: Ki1
    !PetscReal, dimension(MRMT) :: alpha
    PetscReal, pointer :: alpha(:)
    PetscReal, pointer :: th_frac(:)
!    PetscReal :: no3_src
    PetscReal :: stoich_1_doc
    PetscReal :: stoich_1_nh4
    PetscReal :: stoich_1_no3
    PetscReal :: stoich_1_no2
    PetscReal :: stoich_1_co2
    PetscReal :: stoich_1_biomass
    PetscReal :: stoich_2_doc
    PetscReal :: stoich_2_nh4
    PetscReal :: stoich_2_no2
    PetscReal :: stoich_2_n2
    PetscReal :: stoich_2_co2
    PetscReal :: stoich_2_biomass
    PetscReal :: stoich_3_doc
    PetscReal :: stoich_3_nh4
    PetscReal :: stoich_3_o2
    PetscReal :: stoich_3_co2
    PetscReal :: stoich_3_biomass
!SWAT rxn 1
    PetscReal, pointer :: stoich_sw1_o2  !alpha4
    PetscReal, pointer :: stoich_sw1_a !-1.d0
    PetscReal, pointer :: stoich_sw1_orgn !alpha1
    PetscReal, pointer :: stoich_sw1_orgp !alpha2
!SWAT rxn 2
    PetscReal, pointer :: stoich_sw2_dop !alpha2
    PetscReal, pointer :: stoich_sw2_no3 ! (1.d0-f)*alpha1
    PetscReal, pointer :: stoich_sw2_nh4 ! f*alpha1
    PetscReal, pointer :: stoich_sw2_a !-1.d0
    PetscReal, pointer :: stoich_sw2_o2 !alpha3
!SWAT rxn 3
    PetscReal, pointer :: stoich_sw3_o2 !alpha5
    PetscReal, pointer :: stoich_sw3_nh4 !-1.d0
    PetscReal, pointer :: stoich_sw3_no2 !1.d0
!SWAT rxn 4
    PetscReal, pointer :: stoich_sw4_o2 !alpha6
    PetscReal, pointer :: stoich_sw4_no2 !-1.d0
    PetscReal, pointer :: stoich_sw4_no3 !1.d0
!SWAT rxn 5
    PetscReal, pointer :: stoich_sw5_orgp !-1.d0
    PetscReal, pointer :: stoich_sw5_dop !1.d0
!SWAT rxn 6
    PetscReal, pointer :: stoich_sw6_o2 !-1.d0
    PetscReal, pointer :: stoich_sw6_cbod !-1.d0
!SWAT rxn 7
    PetscReal, pointer :: stoich_sw7_cbod !-1.d0
!SWAT rxn 8
    PetscReal, pointer :: stoich_sw8_orgn !-1.d0
    PetscReal, pointer :: stoich_sw8_nh4 !1.d0
!Dynamic rate constants from SWAT
!    PetscReal, pointer :: rate_cs(:)
!
    PetscReal, dimension(MRMT_T) :: stoich_doc
    PetscReal, dimension(MRMT_T) :: stoich_nh4
    PetscReal, dimension(MRMT_T) :: stoich_no3
    PetscReal, dimension(MRMT_T) :: stoich_no2
    PetscReal, dimension(MRMT_T) :: stoich_n2
    PetscReal, dimension(MRMT_T) :: stoich_co2
    PetscReal, dimension(MRMT_T) :: stoich_o2
    PetscReal, dimension(MRMT_T) :: stoich_biomass
!
    PetscReal, dimension(MRMT_T) :: stoich_mrmt_doc
    PetscReal, dimension(MRMT_T) :: stoich_mrmt_nh4
    PetscReal, dimension(MRMT_T) :: stoich_mrmt_no3
    PetscReal, dimension(MRMT_T) :: stoich_mrmt_no2
    PetscReal, dimension(MRMT_T) :: stoich_mrmt_n2
    PetscReal, dimension(MRMT_T) :: stoich_mrmt_co2
    PetscReal, dimension(MRMT_T) :: stoich_mrmt_o2
    PetscReal, dimension(MRMT_T) :: stoich_mrmt_biomass
!storage zone porosity
    PetscReal :: th_st
! volume fraciton of river to storage zone
    PetscReal, pointer :: vol_frac_r2st
    PetscReal :: activation_energy
    PetscReal :: reference_temperature
    PetscInt :: nrxn
    PetscInt :: nrxn_swat
    PetscInt :: nrxn_mt
    PetscInt :: n_stoic
    PetscInt :: n_rates
    PetscInt :: offset_auxiliary
    PetscBool :: store_cumulative_mass
    PetscInt, pointer :: nrow(:)
    PetscInt, pointer :: ncol(:)
    PetscInt, pointer :: irow(:,:)
    PetscInt, pointer :: icol(:,:)
    type(myptr), dimension(:,:), allocatable :: stoich_row
    type(myptr), dimension(:), allocatable :: vol_frac_mrmt
  contains
    procedure, public :: ReadInput => CyberRead
    procedure, public :: Setup => CyberSetup
    procedure, public :: Evaluate => CyberReact
    procedure, public :: UpdateKineticState => CyberUpdateKineticState
    procedure, public :: AuxiliaryPlotVariables => CyberAuxiliaryPlotVariables
    procedure, public :: Destroy => CyberDestroy
  end type reaction_sandbox_cyber_mrmt_type
  
  public :: CyberMrmtCreate

contains

! ************************************************************************** !

function CyberMrmtCreate()
  ! 
  ! Allocates PNNL N reaction object.
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/01/15
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  implicit none
  
  class(reaction_sandbox_cyber_mrmt_type), pointer :: CyberMrmtCreate

  allocate(CyberMrmtCreate)
  CyberMrmtCreate%o2_id = UNINITIALIZED_INTEGER
  CyberMrmtCreate%nh4_id = UNINITIALIZED_INTEGER
  CyberMrmtCreate%no3_id = UNINITIALIZED_INTEGER
  CyberMrmtCreate%no2_id = UNINITIALIZED_INTEGER
  CyberMrmtCreate%n2_id = UNINITIALIZED_INTEGER
  CyberMrmtCreate%doc_id = UNINITIALIZED_INTEGER
  CyberMrmtCreate%biomass_id = UNINITIALIZED_INTEGER
  CyberMrmtCreate%co2_id = UNINITIALIZED_INTEGER
  CyberMrmtCreate%cbod_id = UNINITIALIZED_INTEGER
  CyberMrmtCreate%a_id = UNINITIALIZED_INTEGER
  CyberMrmtCreate%orgn_id = UNINITIALIZED_INTEGER
  CyberMrmtCreate%orgp_id = UNINITIALIZED_INTEGER
  CyberMrmtCreate%carbon_consumption_species_id = UNINITIALIZED_INTEGER
  CyberMrmtCreate%f1 = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%f2 = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%f3 = UNINITIALIZED_DOUBLE  
!  CyberCreate%f_act = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%k_deg = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%k1 = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%k2 = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%k3 = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%Kd1 = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%Ka1 = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%Kd2 = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%Ki1 = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%Ki2 = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%Ka2 = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%Kd3 = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%Ka3 = UNINITIALIZED_DOUBLE 
! multirate 
  CyberMrmtCreate%o2_st_id(:) = UNINITIALIZED_INTEGER
  CyberMrmtCreate%nh4_st_id(:) = UNINITIALIZED_INTEGER
  CyberMrmtCreate%no3_st_id(:) = UNINITIALIZED_INTEGER
  CyberMrmtCreate%no2_st_id(:) = UNINITIALIZED_INTEGER
  CyberMrmtCreate%n2_st_id(:) = UNINITIALIZED_INTEGER
  CyberMrmtCreate%doc_st_id(:) = UNINITIALIZED_INTEGER
  CyberMrmtCreate%biomass_st_id(:) = UNINITIALIZED_INTEGER
  CyberMrmtCreate%co2_st_id(:) = UNINITIALIZED_INTEGER
!  CyberMrmtCreate%alpha(:) = UNINITIALIZED_DOUBLE 
  CyberMrmtCreate%stoich_1_doc = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%stoich_1_nh4 = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%stoich_1_no3 = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%stoich_1_no2 = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%stoich_1_co2 = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%stoich_1_biomass = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%stoich_2_doc = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%stoich_2_nh4 = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%stoich_2_no2 = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%stoich_2_n2 = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%stoich_2_co2 = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%stoich_2_biomass = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%stoich_3_doc = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%stoich_3_nh4 = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%stoich_3_o2 = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%stoich_3_co2 = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%stoich_3_biomass = UNINITIALIZED_DOUBLE 
!mrmt 
  CyberMrmtCreate%stoich_doc(:) = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%stoich_nh4(:) = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%stoich_no3(:) = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%stoich_no2(:) = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%stoich_co2(:) = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%stoich_o2(:) = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%stoich_n2(:) = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%stoich_biomass(:) = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%stoich_mrmt_doc(:) = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%stoich_mrmt_nh4(:) = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%stoich_mrmt_no3(:) = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%stoich_mrmt_no2(:) = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%stoich_mrmt_co2(:) = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%stoich_mrmt_o2(:) = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%stoich_mrmt_n2(:) = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%stoich_mrmt_biomass(:) = UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%th_st = UNINITIALIZED_DOUBLE  
  allocate(CyberMrmtCreate%vol_frac_r2st)
  CyberMrmtCreate%vol_frac_r2st= UNINITIALIZED_DOUBLE  
  CyberMrmtCreate%activation_energy = UNINITIALIZED_DOUBLE
  CyberMrmtCreate%reference_temperature = 298.15d0 ! 25 C
  CyberMrmtCreate%nrxn = UNINITIALIZED_INTEGER
  CyberMrmtCreate%offset_auxiliary = UNINITIALIZED_INTEGER
  CyberMrmtCreate%carbon_consumption_species = ''
  CyberMrmtCreate%store_cumulative_mass = PETSC_FALSE
!  CyberMrmtCreate%no3_src = UNINITIALIZED_INTEGER
!  CyberMrmtCreate%rate_cs(:) = UNINITIALIZED_INTEGER
!SWAT rxn 1
  allocate(CyberMrmtCreate%stoich_sw1_o2)
  allocate(CyberMrmtCreate%stoich_sw1_a)
  allocate(CyberMrmtCreate%stoich_sw1_orgn)
  allocate(CyberMrmtCreate%stoich_sw1_orgp)
  CyberMrmtCreate%stoich_sw1_o2 = UNINITIALIZED_DOUBLE
  CyberMrmtCreate%stoich_sw1_a = UNINITIALIZED_DOUBLE
  CyberMrmtCreate%stoich_sw1_orgn = UNINITIALIZED_DOUBLE
  CyberMrmtCreate%stoich_sw1_orgp = UNINITIALIZED_DOUBLE
!SWAT rxn 2
  allocate(CyberMrmtCreate%stoich_sw2_dop)
  allocate(CyberMrmtCreate%stoich_sw2_no3)
  allocate(CyberMrmtCreate%stoich_sw2_nh4)
  allocate(CyberMrmtCreate%stoich_sw2_a)
  allocate(CyberMrmtCreate%stoich_sw2_o2)
  CyberMrmtCreate%stoich_sw2_dop = UNINITIALIZED_DOUBLE
  CyberMrmtCreate%stoich_sw2_no3 = UNINITIALIZED_DOUBLE
  CyberMrmtCreate%stoich_sw2_nh4 = UNINITIALIZED_DOUBLE
  CyberMrmtCreate%stoich_sw2_a = UNINITIALIZED_DOUBLE
  CyberMrmtCreate%stoich_sw2_o2 = UNINITIALIZED_DOUBLE
!SWAT rxn 3
  allocate(CyberMrmtCreate%stoich_sw3_o2)
  allocate(CyberMrmtCreate%stoich_sw3_nh4)
  allocate(CyberMrmtCreate%stoich_sw3_no2)
  CyberMrmtCreate%stoich_sw3_o2 = UNINITIALIZED_DOUBLE
  CyberMrmtCreate%stoich_sw3_nh4 = UNINITIALIZED_DOUBLE
  CyberMrmtCreate%stoich_sw3_no2 = UNINITIALIZED_DOUBLE
!SWAT rxn 4
  allocate(CyberMrmtCreate%stoich_sw4_o2)
  allocate(CyberMrmtCreate%stoich_sw4_no2)
  allocate(CyberMrmtCreate%stoich_sw4_no3)
  CyberMrmtCreate%stoich_sw4_o2 = UNINITIALIZED_DOUBLE
  CyberMrmtCreate%stoich_sw4_no2 = UNINITIALIZED_DOUBLE
  CyberMrmtCreate%stoich_sw4_no3 = UNINITIALIZED_DOUBLE
!SWAT rxn 5
  allocate(CyberMrmtCreate%stoich_sw5_orgp)
  allocate(CyberMrmtCreate%stoich_sw5_dop)
  CyberMrmtCreate%stoich_sw5_orgp = UNINITIALIZED_DOUBLE
  CyberMrmtCreate%stoich_sw5_dop = UNINITIALIZED_DOUBLE
!SWAT rxn 6
  allocate(CyberMrmtCreate%stoich_sw6_o2)
  allocate(CyberMrmtCreate%stoich_sw6_cbod)
  CyberMrmtCreate%stoich_sw6_o2 = UNINITIALIZED_DOUBLE
  CyberMrmtCreate%stoich_sw6_cbod = UNINITIALIZED_DOUBLE
!SWAT rxn 7
  allocate(CyberMrmtCreate%stoich_sw7_cbod)
  CyberMrmtCreate%stoich_sw7_cbod = UNINITIALIZED_DOUBLE
!SWAT rxn 8
  allocate(CyberMrmtCreate%stoich_sw8_orgn)
  allocate(CyberMrmtCreate%stoich_sw8_nh4)
  CyberMrmtCreate%stoich_sw8_orgn = UNINITIALIZED_DOUBLE
  CyberMrmtCreate%stoich_sw8_nh4 = UNINITIALIZED_DOUBLE
!
  CyberMrmtCreate%n_stoic = UNINITIALIZED_INTEGER
  CyberMrmtCreate%n_rates = UNINITIALIZED_INTEGER
!
  nullify(CyberMrmtCreate%alpha)
  nullify(CyberMrmtCreate%th_frac)
  nullify(CyberMrmtCreate%nrow)
  nullify(CyberMrmtCreate%ncol)
  nullify(CyberMrmtCreate%irow)
  nullify(CyberMrmtCreate%icol)
!  nullify(CyberMrmtCreate%stoich_row)
!  nullify(CyberMrmtCreate%rate_cs)

  nullify(CyberMrmtCreate%next)  
      
end function CyberMrmtCreate

! ************************************************************************** !

subroutine CyberRead(this,input,option)
  ! 
  ! Reads input deck for PNNL N reaction parameters (if any)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/01/15
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  
  implicit none
  
  class(reaction_sandbox_cyber_mrmt_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  PetscInt :: i
  character(len=MAXWORDLENGTH) :: word, internal_units, units
  character(len=MAXSTRINGLENGTH) :: error_string 
  character(len=MAXSTRINGLENGTH) :: string
  PetscReal :: storage_number

  error_string = 'CHEMISTRY,REACTION_SANDBOX,CYBERNETIC Multirate'
  
  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)   

    select case(trim(word))
! reaction stoichiometry
      case('F1')
        call InputReadDouble(input,option,this%f1)  
        call InputErrorMsg(input,option,'f1',error_string)
      case('F2')
        call InputReadDouble(input,option,this%f2)  
        call InputErrorMsg(input,option,'f2',error_string)
      case('F3')
        call InputReadDouble(input,option,this%f3)  
        call InputErrorMsg(input,option,'f3',error_string)
! reaction rate for no3 reduction
      case('K1','K_NO3-')
        call InputReadDouble(input,option,this%k1)  
        call InputErrorMsg(input,option,'k1',error_string)
        call InputReadAndConvertUnits(input,this%k1,'1/sec', &
                                      trim(error_string)//',k1',option)
! reaction rate for no2 reduction
      case('K2','K_NO2-')
        call InputReadDouble(input,option,this%k2)  
        call InputErrorMsg(input,option,'k2',error_string)
        call InputReadAndConvertUnits(input,this%k2,'1/sec', &
                                      trim(error_string)//',k2',option)
! reaction rate for o2 reduction
      case('K3','K_O2(aq)')
        call InputReadDouble(input,option,this%k3)  
        call InputErrorMsg(input,option,'k3',error_string)
        call InputReadAndConvertUnits(input,this%k3,'1/sec', &
                                      trim(error_string)//',k3',option)
! no3 half saturation constant for no3 reduction
      case('KA1','KA_NO3-')
        call InputReadDouble(input,option,this%Ka1)  
        call InputErrorMsg(input,option,'Ka1',error_string)
        call InputReadAndConvertUnits(input,this%Ka1,'M', &
                                      trim(error_string)//',Ka1',option)

! no2 half saturation constant for no2 reduction
      case('KA2','KA_NO2-')
        call InputReadDouble(input,option,this%Ka2)  
        call InputErrorMsg(input,option,'Ka2',error_string)
        call InputReadAndConvertUnits(input,this%Ka2,'M', &
                                      trim(error_string)//',Ka2',option)
! o2 half saturation constant for o2 reduction
      case('KA3','KA_O2(aq)')
        call InputReadDouble(input,option,this%Ka3)  
        call InputErrorMsg(input,option,'Ka3',error_string)
        call InputReadAndConvertUnits(input,this%Ka3,'M', &
                                      trim(error_string)//',Ka3',option)
! donor half saturation for no3 reduction
      case('KD1','KD_NO3-')
        call InputReadDouble(input,option,this%Kd1)  
        call InputErrorMsg(input,option,'Kd1',error_string)
        call InputReadAndConvertUnits(input,this%Kd1,'M', &
                                      trim(error_string)//',Kd1',option)
!o2 inhibition for no3 reduction, not used
      case('KI1','KI_NO3-')
        call InputReadDouble(input,option,this%Ki1)  
        call InputErrorMsg(input,option,'Ki1',error_string)
        call InputReadAndConvertUnits(input,this%Ki1,'M', &
                                      trim(error_string)//',Ki1',option)
! donor half saturation for no2 reduction
      case('KD2','KD_NO2-')
        call InputReadDouble(input,option,this%Kd2)  
        call InputErrorMsg(input,option,'Kd2',error_string)
        call InputReadAndConvertUnits(input,this%Kd2,'M', &
                                      trim(error_string)//',Kd2',option)
! not used
      case('KI2','KI_NO2-')
        call InputReadDouble(input,option,this%Ki2)  
        call InputErrorMsg(input,option,'Ki2',error_string)
        call InputReadAndConvertUnits(input,this%Ki2,'M', &
                                      trim(error_string)//',Ki2',option)
! donor half saturation for o2 reduction
      case('KD3','KD_O2(aq)')
        call InputReadDouble(input,option,this%Kd3)  
        call InputErrorMsg(input,option,'Kd3',error_string)
        call InputReadAndConvertUnits(input,this%Kd3,'M', &
                                      trim(error_string)//',Kd3',option)
! biomass degradation rate
      case('KDEG')
        call InputReadDouble(input,option,this%k_deg)  
        call InputErrorMsg(input,option,'kdeg',error_string)
        call InputReadAndConvertUnits(input,this%k_deg,'1/sec', &
                                      trim(error_string)//',kdeg',option)
!      case('F_ACT')
!        call InputReadDouble(input,option,this%f_act)  
!        call InputErrorMsg(input,option,'f_act',error_string)
      case('ACTIVATION_ENERGY')
        call InputReadDouble(input,option,this%activation_energy)  
        call InputErrorMsg(input,option,'activation energy',error_string)
        call InputReadAndConvertUnits(input,this%activation_energy,'J/mol', &
                              trim(error_string)//',activation energy',option)
      case('REFERENCE_TEMPERATURE')
        call InputReadDouble(input,option,this%reference_temperature)  
        call InputErrorMsg(input,option,'reference temperature [C]', &
                           error_string)
        this%reference_temperature = this%reference_temperature + 273.15d0
      case('CARBON_CONSUMPTION_SPECIES')
        call InputReadWord(input,option, &
                           this%carbon_consumption_species,PETSC_TRUE)
        call InputErrorMsg(input,option,'carbon consumption species', &
                           error_string)
      case('STORE_CONSUMPTION_PRODUCTION')
        this%store_cumulative_mass = PETSC_TRUE
! place holder for reach to hyporheic zone volume ratio, dynamically changing
      case('VOL_FRAC_R2ST')
        call InputReadDouble(input,option,this%vol_frac_r2st)  
        call InputErrorMsg(input,option,'vol_frac_r2st',error_string)
! place holder for hyporheic zone saturation
      case('TH_ST')
        call InputReadDouble(input,option,this%th_st)  
        call InputErrorMsg(input,option,'th_st',error_string)
      case('TH_FRAC')
  ! pore fraction of storage zone
        call UtilityReadArray(this%th_frac,NEG_ONE_INTEGER,string,input, &
                              option) 
      case('ALPHA')
        string = 'RATES for Multirate MassTransfer'
        call UtilityReadArray(this%alpha,NEG_ONE_INTEGER,string,input, &
                              option) 
        MRMT = size(this%alpha)
 
!        do i = 1, MRMT
!          call InputReadDouble(input,option,this%alpha(i))  
!          call InputErrorMsg(input,option,'alpha(i)',error_string)
         call InputReadPflotranString(input,option)
         if (InputError(input)) exit
         if (InputCheckExit(input,option)) exit
         call InputReadWord(input,option,word,PETSC_TRUE)
         call InputErrorMsg(input,option,'keyword',error_string)
         call StringToUpper(word)   
         call InputReadDouble(input,option,storage_number)  
          
        if (size(this%alpha) /= int(storage_number)) then
          write(word,*) size(this%alpha)
          write(string,*) int(storage_number)
          option%io_buffer = 'Number of kinetic mass transfer rates (' // &
            trim(adjustl(word)) // &
            ') does not match the number of storages (' // &
            trim(adjustl(string)) // ').'
          call printErrMsg(option)
        endif
         call InputReadAndConvertUnits(input,storage_number,'1/sec', &
                                      trim(error_string)//',alpha',option)
        do i = 1, MRMT
          this%alpha(i) = this%alpha(i)*storage_number/MRMT
        enddo
      case default
        call InputKeywordUnrecognized(word,error_string,option)
    end select
  enddo
  
end subroutine CyberRead

! ************************************************************************** !

subroutine CyberSetup(this,reaction,option)
  ! 
  ! Sets up the PNNL N reaction either with parameters either
  ! read from the input deck or hardwired.
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/01/15
  ! 

  use Reaction_Aux_module, only : reaction_type, GetPrimarySpeciesIDFromName
  use Reaction_Immobile_Aux_module, only : GetImmobileSpeciesIDFromName
  use Reaction_Mineral_Aux_module, only : GetKineticMineralIDFromName
  use Option_module

  implicit none
  
  class(reaction_sandbox_cyber_mrmt_type) :: this
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: irxn, i, ic, is, j
  PetscInt, dimension(MRMT_T) :: order
  
  PetscReal, parameter :: per_day_to_per_sec = 1.d0 / 24.d0 / 3600.d0
  character(len=4) :: form1
  character(len=16) :: name_
  data form1 / '(I )'/
  word = 'NH4+'
  this%nh4_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'O2(aq)'
  this%o2_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'NO3-'
  this%no3_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'NO2-'
  this%no2_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'OrgN'
  this%orgn_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'Algae'
  this%a_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'OrgP'
  this%orgp_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'DOP'
  this%dop_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'CBOD'
  this%cbod_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'N2(aq)'
  this%n2_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'CH2O(aq)'
  this%doc_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
!biomass not in the river
  !word = 'C5H7O2N(aq)'
  !this%biomass_id = &
  !  GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'CO2(aq)'
  this%co2_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)

! multirate
  order(:) = 0
  do i = 1, MRMT
    ic = 1
    is = i
    is = is/10
    do while (is >=1) 
      ic = ic + 1
      is = is/10
    enddo
    order(i) = ic
  enddo
  name_(:) = ''
  do i = 1, MRMT
    name_(1:3) = 'st('
    write(form1(3:3),'(I1)') order(i)
    write(name_(4:),form1) i  
    word = trim(name_) // ')_NO3-'
!    word = 'st(i)_NO3-'
    this%no3_st_id(i) = &
      GetImmobileSpeciesIDFromName(word,reaction%immobile,option)
  enddo

  do i = 1, MRMT
!    word = 'st(i)_NO2-'
    name_(1:3) = 'st('
    write(form1(3:3),'(I1)') order(i)
    write(name_(4:),form1) i
    word = trim(name_) // ')_NO2-'
    this%no2_st_id(i) = &
      GetImmobileSpeciesIDFromName(word,reaction%immobile,option)
  enddo
  do i = 1, MRMT
!    word = 'st(i)_O2(aq)'
    name_(1:3) = 'st('
    write(form1(3:3),'(I1)') order(i)
    write(name_(4:),form1) i
    word = trim(name_) // ')_O2'
    this%o2_st_id(i) = &
      GetImmobileSpeciesIDFromName(word,reaction%immobile,option)
  enddo
  do i = 1, MRMT
!    word = 'st(i)_NH4+'
    name_(1:3) = 'st('
    write(form1(3:3),'(I1)') order(i)
    write(name_(4:),form1) i
    word = trim(name_) // ')_NH4+'
    this%nh4_st_id(i) = &
      GetImmobileSpeciesIDFromName(word,reaction%immobile,option)
  enddo
  do i = 1, MRMT
!    word = 'st(i)_CO2(aq)'
    name_(1:3) = 'st('
    write(form1(3:3),'(I1)') order(i)
    write(name_(4:),form1) i
    word = trim(name_) // ')_CO2'
    this%co2_st_id(i) = &
      GetImmobileSpeciesIDFromName(word,reaction%immobile,option)
  enddo
  do i = 1, MRMT
!    word = 'st(i)_N2(aq)'
    name_(1:3) = 'st('
    write(form1(3:3),'(I1)') order(i)
    write(name_(4:),form1) i
    word = trim(name_) // ')_N2'
    this%n2_st_id(i) = &
      GetImmobileSpeciesIDFromName(word,reaction%immobile,option)
  enddo

  do i = 1, MRMT
!    word = 'st(i)_C5H7O2N(aq)'
    name_(1:3) = 'st('
    write(form1(3:3),'(I1)') order(i)
    write(name_(4:),form1) i
    word = trim(name_) // ')_C5H7O2N'
    this%biomass_st_id(i) = &
      GetImmobileSpeciesIDFromName(word,reaction%immobile,option)
  enddo
  do i = 1, MRMT
!    word = 'st(i)_CH2O(aq)'
    name_(1:3) = 'st('
    write(form1(3:3),'(I1)') order(i)
    write(name_(4:),form1) i
    word = trim(name_) // ')_CH2O'
    this%doc_st_id(i) = &
      GetImmobileSpeciesIDFromName(word,reaction%immobile,option)
  enddo

!
  
  this%offset_auxiliary = reaction%nauxiliary
    ! source + swat reaction parameters (21 = stoich + rate constants)
!Will be updated during simulation
!SWAT rxn 1
  this%stoich_sw1_o2 = -1.d0 
  this%stoich_sw1_a = -1.d0
  this%stoich_sw1_orgn = 1.d0
  this%stoich_sw1_orgp = 1.d0
!SWAT rxn 2
  this%stoich_sw2_dop = -1.d0
  this%stoich_sw2_no3 = -1.d0
  this%stoich_sw2_nh4 = -1.d0
  this%stoich_sw2_a = 1.d0
  this%stoich_sw2_o2 = 1.d0
!SWAT rxn 3
  this%stoich_sw3_o2 = -1.d0
  this%stoich_sw3_nh4 = -1.d0
  this%stoich_sw3_no2 = 1.d0
!SWAT rxn 4
  this%stoich_sw4_o2 = -1.d0
  this%stoich_sw4_no2 = -1.d0
  this%stoich_sw4_no3 = 1.d0
!SWAT rxn 5
  this%stoich_sw5_orgp = -1.d0
  this%stoich_sw5_dop = 1.d0
!SWAT rxn 6
  this%stoich_sw6_o2 = -1.d0
  this%stoich_sw6_cbod = -1.d0
!SWAT rxn 7
  this%stoich_sw7_cbod = -1.d0
!SWAT rxn 8
  this%stoich_sw8_orgn = -1.d0
  this%stoich_sw8_nh4 = 1.d0
    !
!  reaction%nauxiliary = reaction%nauxiliary + 23 + 4
  reaction%nauxiliary = reaction%nauxiliary + 23 + mrmt*2 

  ! constants based on Hyun's writeup on 10/26/16 entitled "Mini-cybernetic
  ! model of batch denitrification process"
  if (Uninitialized(this%f1)) this%f1 = 0.65d0
  if (Uninitialized(this%f2)) this%f2 = 0.99d0
  if (Uninitialized(this%f3)) this%f3 = 0.2167d0
  if (Uninitialized(this%k1)) this%k1 = 28.26d0 * per_day_to_per_sec
  if (Uninitialized(this%k2)) this%k2 = 23.28d0 * per_day_to_per_sec
  if (Uninitialized(this%k3)) this%k3 = 84.78d0 * per_day_to_per_sec
  if (Uninitialized(this%Kd1)) this%Kd1 = 0.25d-3
  if (Uninitialized(this%Kd2)) this%Kd2 = 0.25d-3
  if (Uninitialized(this%Ki1)) this%Ki1 = 3.13d-2
  if (Uninitialized(this%Ki2)) this%Ki2 = 3.13d-2
  if (Uninitialized(this%Kd3)) this%Kd3 = 0.25d-3
  if (Uninitialized(this%Ka1)) this%Ka1 = 0.001d-3
  if (Uninitialized(this%Ka2)) this%Ka2 = 0.004d-3
  if (Uninitialized(this%Ka3)) this%Ka3 = 0.001d-3
!  if (Uninitialized(this%f_act)) this%f_act = 0.126d0
  if (Uninitialized(this%k_deg)) this%k_deg = 0.242d0 * per_day_to_per_sec
  do i = 1, MRMT
    if (Uninitialized(this%alpha(i))) this%alpha(i) = 0.001d0
  enddo

! uncomment these to zero out reactions
!  this%k1 = 0.d0
!  this%k2 = 0.d0
!  this%k3 = 0.d0
!  this%Kd1 = 0.d0
!  this%Ka1 = 0.d0
!  this%f_act = 1.d40
!  this%k_deg = 0.d0

  ! NOTE: CO2 stiochiometries below factor in on carbon consumption.
  ! Take care to ensure that changes do not adversely affect consumption.
  this%stoich_1_doc = -1.d0
  !this%stoich_1_nh4 = -0.2d0*(1.d0-this%f1)
  this%stoich_1_nh4 = 0.d0
  this%stoich_1_no3 = -2.d0*this%f1
  this%stoich_1_no2 = 2.d0*this%f1
  this%stoich_1_co2 = this%f1
  this%stoich_1_biomass = 0.2d0*(1.d0-this%f1)

  this%stoich_2_doc = -1.d0
  !this%stoich_2_nh4 = -0.2d0*(1.d0-this%f2)
  this%stoich_2_nh4 = 0.d0
  this%stoich_2_no2 = -4.d0/3.d0*this%f2
  this%stoich_2_n2 = 2.d0/3.d0*this%f2
  this%stoich_2_co2 = this%f2
  this%stoich_2_biomass = 0.2d0*(1.d0-this%f2)

  this%stoich_3_doc = -1.d0
  !this%stoich_3_nh4 = -0.2d0*(1.d0-this%f3)
  this%stoich_3_nh4 = 0.d0
  this%stoich_3_o2 = -1.d0*this%f3
  this%stoich_3_co2 = this%f3
  this%stoich_3_biomass = 0.2d0*(1.d0-this%f3)


!multirate 
  do i = 1, MRMT
    this%stoich_doc(i) = -1.d0
! mol/aqueous of stage volume
   ! this%vol_frac_mrmt(i)%ptr = this%vol_frac_r2st
    this%stoich_mrmt_doc(i) = this%vol_frac_r2st/this%th_frac(i)/this%th_st
    this%stoich_nh4(i) = -1.d0
    this%stoich_mrmt_nh4(i) = this%vol_frac_r2st/this%th_frac(i)/this%th_st
    this%stoich_no3(i) = -1.d0
    this%stoich_mrmt_no3(i) = this%vol_frac_r2st/this%th_frac(i)/this%th_st
    this%stoich_no2(i) = -1.d0
    this%stoich_mrmt_no2(i) = this%vol_frac_r2st/this%th_frac(i)/this%th_st
    this%stoich_co2(i) = -1.d0
    this%stoich_mrmt_co2(i) = this%vol_frac_r2st/this%th_frac(i)/this%th_st
    this%stoich_o2(i) = -1.d0
    this%stoich_mrmt_o2(i) = this%vol_frac_r2st/this%th_frac(i)/this%th_st
    this%stoich_n2(i) = -1.d0
    this%stoich_mrmt_n2(i) = this%vol_frac_r2st/this%th_frac(i)/this%th_st
! assuming biomass is not exchanged
!    this%stoich_biomass(i) = -1.d0
!    this%stoich_mrmt_biomass(i) = this%vol_frac_r2st/this%th_frac(i)/this%th_st
  enddo

!
  this%nrxn_swat = 8
! only nrxn_mt species participated in mass transfer
  this%nrxn_mt = 7

  if (this%k3 > 1.d-40) then
  this%nrxn = 3 * MRMT + this%nrxn_swat
  else 
    this%nrxn = 2 * MRMT + this%nrxn_swat
  endif

  this%nrxn = MRMT*this%nrxn_mt + this%nrxn

  allocate(this%nrow(this%nrxn))
  this%nrow = UNINITIALIZED_INTEGER
  allocate(this%ncol(this%nrxn))
  this%ncol = UNINITIALIZED_INTEGER
! 
  allocate(this%irow(6,this%nrxn))
  this%irow = UNINITIALIZED_INTEGER
  allocate(this%icol(4,this%nrxn))
  this%icol = UNINITIALIZED_INTEGER
  allocate(this%stoich_row(6,this%nrxn))
do j=1,this%nrxn
do i=1,6
  allocate(this%stoich_row(i,j)%ptr)
  this%stoich_row(i,j)%ptr = UNINITIALIZED_DOUBLE
enddo
enddo

  irxn = 0

!SWAT rxn 1, each rxn in a column, so species in row
!O2 + A -> Orgn + orgp

  irxn = irxn + 1
  this%nrow(irxn) = 4
  this%irow(1,irxn) = this%o2_id
  this%irow(2,irxn) = this%a_id
  this%irow(3,irxn) = this%orgn_id
  this%irow(4,irxn) = this%orgp_id
  this%stoich_row(1,irxn)%ptr => this%stoich_sw1_o2
  this%stoich_row(2,irxn)%ptr => this%stoich_sw1_a
  this%stoich_row(3,irxn)%ptr => this%stoich_sw1_orgn
  this%stoich_row(4,irxn)%ptr => this%stoich_sw1_orgp
! dependent variables
  this%ncol(irxn) = 4
  this%icol(1,irxn) = this%a_id
  this%icol(2,irxn) = this%nh4_id
  this%icol(3,irxn) = this%no3_id
  this%icol(4,irxn) = this%dop_id


!SWAT rxn 2
! DOP + NH4+ -> A + O2 
  irxn = irxn + 1
  this%nrow(irxn) = 5
  this%irow(1,irxn) = this%dop_id
  this%irow(2,irxn) = this%no3_id
  this%irow(3,irxn) = this%nh4_id
  this%irow(4,irxn) = this%a_id
  this%irow(5,irxn) = this%o2_id
  this%stoich_row(1,irxn)%ptr => this%stoich_sw2_dop
  this%stoich_row(2,irxn)%ptr => this%stoich_sw2_no3
  this%stoich_row(3,irxn)%ptr => this%stoich_sw2_nh4
  this%stoich_row(4,irxn)%ptr => this%stoich_sw2_a
  this%stoich_row(5,irxn)%ptr => this%stoich_sw2_o2
!
  this%ncol(irxn) = 4
  this%icol(1,irxn) = this%a_id
  this%icol(2,irxn) = this%nh4_id
  this%icol(3,irxn) = this%no3_id
  this%icol(4,irxn) = this%dop_id
  
!SWAT rxn 3
! O2 + NH4- -> NO2-
  irxn = irxn + 1
  this%nrow(irxn) = 3
  this%irow(1,irxn) = this%o2_id
  this%irow(2,irxn) = this%nh4_id
  this%irow(3,irxn) = this%no2_id
  this%stoich_row(1,irxn)%ptr => this%stoich_sw3_o2
  this%stoich_row(2,irxn)%ptr => this%stoich_sw3_nh4
  this%stoich_row(3,irxn)%ptr => this%stoich_sw3_no2
!
  this%ncol(irxn) = 1
  this%icol(1,irxn) = this%nh4_id
  
!SWAT rxn 4
! O2 + NO2- => NO3-
  irxn = irxn + 1
  this%nrow(irxn) = 3
  this%irow(1,irxn) = this%o2_id
  this%irow(2,irxn) = this%no2_id
  this%irow(3,irxn) = this%no3_id
  this%stoich_row(1,irxn)%ptr => this%stoich_sw4_o2
  this%stoich_row(2,irxn)%ptr => this%stoich_sw4_no2
  this%stoich_row(3,irxn)%ptr => this%stoich_sw4_no3
!
  this%ncol(irxn) = 1
  this%icol(1,irxn) = this%no2_id
  
!SWAT rxn 5
! OrgP -> DOP
  irxn = irxn + 1
  this%nrow(irxn) = 2
  this%irow(1,irxn) = this%orgp_id
  this%irow(2,irxn) = this%dop_id
  this%stoich_row(1,irxn)%ptr => this%stoich_sw5_orgp
  this%stoich_row(2,irxn)%ptr => this%stoich_sw5_dop
!
  this%ncol(irxn) = 1
  this%icol(1,irxn) = this%orgp_id
  
!SWAT rxn 6
! O2 + CBOD ->
  irxn = irxn + 1
  this%nrow(irxn) = 2
  this%irow(1,irxn) = this%o2_id
  this%irow(2,irxn) = this%cbod_id
  this%stoich_row(1,irxn)%ptr => this%stoich_sw6_o2
  this%stoich_row(2,irxn)%ptr => this%stoich_sw6_cbod
!
  this%ncol(irxn) = 1
  this%icol(1,irxn) = this%cbod_id

!SWAT rxn 7
! CBOD ->
  irxn = irxn + 1
  this%nrow(irxn) = 1
  this%irow(1,irxn) = this%cbod_id
  this%stoich_row(1,irxn)%ptr => this%stoich_sw7_cbod
!
  this%ncol(irxn) = 1
  this%icol(1,irxn) = this%cbod_id
!SWAT rxn 8
! Orgn -> NH4+
  irxn = irxn + 1
  this%nrow(irxn) = 2
  this%irow(1,irxn) = this%orgn_id
  this%irow(2,irxn) = this%nh4_id
  this%stoich_row(1,irxn)%ptr => this%stoich_sw8_orgn
  this%stoich_row(2,irxn)%ptr => this%stoich_sw8_nh4
!
  this%ncol(irxn) = 1
  this%icol(1,irxn) = this%orgn_id
!
  this%n_stoic = 7
  this%n_rates = 16

  ! NO3- -> NO2- in HZ
!  irxn = 0
  do i = 1, MRMT
    irxn = irxn + 1
    this%nrow(irxn) = 6
! species id in the reaction
    this%irow(1,irxn) = this%doc_st_id(i)
    this%irow(2,irxn) = this%nh4_st_id(i)
    this%irow(3,irxn) = this%no3_st_id(i)
    this%irow(4,irxn) = this%no2_st_id(i)
    this%irow(5,irxn) = this%co2_st_id(i)
    this%irow(6,irxn) = this%biomass_st_id(i)
! species stoich in the reaction
    this%stoich_row(1,irxn)%ptr = this%stoich_1_doc
    this%stoich_row(2,irxn)%ptr = this%stoich_1_nh4
    this%stoich_row(3,irxn)%ptr = this%stoich_1_no3
    this%stoich_row(4,irxn)%ptr = this%stoich_1_no2
    this%stoich_row(5,irxn)%ptr = this%stoich_1_co2
    this%stoich_row(6,irxn)%ptr = this%stoich_1_biomass

    this%ncol(irxn) = 4
! species need jacobian evaluation, because of the weight of the sum of rates,
! including the 4
! species
    this%icol(1,irxn) = this%doc_st_id(i)
    this%icol(2,irxn) = this%no3_st_id(i)
    this%icol(3,irxn) = this%no2_st_id(i)
    this%icol(4,irxn) = this%o2_st_id(i)
  end do
  ! NO2- -> N2
  do i = 1, MRMT
    irxn = irxn + 1
    this%nrow(irxn) = 6
    this%irow(1,irxn) = this%doc_st_id(i)
    this%irow(2,irxn) = this%nh4_st_id(i)
    this%irow(3,irxn) = this%no2_st_id(i)
    this%irow(4,irxn) = this%n2_st_id(i)
    this%irow(5,irxn) = this%co2_st_id(i)
    this%irow(6,irxn) = this%biomass_st_id(i)
    this%stoich_row(1,irxn)%ptr = this%stoich_2_doc
    this%stoich_row(2,irxn)%ptr = this%stoich_2_nh4
    this%stoich_row(3,irxn)%ptr = this%stoich_2_no2
    this%stoich_row(4,irxn)%ptr = this%stoich_2_n2
    this%stoich_row(5,irxn)%ptr = this%stoich_2_co2
    this%stoich_row(6,irxn)%ptr = this%stoich_2_biomass
    this%ncol(irxn) = 4
    this%icol(1,irxn) = this%doc_st_id(i)
    this%icol(2,irxn) = this%no3_st_id(i)
    this%icol(3,irxn) = this%no2_st_id(i)
    this%icol(4,irxn) = this%o2_st_id(i)
  enddo
  if (this%k3 > 1.d-40) then
  ! O2 -> H2O + CO2
   do i = 1, MRMT
    irxn = irxn + 1
    this%nrow(irxn) = 5
    this%irow(1,irxn) = this%doc_st_id(i)
    this%irow(2,irxn) = this%nh4_st_id(i)
    this%irow(3,irxn) = this%o2_st_id(i)
    this%irow(4,irxn) = this%co2_st_id(i)
    this%irow(5,irxn) = this%biomass_st_id(i)
    this%stoich_row(1,irxn)%ptr = this%stoich_3_doc
    this%stoich_row(2,irxn)%ptr = this%stoich_3_nh4
    this%stoich_row(3,irxn)%ptr = this%stoich_3_o2
    this%stoich_row(4,irxn)%ptr = this%stoich_3_co2
    this%stoich_row(5,irxn)%ptr = this%stoich_3_biomass
    this%ncol(irxn) = 4
    this%icol(1,irxn) = this%doc_st_id(i)
    this%icol(2,irxn) = this%no3_st_id(i)
    this%icol(3,irxn) = this%no2_st_id(i)
    this%icol(4,irxn) = this%o2_st_id(i)
   enddo
  endif

!multirate mass transfer
! multirate mass transfer
! doc -> doc_st
  allocate( this%vol_frac_mrmt(MRMT))
  do i=1,MRMT
    allocate(this%vol_frac_mrmt(i)%ptr)
    this%vol_frac_mrmt(i)%ptr = UNINITIALIZED_DOUBLE
  enddo
! cbod -> cbod_st
  do i = 1, MRMT
    irxn = irxn + 1
    this%nrow(irxn) = 2
    this%irow(1,irxn) = this%doc_id
    !this%irow(1,irxn) = this%cbod_id
    this%irow(2,irxn) = -this%doc_st_id(i)
    this%stoich_row(1,irxn)%ptr = this%stoich_doc(i)
    this%stoich_row(2,irxn)%ptr = this%stoich_mrmt_doc(i)

    this%ncol(irxn) = 2
    this%icol(1,irxn) = this%doc_id
!    this%icol(1,irxn) = this%cbod_id
    this%icol(2,irxn) = -this%doc_st_id(i)
  enddo

! no3 -> no3_st
  do i = 1, MRMT
    irxn = irxn + 1
    this%nrow(irxn) = 2
    this%irow(1,irxn) = this%no3_id
    this%irow(2,irxn) = -this%no3_st_id(i)
    this%stoich_row(1,irxn)%ptr = this%stoich_no3(i)
    this%stoich_row(2,irxn)%ptr = this%stoich_mrmt_no3(i)

    this%ncol(irxn) = 2
    this%icol(1,irxn) = this%no3_id
    this%icol(2,irxn) = -this%no3_st_id(i)
  enddo
! no2 -> no2_st
  do i = 1, MRMT
    irxn = irxn + 1
    this%nrow(irxn) = 2
    this%irow(1,irxn) = this%no2_id
    this%irow(2,irxn) = -this%no2_st_id(i)
    this%stoich_row(1,irxn)%ptr = this%stoich_no2(i)
    this%stoich_row(2,irxn)%ptr = this%stoich_mrmt_no2(i)

    this%ncol(irxn) = 2
    this%icol(1,irxn) = this%no2_id
    this%icol(2,irxn) = -this%no2_st_id(i)
  enddo
! nh4 -> nh4_st
  do i = 1, MRMT
    irxn = irxn + 1
    this%nrow(irxn) = 2
    this%irow(1,irxn) = this%nh4_id
    this%irow(2,irxn) = -this%nh4_st_id(i)
    this%stoich_row(1,irxn)%ptr = this%stoich_nh4(i)
    this%stoich_row(2,irxn)%ptr = this%stoich_mrmt_nh4(i)

    this%ncol(irxn) = 2
    this%icol(1,irxn) = this%nh4_id
    this%icol(2,irxn) = -this%nh4_st_id(i)
  enddo
! o2 -> o2_st
  do i = 1, MRMT
    irxn = irxn + 1
    this%nrow(irxn) = 2
    this%irow(1,irxn) = this%o2_id
    this%irow(2,irxn) = -this%o2_st_id(i)
    this%stoich_row(1,irxn)%ptr = this%stoich_o2(i)
    this%stoich_row(2,irxn)%ptr = this%stoich_mrmt_o2(i)

    this%ncol(irxn) = 2
    this%icol(1,irxn) = this%o2_id
    this%icol(2,irxn) = -this%o2_st_id(i)
  enddo
! n2 -> n2_st
  do i = 1, MRMT
    irxn = irxn + 1
    this%nrow(irxn) = 2
    this%irow(1,irxn) = this%n2_id
    this%irow(2,irxn) = -this%n2_st_id(i)
    this%stoich_row(1,irxn)%ptr = this%stoich_n2(i)
    this%stoich_row(2,irxn)%ptr = this%stoich_mrmt_n2(i)

    this%ncol(irxn) = 2
    this%icol(1,irxn) = this%n2_id
    this%icol(2,irxn) = -this%n2_st_id(i)
  enddo
! co2 -> co2_st
  do i = 1, MRMT
    irxn = irxn + 1
    this%nrow(irxn) = 2
    this%irow(1,irxn) = this%co2_id
    this%irow(2,irxn) = -this%co2_st_id(i)
    this%stoich_row(1,irxn)%ptr = this%stoich_co2(i)
    this%stoich_row(2,irxn)%ptr = this%stoich_mrmt_co2(i)

    this%ncol(irxn) = 2
    this%icol(1,irxn) = this%co2_id
    this%icol(2,irxn) = -this%co2_st_id(i)
  enddo
! biomass -> biomass_st
!  do i = 1, MRMT
!    irxn = irxn + 1
!    this%nrow(irxn) = 2
!    this%irow(1,irxn) = this%biomass_id
!    this%irow(2,irxn) = -this%biomass_st_id(i)
!    this%stoich_row(1,irxn)%ptr = this%stoich_biomass(i)
!    this%stoich_row(2,irxn)%ptr = this%stoich_mrmt_biomass(i)

!    this%ncol(irxn) = 2
!    this%icol(1,irxn) = this%biomass_id
!    this%icol(2,irxn) = -this%biomass_st_id(i)
!  enddo
  
end subroutine CyberSetup

! ************************************************************************** !

subroutine CyberAuxiliaryPlotVariables(this,list,reaction,option)
  ! 
  ! Adds plot variables to output list
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/21/16
  
  use Option_module
  use Reaction_Aux_module
  use Output_Aux_module
  use Variables_module, only : REACTION_AUXILIARY

  implicit none

  class(reaction_sandbox_cyber_mrmt_type) :: this
  type(output_variable_list_type), pointer :: list
  type(option_type) :: option
  type(reaction_type) :: reaction
  
  character(len=MAXWORDLENGTH) :: names(6)
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: units
  PetscInt :: indices(6)
  PetscInt :: i

  names(1) = 'DOC'
  names(2) = 'NH4'
  names(3) = 'O2'
  names(4) = 'NO3'
  names(5) = 'NO2'
  names(6) = 'CO2'
  indices(1) = DOC_MASS_STORAGE_INDEX
  indices(2) = NH4_MASS_STORAGE_INDEX
  indices(3) = O2_MASS_STORAGE_INDEX
  indices(4) = NO3_MASS_STORAGE_INDEX
  indices(5) = NO2_MASS_STORAGE_INDEX
  indices(6) = CO2_MASS_STORAGE_INDEX
  if (this%store_cumulative_mass) then
    do i = 1, 6
      word = trim(names(i)) // ' Rate'
      units = 'mol/m^3-sec'
      call OutputVariableAddToList(list,word,OUTPUT_RATE,units, &
                                   REACTION_AUXILIARY, &
                                   this%offset_auxiliary+indices(i))
    enddo
    do i = 1, 6
      word = trim(names(i)) // ' Cum. Mass'
      units = 'mol/m^3'
      call OutputVariableAddToList(list,word,OUTPUT_GENERIC,units, &
                                   REACTION_AUXILIARY, &
                                   this%offset_auxiliary+6+indices(i))
    enddo
  endif

end subroutine CyberAuxiliaryPlotVariables

! ************************************************************************** !

subroutine CyberReact(this,Residual,Jacobian,compute_derivative, &
                         rt_auxvar,global_auxvar,material_auxvar,reaction, &
                         option)
  ! 
  ! Evaluates reaction storing residual and/or Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/01/15
  ! 

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_class
  
  implicit none
  
  class(reaction_sandbox_cyber_mrmt_type) :: this  
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscBool :: compute_derivative
  ! the following arrays must be declared after reaction
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1
  PetscReal :: L_water
  PetscReal :: kg_water
  
  PetscInt :: i, j, irxn, non_mrmt, id_, jd_, i_ord, j_ord

  PetscReal :: Co2, Cno3, Cno2, Cn2, Cdoc, Cco2, Cnh4, X
  PetscReal :: r1docmonod(MRMT), r1docmonod_denom(MRMT)
  PetscReal :: r2docmonod(MRMT), r2docmonod_denom(MRMT)
  PetscReal :: r3docmonod(MRMT), r3docmonod_denom(MRMT)
  PetscReal :: r1no3monod(MRMT), r1no3monod_denom(MRMT)
  PetscReal :: r2no2monod(MRMT), r2no2monod_denom(MRMT)
  PetscReal :: r3o2monod(MRMT), r3o2monod_denom(MRMT)
  PetscReal :: r1kin(MRMT), r2kin(MRMT), r3kin(MRMT)
  PetscReal :: sumkin(MRMT), sumkinsq(MRMT)
  PetscReal :: u1(MRMT), u2(MRMT), u3(MRMT)
  PetscReal :: dr1kin_ddoc, dr1kin_dno3
  PetscReal :: dr2kin_ddoc, dr2kin_dno2
  PetscReal :: dr3kin_ddoc, dr3kin_do2
  PetscReal :: du_denom_dr
  PetscReal :: du1_dr1kin, du1_dr2kin, du1_dr3kin
  PetscReal :: du2_dr1kin, du2_dr2kin, du2_dr3kin
  PetscReal :: du3_dr1kin, du3_dr2kin, du3_dr3kin
  PetscReal :: dr1_ddoc, dr1_dno3, dr1_dno2, dr1_do2
  PetscReal :: dr2_ddoc, dr2_dno3, dr2_dno2, dr2_do2
  PetscReal :: dr3_ddoc, dr3_dno3, dr3_dno2, dr3_do2
  PetscReal :: du1_ddoc, du1_dno3, du1_dno2, du1_do2
  PetscReal :: du2_ddoc, du2_dno3, du2_dno2, du2_do2
  PetscReal :: du3_ddoc, du3_dno3, du3_dno2, du3_do2
  PetscReal :: molality_to_molarity
  PetscReal :: temperature_scaling_factor
  PetscReal :: k1_scaled, k2_scaled, k3_scaled, k_deg_scaled
  PetscReal :: volume, rate_scale
  PetscReal :: Ca, Ccbod, Corgp, Corgn, Cdop, sum_rates, cdt
  PetscReal :: f1,p_n
  PetscReal :: Residual_tmp(reaction%ncomp)
  PetscReal, dimension(MRMT_T) :: rkin_no3_mrmt, rkin_no2_mrmt, rkin_nh4_mrmt, &
            rkin_o2_mrmt, rkin_X_mrmt, rkin_co2_mrmt, rkin_n2_mrmt, &
            rkin_doc_mrmt 
!
  PetscReal, dimension(MRMT_T) :: Co2_st, Cno3_st, Cno2_st, Cn2_st, Cdoc_st, &
             Cco2_st, Cnh4_st, X_st

  PetscReal :: rate(8 + 3*MRMT+8*MRMT), derivative_col(6,3*MRMT+8*MRMT)
  
  volume = material_auxvar%volume
  L_water = material_auxvar%porosity*global_auxvar%sat(iphase)* &
            volume*1.d3 ! m^3 -> L
  kg_water = material_auxvar%porosity*global_auxvar%sat(iphase)* &
             global_auxvar%den_kg(iphase)*volume
! molarlity - mol/kg.  molarity - mol/volume => convert to mol/L
  molality_to_molarity = global_auxvar%den_kg(iphase)*1.d-3
    
  if (reaction%act_coef_update_frequency /= ACT_COEF_FREQUENCY_OFF) then
    option%io_buffer = 'Activity coefficients not currently supported in &
      &CyberReact().'
    call printErrMsg(option)
  endif
  
  temperature_scaling_factor = 1.d0
  if (Initialized(this%activation_energy)) then
    temperature_scaling_factor = &
      exp(this%activation_energy/IDEAL_GAS_CONSTANT* &
          (1.d0/this%reference_temperature-1.d0/(global_auxvar%temp+273.15d0)))
  endif
  
  ! concentrations are molarities [M], convert to mol/L
  Co2 = rt_auxvar%pri_molal(this%o2_id)* &
        rt_auxvar%pri_act_coef(this%o2_id)*molality_to_molarity
  Cnh4 = rt_auxvar%pri_molal(this%nh4_id)* &
         rt_auxvar%pri_act_coef(this%nh4_id)*molality_to_molarity
  Cno3 = rt_auxvar%pri_molal(this%no3_id)* &
         rt_auxvar%pri_act_coef(this%no3_id)*molality_to_molarity
  Cno2 = rt_auxvar%pri_molal(this%no2_id)* &
         rt_auxvar%pri_act_coef(this%no2_id)*molality_to_molarity
  Cn2 = rt_auxvar%pri_molal(this%n2_id)* &
        rt_auxvar%pri_act_coef(this%n2_id)*molality_to_molarity
  Cdoc = rt_auxvar%pri_molal(this%doc_id)* &
         rt_auxvar%pri_act_coef(this%doc_id)*molality_to_molarity
!Fang - fixed doc concentration in the river
!  Cdoc = 0.11e-3 !average river chem
  Ccbod = rt_auxvar%pri_molal(this%cbod_id)* &
         rt_auxvar%pri_act_coef(this%cbod_id)*molality_to_molarity
!assuming doc is cbod
! in mg/L
  Ca = rt_auxvar%pri_molal(this%a_id)
  Corgp = rt_auxvar%pri_molal(this%orgp_id)* &
         rt_auxvar%pri_act_coef(this%orgp_id)*molality_to_molarity
  Corgn = rt_auxvar%pri_molal(this%orgn_id)* &
         rt_auxvar%pri_act_coef(this%orgn_id)*molality_to_molarity
  Cco2 = rt_auxvar%pri_molal(this%co2_id)* &
         rt_auxvar%pri_act_coef(this%co2_id)*molality_to_molarity
  Cdop = rt_auxvar%pri_molal(this%dop_id)* &
         rt_auxvar%pri_act_coef(this%dop_id)*molality_to_molarity

!  X = rt_auxvar%pri_molal(this%biomass_id)* &
!      rt_auxvar%pri_act_coef(this%biomass_id)*molality_to_molarity
  
!
  i = this%offset_auxiliary+1 !stoichiometric coefficients
!SWAT rxn 1
  this%stoich_sw1_o2 = sign(rt_auxvar%auxiliary_data(i),this%stoich_sw1_o2) !alpha4
  this%stoich_sw1_orgn = sign(rt_auxvar%auxiliary_data(i+1),this%stoich_sw1_orgn) !alpha1
  this%stoich_sw1_orgp = sign(rt_auxvar%auxiliary_data(i+2),this%stoich_sw1_orgp) !alpha2
!SWAT rxn 2
  this%stoich_sw2_dop = sign(rt_auxvar%auxiliary_data(i+2),this%stoich_sw2_dop) !alpha2
  p_n = rt_auxvar%auxiliary_data(i+3)
  f1 = p_n * Cnh4 / (p_n * Cnh4 + (1. - p_n) * Cno3 + 1.e-6/14.d0*1.d-3)
  p_n = rt_auxvar%auxiliary_data(i+1)
  this%stoich_sw2_no3 = sign((1.d0 - f1)*p_n, this%stoich_sw2_no3)  !1-falpha1
  this%stoich_sw2_nh4 = sign(f1*p_n,this%stoich_sw2_nh4) !falpha1
  this%stoich_sw2_o2 = sign(rt_auxvar%auxiliary_data(i+4),this%stoich_sw2_o2) !alpha3
!SWAT rxn 3
  this%stoich_sw3_o2 = sign(rt_auxvar%auxiliary_data(i+5),this%stoich_sw3_o2) !alpha5
!SWAT rxn 4
  this%stoich_sw4_o2 = sign(rt_auxvar%auxiliary_data(i+6),this%stoich_sw4_o2) !alpha6

! rate constants
  rate(1) = rt_auxvar%auxiliary_data(i+7) * Ca !rho
  rate(2) = rt_auxvar%auxiliary_data(i+8) * Ca !mu 
  rate(3) = rt_auxvar%auxiliary_data(i+9) * Cnh4 !beta1
  rate(4) = rt_auxvar%auxiliary_data(i+10) * Cno2 !beta2
  rate(5) = rt_auxvar%auxiliary_data(i+11) * Corgp !beta4
  rate(6) = rt_auxvar%auxiliary_data(i+12) * Ccbod !k1
  rate(7) = rt_auxvar%auxiliary_data(i+13) * Ccbod !k3
  rate(8) = rt_auxvar%auxiliary_data(i+21) * Corgn !beta3

!update multirate stoichiometry
  if(this%k3 > 1.d-40) then
    irxn = 8+3*mrmt 
  else
    irxn = 8+2*mrmt
  endif
  j = this%offset_auxiliary + 23
  do i = 1, MRMT
    this%vol_frac_mrmt(i)%ptr = rt_auxvar%auxiliary_data(i+j+mrmt) !
    this%alpha(i) = rt_auxvar%auxiliary_data(i+j)
!print *,'i-',i,this%alpha(i),this%vol_frac_mrmt(i)%ptr
  enddo
  do j=1,7
  do i = 1, MRMT
    irxn = irxn + 1
    this%stoich_row(2,irxn)%ptr = this%vol_frac_mrmt(i)%ptr !update multirate reaction stoichometry
  enddo
  enddo
! 
!
! multirate, storage zone, immobile aqueous units: mol/L aqueous
  do i = 1, MRMT
    Co2_st(i) = rt_auxvar%immobile(this%o2_st_id(i))
    Cnh4_st(i) = rt_auxvar%immobile(this%nh4_st_id(i))
    Cno3_st(i) = rt_auxvar%immobile(this%no3_st_id(i))
    Cno2_st(i) = rt_auxvar%immobile(this%no2_st_id(i))
    Cco2_st(i) = rt_auxvar%immobile(this%co2_st_id(i))
    Cdoc_st(i) = rt_auxvar%immobile(this%doc_st_id(i))
    Cn2_st(i) = rt_auxvar%immobile(this%n2_st_id(i))
    X_st(i) = rt_auxvar%immobile(this%biomass_st_id(i))
  enddo
!if(option%time <= 3600*365*24.0) then
!  k1_scaled = 0.
!  k2_scaled = 0.
!  k3_scaled = 0.
!  k_deg_scaled = 0.
!else
  k1_scaled = this%k1 * temperature_scaling_factor
  k2_scaled = this%k2 * temperature_scaling_factor
  k3_scaled = this%k3 * temperature_scaling_factor
  k_deg_scaled = this%k_deg * temperature_scaling_factor
!endif
  ! NO3- -> NO2-
  do i = 1, MRMT
    r1docmonod_denom(i) = Cdoc_st(i)+this%Kd1
!Fang test
    r1docmonod(i) = Cdoc_st(i)/r1docmonod_denom(i)
!    r1docmonod(i) = 1.d0
    r1no3monod_denom(i) = Cno3_st(i)+this%Ka1
  r1no3monod(i) = Cno3_st(i)/r1no3monod_denom(i)
  r1kin(i) = k1_scaled*r1docmonod(i)*r1no3monod(i) & 
               *this%Ki1/(Co2_st(i)+this%Ki1)
!&

  ! NO2- -> N2
  r2docmonod_denom(i) = Cdoc_st(i)+this%Kd2
  r2docmonod(i) = Cdoc_st(i)/r2docmonod_denom(i)
!Fang test, no doc dependency
  !r2docmonod(i) = 1.d0
  r2no2monod_denom(i) = Cno2_st(i)+this%Ka2
  r2no2monod(i) = Cno2_st(i)/r2no2monod_denom(i)
  r2kin(i) = k2_scaled*r2docmonod(i)*r2no2monod(i) 
!&
!               *this%Ki2/(Co2_st(i)+this%Ki2)

  ! O2 ->
  r3docmonod_denom(i) = Cdoc_st(i)+this%Kd3
!Fang test
  r3docmonod(i) = Cdoc_st(i)/r3docmonod_denom(i)
  !r3docmonod(i) = 1.d0
  r3o2monod_denom(i) = Co2_st(i)+this%Ka3
  r3o2monod(i) = Co2_st(i)/r3o2monod_denom(i)
  r3kin(i) = k3_scaled*r3docmonod(i)*r3o2monod(i)

  sumkin(i) = r1kin(i) + r2kin(i) + r3kin(i)
  sumkinsq(i) = sumkin(i) * sumkin(i)

  u1(i) = 0.d0
  if (r1kin(i) > 0.d0) u1(i) = r1kin(i)/sumkin(i)
  u2(i) = 0.d0
  if (r2kin(i) > 0.d0) u2(i) = r2kin(i)/sumkin(i)
  u3(i) = 0.d0
  if (r3kin(i) > 0.d0) u3(i) = r3kin(i)/sumkin(i)

  rate(this%nrxn_swat+i) = u1(i)*r1kin(i)  ! mol/mol biomass/sec
  rate(this%nrxn_swat+i+MRMT) = u2(i)*r2kin(i)
  rate(this%nrxn_swat+i+2*MRMT) = u3(i)*r3kin(i)
 enddo
! mass transfer, rate unit in mol/s/L
! assuming porosity in storage zone is 1
  ! no3 -
  do i = 1, MRMT
    rkin_no3_mrmt(i) = this%alpha(i)*(Cno3 - Cno3_st(i))
    rkin_no2_mrmt(i) = this%alpha(i)*(Cno2 - Cno2_st(i))
    rkin_nh4_mrmt(i) = this%alpha(i)*(Cnh4 - Cnh4_st(i))
    rkin_n2_mrmt(i) = this%alpha(i)*(Cn2 - Cn2_st(i))
    rkin_o2_mrmt(i) = this%alpha(i)*(Co2 - Co2_st(i))
    rkin_co2_mrmt(i) = this%alpha(i)*(Cco2 - Cco2_st(i))
!Fang
    rkin_doc_mrmt(i) = this%alpha(i)*(Cdoc - Cdoc_st(i))
!    rkin_doc_mrmt(i) = this%alpha(i)*(Ccbod - Cdoc_st(i))
!    rkin_X_mrmt(i) = this%alpha(i)*(X - X_st(i))
  enddo
  non_mrmt = this%nrxn - MRMT*this%nrxn_mt
  do i = 1, MRMT
    non_mrmt = non_mrmt + 1
    rate(non_mrmt) = rkin_doc_mrmt(i)
  enddo
  do i = 1, MRMT
    non_mrmt = non_mrmt + 1
    rate(non_mrmt) = rkin_no3_mrmt(i)
  enddo
  do i = 1, MRMT
    non_mrmt = non_mrmt + 1
    rate(non_mrmt) = rkin_no2_mrmt(i)
  enddo
  do i = 1, MRMT
    non_mrmt = non_mrmt + 1
    rate(non_mrmt) = rkin_nh4_mrmt(i)
  enddo
  do i = 1, MRMT
    non_mrmt = non_mrmt + 1
    rate(non_mrmt) = rkin_o2_mrmt(i)
  enddo
  do i = 1, MRMT
    non_mrmt = non_mrmt + 1
    rate(non_mrmt) = rkin_n2_mrmt(i)
  enddo
  do i = 1, MRMT
    non_mrmt = non_mrmt + 1
    rate(non_mrmt) = rkin_co2_mrmt(i)
  enddo
!  do i = 1, MRMT
!    non_mrmt = non_mrmt + 1
!    rate(non_mrmt) = rkin_X_mrmt(i)
!  enddo
! Before subtracting sum of rates
  Residual_tmp(:) = 0.d0 
  Residual_tmp(:) = Residual(:)
!
! Calculate residual
!
! swat reactions
  do irxn = 1, this%nrxn_swat
    do i = 1, this%nrow(irxn)
! bulk volume
      Residual(this%irow(i,irxn)) = Residual(this%irow(i,irxn)) - &
        this%stoich_row(i,irxn)%ptr * rate(irxn) * L_water
    enddo
  enddo
!
! cybernetic reactions, only in HZ, all immobile
!
  do irxn = this%nrxn_swat+1, this%nrxn - MRMT*this%nrxn_mt
    do i = 1, this%nrow(irxn)
      ! mol/sec
      ! X is in [M]
      j = reaction%offset_immobile + abs(this%irow(i,irxn))
      if(irxn > 2*MRMT + this%nrxn_swat) then
        i_ord = irxn - 2*MRMT - this%nrxn_swat
      elseif(irxn > MRMT + this%nrxn_swat) then
        i_ord = irxn - MRMT - this%nrxn_swat
      else
        i_ord = irxn - this%nrxn_swat
      endif
      Residual(j) = Residual(j) - &
        this%stoich_row(i,irxn)%ptr * rate(irxn) * X_st(i_ord) * &
        volume
        ! if biomass is aqueous multiply by L_water
        ! if biomass is immobile multiply by volume
!        L_water
    enddo
  enddo
!
!  i = this%offset_auxiliary !stoichiometric coefficients
!  Residual(this%a_id) = Residual(this%a_id) + &
!         rt_auxvar%auxiliary_data(i+15) * Ca * L_water !sig1OVd
!  Residual(this%orgp_id) = Residual(this%orgp_id) + &
!         rt_auxvar%auxiliary_data(i+16) * Corgp * L_water !sig5
!  Residual(this%orgn_id) = Residual(this%orgn_id) + &
!         rt_auxvar%auxiliary_data(i+17) * Corgn * L_water !sig4
!  Residual(this%nh4_id) = Residual(this%nh4_id) - &
!         rt_auxvar%auxiliary_data(i+23) * L_water  !sig3OVd
!  Residual(this%dop_id) = Residual(this%dop_id) - &
!         rt_auxvar%auxiliary_data(i+18) * L_water  !sig2OVd
!  Residual(this%o2_id) = Residual(this%o2_id) - &
!         rt_auxvar%auxiliary_data(i+19)*(rt_auxvar%auxiliary_data(i+20) - &
!         Co2) * L_water !k2,eqO2
!  if(Co2/option%tran_dt >= rt_auxvar%auxiliary_data(i+21)) then
!     Residual(this%o2_id) = Residual(this%o2_id) + rt_auxvar%auxiliary_data(i+21) * L_water ! k4/d, in case it's consumed more than present
!  else
!     Residual(this%o2_id) = Residual(this%o2_id) + Co2/option%tran_dt * L_water
!  end if
! Fang fix doc in stream
!  Residual(this%doc_id) = 0.0
!
  i = this%offset_auxiliary !stoichiometric coefficients
  Residual(this%a_id) = Residual(this%a_id) + &
         rt_auxvar%auxiliary_data(i+15) * Ca * L_water !sig1OVd
! In case not enough to react
  sum_rates = Residual_tmp(this%a_id) - Residual(this%a_id)
  cdt = Ca/option%tran_dt
  if(Ca > 1.d-10 .and. sum_rates < 0.d0) then
!    if(cdt < abs(sum_rates)) &
!    Residual(this%a_id) = Residual_tmp(this%a_id) + cdt * L_water
  endif
  Residual(this%orgp_id) = Residual(this%orgp_id) + &
         rt_auxvar%auxiliary_data(i+16) * Corgp * L_water !sig5
  sum_rates = Residual_tmp(this%orgp_id) - Residual(this%orgp_id)
  cdt = Corgp/option%tran_dt
  if(Corgp > 1.d-10 .and. sum_rates < 0.d0) then
!    if(cdt < abs(sum_rates)) &
!    Residual(this%orgp_id) = Residual_tmp(this%orgp_id) + cdt * L_water
  endif
  Residual(this%orgn_id) = Residual(this%orgn_id) + &
         rt_auxvar%auxiliary_data(i+17) * Corgn * L_water !sig4
  sum_rates = Residual_tmp(this%orgn_id) - Residual(this%orgn_id)
  cdt = Corgn/option%tran_dt
  if(Corgn > 1.d-10 .and. sum_rates < 0.d0) then
!    if(cdt < abs(sum_rates)) &
!    Residual(this%orgn_id) = Residual_tmp(this%orgn_id) + cdt * L_water
  endif
  Residual(this%nh4_id) = Residual(this%nh4_id) - &
         rt_auxvar%auxiliary_data(i+23) * L_water  !sig3OVd
  sum_rates = Residual_tmp(this%nh4_id) - Residual(this%nh4_id)
  cdt = Cnh4/option%tran_dt
  if(Cnh4 > 1.d-10 .and. sum_rates < 0.d0) then
!    if(cdt < abs(sum_rates)) &
!    Residual(this%nh4_id) = Residual_tmp(this%nh4_id) + cdt * L_water
  endif
  Residual(this%dop_id) = Residual(this%dop_id) - &
         rt_auxvar%auxiliary_data(i+18) * L_water  !sig2OVd
  sum_rates = Residual_tmp(this%dop_id) - Residual(this%dop_id)
  cdt = Cdop/option%tran_dt
  if(Cdop > 1.d-10 .and. sum_rates < 0.d0) then
!    if(cdt < abs(sum_rates)) &
!    Residual(this%dop_id) = Residual_tmp(this%dop_id) + cdt * L_water
  endif
  Residual(this%o2_id) = Residual(this%o2_id) - &
         rt_auxvar%auxiliary_data(i+19)*(rt_auxvar%auxiliary_data(i+20) - &
         Co2) * L_water !k2,eqO2
  Residual(this%o2_id) = Residual(this%o2_id) + rt_auxvar%auxiliary_data(i+21) * L_water ! k4/d, in case it's consumed more than present
  sum_rates = Residual_tmp(this%o2_id) - Residual(this%o2_id)
  cdt = Co2/option%tran_dt
  if(sum_rates < 0.d0) then
    if(cdt < abs(sum_rates/L_water)) &
    Residual(this%o2_id) = Residual_tmp(this%o2_id) + cdt * L_water
  endif

  ! decay of biomass
  ! if biomass is aqueous multiply by L_water
  ! if biomass is immobile multiply by volume
  do i = 1, MRMT
    j = reaction%offset_immobile + abs(this%biomass_st_id(i))
    Residual(j) = Residual(j) + &
!                              k_deg_scaled * X * L_water
                              k_deg_scaled * X_st(i) * volume
  ! production of doc by biomass decay
  ! note the addition
  ! mol/sec
    j = reaction%offset_immobile + abs(this%doc_st_id(i))
    Residual(j) = Residual(j) - &
!                          k_deg_scaled/this%f_act * X * L_water
                          5.d0*k_deg_scaled * X_st(i) * volume
  enddo
  ! uptake of nitrate
  do i = 1, MRMT
    j = reaction%offset_immobile + abs(this%no3_st_id(i))
    !Residual(j) = Residual(j) + &
    !                        Cno3_st(i) * 11.e-6/3600.d0/(52e-6+Cno3_st(i)) * volume !11umol/L/hr,ksat = 52 umol/L
!    Residual(j) = Residual(j) + &
!                            Cno3_st(i) * 1e-8/3600.d0/(20e-6+Cno3_st(i)) * volume !11umol/L/hr,ksat = 52 umol/L
!Fang - commented out
!    Residual(j) = Residual(j) + 3.13e-7/(3.13e-5+Co2_st(i))* &
!                            Cno3_st(i) * 8.64e-5/24/3600.d0/(32.3e-6+Cno3_st(i)) * volume !8.64e-2 mM/d, two orders lower used by me, for all the monthly simulations
!    Residual(j) = Residual(j) + 3.13e-6/(3.13e-5+Co2_st(i))* &
!                            Cno3_st(i) * 8.64e-5/24/3600.d0/(32.3e-6+Cno3_st(i)) * volume !8.64e-2 mM/d, two orders lower used by me
  enddo
  ! uptake of o2
!Fang - commented out
!  do i = 1, MRMT
!    j = reaction%offset_immobile + abs(this%o2_st_id(i))
!    Residual(j) = Residual(j) + &
!                            Co2_st(i) * 4.78e-6/24/3600.d0/(6.25e-6+Co2_st(i)) * volume !4.78e-1mM/d
!  enddo
! multirate
  do irxn = this%nrxn - MRMT*this%nrxn_mt + 1, this%nrxn
    do i = 1, this%nrow(irxn)
      ! mol/sec
      ! X is in [M]
      ! immobile aqueous, multiply * volume for book keeping
      if( this%irow(i,irxn) <  0) then
        j = reaction%offset_immobile + abs(this%irow(i,irxn))
        Residual(j) = Residual(j) - &
          this%stoich_row(i,irxn)%ptr * rate(irxn) *volume
      else
        j = this%irow(i,irxn)
        Residual(j) = Residual(j) - &
          this%stoich_row(i,irxn)%ptr * rate(irxn) * L_water
      endif
        ! if biomass is immobile multiply by volume
    enddo
  enddo
! assign cbod to doc
  
!  residual(this%cbod_id) = residual(this%doc_id)

! multirate, storage zone, immobile aqueous units: mol/L aqueous
!  do i = 1, MRMT
!    Cno3_st(i) = rt_auxvar%immobile(this%no3_st_id(i)) 
!  enddo  

! mass transfer, rate unit in mol/s/L  
! assuming porosity in storage zone is 1
  ! no3 - 
!  do i = 1, MRMT
!    rkin_no3_mrmt(i) = this%alpha(i)*(Cno3 - Cno3_st(i))
!  enddo
!  non_mrmt = 0
!  do i = 1, MRMT
!    non_mrmt = non_mrmt + 1
!    rate(non_mrmt) = rkin_no3_mrmt(i)  
!  enddo

! multirate
!  do irxn = this%nrxn+1, this%nrxn+MRMT
!    do i = 1, this%nrow(irxn)
      ! mol/sec
      ! X is in [M]
      ! immobile aqoues, multiply * volume for book keeping
!      if( this%irow(i,irxn) <  0) then
!        j = reaction%offset_immobile + abs(this%irow(i,irxn))
!        Residual(j) = Residual(j) - &
!          this%stoich_row(i,irxn)%ptr * rate(irxn) *volume
!      else
!        j = this%irow(i,irxn)
!        Residual(j) = Residual(j) - &
!          this%stoich_row(i,irxn)%ptr * rate(irxn) * L_water
!      endif
        ! if biomass is immobile multiply by volume
!    enddo
!  enddo
! external source
!  i = reaction%nauxiliary
!  this%no3_src = rt_auxvar%auxiliary_data(i)
!  residual(this%no3_id) = residual(this%no3_id) - this%no3_src * volume
                 
  
end subroutine CyberReact

! ************************************************************************** !

subroutine CyberUpdateKineticState(this,rt_auxvar,global_auxvar, &
                                   material_auxvar,reaction,option)
  ! 
  ! Updates kinetic state (e.g. increments mass balance based on current rate
  ! and time step size)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/21/16
  ! 
  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_class

  implicit none

  class(reaction_sandbox_cyber_mrmt_type) :: this
  type(option_type) :: option
  type(reaction_type) :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar

  PetscInt :: i, irate, icum 

  if (this%store_cumulative_mass) then
    irate = this%offset_auxiliary
    icum = this%offset_auxiliary + 6
    do i = 1, 6
      ! rate is in mol/sec
      icum = icum + 1
      rt_auxvar%auxiliary_data(icum) = rt_auxvar%auxiliary_data(icum) + &
        rt_auxvar%auxiliary_data(irate+i) * &  ! mole/m^3-second
        option%tran_dt
    enddo
  endif

end subroutine CyberUpdateKineticState

! ************************************************************************** !

subroutine CyberDestroy(this)
  ! 
  ! Destroys allocatable or pointer objects created in this
  ! module
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/01/15
  ! 
  use Utility_module
  
  implicit none
  
  class(reaction_sandbox_cyber_mrmt_type) :: this  

  call DeallocateArray(this%alpha)
  call DeallocateArray(this%th_frac)
  call DeallocateArray(this%nrow)
  call DeallocateArray(this%irow)
  call DeallocateArray(this%icol)
  !nullify(this%stoich_row)

end subroutine CyberDestroy

end module Reaction_Sandbox_Cyber_Mrmt_class
