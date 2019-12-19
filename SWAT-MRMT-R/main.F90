      !include 'modparm.f'
      program main
!!    this is the main program that reads input, calls the main simulation
!!    model, and writes output.
!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!         ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    date        |NA            |date simulation is performed where leftmost
!!                               |eight characters are set to a value of
!!                               |yyyymmdd, where yyyy is the year, mm is the 
!!                               |month and dd is the day
!!    isproj      |none          |special project code:
!!                               |1 test rewind (run simulation twice)
!!    time        |NA            |time simulation is performed where leftmost
!!                               |ten characters are set to a value of
!!                               |hhmmss.sss, where hh is the hour, mm is the 
!!                               |minutes and ss.sss is the seconds and
!!                               |milliseconds
!!    values(1)   |year          |year simulation is performed
!!    values(2)   |month         |month simulation is performed
!!    values(3)   |day           |day in month simulation is performed
!!    values(4)   |minutes       |time difference with respect to Coordinated
!!                               |Universal Time (ie Greenwich Mean Time)
!!    values(5)   |hour          |hour simulation is performed
!!    values(6)   |minutes       |minute simulation is performed
!!    values(7)   |seconds       |second simulation is performed
!!    values(8)   |milliseconds  |millisecond simulation is performed
!!    zone        |NA            |time difference with respect to Coordinated
!!                               |Universal Time (ie Greenwich Mean Time)
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    prog        |NA            |program name and version
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    i           |none          |counter
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: date_and_time
!!    SWAT: getallo, allocate_parms, readfile, readfig
!!    SWAT: readbsn, std1, readwwq, readinpt, std2, storeinitial
!!    SWAT: openwth, headout, simulate, finalbal, writeaa, pestw 
!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use parm
!      use Reaction_module
!      use Reaction_Aux_module
!      use Reactive_Transport_Aux_module
!      use Global_Aux_module
!      use Material_Aux_class
!      use Reaction_Database_module
!      use Option_module
      use pflotran_rxn, only: pflotran_rxn_init
      implicit none

!      interface
!        subroutine pflotran_rxn_init(option, reaction, global_auxvars, &
!          rt_auxvars,material_auxvars)
  
!#include "petsc/finclude/petscsys.h"
!          use petscsys
!          use Reaction_module
!          use Reaction_Aux_module
!          use Reactive_Transport_Aux_module
!          use Global_Aux_module
!          use Material_Aux_class
!          use Reaction_Database_module
!          use Option_module
!          use Input_Aux_module
!          use String_module
  
!          use Transport_Constraint_module
!          use PFLOTRAN_Constants_module

!          use BatchChem
!          use parm

!          implicit none
!
!
!          PetscErrorCode :: ierr
!          PetscBool :: option_found  
!          character(len=MAXSTRINGLENGTH) :: string
!          character(len=MAXSTRINGLENGTH) :: filename_out
!          type(reaction_type), pointer :: reaction
!          type(option_type), pointer :: option
!          type(input_type), pointer :: input

!          type(global_auxvar_type), pointer :: global_auxvars
!          type(reactive_transport_auxvar_type), pointer :: rt_auxvars
!          class(material_auxvar_type), pointer :: material_auxvars

!          character(len=MAXWORDLENGTH) :: card
!          character(len=MAXWORDLENGTH) :: word
!          type(tran_constraint_list_type), pointer :: transport_constraints
!          type(tran_constraint_coupler_type), pointer :: constraint_coupler
! 
!          PetscInt :: istart, iend, iendaq
!          PetscInt :: immobile_start, immobile_end
!          PetscInt :: step, max_steps, num_iterations
!          PetscReal, dimension(:), allocatable :: tran_xx_p
!          PetscReal :: t_max, dt
!          PetscInt :: jrch
!        end subroutine pflotran_rxn_init
!      end interface

!      type(reaction_type), pointer :: reaction
!      type(option_type), pointer :: option

!      type(global_auxvar_type), pointer :: global_auxvars
!      type(reactive_transport_auxvar_type), pointer :: rt_auxvars
!      class(material_auxvar_type), pointer :: material_auxvars
      prog = "SWAT Dec 23 2016    VER 2016/Rev 664"
      write (*,1000)
 1000 format(1x,"               SWAT2016               ",/,            & 
               "               Rev. 664               ",/,            & 
               "      Soil & Water Assessment Tool    ",/,            & 
               "               PC Version             ",/,            & 
               " Program reading from file.cio . . . executing",/)
		
      open(9003,file='hz_conc.dat',status='unknown')
      open(9004,file='river_doc93.dat',status='unknown')
      open(9005,file='river_doc101.dat',status='unknown')
      call getallo
      call allocate_parms
      call readfile
      call readbsn
      call readwwq
! initialize pflotran chem
!     option => OptionCreate()
      call pflotran_rxn_init()
!       print *,'reaction%ncomp---',reaction%ncomp
!! process input
      if (fcstyr > 0 .and. fcstday > 0) call readfcst
      call readplant             !! read in the landuse/landcover database
      call readtill              !! read in the tillage database
      call readpest              !! read in the pesticide database
      call readfert              !! read in the fertilizer/nutrient database
      call readurban             !! read in the urban land types database
      call readseptwq            !! read in the septic types database
      call readlup
      call readfig
      call readatmodep
      call readinpt
      call std1
      call std2
      call openwth
      call headout

      !! convert integer to string for output.mgt file
      subnum = ""
      hruno = ""
      do i = 1, mhru
        write (subnum(i),fmt=' (i5.5)') hru_sub(i)
        write (hruno(i),fmt=' (i4.4)') hru_seq(i)  
      end do

      if (isproj == 2) then 
        hi_targ = 0.0
      end if

!! save initial values
      if (isproj == 1) then
        scenario = 2
        call storeinitial
      else if (fcstcycles > 1) then
        scenario =  fcstcycles
        call storeinitial
      else
        scenario = 1
      endif
        if (iclb /= 4) then
      do iscen = 1, scenario

     
        !! simulate watershed processes
        call simulate

        !! perform summary calculations
        call finalbal
        call writeaa
        call pestw

        !!reinitialize for new scenario
        if (scenario > iscen) call rewind_init
      end do
         end if
      do i = 101, 109       !Claire 12/2/09: change 1, 9  to 101, 109.
        close (i)
      end do
      close(124)
      write (*,1001)
 1001 format (/," Execution successfully completed ")
	
        iscen=1
!! file for Mike White to review to ensure simulation executed normally
      open (9999,file='fin.fin')
      write (9999,*) 'Execution successful'
      close (9999)
      
	stop
      end
