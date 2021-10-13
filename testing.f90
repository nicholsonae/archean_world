MODULE random_n
  IMPLICIT NONE
  REAL, PRIVATE      :: zero = 0.0, half = 0.5, one = 1.0, two = 2.0, vsmall = TINY(1.0), vlarge = HUGE(1.0)
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)


CONTAINS


FUNCTION random_normal() RESULT(fn_val)

! Adapted from the following Fortran 77 code
!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

!  The function random_normal() returns a normally distributed pseudo-random
!  number with zero mean and unit variance.

!  The algorithm uses the ratio of uniforms method of A.J. Kinderman
!  and J.F. Monahan augmented with quadratic bounding curves.

REAL :: fn_val

!     Local variables
REAL     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472
REAL     :: r1 = 0.27597, r2 = 0.27846, u, v, x, y, q

!     Generate P = (u,v) uniform in rectangle enclosing acceptance region

DO
  CALL RANDOM_NUMBER(u)
  CALL RANDOM_NUMBER(v)
  v = 1.7156 * (v - 0.5)

!     Evaluate the quadratic form
  x = u - s
  y = ABS(v) - t
  q = x**2 + y*(a*y - b*x)

!     Accept P if inside inner ellipse
  IF (q < r1) EXIT
!     Reject P if outside outer ellipse
  IF (q > r2) CYCLE
!     Reject P if outside acceptance region
  IF (v**2 < -4.0*LOG(u)*u**2) EXIT
END DO

!     Return ratio of P's coordinates as the normal deviate
fn_val = v/u
RETURN

END FUNCTION random_normal

END MODULE random_n



program testing
  
  use random_n
  implicit none

  
         !! real, parameter                   :: dp = selected_real_kind(15, 307)
  
          TYPE :: Microbe
             integer                        :: id
             real(16)                       :: population
             real(16)                       :: ATP
             integer                        :: metabolism
          END TYPE Microbe

          TYPE :: Environment
             real(16)                       :: H2
             real(16)                       :: CO2
             real(16)                       :: CH4
             real                           :: current_T
             real                           :: eq_T
             real                           :: T_increment
          END TYPE Environment
          
          
          real, parameter                   :: pi = 3.1415927
          real, parameter                   :: av_constant = 6.02214076E23 !! avogadros constant
          real, parameter                   :: H2_Cmax = (60*60) * 3.76E-17        !! mol_H2 cell^-1 s^-1
          real(16), parameter               :: cell_main = (60*60) * 2.16E-19      !! mol_ATP cell^-1 s^-1
          real(16), parameter               :: cell_growth  = 4.237E-14    !! mole ATP to make a cell
          real, parameter                   :: death_starve = 2.5E-7       !! s^-1
          real, parameter                   :: death_other  = 1E-15        !! s^-1
          real, parameter                   :: protein_cell = 7.4E-15      !! mol_CH2O cell^-1
          real, parameter                   :: T_ideal      = 283          !! microbe ideal T
          real, parameter                   :: T_sens       = 0.0          !! microbe sensivity to temperature
          real, parameter                   :: ATP_CH4      = 0.6          !! moles of ATP per moles of CH4 produced
          real, parameter                   :: ATP_H2       = 0.15
          real, parameter                   :: H2_lim       = 0.0          !! lower bound of H2 consumption
          real, dimension(90)               :: lat                         !! lattitude points on globe
          real, dimension(144)              :: lon                         !! longitude points on globe
          integer, parameter                :: t_array_length = 10001      !! length of temeperature array, was 10001
          
          real, parameter                   :: planet_r = 6051.8E3         !! radius of planet
          real, parameter                   :: atmosphere_mass = 5.15E18   !! mass of the atmosphere
          real(16), parameter               :: moles_air = 1.73E20 * 0.5   !! moles of air in atmosphere
          real(16), parameter               :: ocean_volume = 9.2E14       !! volume of 2m deep ocean in m^3
          real(16), parameter               :: ocean_surf_area = 4602.3E11 !! ocean surface area m^2
          real, parameter                   :: atmo_pressure = 1.0         !! atmopsheric pressure in bar in archean
          
          real, parameter                   :: piston_vel_CO2 = 7.3E-6     !! piston velocity CO2 at 25 degrees C
          real, parameter                   :: solubility_CO2 = 0.035      !! solubility of CO2
          real, parameter                   :: piston_vel_CH4 = 4.5E-5     !! piston velocity of CH4 at 25 degrees C
          real, parameter                   :: solubility_CH4 = 1.4E-3     !! solubility of CH4
          real, parameter                   :: piston_vel_H2  = 1.3E-4     !! piston velocity for H2 at 25 degrees
          real, parameter                   :: solubility_H2  = 7.8E-4     !! solubility of H2 mol L^-1 bar^-1

          real, parameter                   :: CO2_burial       = 0.001    !! % of CO2 removed from atmosphere per year
          real, parameter                   :: CH4_burial       = 0.001    !! % of CH4 removed from atmosphere per year
          real, parameter                   :: H2_burial        = 0.001    !! % of H2 removed from atmosphere per year
          real, parameter                   :: CO2_mole_weight  = 44.01    !! molar mass of CO2
          real, parameter                   :: CH4_mole_weight  = 16.043   !! molar mass of CH4
          real, parameter                   :: H2_mole_weight   = 2.01588  !! molar mass of H2
          
          real                              :: CO2_outflux_year = 10E15   !! outflux in moles CO2 20x as much as current was 5.5E  
          real                              :: H2_outflux_year  = 10E15     !! outflux in moles H2  was 5E13  !! 10E15
          
          real                              :: H2_pp, CH4_pp, CO2_pp, T    !! environment
          real                              :: c, d, mu, omega
          real                              :: new_bugs, bug_starve, bugs_death, bugs_culled, tot_died
          real                              :: fit_level = -1
          real                              :: a_r, ATP_used, ATP_made
          real                              :: ATP_maintain
          real                              :: max_H2
          real                              :: ATP_starve
          real(16)                          :: CO2_MMR, CH4_MMR
          real                              :: CO2_growth, H2_growth
          real                              :: towrite_CO2, towrite_CH4
          real                              :: starve_random, repro_random, death_random
          real                              :: CO2_atmo_2_ocean
          real                              :: H2_atmo_2_ocean
          real                              :: CH4_atmo_2_ocean
          real                              :: biotic_CH4_output

          real                              :: pp_CO2
          real                              :: pp_CH4
          real                              :: pp_H2
         
          TYPE(Microbe)                     :: species1
          TYPE(Microbe), dimension(1)       :: species_list
          
          TYPE(Environment)                 :: atmosphere
          TYPE(Environment)                 :: ocean
          
          real, dimension(t_array_length,t_array_length)      :: temp_array
          real, dimension(t_array_length)            :: CO2_array, CH4_array
          integer                           :: i, j, n, t_step, IOStatus
          integer                           :: biotic_step
          integer                           :: d_step

          integer                           :: temp_value_CH4
          integer                           :: temp_value_CO2
          integer                           :: seeded = 0
          integer                           :: habitable = 1 !! CHANGE THIS
          integer                           :: step_length = 365*24

          real                              :: H2_to_add, CO2_to_add, CH4_to_add
          real                              :: H2_ocean_to_add, CO2_ocean_to_add, CH4_ocean_to_add
          !!real                              :: t_bugs_starved, t_bugs_born
          integer                           :: seed_input
          character(100)                    :: num1char, num2char, file_number

          CALL GET_COMMAND_ARGUMENT(1,num1char)   !first, read in the two values
          !!print *, num1char
          CALL GET_COMMAND_ARGUMENT(2,file_number)
          !!print *, file_number
          READ(num1char,*)seed_input                    !then, convert them to REALs
          !!print *, seed_input
          !!READ(num2char,*)file_number

          call random_seed(put=[seed_input,seed_input,seed_input,seed_input,seed_input,seed_input,seed_input,seed_input])


          open(10, file="temperature_file.txt", status='old', action='read')
          DO i = 1, t_array_length
             read(10,*) (temp_array(j,i), j = 1,t_array_length)
          END DO
          close(10)

          open(20, file="data_dummy_"//trim(file_number)//".txt", status='replace', action='write')

          open(30, file="data_random_"//trim(file_number)//".txt", status='replace', action='write')

          DO i = 1, t_array_length
             CO2_array(i) = 0.005 + 0.0000095  * (i-1)
             CH4_array(i) = 0.000001 * (i - 1)
          END DO
          !!print *, "CO2_array", CO2_array
          !!print *, "CH4_array", CH4_array

          species1%id           = 0
          species1%population   = 0
          species1%ATP          = 0
          species1%metabolism   = 0

          species_list(1)       = species1

          !!species_list(1)%population = species_list(1)%population

          CO2_MMR               = 0.0 !! atmosphere%CO2 / moles_air
          CH4_MMR               = 0.0 !! atmosphere%CH4 / moles_air

          temp_value_CO2 = CO2_array_func(CO2_MMR, CO2_array, t_array_length)
          temp_value_CH4 = CH4_array_func(CH4_MMR, CH4_array, t_array_length)
          
          !!print *, "temp_value_CO2", temp_value_CO2, "temp_value_CH4", temp_value_CH4

          !! set up the atmosphere

          atmosphere%CO2          = moles_air * 0.01
          atmosphere%CH4          = 0.0
          atmosphere%H2           = moles_air * 0.0001
          atmosphere%current_T    = temp_array(temp_value_CO2, temp_value_CH4)
          atmosphere%eq_T         = temp_array(temp_value_CO2, temp_value_CH4)
          atmosphere%T_increment  = 0

          !!DO i = 1, 1000
          !!   atmosphere%CO2 = atmosphere%CO2 * (1 - CO2_burial) + CO2_outflux_year
          !!   atmosphere%H2  = atmosphere%H2  * (1 - H2_burial)  + H2_outflux_year
          !!END DO

          ocean%CO2               = 0.0 !! moles_air * 0.01
          ocean%CH4               = 0.0
          ocean%H2                = 0.0 !!moles_air * 0.0001
          ocean%current_T         = temp_array(temp_value_CO2, temp_value_CH4)
          ocean%eq_T              = temp_array(temp_value_CO2, temp_value_CH4)
          ocean%T_increment       = 0

          biotic_CH4_output = 0


          lat  = (/ (I, I=-89, 89, 2) /)
          DO i = 1, 144
             lon(i) = 1.25 + 2.5*(i-1)
          END DO

          DO t_step = 0, 41000  !! run the simulation for 1 year then update the temperature

             write(20, *) t_step, atmosphere%current_T, species1%population, species1%ATP, &
                  H2_to_add, atmosphere%H2, atmosphere%CO2, atmosphere%CH4, ocean%H2, ocean%CO2, ocean%CH4, &
                  biotic_CH4_output, fit_level

             CO2_to_add = -1.0 * CO2_burial * atmosphere%CO2 + (CH4_burial *  atmosphere%CH4) + CO2_outflux_year
             H2_to_add  = -1.0 * H2_burial  * atmosphere%H2  + (4 * CH4_burial *  atmosphere%CH4) + H2_outflux_year
             CH4_to_add = -1.0 * CH4_burial * atmosphere%CH4

             CO2_MMR        = (atmosphere%CO2 + CO2_to_add) / moles_air
             CH4_MMR        = (atmosphere%CH4 + CH4_to_add) / moles_air

             temp_value_CO2 = CO2_array_func(CO2_MMR, CO2_array, t_array_length)
             temp_value_CH4 = CH4_array_func(CH4_MMR, CH4_array, t_array_length)
             IF (temp_value_CO2 .LT. 1) THEN
                PRINT *, "co2 value error"
                PRINT *, atmosphere%CO2
                PRINT *, CO2_to_add
                PRINT *,  moles_air
                PRINT *, "ocean CO2"
                PRINT *, ocean%CO2
             END IF
              IF (temp_value_CO2 .LT. 1) THEN
                PRINT *, "ch4 value error"
                PRINT *, atmosphere%CH4
                PRINT *, CH4_to_add
                PRINT *,  moles_air
                PRINT *, "Ocean CH4"
                PRINT *, ocean%CH4
                PRINT *, "atmosphere H2"
                PRINT *, atmosphere%H2
                PRINT *, "ocean H2"
                PRINT *, ocean%H2
                PRINT *, "listing vairables in order for CO2_to_add"
                PRINT *, CO2_burial
                PRINT *, atmosphere%CO2
                PRINT *, CH4_burial
                PRINT *, atmosphere%CH4
                PRINT *, CO2_outflux_year
                PRINT *, "t_step"
                PRINT *, t_step
                PRINT *, "population"
                PRINT *, species1%population
             END IF
         

             atmosphere%eq_T = temp_array(temp_value_CO2, temp_value_CH4)
             atmosphere%T_increment = (atmosphere%eq_T - atmosphere%current_T) / step_length

             IF ((seeded .EQ. 0) .AND. (T_ideal - 3 < atmosphere%current_T) .AND. (atmosphere%current_T < T_ideal + 3)) THEN
                habitable = 1
             END IF
             towrite_CO2 = atmosphere%CO2 / moles_air
             towrite_CH4 = atmosphere%CH4 / moles_air
             
             IF (habitable .EQ. 1 .AND. t_step > 20000 .AND. seeded .EQ. 0) THEN !! seed if appropiate
                !!IF (seeded .EQ. 0) THEN   !! seed the planet with life
                   species1%id           = 0
                   species1%population   = 1E2
                   species1%ATP          = 2.0 * cell_main * 1E2 !!* av_constant
                   species1%metabolism   = 0
                   seeded                = 1
             END IF
             
             DO biotic_step = 1, step_length !! years in hours
                !!print *, "biotic_step", biotic_step, "population", species1%population

                !!temp_value_CO2 = CO2_array_func(CO2_MMR, CO2_array)
                !!temp_value_CH4 = CH4_array_func(CH4_MMR, CH4_array)
                !!atmosphere%eq_T = temp_array(temp_value_CO2, temp_value_CH4)
                !!atmosphere%current_T = atmosphere%eq_T
                
                atmosphere%current_T = atmosphere%current_T + atmosphere%T_increment
                atmosphere%CO2 = atmosphere%CO2   + CO2_to_add / step_length
                atmosphere%H2  = atmosphere%H2    + H2_to_add  / step_length
                atmosphere%CH4 = atmosphere%CH4   + CH4_to_add / step_length

                CO2_atmo_2_ocean =  gas_flux_atmo_ocean(atmosphere%CO2, piston_vel_CO2, solubility_CO2, ocean%CO2)
                CH4_atmo_2_ocean =  gas_flux_atmo_ocean(atmosphere%CH4, piston_vel_CH4, solubility_CH4, ocean%CH4)
                H2_atmo_2_ocean  =  gas_flux_atmo_ocean(atmosphere%H2,  piston_vel_H2, solubility_H2,   ocean%H2 )

                CO2_ocean_to_add = 60*60* CO2_atmo_2_ocean * ocean_surf_area
                CH4_ocean_to_add = 60*60* CH4_atmo_2_ocean * ocean_surf_area
                H2_ocean_to_add  = 60*60* H2_atmo_2_ocean  * ocean_surf_area

                atmosphere%CO2 = atmosphere%CO2 - CO2_ocean_to_add
                atmosphere%CH4 = atmosphere%CH4 - CH4_ocean_to_add
                atmosphere%H2  = atmosphere%H2  - H2_ocean_to_add

                ocean%CO2 = ocean%CO2 + CO2_ocean_to_add
                ocean%CH4 = ocean%CH4 + CH4_ocean_to_add
                ocean%H2  = ocean%H2  + H2_ocean_to_add
                   
                new_bugs = 0
                bug_starve = 0

                IF (t_step .GT. 30000) THEN
                   species1%population = 0
                END IF

                IF (species1%population .GT. 0) THEN

                   death_random = 1 + 0.02 * random_normal()
                   !!IF (death_random > 1) THEN
                   !!   death_random = 1
                   !!END IF
                   
                   species1%population = species1%population * 0.98 * death_random
                   species1%ATP        = species1%ATP        * 0.98 * death_random

                   IF (species1%population < 1) THEN
                      species1%population = 0
                      species1%ATP        = 0

                   ELSE

                      starve_random = 1 + 0.02 * random_normal() !! adds noise to runs
                      !!write(30, *) starve_random
                      
                      bug_starve = starve_random * num_bugs_starved(species1%ATP, species1%population, cell_main)

                      !!print *, "starved pop", bug_starve
                      IF (bug_starve >= species1%population) THEN
                         ATP_starve = species1%ATP
                         bug_starve = species1%population
                      ELSE IF (bug_starve > 0) THEN
                         ATP_starve = ATP_of_starved(species1%ATP, species1%population, cell_main)
                      ELSE
                         ATP_starve = 0
                         bug_starve = 0
                      END IF

                      !!IF (species1%population < 1) THEN
                      !!   species1%population = 0
                      !!   species1%ATP        = 0

                      !!ELSE

                         ATP_maintain = cell_main * (species1%population - bug_starve)
                         repro_random = 1 + 0.02 * random_normal()
                         !!write(30, *) repro_random
                         
                         new_bugs = repro_random * bugs_made(species1%ATP, species1%population, cell_growth, cell_main)
                         IF (new_bugs < 1) THEN
                            new_bugs = 0
                         END IF

                         fit_level = fitness(atmosphere%current_T, ocean%H2, H2_Cmax, H2_lim, T_ideal, T_sens)

                         max_H2 = H2_Cmax * fit_level * (species1%population - bug_starve)


                         IF (max_H2 .GT. ocean%H2) THEN
                            max_H2 =  ocean%H2
                         END IF
 
                         H2_growth =  2 * protein_cell * new_bugs


                         IF (H2_growth .GT. max_H2) THEN
                            !!print *, "yes too many new bugs for amount of H2 we can eat"
                            !!new_bugs  = max_H2 / (2 * protein_cell * new_bugs * av_constant)
                            H2_growth =  max_H2
                            new_bugs  = (H2_growth / 2.0) / protein_cell
                            !!print *, "updated new bugs", new_bugs
                            !!H2_growth =  2 * protein_cell * new_bugs * av_constant
                            !!2 * protein_cell * new_bugs
                         END IF

!!$                         IF (H2_growth .GT. 2*ocean%CO2) THEN  !! do we have sufficient carbon dioxide for cells?
!!$                            H2_growth = 2*ocean%CO2
!!$                            new_bugs  = H2_growth / (2 * protein_cell * new_bugs)
!!$                         END IF
!!$
!!$                         IF (4*(max_H2 - H2_growth) .GT. (ocean%CO2 - H2_growth/2.0)) THEN
!!$                            !! function
!!$                         END IF
                         

                         ATP_used = new_bugs * cell_growth!! * av_constant !! atp used in creating new bugs
 
                         CO2_growth = protein_cell * new_bugs

                         IF (CO2_growth + ((max_H2 - H2_growth) / 4.0) .GT. ocean%CO2) THEN
                            PRINT *, "ERROR IN CO2 OCEAN"
                         END IF
                         

                         !!print *, new_bugs
                         ATP_made = (max_H2 - H2_growth) * ATP_H2 !! use rest of H2 for ATP production

                         ocean%H2              = ocean%H2       - max_H2
                         ocean%CO2             = ocean%CO2      - CO2_growth
                         ocean%CO2             = ocean%CO2      - ((max_H2 - H2_growth) / 4.0)
                         ocean%CH4             = ocean%CH4      + ((max_H2 - H2_growth) / 4.0)
                         biotic_CH4_output = (max_H2 - H2_growth) / 4.0

                         species1%ATP        = species1%ATP        + ATP_made - ATP_starve - ATP_used - ATP_maintain

                         species1%population = species1%population + new_bugs - bug_starve

                         IF (species1%population < 1 .OR. species1%ATP < 0) THEN
                            !!print *, "pop before", species1%population
                            species1%population = 0
                            species1%ATP        = 0
                            !!print *, "pop after", species1%population
                         END IF
                   END IF
    
                   !!write(20, *) t_step, atmosphere%current_T, species1%population, H2_to_add, &
                   !!atmosphere%H2, atmosphere%CO2, atmosphere%CH4, species1%ATP, ocean%H2, ocean%CO2, ocean%CH4
                END IF
             END DO
             !!END IF

             !!print *, "t_step", t_step
          END DO

          close(20)
          close(30)
          
          !!print *, "largest num poss", huge(0.0d0)
          !!print *, "largest num poss", huge(0.0)
          !!print *, "largest H2 poss",  huge(max_H2)
       contains

       function grid_area(lo,la, planet_r)
          real :: grid_area, lo, la, alpha, u_angle, l_angle
          real, parameter  :: pi = 3.1415927
          !!real, parameter  :: planet_r = 6051.8E3
          real :: planet_r
          real, parameter  :: num_longs = 144.0
          alpha = 2*pi*planet_r*planet_r
          u_angle = (pi/180)*(abs(la)+1)
          l_angle = (pi/180)*(abs(la)-1)
          grid_area = (alpha/num_longs) * (sin(u_angle) - sin(l_angle))
        end function grid_area

        function bugs_made(ATP, population, cell_growth, cell_main) result(new_bugs)
          real(16)            :: ATP, population, new_bugs, mu, sigma, d, a_r
          real(16)            :: cell_main, cell_growth
          a_r = (cell_growth + cell_main)
          !!print *, "in function pop", population
          !!print *, "in function ATP", ATP
          !!print *, "in function a_r", a_r
          mu = ATP / population
          !!print *, "in function mu", mu
          sigma = 0.1*cell_main*((mu/cell_main)**0.5) !!0.35 * mu!!0.3*cell_main !!1.0
          d = (a_r - mu) / sigma
          !!print *, "d", d
          new_bugs = 1 - 0.5 * (1 + ERF((a_r - mu)/(sigma * SQRT(2.0))))
          !!print *, "in function new bugs 111", new_bugs
          new_bugs = population * new_bugs
          !!print *, "in function new bugs", new_bugs
        end function bugs_made



        function ATP_of_starved(ATP, population, cell_main) result(ATP_starve)
          real(16) :: ATP, population, ATP_starve, mu, sigma, c
          real(16) :: cell_main
          real, parameter :: pi = 3.1415927
          mu = ATP / population
          !!print *, "mu", mu
          !!print *, "cell_main", cell_main
          !!print *, "in function ATP", ATP, "in function population", population, mu
          sigma = 0.1*cell_main*((mu/cell_main)**0.5) !!0.35 * mu!!0.3*cell_main !!1.0 !!mu / cell_main
          !!print *, "in function sigma", sigma
          c = (cell_main - mu) / sigma
          !!print *, "in function c", c
          ATP_starve = -1.0*(sigma / SQRT(2*pi) )* EXP(-0.5*c**2)
          !!print *, "in function ATP starve first bit", ATP_starve
          ATP_starve = ATP_starve + (mu/2.0)*(1+ERF(c*SQRT(0.5)))
          ATP_starve = ATP_starve * population
          IF (ATP_starve .LT. 0) THEN
             ATP_starve = 0
          ELSE IF (ATP_starve .GT. ATP) THEN
             ATP_starve = ATP
          END IF
          !!print *, "in function ATP starve", ATP_starve
          !!print *, "bug starve: ", bug_starve, "ATP starve: ", ATP_starve
        end function ATP_of_starved

        function num_bugs_starved(ATP, population, cell_main) result(bug_starve_out)
          real(16)          :: ATP, population, mu, sigma
          real(16)          :: bug_starve_out, cell_main
          mu = ATP / population
          !!print *, "mu", mu
          !!print *, "maintaining", cell_main
          sigma = 0.1*cell_main*((mu/cell_main)**0.5) !!0.35 * mu!!0.3*cell_main !!1.0 !!0.35 * mu
          bug_starve_out = population * 0.5 * (1 + ERF((cell_main - mu) / (SQRT(2.0) * sigma)))
          !!print *, "bug starve bit first", filler
          !!filler = 0.5 * (1 + filler)
          !!print *, "bug starve bit second", filler
          !!print *, "in function population", population
          !!filler = population * filler
          !!print *, "bug starve bit third", filler
          !!bug_starve_out = filler
          !!print *, "bug starve bit final", bug_starve_out
          
          IF (bug_starve_out < 0) THEN
             bug_starve_out = 0
          ELSE IF (bug_starve_out > population) THEN
             bug_starve_out = population
          END IF
          !!print *, "in function mu", mu, "cell main", cell_main*av_constant, "bugs starve", bug_starve_out
        end function num_bugs_starved
        
          
        function fitness(T, H2_pp, H2_Cmax, H2_lim, T_ideal, T_sens) result(fitness_val)
          real(16)   :: H2_pp
          real :: T, H2, H2_fit, T_fit, factor_i, fitness_val
          real ::  H2_lim, T_ideal, T_sens, H2_Cmax
          factor_i = T_sens * sqrt((T - T_ideal)**2)
          !!print *, "factor_i", factor_i
          T_fit = exp (-1.0 * (factor_i ** 2))
          !!print *, "T_fit", T_fit
          H2_fit = tanh(H2_pp - H2_lim)  !! in number of moles
          !!print *, "H2_fit", H2_fit
          IF (H2_fit < 0) THEN
             H2_fit = 0
          END IF
          fitness_val = H2_fit * T_fit !!T_fit 
          !!print *, "fitness_val", fitness_val
        end function fitness

        function CO2_array_func(CO2_MMR, CO2_array,t_array_length) result(temp_value_CO2)
          real(16)               :: CO2_MMR
          real, dimension(t_array_length) :: CO2_array
          real                   :: temp_value_CO2
          integer                :: i, t_array_length
          IF (CO2_MMR < CO2_array(1)) THEN
                temp_value_CO2 = 1
          ELSE IF (CO2_MMR > CO2_array(t_array_length)) THEN
                temp_value_CO2 = t_array_length
          END IF
          DO i = 1,  t_array_length - 1
             IF (CO2_array(i)  <= CO2_MMR .AND. CO2_MMR < CO2_array(i + 1)) THEN
                temp_value_CO2 = i
             END IF
          END DO
          If (temp_value_CO2 .LT. 1) THEN
             PRINT *, "CO2 in function error"
             PRINT *, CO2_MMR
             PRINT *, temp_value_CO2
          END IF
        end function CO2_array_func

        function CH4_array_func(CH4_MMR, CH4_array, t_array_length) result(temp_value_CH4)
          real(16)                 :: CH4_MMR
          !!integer                  :: max_val = 10001
          real, dimension(t_array_length)   :: CH4_array
          real                     :: temp_value_CH4
          integer                  :: i, t_array_length
          IF (CH4_MMR < CH4_array(1)) THEN
                temp_value_CH4 = 1
          ELSE IF (CH4_MMR > CH4_array(t_array_length)) THEN
                temp_value_CH4 = t_array_length
          END IF
          DO i = 1,  t_array_length-1
             IF (CH4_array(i)  <= CH4_MMR .AND. CH4_MMR < CH4_array(i + 1)) THEN
                temp_value_CH4 = i
             END IF
          END DO
          IF (temp_value_CH4 .LT. 1) THEN
             PRINT *, "CH4 in function error"
             PRINT *, CH4_MMR
             PRINT *, temp_value_CH4
          END IF
        end function CH4_array_func

        function gas_flux_atmo_ocean(num_moles, piston_velocity, solubility, &
          ocean_moles) result(gas_flux)
          real(16)                 :: num_moles
          real(16)                 :: ocean_moles
          real                     :: atmo_pressure = 1.0         !! atmospheric pressure in bar
          real                     :: piston_velocity
          real                     :: solubility
          real                     :: partial_pressure
          real                     :: dissolved_concentration
          real                     :: gas_flux
          real(16), parameter      :: atmo_moles = 1.73E20   !! total moles of air in atmosphere ( current)
          real(16), parameter      :: ocean_volume = 9.2E14  !! m^3
          dissolved_concentration = ocean_moles / ocean_volume
          partial_pressure = (num_moles * atmo_pressure) / atmo_moles
          gas_flux = piston_velocity * (solubility * partial_pressure - dissolved_concentration)
          !!print *, "dissolved concentration", dissolved_concentration, "partial pressure", partial_pressure, &
          !!     "gas flux", gas_flux
        end function gas_flux_atmo_ocean
        

        
!!$      function metabolism(pop,ATP,met,CO2,H2,CH4,T,timestep)
!!$          real :: fit, H_cons, ATP_req, ATP_def, ATP_current, pop_starve
!!$          real :: growth_possible, H2_for_growth, pop, CO2, H2, CH4, T, timestep
!!$          real, parameter    :: H2_Cmax = 3.76E-17   * (60*60)             !! max cons moles H2 possible per cell per second
!!$          real, parameter    :: cell_main = 2.16E-19 * (60*60)             !! mol_ATP cell^-1 s^-1
!!$          real, parameter    :: death_starve = 2.5E-7             !! s^-1
!!$          real, parameter    :: cell_growth = 2.36E13             !! cells mol_ATP^-1
!!$          real, parameter    :: cells_per_CH4 = 1416E10           !! cells per mole of CH4 produced
!!$          !!real, parameter  :: CH2O_per_cell = 7.4E-15           !! moles of CH2O per cell
!!$          real, parameter :: CH2O_per_cell = 8.912768470E9        !! atoms of CH20 per cell
!!$          real, parameter    :: a_constant = 6.02214076E23
!!$          !! first use ATP to maintain population
!!$          ATP_req = pop * cell_main
!!$          IF (ATP_current >= ATP_req) THEN
!!$             ATP_current = ATP_current - ATP_req
!!$          ELSE
!!$             ATP_def = ATP_req - ATP_current
!!$             pop_starve = ATP_def / cell_main
!!$             pop = pop - pop_starve * death_starve
!!$          END IF
!!$       
!!$          !! HOW MUCH H2 IS CONSUMED
!!$       
!!$          fit = fitness(T,H2_pp)
!!$          H_cons = fit * H2_Cmax * timestep * pop
!!$          IF (H_cons > H2) THEN
!!$             H_cons = H2
!!$          END IF
!!$          growth_possible = ATP_current * cell_growth
!!$          H2_for_growth = growth_possible * (CH2O_per_cell / a_constant) * 2
!!$          IF (H2_for_growth > H_cons) THEN
!!$             H2_for_growth = H_cons
!!$          END IF
!!$       
!!$           !! USE H2 FOR CELL GROWTH
!!$       
!!$          H_cons = H_cons - H2_for_growth
!!$          H2 = H2 - H2_for_growth
!!$          CO2 = CO2 - H2_for_growth/2.0
!!$          pop = pop + (H2_for_growth * a_constant / 2.0) / CH2O_per_cell
!!$          ATP_current = ATP_current - H2_for_growth/(cell_growth * (CH2O_per_cell / a_constant) * 2)
!!$       
!!$          !! BUG CHECK
!!$          IF (ATP_current < 0) THEN
!!$             print *, "ATP BUG!"
!!$             ATP_current = 0
!!$          END IF
!!$           
!!$          !! now use rest of H2 to create ATP
!!$          H2 = H2 - H_cons
!!$          CO2 = CO2 - H_cons / 4.0
!!$          CH4 = CH4 + H_cons / 4.0
!!$          ATP_current = ATP_current + (0.6 * H_cons) / 4.0 !! 0.6 moles of ATP per 4 moles of H2 used
!!$       
!!$      end function metabolism
!!$               

 end program testing
      
