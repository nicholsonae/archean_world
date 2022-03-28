
program archean_world

  implicit none
  
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
             real(16)                       :: current_T
             real(16)                       :: eq_T
             real(16)                       :: T_increment
          END TYPE Environment
          
          
          real, parameter                   :: pi = 3.1415927
          real, parameter                   :: av_constant = 6.02214076E23 !! avogadros constant
          real, parameter                   :: molar_mass_CH2O = 30.0      !! molar mass of CH2O

          real(16), parameter               :: H2_Cmax = (60*60) * 3.76E-17        !! mol_H2 cell^-1 s^-1
          real(16), parameter               :: CO2_Cmax = (60*60) * 0.5 * 3.76E-17        !! mol_H2 cell^-1 s^-1
          real, parameter                   :: energy_to_ATP = 32.5        !! kJ to make a mol of ATP

          real(16), parameter               :: cell_main = 1.0 * (60*60) * 2.16E-19  !! mol_ATP cell^-1 s^-1 -19
          real(16), parameter               :: protein_cell = 1.0 *7.4E-15           !! mol_CH2O cell^-1 
          real(16), parameter               :: cell_growth  = 1.0 * 4.237E-14        !! mole ATP to make a cell
          
          real, parameter                   :: T_ideal      = 283.0        !! microbe ideal T
          real, parameter                   :: T_sens       = 0.0          !! microbe sensivity to temperature
          real, parameter                   :: ATP_CH4      = 0.6          !! moles of ATP per moles of CH4 produced
          real, parameter                   :: ATP_H2       = 0.15
          real, parameter                   :: H2_lim       = 0.0          !! lower bound of H2 consumption
          real, parameter                   :: rand_dr      = 0.01         !! random death rate 1% currently
          real, parameter                   :: add_rand     = 0.0          !! add randomness to the model
          real, dimension(90)               :: lat                         !! lattitude points on globe
          real, dimension(144)              :: lon                         !! longitude points on globe
          integer, parameter                :: t_array_length = 101        !! length of temeperature array, was 10001
          
          real, parameter                   :: planet_r = 6051.8E3         !! radius of planet
          real, parameter                   :: ocean_depth = 100           !! ocean depth in m
          real, parameter                   :: atmosphere_mass = 5.15E18   !! mass of the atmosphere
          real(16), parameter               :: moles_air = 1.73E20 !!* 0.5 !! moles of air in atmosphere
          real, parameter                   :: atmo_pressure = 1.013       !! atmopsheric pressure in bar
          real(16)                          :: ocean_volume                !! volume ocean in m^3
          real(16)                          :: ocean_surf_area             !! ocean surface area m^2
          
          real, parameter                   :: piston_vel_CO2 = 6.7E-4     !! piston velocity CO2 at 25 degrees C -4
          real, parameter                   :: solubility_CO2 = 0.035      !! solubility of CO2   mol L−1 bar−1
          real, parameter                   :: piston_vel_CH4 = 4.5E-5     !! piston velocity of CH4 at 25 degrees C
          real, parameter                   :: solubility_CH4 = 1.4E-3     !! solubility of CH4  mol L−1 bar−1
          real, parameter                   :: piston_vel_H2  = 1.3E-4     !! piston velocity for H2 at 25 degrees
          real, parameter                   :: solubility_H2  = 7.8E-4     !! solubility of H2 mol L^-1 bar^-1
          real, parameter                   :: diffusion_H2   = 4.8E-9     !! m^2 s^-1
          real, parameter                   :: diffusion_CO2  = 1E-9       !! for DIC m^2 s^-1

          real, parameter                   :: CO2_burial       = 0.001    !! % of CO2 removed from atmosphere per year
          real, parameter                   :: CH4_burial       = 0.001    !! % of CH4 removed from atmosphere per year
          real, parameter                   :: H2_burial        = 0.001    !! % of H2 removed from atmosphere per year
          real, parameter                   :: CO2_mole_weight  = 44.0095  !! molar mass of CO2
          real, parameter                   :: CH4_mole_weight  = 16.043   !! molar mass of CH4
          real, parameter                   :: H2_mole_weight   = 2.01588  !! molar mass of H2
          real, parameter                   :: N2_mole_weight   = 28.0134  !! molar mass of N2
          !!real, parameter                   :: pN2              = 100000.0 !! pressure of N2
          
          real                              :: CO2_outflux_year = 10E15    !! outflux in moles CO2 20x as much as current was 5.5E 
          real                              :: H2_outflux_year  = 10E13     !! outflux in moles H2  was 10E13
          
          real(16)                          :: H2_pp, CH4_pp, CO2_pp, T    !! environment
          real(16)                          :: c, d, mu, omega
          real(16)                          :: new_bugs, bug_starve, bugs_death, bugs_culled, tot_died
          real(16)                          :: fit_level = -1
          real(16)                          :: a_r, ATP_used, ATP_made
          real(16)                          :: ATP_maintain
          real(16)                          :: max_H2, max_CO2, CO2_maint, H2_maint
          real(16)                          :: ATP_starve
          real(16)                          :: CO2_MMR, CH4_MMR
          real(16)                          :: CO2_growth, H2_growth
          real(16)                          :: towrite_CO2, towrite_CH4
          real(16)                          :: starve_random, repro_random, death_random
          real(16)                          :: CO2_atmo_2_ocean
          real(16)                          :: H2_atmo_2_ocean
          real(16)                          :: CH4_atmo_2_ocean
          real(16)                          :: biotic_CH4_output
          real(16)                          :: born_timestep, died_timestep
          real(16)                          :: H2_out, total_C_eat, C_growth, growth_eff, current_G

          real(16)                          :: pp_CO2
          real(16)                          :: pp_CH4
          real(16)                          :: pp_H2
          real(16)                          :: mmCO2, mmCH4, mmN2
         
          TYPE(Microbe)                     :: species1
          TYPE(Microbe), dimension(1)       :: species_list
          
          TYPE(Environment)                 :: atmosphere
          TYPE(Environment)                 :: ocean
          
          real, dimension(t_array_length,t_array_length)      :: temp_array
          real, dimension(t_array_length)                     :: CO2_array, CH4_array
          integer                           :: i, j, n, t_step, IOStatus
          integer                           :: biotic_step
          integer                           :: d_step

          integer                           :: temp_value_CH4
          integer                           :: temp_value_CO2
          integer                           :: seeded = 0
          integer                           :: habitable = 1                 !! CHANGE THIS
          integer                           :: step_length = 365*24          !! hours in a year

          real(16)                          :: H2_to_add, CO2_to_add, CH4_to_add
          real(16)                          :: H2_ocean_to_add, CO2_ocean_to_add, CH4_ocean_to_add
          !!real                              :: t_bugs_starved, t_bugs_born
          integer                           :: seed_input
          integer                           :: fix_reseed
          character(100)                    :: num1char, num2char, file_number

          CALL GET_COMMAND_ARGUMENT(1,num1char)   !first, read in the two values
          !!print *, num1char
          CALL GET_COMMAND_ARGUMENT(2,file_number)
          !!print *, file_number
          READ(num1char,*)seed_input                    !then, convert them to REALs
          !!print *, seed_input
          !!READ(num2char,*)file_number


          call random_seed(put=[seed_input,seed_input,seed_input,seed_input,seed_input,seed_input,seed_input,seed_input])


          open(10, file="temperature_file2.txt", status='old', action='read')
          DO i = 1, t_array_length
             read(10,*) (temp_array(j,i), j = 1,t_array_length)
          END DO
          close(10)

          open(20, file="data_abstract_"//trim(file_number)//".txt", status='replace', action='write')

          !!open(30, file="data_random_"//trim(file_number)//".txt", status='replace', action='write')

          DO i = 1, t_array_length
             CO2_array(i) = 0.005 + 0.00095  * (i-1)
             CH4_array(i) = 0.0001 * (i - 1)
          END DO

          ocean_volume = (4/3.0) * pi * ((planet_r + ocean_depth)**3 - planet_r**3 )
          ocean_surf_area = 4.0 * pi * (planet_r + ocean_depth)**2


          
          species1%id           = 0
          species1%population   = 0
          species1%ATP          = 0
          species1%metabolism   = 0

          born_timestep = 0
          died_timestep = 0

          species_list(1)       = species1

          fix_reseed            = 0
          C_growth = 0.0
          total_C_eat = 0.0

          CO2_MMR               = 0.0
          CH4_MMR               = 0.0

          atmosphere%CO2          = 100
          atmosphere%CH4          = 0.0
          atmosphere%H2           = 100
          atmosphere%current_T    = interpolate_temp(CO2_MMR, CH4_MMR, t_array_length, temp_array,&
               CH4_array, CO2_array)
          atmosphere%eq_T         = interpolate_temp(CO2_MMR, CH4_MMR, t_array_length, temp_array,&
               CH4_array, CO2_array)
          atmosphere%T_increment  = 0


          ocean%CO2               = 0.0
          ocean%CH4               = 0.0
          ocean%H2                = 0.0
          ocean%current_T         = interpolate_temp(CO2_MMR, CH4_MMR, t_array_length, temp_array,&
               CH4_array, CO2_array)
          ocean%eq_T              = interpolate_temp(CO2_MMR, CH4_MMR, t_array_length, temp_array,&
               CH4_array, CO2_array)
          ocean%T_increment       = 0

          biotic_CH4_output = 0
          H2_out = 0.0

          lat  = (/ (I, I=-89, 89, 2) /)
          DO i = 1, 144
             lon(i) = 1.25 + 2.5*(i-1)
          END DO

          current_G =  delta_G(atmosphere%current_T, ocean%H2, ocean%CO2, &
                              ocean%CH4, solubility_H2, solubility_CO2, solubility_CH4, ocean_volume)

          DO t_step = 0, 40000  !! run the simulation for 1 year then update the temperature

             IF (total_C_eat .LE. 0) THEN
                growth_eff = -1
             ELSE
                growth_eff = C_growth / total_C_eat
             END IF

             write(20, *) t_step, atmosphere%current_T, species1%population, species1%ATP, &
                  H2_to_add, atmosphere%H2/moles_air, atmosphere%CO2/moles_air, &
                  atmosphere%CH4/moles_air, ocean%H2/ocean_volume, &
                  ocean%CO2/ocean_volume, ocean%CH4/ocean_volume, biotic_CH4_output, fit_level, &
                  new_bugs, died_timestep, H2_out, growth_eff, current_G


             born_timestep = 0
             died_timestep = 0
             biotic_CH4_output = 0
             H2_out = 0
             C_growth = 0.0
             total_C_eat = 0.0

             atmosphere%CO2 = (1 - CO2_burial) * (atmosphere%CO2 + CH4_burial *  atmosphere%CH4 + CO2_outflux_year)
             atmosphere%H2 = atmosphere%H2 + 4 * CH4_burial * atmosphere%CH4 + H2_outflux_year
             atmosphere%CH4 = (1 - CH4_burial) * atmosphere%CH4
             H2_out = H2_burial * (atmosphere%H2 + 2 * atmosphere%CH4)
             atmosphere%H2 =  atmosphere%H2 - H2_burial * (atmosphere%H2 + 2 * atmosphere%CH4)

             mmCO2 = atmosphere%CO2 * CO2_mole_weight
             mmCH4 = atmosphere%CH4 * CH4_mole_weight
             mmN2  = (moles_air - atmosphere%CO2 - atmosphere%CH4)*N2_mole_weight

             CO2_MMR        =  mmCO2 / (mmCO2 + mmCH4 + mmN2)
             CH4_MMR        =  mmCH4 / (mmCO2 + mmCH4 + mmN2)

             atmosphere%eq_T = interpolate_temp(CO2_MMR, CH4_MMR, t_array_length, temp_array, &
               CH4_array, CO2_array)
             atmosphere%T_increment = (atmosphere%eq_T - atmosphere%current_T) / step_length

             atmosphere%current_T =  interpolate_temp(CO2_MMR, CH4_MMR, t_array_length, temp_array, &
               CH4_array, CO2_array)

             IF ((seeded .EQ. 0) .AND. (T_ideal - 3 < atmosphere%current_T) .AND. (atmosphere%current_T < T_ideal + 3)) THEN
                habitable = 1
             END IF
             towrite_CO2 = atmosphere%CO2 / moles_air
             towrite_CH4 = atmosphere%CH4 / moles_air
             
             IF (habitable .EQ. 1 .AND. t_step > 20000 .AND. seeded .EQ. 0) THEN !! seed if appropiate
                species1%id           = 0
                species1%population   = 1E2
                species1%ATP          = 2.0 * cell_main * 1E2 + 2.16E-19
                species1%metabolism   = 0
                seeded                = 1
                fix_reseed            = 1
             END IF

             IF (fix_reseed .EQ. 1 .AND. species1%population .EQ. 0) THEN
                species1%id           = 0
                species1%population   = 1E2
                species1%ATP          = 2.0 * cell_main * 1E2 + 2.16E-19
                species1%metabolism   = 0
                fix_reseed = 0
             END IF
                
             
             DO biotic_step = 1, step_length !! update biology for 1 year in hour steps

                CO2_atmo_2_ocean =  gas_flux_atmo_ocean(atmosphere%CO2, piston_vel_CO2, solubility_CO2, ocean%CO2, &
                     atmo_pressure, moles_air, ocean_volume)
                CH4_atmo_2_ocean =  gas_flux_atmo_ocean(atmosphere%CH4, piston_vel_CH4, solubility_CH4, ocean%CH4, &
                     atmo_pressure, moles_air, ocean_volume)
                H2_atmo_2_ocean  =  gas_flux_atmo_ocean(atmosphere%H2,  piston_vel_H2, solubility_H2,   ocean%H2, &
                     atmo_pressure, moles_air, ocean_volume)

                CO2_ocean_to_add = 60*60* CO2_atmo_2_ocean * ocean_surf_area  !! the amount that enters the ocean each hour
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

                IF (species1%population .GT. 0) THEN
                   
                   H2_out = H2_out + 2.0 *  protein_cell * species1%population * rand_dr !! removal of H2 via microbe death
                   species1%population = species1%population * (1-rand_dr) !!* death_random
                   species1%ATP        = species1%ATP        * (1-rand_dr) !!* death_random 

                   IF (species1%population < 1.0) THEN
                      species1%population = 0
                      species1%ATP        = 0
                   ELSE

                      bug_starve = num_bugs_starved(species1%ATP, species1%population, cell_main)
                      died_timestep = died_timestep + bug_starve

                      ATP_starve = ATP_of_starved(species1%ATP, species1%population, cell_main)

                      ATP_maintain = cell_main * (species1%population - bug_starve)

                      new_bugs = bugs_made(species1%ATP, species1%population, cell_growth, cell_main)
                      
                      fit_level = fitness(atmosphere%current_T, ocean%H2, H2_lim, T_ideal, T_sens)

                      max_H2 = H2_Cmax * fit_level * (species1%population - bug_starve)

                      max_CO2 = CO2_Cmax * fit_level * (species1%population - bug_starve)

                      IF (max_H2 .GT. ocean%H2) THEN
                         max_H2 =  ocean%H2
                      END IF

                      IF (max_CO2 .GT. ocean%CO2) THEN
                         max_CO2 =  ocean%CO2
                      END IF
 
                      H2_growth  =  2 * protein_cell * new_bugs
                      CO2_growth =  protein_cell * new_bugs

                      IF (H2_growth .GT. max_H2) THEN
                         H2_growth =  max_H2
                         new_bugs  = (H2_growth / 2.0) / protein_cell
                         CO2_growth =  protein_cell * new_bugs
                      END IF

                      IF (CO2_growth .GT. max_CO2) THEN
                         CO2_growth =  max_CO2
                         new_bugs  = CO2_growth / protein_cell
                         H2_growth = 2.0 * protein_cell * new_bugs
                      END IF
                         
                      born_timestep = born_timestep + new_bugs

                      ATP_used = new_bugs * cell_growth !! atp used in creating new bugs

                      CO2_maint = max_CO2 - CO2_growth
                      H2_maint  = max_H2  - H2_growth

                      IF (H2_maint .GE. 4.0 * CO2_maint) THEN
                         H2_maint = CO2_maint * 4.0
                      ELSE
                         CO2_maint = H2_maint / 4.0
                      END IF

                      ATP_made = ATP_H2 * H2_maint

                      ocean%CO2             = ocean%CO2         - CO2_growth - CO2_maint
                      ocean%CH4             = ocean%CH4         + CO2_maint
                      biotic_CH4_output     = biotic_CH4_output + CO2_maint
                      ocean%H2              = ocean%H2 -  H2_growth - H2_maint

                      C_growth = C_growth + CO2_growth
                      total_C_eat = total_C_eat + CO2_growth + CO2_maint
                         
                      IF (ocean%H2 < 0.0) THEN
                         !!ocean%H2 = 0
                         PRINT *, "ocean H2 error in bio loop"
                      END IF

                      species1%ATP          = species1%ATP + ATP_made - ATP_starve - ATP_used &
                           - ATP_maintain

                      species1%population   = species1%population + new_bugs - bug_starve

                      H2_out = H2_out + 2.0 * bug_starve * protein_cell

                      IF (species1%population < 1 .OR. species1%ATP < 0) THEN
                         species1%population = 0
                         species1%ATP        = 0
                      END IF
                   END IF
    
                END IF
             END DO

          END DO

          close(20)
          !!close(30)

       contains

       function grid_area(lo,la, planet_r)
          real :: grid_area, lo, la, alpha, u_angle, l_angle
          real, parameter  :: pi = 3.1415927
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
          mu = ATP / population
          sigma = 0.1*mu 
          d = (a_r - mu) / sigma
          new_bugs = 1 - 0.5 * (1 + ERF((a_r - mu)/(sigma * SQRT(2.0))))
          new_bugs = population * new_bugs
          IF (new_bugs .LE. 1) THEN
             new_bugs = 0
          END IF
        end function bugs_made



        function ATP_of_starved(ATP, population, cell_main) result(ATP_starve)
          real(16) :: ATP, population, ATP_starve, mu, sigma, c
          real(16) :: cell_main
          real, parameter :: pi = 3.1415927
          mu = ATP / population
          sigma = 0.1*mu 
          c = (cell_main - mu) / sigma
          ATP_starve = -1.0*(sigma / SQRT(2*pi) )* EXP(-0.5*c**2)
          ATP_starve = ATP_starve + (mu/2.0)*(1+ERF(c*SQRT(0.5)))
          ATP_starve = ATP_starve * population
          IF (ATP_starve .LT. 0) THEN
             ATP_starve = 0
          ELSE IF (ATP_starve .GT. ATP) THEN
             ATP_starve = ATP
          END IF
        end function ATP_of_starved


        function num_bugs_starved(ATP, population, cell_main) result(bug_starve_out)
          real(16)          :: ATP, population, mu, sigma
          real(16)          :: bug_starve_out, cell_main
          mu = ATP / population
          sigma = 0.1*mu
          bug_starve_out = population * 0.5 * (1 + ERF((cell_main - mu) / (SQRT(2.0) * sigma)))
          IF (bug_starve_out < 0) THEN
             bug_starve_out = 0
          ELSE IF (bug_starve_out > population) THEN
             bug_starve_out = population
          END IF
        end function num_bugs_starved
        
          
        function fitness(T, H2_pp, H2_lim, T_ideal, T_sens) result(fitness_val)
          real(16)   :: H2_pp
          real(16) :: H2, H2_fit, T_fit, factor_i, fitness_val
          real(16) :: T, H2_Cmax
          real ::  H2_lim, T_ideal, T_sens
          factor_i = T_sens * sqrt((T - T_ideal)**2)
          T_fit = exp (-1.0 * (factor_i ** 2))
          H2_fit = tanh(H2_pp - H2_lim)  !! in number of moles
          IF (H2_fit < 0) THEN
             H2_fit = 0
          END IF
          fitness_val = H2_fit * T_fit  
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
             print *, "at CO2 limit"
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
             print *, "at CH4 limit"
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
             ocean_moles, atmo_pressure, atmo_moles, ocean_volume) result(gas_flux)
          real, parameter          :: av_constant = 6.02214076E23
          real(16)                 :: num_moles
          real(16)                 :: ocean_moles
          real                     :: atmo_pressure         !! atmospheric pressure in bar
          real                     :: piston_velocity
          real                     :: solubility
          real(16)                 :: partial_pressure
          real(16)                 :: dissolved_concentration
          real(16)                 :: gas_flux
          real(16)                 :: atmo_moles   !! total moles of air in atmosphere ( current)
          real(16)                 :: ocean_volume  !! m^3
          dissolved_concentration = ocean_moles / ocean_volume  !! dissolved L-1
          partial_pressure = (num_moles / atmo_moles) * atmo_pressure !! in bar
          gas_flux = piston_velocity * (1000 * solubility * partial_pressure - dissolved_concentration)

        end function gas_flux_atmo_ocean

        function delta_G(T, ocean_H2, ocean_CO2, ocean_CH4, solubility_H2, solubility_CO2, solubility_CH4, &
             ocean_volume) result(G)
          real, parameter  :: R = 0.008314   !! gas constant in kJ mol^-1 K-1
          real(16)         :: T
          real(16)         :: ocean_H2
          real(16)         :: ocean_CO2
          real(16)         :: ocean_CH4
          real             :: solubility_H2
          real             :: solubility_CO2
          real             :: solubility_CH4
          real(16)         :: ocean_volume
          real(16)         :: Q
          real(16)         :: G
          real(16)         :: q_CH4
          real(16)         :: q_H2
          real(16)         :: q_CO2

          q_CH4 = (ocean_CH4 / ocean_volume) / (1000*solubility_CH4)
          q_H2  = (ocean_H2  / ocean_volume) / (1000*solubility_H2 )
          q_CO2 = (ocean_CO2 / ocean_volume) / (1000*solubility_CO2)

          IF ((q_H2 .EQ. 0) .OR. (q_CO2 .EQ. 0)) THEN
             Q = 1E30
          ELSE IF (q_CH4 .EQ. 0) THEN
             Q = 1E-30
          ELSE
             Q = q_CH4 / (q_CO2 * q_H2**4 )
          END IF

          !!T = 298
          
          G = -253 + 0.41*T + R*T*LOG(Q)
          G = -1.0 * G
          
        end function delta_G   

        function interpolate_temp(CO2MMR, CH4MMR, t_array_length, temp_array, CH4_array, CO2_array) &
             result(temperature)
          real, dimension(t_array_length)   :: CH4_array, CO2_array
          real(16)     :: CO2MMR, CH4MMR, temperature
          integer      :: t_array_length, i, CO2_low_i, CO2_high_i, CH4_low_i, CH4_high_i
          real(16)     :: CO2_low, CO2_high, CH4_low, CH4_high
          real, dimension(t_array_length,t_array_length)      :: temp_array
          real(16)     :: Q11, Q12, Q21, Q22  !! sides of our square
          real(16)     :: fxy1, fxy2      !! used for the bilinear interpolation

          IF (CH4_MMR < CH4_array(1)) THEN
             CH4_low_i  = 1
             CH4_high_i = 1
             CH4_low = CH4_array(1)
             CH4_high = CH4_array(1)
          ELSE IF (CH4_MMR > CH4_array(t_array_length)) THEN

             CH4_low_i  = t_array_length
             CH4_high_i = t_array_length
             CH4_low = CH4_array(t_array_length)
             CH4_high = CH4_array(t_array_length)
             print *, "at CH4 limit"
          ELSE
             DO i = 1,  t_array_length-1
                IF (CH4_array(i)  <= CH4MMR .AND. CH4MMR < CH4_array(i + 1)) THEN

                   CH4_low_i = i
                   CH4_high_i = i+1
                   CH4_low = CH4_array(i)
                   CH4_high = CH4_array(i+1)
                END IF
             END DO
          END IF

          IF (CO2_MMR < CO2_array(1)) THEN
             CO2_low_i  = 1
             CO2_high_i = 1
             CO2_low = CO2_array(1)
             CO2_high = CO2_array(1)
          ELSE IF (CO2_MMR > CO2_array(t_array_length)) THEN
             CO2_low_i  = t_array_length
             CO2_high_i = t_array_length
             CO2_low = CO2_array(t_array_length)
             CO2_high = CO2_array(t_array_length)
             print *, "at CO2 limit"
          ELSE
             DO i = 1,  t_array_length-1
                IF (CO2_array(i)  <= CO2MMR .AND. CO2MMR < CO2_array(i + 1)) THEN
                   CO2_low_i = i
                   CO2_high_i = i+1
                   CO2_low = CO2_array(i)
                   CO2_high = CO2_array(i+1)
                END IF
             END DO
          END IF

          Q11 = temp_array(CO2_low_i, CH4_low_i)
          Q12 = temp_array(CO2_low_i, CH4_high_i)
          Q21 = temp_array(CO2_high_i, CH4_low_i)
          Q22 = temp_array(CO2_high_i, CH4_high_i)

          
          IF (CO2_low_i .EQ. CO2_high_i) THEN
              temperature = Q11 + Q12*(CH4MMR - CH4_low)/(CH4_high - CH4_low)
          ELSE IF (CH4_low_i .EQ. CH4_high_i) THEN
              temperature = Q11 + Q21*(CO2MMR - CO2_low)/(CO2_high - CO2_low)
          ELSE
             fxy1 = Q11*(CO2_high - CO2MMR)/(CO2_high - CO2_low) &
                  + Q21*(CO2MMR - CO2_low)/(CO2_high - CO2_low)
             fxy2 = Q12*(CO2_high - CO2MMR)/(CO2_high - CO2_low) &
                  + Q22*(CO2MMR - CO2_low)/(CO2_high - CO2_low)
             temperature =  fxy1*(CH4_high -CH4MMR)/(CH4_high - CH4_low) &
                  + fxy2*(CH4MMR - CH4_low)/(CH4_high - CH4_low)
          END IF

        end function interpolate_temp
        
  
 end program archean_world
      
