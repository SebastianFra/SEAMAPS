using JuMP, Gurobi, Ipopt, CSV, DataFrames, DelimitedFiles

#Scenarios
#define Market Based Measure
MBM = "CO2_tax" #CO2_tax, CO2_tax_flat, Global_cap, Baseline
#Emission Accounting standard
Settings = "Infra"#WTW,Infra,TTW, GWP20 (Global Warming Potential 20 Years)
#Taxation of Emissions
Taxed = "WTW"#WTW,TTW
#Sensitivity Parameters
#green push - decrease green fuel cost by 25%
green_push = "FALSE"
#green decay - increase green fuel cost by 25%
green_decay = "FALSE"
#green push - decrease green fuel WTTemissions by 25%
green_push_em = "FALSE"
#green decay - increase green fuel WTT emissions by 25%
green_decay_em = "FALSE"
#input data file
Data_file_name = "Input_data_SEAMAPS.xls"
#define your local folder
Main_folder = "C:/Users/semfr/Documents/SEAMAPS/SEAMAPS/data"
#define your scenario data
Multiple_Runs = "TRUE" #FALSE, TRUE
if Multiple_Runs == "TRUE"
    Scenarios_file_name = "Scenarios.xls"
    Scenarios_set = "All_scenarios"
elseif Multiple_Runs != "FALSE"
    Scenarios_file_name = "Scenarios_multiple_runs.xls"
    Scenarios_set = "All_scenarios"
end
include("Import_scenarios.jl");
#define ambition towards climate mitigation actions (always stay at Realisic unless different implementation)
Ambition = "Realistic"
#define the name of your output folder in your local folder
Bottleneck_tech = ["Electrolyser", "CCS"] ; BT = length(Bottleneck_tech)
scale ="Conventional" #Conventional, Unconventional
dynamic_learning="TRUE"
if scale == "Unconventional"
Growth_rate = [2 2]
elseif scale == "Conventional"
Growth_rate = [1.5 1.5]
end
Cap_Init = [5 1.8]
Cap_0 = [0 0]
All_results_folder = "Results_" * Settings * "_" * MBM * "_" * Ambition * "_" * Taxed * "_" * scale
#define limitations on ramping up production facilities to be enabled or not
if Ambition == "Realistic"
    Option_ramping_fuels = "TRUE"
    Option_ramping_ships = "TRUE"
elseif Ambition != "Realistic"
    Option_ramping_fuels = "FALSE"
    Option_ramping_ships = "FALSE"
end

#prepare scenarios
#MBM differentiation
decommission_value = 0.5
Y = 2050 - 2020 + 1

# Scenario under study (all between N_scen_0 and N_scen_end)
N_scen_0 = 1;
N_scen_end = N_scenarios; # N_scenarios = total number of scenarios

NScen = N_scen_0

while NScen < N_scen_end + 1 #Run the optimization model for all scenarios

include("SEAMAPS_import.jl")
#To add in the excel sheet and import file




#Model
Shipping_stock = Model(Gurobi.Optimizer)


#variables
@variable(Shipping_stock, x[1:S, 1:Y] >= 0) #number of ships bought per year [# of av ships]
@variable(Shipping_stock, q[1:S, 1:Y] >= 0) #Total ship stock at year Y [# of av ships]
@variable(Shipping_stock, z[1:F, 1:S, 1:Y] >= 0) #Fuel demand per fueltype, ship, and year [PJ]
@variable(Shipping_stock, TotFuel[1:F,1:Y] >= 0) #Total fuel demand per fueltype per year [PJ]
@constraint(Shipping_stock, [f=1:F,y=1:Y], TotFuel[f,y] == sum(z[f, s, y] for s=1:S))
@variable(Shipping_stock, d[1:S, 1:Y] >= 0) #ship stock to decomission new ships [# of av ships]
if MBM != "ETS"
elseif MBM == "ETS"
    @variable(Shipping_stock, st[1:Y] >= 0) #stored CO2 allowances each year (in EUR/tCO2eq)
end
@variable(Shipping_stock, EF[1:S, 1:Y] >= 0) # Stock of ships from the existing fleet [# of av ships]
@variable(Shipping_stock, NB[1:S, 1:Y] >= 0) # Stock of new build ships [# of av ships]
@variable(Shipping_stock, d_EF[1:S, 1:Y] >= 0) #ship stock to decomission existing ships [# of av ships]
@variable(Shipping_stock, Bio_used_MeOH[1:Y] >= 0) #Biomass used to produce e_methanol [PJ]
@variable(Shipping_stock, Bio_used_PO[1:Y] >= 0) #Biomass used to produce PO [PJ]
@variable(Shipping_stock, Cap_bt[1:BT,1:Y]) # Production capacities of fuel bottleneck technologies [Any unit is possible MW installed, Mt/y, ...]
#@variable(Shipping_stock, Inv_dec_f_1[1:F, 1:S, 1:Y],Bin)# Investment decision for dynamic learning curves
#@variable(Shipping_stock, Inv_dec_f_2[1:F, 1:S, 1:Y],Bin)# Investment decision for dynamic learning curves
@variable(Shipping_stock, Inv_dec[1:BT,1:Y],Bin)
#@variable(Shipping_stock, Inv_dec_f_0[1:F, 1:S, 1:Y],Bin) # Investment decision for production capacities of bottleneck technology (driven by costs and demand)
@variable(Shipping_stock, fuel_cost_dynamic[f=1:F,y=1:Y]>= 0) #Fuel cost, ship, and year [PJ]


if MBM != "ETS"
    @objective(
        Shipping_stock,
        Min,
        sum(
            ship_inv[s] * x[s, y] + ship_OM[s] * q[s, y] for s = 1:S,
            y = 1:Y
        ) +
        sum(
            z[f, s, y] * (fuel_cost_dynamic[f, y] + fuel_tax[f, y]) for s = 1:S,
            f = 1:F, y = 1:Y
        ) +
        sum(
            d_EF[s, y] * ship_inv[s] * decommission_value for s = 1:S,
            y = 1:Y
        )
    )

elseif MBM == "ETS"
    @objective(
        Shipping_stock,
        Min,
        sum(
            ship_inv[s] * x[s, y] + ship_OM[s] * q[s, y] for s = 1:S,
            y = 1:Y
        ) +
        sum(
            z[f, s, y] * (fuel_cost[f, y] + fuel_tax[f, y]) for s = 1:S,
            f = 1:F, y = 1:Y
        ) +
        sum(
            d_EF[s, y] * ship_inv[s] * decommission_value for s = 1:S,
            y = 1:Y
        ) - sum(st[y] for y = 1:Y)
    )
end

#----------------------Ships related constraints---------------------

#Maximum amount of ship that can be produced per year
@constraint(
    Shipping_stock,
    [s = 1:S, y = 1:Y],
    x[s, y] <= ship_prod_cap[s, y]
)

#Stock of existing fleet
@constraint(
    Shipping_stock,
    [s = 1:S, y = 1:Y],
    EF[s, y] == (
        y > 1 ? EF[s, y-1] - d_EF[s, y] :
        preexisting_fleet[y, s] - d_EF[s, 1]
    )
)
@constraint(
    Shipping_stock,
    [s = 1:S, y = 1:Y],
    EF[s, y] <= preexisting_fleet[y, s]
)
#Stock of new build fleet
@constraint(
    Shipping_stock,
    [s = 1:S, y = 1:Y],
    NB[s, y] ==
    x[s, y] + (y > 1 ? sum(x[s, y] for y = 1:y-1) : 0) -
    (y > ship_lifetime[s] ? sum(x[s, y] for y = 1:y-ship_lifetime[s]) : 0) -
    (y > 1 ? sum(d[s, y] for y = 1:y-1) : d[s, 1])
)
#Total stock of ships
@constraint(
    Shipping_stock,
    [s = 1:S, y = 1:Y],
    q[s, y] == EF[s, y] + NB[s, y]
)

if dynamic_learning="TRUE":
#------------------------Learning Curves constraint--------------------------
# Set the initial parameter values and thresholds
thresholds = [0, 100, 250, 500, 1000, 2000, 3500, 5000, 7000, 10000]
param_values_exogenous = [1,0.8,0.6,0.5,0.45,0.425,0.4,0.375,0.35,0.3]
#set solver to solve NonConvex problems
set_optimizer_attribute(Shipping_stock, "NonConvex", 2)
#set solver MIPGap 
set_optimizer_attribute(Shipping_stock,"MIPGap", 0.015)
# Set the initial cost
fuel_cost_0 = [10.17 10.91 10.91 8.77 8.09 10.50 950.51 41.31 51.99 61.77 10.66 29.02 45.90 24.00 40.00 400.00]
#add variable for the learning parameter, which depends on the respective fuel used
@variable(Shipping_stock, learning_parameter[1:F,1:Y])
J = length(thresholds)
# Create binary variables to indicate which threshold has been reached in each year
@variable(Shipping_stock, b[1:F,1:Y,1:J],Bin)
# Add constraints to active the respective binary for the threshold reached
@constraint(Shipping_stock, [f = 1:F,j=1:J,y=1:Y-1],TotFuel[f,y+1]>=thresholds[j]*b[f,y,j])
# Add constraints to ensure only one threshold is triggered
@constraint(Shipping_stock, [f = 1:F,y=1:Y], sum(b[f,y,j] for j in 1:J)== 1)
# Add a constraint to update the parameter value based on the threshold that has been reached in each year
@constraint(Shipping_stock, [f = 1:F,y=1:Y],learning_parameter[f,y] == sum(b[f,y,j]* param_values_exogenous[j] for j in 1:J))  
# Add a constraint to calibrate the initial fuel cost to our 2020 estimates
@constraint(Shipping_stock,[f = 1:F,y=1], fuel_cost_dynamic[f, y] == fuel_cost_0[f])  
@constraint(Shipping_stock,[f = 1:F,y=2:Y], fuel_cost_dynamic[f, y] == fuel_cost_0[f]-(fuel_cost_0[f]*(1-learning_parameter[f,y])))  

else if dynamic_learning != "TRUE":

#------------------------Demand constraint--------------------------
#Transport demand must be satisfied
@constraint(
    Shipping_stock,
    [t = 1:Ty, y = 1:Y],
    sum(
        q[s, y] * ship_type_relation[t, s] * average_transport_work[y, s]
        for s = 1:S
    ) == Ship_Demands[y, t]
)

#---------------------Fuel related constraints------------------
#Fuel consumption associated with fulfilling demand
@constraint(
    Shipping_stock,
    [s = 1:S, y = 1:Y],
    sum(z[f, s, y] * ship_eff[f, s] for f = 1:F) ==
    q[s, y] * average_transport_work[y, s]
)
@constraint(
    Shipping_stock,
    [f = 1:F, y = 1:Y],
    TotFuel[f,y] <= fuel_available[y, f]
)
#*********Availability constraint for fuels sharing the same resource with exogenous global availability like biomass (Bioemethanol and refined PO)
@constraint(
    Shipping_stock,
    [y = 1:Y],
    Bio_used_MeOH[y] + Bio_used_PO[y] <= Shared_biomass_available[y]
)
#Quantity of fuel produced using biomass
@constraint(
    Shipping_stock,
    [f in BioMeOH, y = 1:Y],
    TotFuel[f,y] <= Bio_used_MeOH[y] * ConvBio_MeOH[y]
)
@constraint(
    Shipping_stock,
    [f in RefPO, y = 1:Y],
    TotFuel[f,y] <= Bio_used_PO[y] * ConvBio_PO[y]
)

#*************Bottleneck technologies ramping up constraints**************

if Option_ramping_fuels == "TRUE"

    @constraint(Shipping_stock,[bt=1:BT], Cap_bt[bt,1] == Cap_0[bt])
    @constraint(Shipping_stock,[bt=1:BT,y=1:Y-1], Cap_bt[bt,y+1] >= Cap_Init[bt]*Inv_dec[bt,y]) #Start having capacities after 1 year of investment decision
    @constraint(Shipping_stock,[bt=1:BT,y=1:Y-1], Cap_bt[bt,y+1] <= Cap_bt[bt,y]*Growth_rate[bt] + Cap_Init[bt]*Inv_dec[bt,y]) #New capacities limited by growth rate
    @constraint(Shipping_stock,[bt=1:BT],sum(Inv_dec[bt,y] for y=1:Y) <= 1) # Only one initial investment
    #Lock-in effect: once you invest you have to keep the production capacities
    #@constraint(Shipping_stock,[bt=1:BT,y=1:Y-1], Cap_bt[bt,y+1] >= Cap_bt[bt,y])


#------------------ From bottleneck technologies to fuels (e.g. H2 to eMeOH)
# Electrolyser is the bottleneck technology for H2 (max x PJ (or kg) prod per year)
# And H2 is a common resource for e-meoh, e-nh3 and e-bioemeoh

#ConveH2toFuel = ones(16) # Should be replaced on excel by quantity of H2 prod capacity needed (units consistent with Cap_bt) to produce 1 PJ of fuel
#Now all the same and equal to 1

#All fuel must be used (equality constraint)
@constraint(Shipping_stock,[y = 1:Y],
sum(TotFuel[f,y]/ConveH2toFuel[f,y] for f in eH2based) == Cap_bt[1,y]) #Cap_bt[1,y] is for electrolyser technology


@constraint(Shipping_stock,[y = 1:Y],
sum(TotFuel[f,y]/ConveCCStoFuel[f,y] for f in CCSbased) == Cap_bt[2,y]) #Cap_bt[2,y] is for CCS technology
end

#-----------------Constraints related to emissions------------

#emission constraint for year 2050 and onwards
#TTW Set settings
if Settings == "TTW"
    for f = 1:F, y = 1:Y
        fuel_emissions_WTT[f, y] = fuel_emissions_WTT[f, y] * 0
    end
elseif Settings != "TTW"
end
if MBM == "Ship_cap"
    if CO2_limit == "WTW_99"
        @constraint(
            Shipping_stock,
            ConCO2_4[t = 1:Ty, y = 1:Y],
            (sum(
                z[f, s, y] *
                (fuel_emissions_WTT[f, y] + fuel_emissions_TTW[f, y]) *
                ship_type_relation[t, s] *
                Ship_Demands_share[y, t] for f = 1:F, s = 1:S
            )) / (Ship_Demands[y, t] / average_transport_work_adj[y, t]) <=
            EL_WTW_99_adj[y, t]
        )
    end
elseif MBM != "Ship_cap"
    if CO2_limit == "WTW_99"
        @constraint(
            Shipping_stock,
            ConCO2_4[y = 1:Y],
            sum(
                z[f, s, y] *
                (fuel_emissions_WTT[f, y] + fuel_emissions_TTW[f, y]) for
                s = 1:S, f = 1:F
            ) <= EL_WTW_99[y]
        )
    end
end
# solve
optimize!(Shipping_stock)
CO2_margcost = zeros(Y)
for y = 1:Y
    if CO2_limit == "WTW_99"
        #CO2_margcost[y] = -dual(ConCO2_4[y]) * 1000 # In €/t CO2
    end
end



#--------------------
#OUTPUTS
Results_folder =
    joinpath(Main_folder, All_results_folder, Results_folder_name)
mkpath(Results_folder)
Year = round.(Int, zeros(Y))
for y = 1:Y
    Year[y] = 2020 + y - 1
end
fuel_fueltype_year = zeros(F, Y)
for f = 1:F
    for y = 1:Y
        fuel_fueltype_year[f, y] = sum(JuMP.value.(TotFuel[f, y]))
    end
end

fuel_ship_year = zeros(S, Y)
for s = 1:S
    for y = 1:Y
        fuel_ship_year[s, y] = sum(JuMP.value.(z[f, s, y]) for f = 1:F)
    end
end

bought_ships = zeros(S, Y)
for s = 1:S
    for y = 1:Y
        bought_ships[s, y] = JuMP.value.(x[s, y])
    end
end
ship_stock = zeros(S, Y)
for s = 1:S
    for y = 1:Y
        ship_stock[s, y] = JuMP.value.(q[s, y])
    end
end
decomission = zeros(S, Y)
for s = 1:S
    for y = 1:Y
        decomission[s, y] = JuMP.value.(d_EF[s, y])
    end
end
EmissionsTTW_f_y = zeros(F, Y)
EmissionsWTT_f_y = zeros(F, Y)
EmissionstotTTW = zeros(Y)
EmissionstotWTT = zeros(Y)
Emissionstot_f_y = zeros(F, Y)
Emissionstot = zeros(Y)
Emissionstot_total = zeros(Y)
Carbon_levy_data = zeros(F, Y)

for f = 1:F, y = 1:Y
    EmissionsTTW_f_y[f, y] = sum(
        JuMP.value.(z[f, s, y]) *
        fuel_emissions_TTW[f, y]  for s = 1:S
    )
    EmissionsWTT_f_y[f, y] = sum(
        JuMP.value.(z[f, s, y]) *
        fuel_emissions_WTT[f, y]  for s = 1:S
    )
    Emissionstot_f_y[f, y] = EmissionsTTW_f_y[f, y] + EmissionsWTT_f_y[f, y]
end
for y = 1:Y
    EmissionstotTTW[y] = sum(EmissionsTTW_f_y[f, y] for f = 1:F)
    EmissionstotWTT[y] = sum(EmissionsWTT_f_y[f, y] for f = 1:F)
    Emissionstot[y] = sum(Emissionstot_f_y[f, y] for f = 1:F)
end
    Emissionstot_total = sum(Emissionstot[y] for y = 1:Y)

for f = 1:F, y = 1:Y
    Carbon_levy_data[f, y] =
        sum(JuMP.value.(z[f, s, y]) * fuel_tax[f, y] for s = 1:S)
end

Total_cost = zeros(Y)
for y = 1:Y
    Total_cost[y] = JuMP.objective_value(Shipping_stock)
end

binary_data = zeros(F,Y)
for f=1:F
    for y = 1:Y
        binary_data[f,y] = sum(JuMP.value.(b[f,y,j]) for j = 1:J)
    end
end

fuel_cost_dynamic_data = zeros(F, Y)
for f = 1:F
    for y = 1:Y
    fuel_cost_dynamic_data[f, y] = JuMP.value.(fuel_cost_dynamic[f, y])
    end
end

Total_cost_year = zeros(Y)
for y = 1:Y
    Total_cost_year[y] = sum(
        ship_inv[s] * bought_ships[s, y] + ship_OM[s] * ship_stock[s, y] for s = 1:S
    )
    +sum(
        fuel_ship_year[s, y] * (fuel_cost[f, y] + fuel_tax[f, y]) for
        s = 1:S, f = 1:F
    )
    +sum(d_EF[s, y] * ship_inv[s] * decommission_value for s = 1:S)
end
Total_fuel_tax = zeros(F, Y)
for f = 1:F, y = 1:Y
    Total_fuel_tax[f, y] = fuel_tax[f, y]
end
#creates dataframes
Carbon_levy = DataFrame([Year transpose(Carbon_levy_data)],:auto)
fuel_cost_dynamic = DataFrame([Year transpose(fuel_cost_dynamic_data)],:auto)
binary = DataFrame([Year transpose(binary_data)],:auto)
ships_bought = DataFrame([Year transpose(JuMP.value.(x))],:auto)
Existing_fleet = DataFrame([Year transpose(JuMP.value.(EF))],:auto)
stock = DataFrame([Year transpose(JuMP.value.(q))],:auto)
decommission_new = DataFrame([Year transpose(JuMP.value.(d))],:auto)
decommission_EF = DataFrame([Year transpose(JuMP.value.(d_EF))],:auto)
fuel_f_y = DataFrame([Year transpose(fuel_fueltype_year)],:auto)
fuel_s_y = DataFrame([Year transpose(fuel_ship_year)],:auto)
Total_fueltax_data = DataFrame([Year transpose(Total_fuel_tax)],:auto)
emissions = DataFrame(
    [Year transpose(EmissionsTTW_f_y) transpose(EmissionsWTT_f_y) transpose(
        Emissionstot_f_y,
    )],:auto
)
emissionstot =
    DataFrame([Year EmissionstotWTT EmissionstotTTW Emissionstot],:auto)
Results_CO2_marcost = DataFrame([Year CO2_margcost],:auto)
obj_value = DataFrame([Year Total_cost],:auto)
obj_value_year = DataFrame([Year Total_cost_year],:auto)
if MBM != "ETS"
elseif MBM == "ETS"
    emission_allowances = DataFrame([Year (JuMP.value.(st))],:auto)
end



#adds headers
rename!(ships_bought, ["Year"; Ships])
if MBM != "ETS"
elseif MBM == "ETS"
    rename!(emission_allowances, ["Year"; "Revenue"])
end
rename!(Existing_fleet, ["Year"; Ships])
rename!(stock, ["Year"; Ships])
rename!(decommission_new, ["Year"; Ships])
rename!(decommission_EF, ["Year"; Ships])
rename!(fuel_cost_dynamic, ["Year"; Fuels])
rename!(fuel_f_y, ["Year"; Fuels])
rename!(fuel_s_y, ["Year"; Ships])
rename!(emissionstot, ["Year", "WTT", "TTW", "WTW"])
rename!(Results_CO2_marcost, ["Year", "CO2 marginal cost(€/t CO2)"])
rename!(
    Carbon_levy,
    [
        "Year",
        "VLSFO/HFOsc",
        "MDO",
        "MGO",
        "LNG",
        "LPG",
        "MeOH-grey",
        "MeOH-blue",
        "MeOH-ebio",
        "MeOH-CCU",
        "MeOH-DAC",
        "NH3-grey",
        "NH3-blue",
        "NH3-green",
        "Refined Pyrolysis Oil",
        "LBG",
        "Gap",
    ],
)

# write DataFrame out to CSV file
CSV.write(joinpath(Results_folder, "carbon_levy.csv"), Carbon_levy)
CSV.write(joinpath(Results_folder, "Emissions.csv"), emissions)
CSV.write(joinpath(Results_folder, "Fueltax.csv"), Total_fueltax_data)
CSV.write(joinpath(Results_folder, "Emissions_tot.csv"), emissionstot)
CSV.write(
    joinpath(Results_folder, "Results_ships_bought.csv"),
    ships_bought,
)
CSV.write(joinpath(Results_folder, "Results_stock.csv"), stock)
CSV.write(joinpath(Results_folder, "Results_EF.csv"), Existing_fleet)
CSV.write(
    joinpath(Results_folder, "Results_decommission_EF.csv"),
    decommission_EF,
)
CSV.write(
    joinpath(Results_folder, "Results_decommission_new.csv"),
    decommission_new,
)
Emissionstot = zeros(Y)
for  y = 1:Y
    Emissionstot[y] = Emissionstot_total
end
    Emissionstot_total_adj = DataFrame([Year Emissionstot],:auto)
CSV.write(joinpath(Results_folder, "Cumulated_Emissions.csv"), Emissionstot_total_adj)
CSV.write(joinpath(Results_folder, "Results_fuels_f_y_PJ.csv"), fuel_f_y)
CSV.write(joinpath(Results_folder, "Fuel_cost_dynamic.csv"), fuel_cost_dynamic)
CSV.write(joinpath(Results_folder, "binary.csv"), binary)
CSV.write(joinpath(Results_folder, "Results_fuels_s_y_PJ.csv"), fuel_s_y)
CSV.write(
    joinpath(Results_folder, "Results_CO2magcost.csv"),
    Results_CO2_marcost,
)
CSV.write(joinpath(Results_folder, "obj_value.csv"), obj_value)
CSV.write(joinpath(Results_folder, "obj_value_year.csv"), obj_value_year)
if MBM != "ETS"
elseif MBM == "ETS"
    CSV.write(
        joinpath(Results_folder, "emission_allowances.csv"),
        emission_allowances,
    )
end
println("Objective value: ---> ", JuMP.objective_value(Shipping_stock))



global NScen += 1 #Increment and run the script for another scenario

end
