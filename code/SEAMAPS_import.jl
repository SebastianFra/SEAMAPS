using ExcelReaders

Datafile = joinpath(Main_folder, Data_file_name)

Biomass_av = All_bio_av[NScen] #HBLD, MBMD, LBLD
CO2_limit = All_CO2_lim[NScen] #TTW_50, WTW_50, WTW_70, WTW_99
carbon_price = All_Carbon_price[NScen]
Power_supply = All_Power_supply[NScen] #Mixed, Off-grid
Results_folder_name = All_results_folder_names[NScen]
carbon_multiplier = All_carbon_multiplier[NScen]
fuel_savings = All_fuel_savings[NScen]

#Biomass_av*"_"*CO2_limit*"_"*Power_supply


#-------------------------Input data-----------------------------------

Data_emissions_cap = readxlsheet(Datafile, "Emission_cap")
Data_emissions_cap_WTW99 = readxlsheet(Datafile, "99%_WTW")
Data_emissions_cap_TTW = readxlsheet(Datafile, "TTW")
Data_fuel_emissions_WTT_infra =
    readxlsheet(Datafile, "i_Fuel_emissions_WTT_" * Power_supply)
Data_fuel_emissions_WTT =
    readxlsheet(Datafile, "Fuel_emissions_WTT_" * Power_supply)
Data_Gap_fuel_emissions_WTT =
    readxlsheet(Datafile, "Gap_Fuel_emissions_WTT_" * Power_supply)
Data_fuel_emissions_TTW = readxlsheet(Datafile, "Fuel_emissions_TTW")
Data_Gap_fuel_emissions_TTW = readxlsheet(Datafile, "Gap_Fuel_emissions_TTW")
Data_fuel_cost = readxlsheet(Datafile, "Fuel_cost_" * Power_supply)
Data_availabity = readxlsheet(Datafile, "Availability_" * Biomass_av)
Data_biomass_availabity = readxlsheet(Datafile, "Shared_biomass_available")
Data_ships = readxlsheet(Datafile, "Ship_data")
Data_prod_cap = readxlsheet(Datafile, "Ship_prod_cap")
if fuel_savings == 90
    Data_Fuel_savings = readxlsheet(Datafile, "Fuel Savings_90")
elseif fuel_savings == 80
    Data_Fuel_savings = readxlsheet(Datafile, "Fuel Savings_80")
elseif fuel_savings == 70
    Data_Fuel_savings = readxlsheet(Datafile, "Fuel Savings_70")
elseif fuel_savings == 60
    Data_Fuel_savings = readxlsheet(Datafile, "Fuel Savings_60")
elseif fuel_savings == 50
    Data_Fuel_savings = readxlsheet(Datafile, "Fuel Savings_50")
elseif fuel_savings == 40
    Data_Fuel_savings = readxlsheet(Datafile, "Fuel Savings_40")
elseif fuel_savings == 30
    Data_Fuel_savings = readxlsheet(Datafile, "Fuel Savings_30")
elseif fuel_savings == 20
    Data_Fuel_savings = readxlsheet(Datafile, "Fuel Savings_20")
elseif fuel_savings == 10
    Data_Fuel_savings = readxlsheet(Datafile, "Fuel Savings_10")
elseif fuel_savings == 0
    Data_Fuel_savings = readxlsheet(Datafile, "Fuel Savings_0")
end
Data_ship_eff = readxlsheet(Datafile, "Ship_fuel_eff")
Data_existing_fleet = readxlsheet(Datafile, "Existing_fleet")
Data_av_trans_work = readxlsheet(Datafile, "Av_trans_work")
Data_av_trans_work_adj = readxlsheet(Datafile, "Av_trans_work_adj")
Data_ship_relation = readxlsheet(Datafile, "Ship_type_relation")
Data_demand = readxlsheet(Datafile, "Ship_demand")
Data_demand_share = readxlsheet(Datafile, "Ship_demand_share")
Data_CO2_price = readxlsheet(Datafile, "CO2-Price")
Data_Conversion_1 = readxlsheet(Datafile, "Conversion_1")
Data_Conversion_2 = readxlsheet(Datafile, "Conversion_2")
Data_up_ramping_factor = readxlsheet(Datafile, "ramping_up_factor")
Data_up_ramping_factor_total = readxlsheet(Datafile, "ramping_up_factor_total")



Conversion_1 = Data_Conversion_1[2:end, 2:end]
Conversion_2 = Data_Conversion_2[2:end, 2:end]
CO2_price = Data_CO2_price[2:end, 2:end]
Fuel_savings = Data_Fuel_savings[2:end, 2:end]
ramping_up_factor = Data_up_ramping_factor[2:end, 2:end]
ramping_factor_total = Data_up_ramping_factor_total[2:end, 2:end]
EL_WTW_99 = Data_emissions_cap[2:end, 5]
CO2_price = Data_CO2_price[2:end, 2:end] * 1
Fuels = Array{String}(Data_fuel_cost[2:end, 1]);
F = length(Fuels);
BioMeOH = findfirst(x -> x == "MET-ebio", Data_fuel_cost[2:end, 1])
RefPO = findfirst(x -> x == "Ref-PO", Data_fuel_cost[2:end, 1])
PSMeOH = findfirst(x -> x == "MET-CCU", Data_fuel_cost[2:end, 1])
NH3_green = findfirst(x -> x == "AMM-green", Data_fuel_cost[2:end, 1])
DACMeOH = findfirst(x -> x == "MET-DAC", Data_fuel_cost[2:end, 1])
NH3_blue = findfirst(x -> x == "AMM-blue", Data_fuel_cost[2:end, 1])
LBG = findfirst(x -> x == "LBG", Data_fuel_cost[2:end, 1])

if Set == "Infra"
    fuel_emissions_WTT = Data_fuel_emissions_WTT_infra[2:end, 2:end]
elseif Set != "Infra"
    fuel_emissions_WTT = Data_fuel_emissions_WTT[2:end, 2:end]
end
if Set == "GWP20"
    fuel_emissions_TTW = Data_fuel_emissions_TTW_GWP20[2:end, 2:end]
elseif Set != "GWP20"
    fuel_emissions_TTW = Data_fuel_emissions_TTW[2:end, 2:end]
end
fuel_emissions_WTT = fuel_emissions_WTT
Gap_fuel_emissions_WTT = Data_Gap_fuel_emissions_WTT[2:end, 2:end]
Gap_fuel_emissions_TTW = Data_Gap_fuel_emissions_TTW[2:end, 2:end]
fuel_available = Data_availabity[2:end, 2:end]


if Biomass_av == "HBLD"
    Shared_biomass_available = Data_biomass_availabity[2:end, 2]
elseif Biomass_av == "MBMD"
    Shared_biomass_available = Data_biomass_availabity[2:end, 3]
elseif Biomass_av == "LBLD"
    Shared_biomass_available = Data_biomass_availabity[2:end, 4]
end
ConvBio_PO = Data_biomass_availabity[2:end, 6]
ConvBio_MeOH = Data_biomass_availabity[2:end, 7]
fuel_cost = Data_fuel_cost[2:end, 2:end]
Ships = Array{String}(Data_ships[2:end, 1]);
S = length(Ships);
ship_lifetime = round.(Int, Data_ships[2:end, 2])
ship_inv = Data_ships[2:end, 3]
ship_OM = Data_ships[2:end, 4]
ship_prod_cap = Data_prod_cap[2:end, 2:end]
if Option_ramping_ships == "FALSE"
    for s = 1:S, y = 1:Y
        ship_prod_cap[s, y] = ship_prod_cap[s, y] * 999999999999
    end
elseif Option_ramping_ships == "TRUE"
end
ship_eff = Data_ship_eff[3:end, 2:end]
ship_type = Data_demand[1, 2:end];
Ty = length(ship_type);
preexisting_fleet = round.(Int, Data_existing_fleet[2:end, 2:end])

average_transport_work = Data_av_trans_work[2:end, 2:end]
average_transport_work_adj = Data_av_trans_work_adj[2:end, 2:end]
ship_type_relation = Data_ship_relation[2:end, 2:end]
Ship_Demands = Data_demand[2:end, 2:end]
for y = 1:Y, t = 1:Ty
    Ship_Demands[y, t] = Ship_Demands[y, t] *Fuel_savings[y]
end

Ship_Demands_share = Data_demand_share[2:end, 2:end]


if Ambition == "Realistic"
    CO2_price = Data_CO2_price[2:end, 2:end]*carbon_multiplier
    EL_WTW_99_adj = Data_emissions_cap_WTW99[2:end, 2:end]
end
#create fuel tax
if Taxed == "WTW"
    if MBM != "CO2_tax_flat"
        for y = 1:Y
            CO2_price[y] = CO2_price[y] *carbon_multiplier
        end
        fuel_tax = zeros(F, Y)
        for f = 1:F, y = 1:Y
            fuel_tax[f, y] =
                (
                    (
                        (
                            (
                                (
                                    fuel_emissions_WTT[f, y] +
                                    fuel_emissions_TTW[f, y]
                                ) * 3.6
                            ) * Conversion_2[f, y]
                        ) / 1000
                    ) * CO2_price[y]
                ) / Conversion_1[f, y]
        end
    end
elseif Taxed == "TTW"
    if MBM != "CO2_tax_flat"
        for y = 1:Y
            CO2_price[y] = CO2_price[y] *carbon_multiplier
        end
        fuel_tax = zeros(F, Y)
        for f = 1:F, y = 1:Y
            fuel_tax[f, y] =
                (
                    (
                        (
                            (fuel_emissions_TTW[f, y] * 3.6) *
                            Conversion_2[f, y]
                        ) / 1000
                    ) * CO2_price[y]
                ) / Conversion_1[f, y]
        end
    end
end

if MBM == "CO2_tax"
    for y = 1:Y
        EL_WTW_99[y] = EL_WTW_99[y] * 10000000000
    end
    for y = 1:Y
        EL_WTW_99_adj[y] = EL_WTW_99_adj[y] * 10000000000
    end
    for f = 1:F, y = 1:Y
        fuel_tax[f, y] = fuel_tax[f, y]
    end
elseif MBM == "CO2_tax_flat"
    for y = 1:Y
        EL_WTW_99[y] = EL_WTW_99[y] * 10000000000
    end
    for y = 1:Y
        EL_WTW_99_adj[y] = EL_WTW_99_adj[y] * 10000000000
    end
    CO2_price = Data_CO2_price[2:end, 2:end] * 1
    for y = 1:Y
        CO2_price[y] = carbon_price
    end
    fuel_tax = zeros(F, Y)
    for f = 1:F, y = 1:Y
        fuel_tax[f, y] =
            (
                (
                    (
                        (
                            (
                                fuel_emissions_WTT[f, y] +
                                fuel_emissions_TTW[f, y]
                            ) * 3.6
                        ) * Conversion_2[f, y]
                    ) / 1000
                ) * CO2_price[y]
            ) / Conversion_1[f, y]
    end
elseif MBM == "Baseline"
    for y = 1:Y
        EL_WTW_99[y] = EL_WTW_99[y] * 10000000000
    end
    for y = 1:Y
        EL_WTW_99_adj[y] = EL_WTW_99_adj[y] * 10000000000
    end
    for f = 1:F, y = 1:Y
        fuel_tax[f, y] = fuel_tax[f, y] * 0
    end
elseif MBM == "Global_cap"
    for y = 1:Y
        EL_WTW_99_adj[y] = EL_WTW_99_adj[y] * 10000000000
    end
    for f = 1:F, y = 1:Y
        fuel_tax[f, y] = fuel_tax[f, y] * 0
    end
    for y = 1:Y
        EL_TTW_50[y] = EL_TTW_50[y]
    end
    for y = 1:Y
        EL_WTW_99[y] = EL_WTW_99[y]
    end
    for f = 1:F, y = 1:Y
        fuel_tax[f, y] = fuel_tax[f, y] * 0
    end
end

if green_push == "TRUE"
    for f in NH3_green, y = 1:Y, fuel_cost[f, y] in fuel_cost[f, y] * 0.75
    end
    for f in PSMeOH, y = 1:Y, fuel_cost[f, y] in fuel_cost[f, y] * 0.75
    end
    for f in DACMeOH, y = 1:Y, fuel_cost[f, y] in fuel_cost[f, y] * 0.75
    end
    for f in BioMeOH, y = 1:Y, fuel_cost[f, y] in fuel_cost[f, y] * 0.75
    end
elseif green_push == "FALSE"
end
if green_decay == "TRUE"
    for f in NH3_green, y = 1:Y, fuel_cost[f, y] in fuel_cost[f, y] * 1.25
    end
    for f in PSMeOH, y = 1:Y, fuel_cost[f, y] in fuel_cost[f, y] * 1.25
    end
    for f in DACMeOH, y = 1:Y, fuel_cost[f, y] in fuel_cost[f, y] * 1.25
    end
    for f in BioMeOH, y = 1:Y, fuel_cost[f, y] in fuel_cost[f, y] * 1.25
    end
elseif green_decay == "FALSE"
end

if green_push_em == "TRUE"
    for f in NH3_green, y = 1:Y, fuel_emissions_WTT[f, y] in fuel_emissions_WTT[f, y] * 0.75
    end
    for f in PSMeOH, y = 1:Y, fuel_emissions_WTT[f, y] in fuel_emissions_WTT[f, y] * 0.75
    end
    for f in DACMeOH, y = 1:Y, fuel_emissions_WTT[f, y] in fuel_emissions_WTT[f, y] * 0.75
    end
    for f in BioMeOH, y = 1:Y, fuel_emissions_WTT[f, y] in fuel_emissions_WTT[f, y] * 0.75
    end
elseif green_push_em == "FALSE"
end
if green_decay_em == "TRUE"
    for f in NH3_green, y = 1:Y, fuel_emissions_WTT[f, y] in fuel_emissions_WTT[f, y] * 1.25
    end
    for f in PSMeOH, y = 1:Y, fuel_emissions_WTT[f, y] in fuel_emissions_WTT[f, y] * 1.25
    end
    for f in DACMeOH, y = 1:Y, fuel_emissions_WTT[f, y] in fuel_emissions_WTT[f, y] * 1.25
    end
    for f in BioMeOH, y = 1:Y, fuel_emissions_WTT[f, y] in fuel_emissions_WTT[f, y] * 1.25
    end
elseif green_decay_em == "FALSE"
end
