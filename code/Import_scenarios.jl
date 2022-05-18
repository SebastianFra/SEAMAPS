using ExcelReaders

Data_scenarios = readxlsheet(joinpath(Main_folder,Scenarios_file_name),Scenarios_set)

L1_scenario = findfirst(x -> x == "Scenario number", Data_scenarios)[1] + 1
C_Power_supply = findfirst(x -> x == "Power supply type", Data_scenarios)[2]
C_bio_av = findfirst(x -> x == "Biomass available", Data_scenarios)[2]
C_CO2_lim = findfirst(x -> x == "CO2 cap", Data_scenarios)[2]
C_carbon_price = findfirst(x -> x == "Carbon Price", Data_scenarios)[2]
C_carbon_multiplier = findfirst(x -> x == "Carbon Multiplier", Data_scenarios)[2]
C_fuel_savings = findfirst(x -> x == "Fuel Savings", Data_scenarios)[2]
#C_input_data = findfirst(x -> x == "Input data sheet", Data_scenarios)[2]
C_results_folder = findfirst(x -> x == "Result folder name", Data_scenarios)[2]

All_Power_supply = Array{String}(Data_scenarios[L1_scenario:end , C_Power_supply]) ; N_scenarios = length(All_Power_supply)
All_bio_av = Array{String}(Data_scenarios[L1_scenario:end , C_bio_av])
All_CO2_lim = Array{String}(Data_scenarios[L1_scenario:end , C_CO2_lim])
All_fuel_savings = Array{Int64}(Data_scenarios[L1_scenario:end , C_fuel_savings])
All_Carbon_price = Array{Int64}(Data_scenarios[L1_scenario:end , C_carbon_price])
All_carbon_multiplier = Array{Float64}(Data_scenarios[L1_scenario:end , C_carbon_multiplier])
All_results_folder_names = Array{String}(Data_scenarios[L1_scenario:end , C_results_folder])
#All_Input_data = Array{String}(Data_scenarios[L1_scenario:end , C_input_data])
