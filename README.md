# SEAMAPS
The SEAMAPS model is a least-cost optimization model dedicated to anaylse future fuel mix, emissions and underlying cost of the maritime industry.
SEAMAPS runs in Julia. The Data-format is an excel document. The analysis of results works via an R-Markdown file.

How to use it:
1. Install Julia and the following packages: JuMP, Gurobi, CSV, DataFrames, DelimitedFiles'
2. Install R & R-Studio
3. Download this repository to your local machine
4. Open the folder and go to /code - here open the file "SEAMAPS.jl"
5. Inside "SEAMAPS.jl" -Insert your local path to the SEAMAPS/data folder
6. Run SEAMAPS for all the relevant scenarios you want to analyse (e.g. CO2-tax, Best practice, CO2-tax-flat)
7. Open SEAMAPS/code and open "Carbon_Price_Paper_Plots.RMD" in R-Studio
8. In the R-Markdown file - Insert your result folders paths into the file
9. Knit the file!
10. Analyse your results!

For further questions please contact semfr@dtu.dk
