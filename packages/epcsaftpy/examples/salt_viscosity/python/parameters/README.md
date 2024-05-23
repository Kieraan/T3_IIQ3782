# kijs.json 
-  Sources: 
    - Ions - water: Table 2 from [Fluid Phase Equilibria 535 (2021) 112967](https://doi.org/10.1016/j.ﬂuid.2021.112967).
    - Ions - water:  Table 3 from [Fluid Phase Equilibria 537 (2021) 112989](https://doi.org/10.1016/j.ﬂuid.2021.112989).
        * The exponential has to be 1e-3 for the thermal contribution instead of 1e-4 in order to be consistent with the kij at 298.15 K in the previous table.
    - Cation - anions: Table 4 from [Fluid Phase Equilibria 537 (2021) 112989](https://doi.org/10.1016/j.ﬂuid.2021.112989).
- The parameters are saved as a dictionary, where the key is the combination between the name of component i and component j, and the value can be a float or a list (for thermal dependency).

# solvents.json
- Sources:
    - Parameters: Table 1 from  [Fluid Phase Equilibria 535 (2021) 112967](https://doi.org/10.1016/j.ﬂuid.2021.112967).
    - Thermal dielectric constant: Table 5 from [Fluid Phase Equilibria 537 (2021) 112989](https://doi.org/10.1016/j.ﬂuid.2021.112989).

# iones.json
- Sources:
    - Parameters: Table 3 from [Fluid Phase Equilibria 537 (2021) 112989](https://doi.org/10.1016/j.ﬂuid.2021.112989).