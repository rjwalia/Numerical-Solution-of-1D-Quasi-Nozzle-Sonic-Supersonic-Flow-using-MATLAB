# Numerical-Solution-of-1D-Quasi-Nozzle-Sonic-Supersonic-Flow-using-MATLAB

# Objective 
To solve Quasi 1D Converging-Diverging Nozzle numerically for both conservative & non-conservative forms of governing equations by implementing MacCormack's Method using MATLAB.

# Problem Description
We consider a steady, isentropic flow through the converging-diverging nozzle. The inlet of the nozzle is connected with the reservoir. The cross-sectional area of the reservoir is large (tending to infinity), and hence the velocity inside the reservoir is very small(tending to 0). The flow is locally subsonic in the converging section of the nozzle, sonic at the nozzle throat, and expands isentropically to the supersonic speed at the nozzle exit. 

The area at the throat is A*. We assume that at a given section, the flow properties are uniform or constant across that section. Hence, Although the area of the nozzle changes as a function of nozzle’s distance x, and therefore in reality the problem is 2D at minimum in reality, We make an assumption that the flow properties vary only along x-direction of the nozzle. Such flow are defined as quasi-one-dimensional flow.	

Non-dimensional governing equations for the Quasi 1D Converging-Diverging Nozzle in the conservative and non-conservative form are discretized & solved with MacCormack's Method using MATLAB. Then the following numerical results obtained from both forms of equations are compared:-

1. Non-Dimensional mass flow rate at the nozzle’s throat location for different time steps.
2. Non-Dimensional flow variables with respect to time at nozzle's throat location.
3. Steady-State flow field along the nozzle X-direction.
4. Comparison of steady-state mass flow rate along the Nozzle X-Direction for the Conservative & Non-Conservative form.
5. Comparison of steady-state flow variable at the throat between Conservative, Non-Conservative and Analytical Results.
6. Grid Sensitivity study.

# Procedure of running the script
1. Download the main_code, non_conservative_form & conservative_form file.
2. Copy and Save all three files in a same folder.
3. Open all three scripts using MATLAB. Do not change any value in any of the code. 
4. Run the main_code. It will start the simulation.
5. Wait for sometime (5-10 minutes), You would be able to see various numerical results plotted on the screen.
