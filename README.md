# Physical Model of Brinicle Growth Rate
<div align ="center" style="font-size:20px;color:white;background-color:black;"> 
For Code Documentation:

Please Visit the `Simulation` Folder Above.
</div>

# Introduction
A brinicle (brine icicle, also known as ice stalactite) is a downward-growing hollow tube of ice enclosing a plume of descending brine that is formed beneath developing sea ice. The main goal of this research is to build a simulation of Brinicles formation and analyze the length of the brinicle as a function of time.
<div align="center">
    <img src="3D Models/brinicle.jpg" alt="drawing" width="300"/>
    <div> Figure 1: Image of brinicle from BBC documentary </div>
</div>

# Growth Model
As found by Paul K Dayton and Seelye Martin in their research on [Observations of ice stalactites in McMurdo Sound, Antarctica](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/JC076i006p01595), the length of the brinicle as a function of time exhibits a 
<img src="https://latex.codecogs.com/gif.latex?\dpi{70}&space;\bg_black&space;\fn_jvn&space;\sqrt{t}" 
/> relationship such that: 
<!--You must use &space instead of " " or it will break-->
<!-----------LATEX IN HTML----------->
<div align ="center"> 
    <img src="https://latex.codecogs.com/gif.latex?\dpi{150}&space;\bg_black&space;\fn_jvn&space;\boxed{L(t)=\alpha\sqrt{t}}"/>
</div>
<!--------------------------------------->
where <img src="https://latex.codecogs.com/gif.latex?\dpi{100}&space;\bg_black&space;\fn_jvn&space;\alpha" 
/> is some constant. Our goal is to create a "brinicle growth tank" to analyze this <img src="https://latex.codecogs.com/gif.latex?\dpi{70}&space;\bg_black&space;\fn_jvn&space;\sqrt{t}" 
/> 

relationship and determine the alpha growth coefficient. We will then compare these results to those found from python simulations created by `Caleb Powell, Fardin Hoque, and Professor Dr. Nathan Tompkins` of Wabash College in the summer of 2020.
<div align="center">
    <img src="3D Models/results.gif" alt="drawing" width="300"/>
    <div> Figure 2: Python Simulation Video Output </div>
</div>
As seen from the simulation above, a value of the coefficient was discovered to be 
<img src="https://latex.codecogs.com/gif.latex?\dpi{100}&space;\bg_black&space;\fn_jvn&space;\alpha\approx0.139" 
/>
. We hope to improve the results of this simulation using observed values for the Advection and Diffusion rate of Temperature, and Salt concentrations.


# Lab Work
## Setup
<div align="center">
    <img src="3D Models/Lab_Setup2.jpg" alt="drawing" width="250"/>
    <img src="3D Models/3D View 1.jpg" alt="drawing" width="250"/>
    <div> Figure 3: Image of lab setup (left) and 3D model of setup (right) </div>
</div>

### Required Materials
- Small Fish Tank
- Two 5 Gallon Buckets
- 50' 1/4" stainless steele tubing.
- Two [Water Cooling Pumps](https://www.amazon.com/Water-Cooling-Integrated-Support-System/dp/B07PB1WJ6Z/ref=sr_1_9?keywords=cpu+pump&qid=1639951705&sr=8-9)
- Two low voltage pumps
- Four circuit switches (1 for each pump)

### Procedure 
1. Fill the left bucket and the fish tank with a combination of ice and water to achieve the lowest temperature possible.
2. Mix the right bucket with rock salt to and ice to allow the ice to melt and dissolve the salt (effectively creating brine)
3. Note the temperature of all the containers of water
4. Turn on the front two pumps (clear ones) that are connected to the steele pipes to maintain a steady temperature.
5. Once ice is melted and a stable temperature is achieved, turn on the pumps above the brine bucket.
6. Record the rate of diffusion and advection for salt and temperature
7. Measure the length of the brinicle over time.
<div align="center">
    <img src="3D Models/sketch.jpg" alt="drawing" width="400"/>
    <div> Figure 4: Sketch of Brine Dispense and Disposal </div>
</div>
The functions of the two pumps above the brine bucket are to dispense brine into the system (create the brinicle) and remove excess brine from the bottom of the tank and dispose it down the sink drain.  The functions of the pumps on the wood blocks are to maintain constant temperatures within each of the buckets and the tank itself.

# Results
<div align ="center" style="font-size:25px;color:white;background-color:black;"> 
Due to delays in shipping and difficulties in the manufacturing of the brinicle tanks, no trials we able to be run in the Fall 2021 Semester. However, tank construction is complete and trials will be run at the start of the Spring 2022 Semester.
</div>

<!--
## Notes
switching between branches
```
git reset --hard origin/<branch_name>
```
-->