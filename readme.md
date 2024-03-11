![alt text](https://raw.githubusercontent.com/Pwhsky/active-matter-thesis/main/thermophoresis/figures/display/display.png?raw=true)
![alt text](https://raw.githubusercontent.com/Pwhsky/active-matter-thesis/main/thermophoresis/figures/display/display2.png?raw=true)
![alt text](https://raw.githubusercontent.com/Pwhsky/active-matter-thesis/main/thermophoresis/figures/display/tangential.png?raw=true)


# Requirements:
G++, Pandas, Matplotlib, pyarrow, and whatever other packages I forgot to list here that python will complain about...


 # Arguments:
 The programs are located in the /thermophoresis directory, and take in 4 arguments, you must specify them or the program will not run.
 The arguments are as follows:
 
```console
   $ python3 temperature.py    <nDeposits> <generateData> <bounds> <periodicity>
   $ python3 gradient.py       <nDeposits> <generateData> <bounds> <periodicity>
   $ python3 time_evolution.py <nDeposits> <nTimeSteps>   <bounds> <periodicity>
```
- nDeposits    (int):   Specifies the number of deposits to be created in the particle, it can go anywhere from 1 to 10000.
- generateData (str):   Specifies whether or not to create new datapoints and overwrite the old data, true or false.    
- bounds       (float): Specifies the size of the simulation box, setting it to 4 gives 4x4 micrometer sized box.
- periodicity  (float): Specifies the laser periodicity in nanometers, I commonly use 80000 (80 microns).
- nTimeSteps   (int):   Number of total time steps to simulate.





# Examples: 
```console
  $ python3 gradient_quiver.py 600  true  3  80000
  $ python3 temperature.py     2000 true  4  80000
  $ python3 gradient.py        600  false 3  80000
  $ python3 time_evolution.py  600  50    5  80000
```


functions.cpp contain:
- initializing particles
- linear algebra 
- geometric functions


particle.cpp contain:
- constants
- main integral
- kinematics
