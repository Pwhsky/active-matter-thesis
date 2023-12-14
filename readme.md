![alt text](https://raw.githubusercontent.com/Pwhsky/active-matter-thesis/main/thermophoresis/figures/temperatureIncrease.png?raw=true)
![alt text](https://raw.githubusercontent.com/Pwhsky/active-matter-thesis/main/thermophoresis/figures/tangential.png?raw=true)
![alt text](https://raw.githubusercontent.com/Pwhsky/active-matter-thesis/main/thermophoresis/figures/quiver.png?raw=true)


# Requirements:
G++, Pandas, Matplotlib, pyarrow, and whatever other packages I forgot to list here that python will complain about...




 # Arguments:
 The program takes in 4 arguments as of now, you must specify them or the program will not run.
 The arguments are as follows:
```console
   $ python3 temperature.py <nDeposits> <generateData> <bounds> <periodicity>
```
- nDeposits    (int): Specifies the number of deposits to be createdin the particle, it can go anywhere from 1 to 10000.
- generateData (str): Specifies whether or not to create new datapoints and overwrite the old data, true or false.    
- bounds       (float): Specifies the size of the simulation box, setting it to 4 gives 4x4 micrometer sized box.
- periodicity  (float): Specifies the laser periodicity in nanometers, I commonly use 80000



# Examples: 
```console
  $ python3 gradient_quiver.py 600 true 3 80000
  $ python3 temperature.py 2000 true 4 80000
  $ python3 gradient.py 600 false 3 80000
```

Note: even if you put in false, you still have to specify the bounds as this will change the plotting properties as well.

functions used to initiate deposits are in:
functions.cpp

Constants and class declarations are in:
functions.h
