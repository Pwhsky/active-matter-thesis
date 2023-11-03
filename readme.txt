![](<https://raw.githubusercontent.com/Pwhsky/active-matter-thesis/main/janus-gradient/figures/quiver3.png>)

#Requirements:
G++, Pandas, Matplotlib, FFmpeg for movie creation, pyarrow, and whatever other packages i forgot to list here that python will complain about


The program takes in 5 arguments as of now, you must specify them or the program will not run.

  $ python3 temperature.py nDeposits newData bounds periodicity

  nDeposits   (integer): Specifies the number of deposits to be createdin the particle, it can go anywhere from 1 to 10000.
  newData     (boolean): Specifies whether or not to create new datapoints and overwrite the old data, true or false.    
  bounds      (double) : Specifies the size of the simulation box, setting it to 4 gives 4 micrometer length per side
  periodicity (double) : Specifies the laser periodicity in nanometers, I commonly use 80000





Examples (be in janus-gradient directory):

  $ python3 gradient_quiver.py 600 true 3 80000
  $ python3 temperature.py 2000 true 4 80000
  $ python3 gradient.py 600 false 3 80000

Note: even if you put in false, you still have to specify the bounds as this will change the plotting properties as well.

functions used to initiate deposits are in:
functions.cpp

Constants are in:
functions.h


The step-size can be computed by dividing the bounds parameter with 300.

