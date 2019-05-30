# Body-size-distribution

The code in the file body-size-distribution.jl has the following functions: 

* t_repr, which calculates the reproduction time of an animal or a set of animals considering their lengths, 

* vecinos, that calculate which are the neighbords a cell in a lattice.

* reproduccion_comida, which literally means food reproduction and calculate the quantity of food in each cell of the lattice after a unity of time, using a Logistic-like model, and considering the volume of animals

* condiciones_iniciales, which just produce random initial conditions for the animals, 

* movimiento, which reproduce the animals and move the babies to a new cell

* dispersion, which use the other functions to simulate during a time the proces of reproducing, migrate and dying. 
