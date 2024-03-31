Application of FEM on truss element using the python codes 

I.	The nodal co-ordinates of the Truss are inserted in the code in the form of an array.
II.	Similarly, an Element array is also inserted which contains the node numbers which are forming the elements, in element wise.
III. Now the no. of nodes, no. of elements, area of cross section, Elastic Modulus values are taken as input while running the code.
IV.	The load acting on the Truss is inserted in the form of an array. The load acting is oriented in x-axis and y-axis acting at each node sequentially.
    Positive value of load indicates the load is acting in the direction of positive x-axis or positive y-axis. In the load vector x-axis is considered first. 
V.	As, here planar truss is considered. There are Two degrees of freedom at each node (One along x-axis and other along y-axis).
    Now, the boundary conditions are applied in the code. The Degrees of freedom where displacement is Zero is taken as input while running the code. 
VI.	For plotting the shape of the structure the nodes are connected element wise.
VII.	For plotting the deformed shape of the structure.
      The displacement at each degree of freedom is added to its respective nodal co-ordinates and then the nodes are connected element wise.
