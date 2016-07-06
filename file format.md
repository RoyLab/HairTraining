# Guide Hairs, \*.guide, \*.info

int a: number of guide strand  
int a1: particle per strand   
int b: number of frame  
int \* a: guide id  
int: frame id  
float \* a \* a1 \* (6+3): rot, translate  

---

# Weights, \*.weights

int b: number of strand  

For Loop \* b  
* int: number of guides  
* int:   guide Id
* float:   guide weight

---

# Ground Truth, \*.anim

int a: frame number  
int b: particle number  

For Loop \* a  
int: frame id
* **For Loop** \* b  
* float \* 3: position

---

# Ground Truth, \*.anim2

* INT b: particle number  
* **For Loop** \* a  
  * INT: frame id
  * FLOAT * 16: rigid motion
  * **For Loop** \* b  
    * FLOAT \* 3: position
  * **For Loop** \* b  
    * FLOAT \* 3: direction

---

# Neighbour Map, \*.neigh

* INT a: group number  
* **For Loop** \* a  
  * INT b: number of neigh, including itself
  * int \* b: neigh group id
