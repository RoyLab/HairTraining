# \*.weights

* **Int** a: particle number
* **Int** b: strand number  
* **Int** k: particle per strand
* **Int** c: guide strand number
* 64Byte: reserve for future, fill with 0
* **For Loop** x c
  * **Int** guide Ids
* **For Loop** x b
  * **Int p:** number of guides  
  * **For Loop** x p
    * **Int:** guide ID (global overall ID)
  * **For Loop** x k
    * **For Loop** x p
      * **Float:** guide weight


# \*.group

* **Int** a: particle number
* **Int** b: strand number  
* **Int** k: particle per strand
* **Int** n: group number
* 64Byte: reserve for future, fill with 0
* **For Loop** x a
  * **Int** group ids


# \*.anim2

* **Int** b: particle number  
* **Int** c: strand number  
* **Int** k: particle per strand
* **Int** a: frame number  
* **For Loop** \* a  
  * INT: frame id
  * FLOAT * 16: rigid motion
  * **For Loop** \* b  
    * FLOAT \* 3: position
  * **For Loop** \* b  
    * FLOAT \* 3: direction
