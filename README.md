
This is the data repository of the article *Python code for modeling multi-layer structures with controlled cracking and delamination* submitted to Software Impacts.

### Contents

The Abaqus Python code `RVE_soft_islands.py`  generates representative volume elements of multi-layer composites structured with embedded islands and stable cracks. It provides the framework to study the composites's stretchability and deformation mechanisms using finite element modelling.

The script is ready to use in Abaqus CAE finite elements software and generates a fully parameterized representative volume element of a soft islands-patterned microstructure. Geometry, material and discretization parameters can be changed by the user in the `Variables` section of the script.

The code automatically generates RVEs for all input parameter combinations and handles meshing and application of the boundary conditions as specified.

### Further information

For more information on the code, the reader is referred to the article [*Python code for modeling multi-layer structures with controlled cracking and delamination*]().

The code has been used to investigate structure-property relationships in soft islands-stable crack-patterned composites and their delamination behaviour in:
* P. Kowol, S. Bargmann, P. Görrn, and J. Wilmers. *Strain relief by controlled cracking in highly stretchable multi-layer composites*. Extreme Mechanics Letters 54:101724, 2022. [https://doi.org/10.1016/j.eml.2022.101724](https://doi.org/10.1016/j.eml.2022.101724)
* P. Kowol, S. Bargmann, P. Görrn, and J. Wilmers. *Delamination behavior of highly stretchable soft islands multi-layer materials*. Applied Mechanics 4(2), 514-527, 2023. [https://doi.org/10.3390/applmech4020029](https://doi.org/10.3390/applmech4020029)



### Licence

This project is licensed under the terms of the MIT license.
