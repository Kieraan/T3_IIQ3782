### IIQ3782 - Termodinámica Avanzada


El objetivo de este curso es presentar a los estudiantes del departamento ingeniería química y bioprocesos, un conjunto interconectado de métodos computacionales necesarios en para la resolución de problemas de termodinámica avanzada. En particular, métodos que ayudan a los estudiantes a resolver problemas en las áreas equilibrios de fase utilizando ecuaciones de estado y modelos de coeficientes de actividad. 

El curso utilizará Jupyter notebooks para las experiencias (numerados del 00 al 06) para que los estudiantes practiquen sus habilidades en termodinámica computacional, usando la resolución de problemas (Ver `notebooks/`). Los notebooks están preparadas como una combinación de métodos computacionales y programación informática (en lenguaje Python) destinado a ayudar a los estudiantes para la realización de tareas y el proyecto del curso.


Los comentarios y la colaboración para mejorar este curso son bienvenidos en GitHub con `pull requests` o `issues` por correo electrónico.

Este curso utiliza Jupyter Notebooks en el lenguaje de programación Python. Se puede acceder al contenido en
las siguientes maneras:

+ La versión HTML estática de los notebooks se mostrará en el navegador actual si el notebook figura en el repositorio de código, en `notebook/`. Esto no permitirá representar siempre las fórmulas matemáticas o interactuar con el código. Alternativamente, puede renderizar los cuadernos en [NBViewer](http://nbviewer.jupyter.org/).
+ De manera local para el cual debe seguir la siguiente guía

Gracias de antemano por los aportes para mejorar este curso.\
Saludos,\
Equipo docente IIQ3782

https://docs.anaconda.com/free/miniconda/
https://visualstudio.microsoft.com/visual-cpp-build-tools/

conda env create --name iiq3782 --file=environment.yml
conda activate iiq3782
jupyter-notebook 01_Introduccion.ipynb
git clone https://github.com/gustavochm/sgtpy
cd sgtpy
pip install .
cd ../
cd epcsaftpy
pip install .


