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