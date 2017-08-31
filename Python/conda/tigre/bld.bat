echo copy distutils.cfg %PREFIX%\LIB\distutils\distutils.cfg
copy distutils.cfg %PREFIX%\LIB\distutils\distutils.cfg
"%PYTHON%" python setup.py build_ext --inplace
if errorlevel 1 exit 1
