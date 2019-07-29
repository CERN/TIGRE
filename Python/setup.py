from __future__ import absolute_import
import os
from os.path import join as pjoin
from setuptools import setup, find_packages, Extension
# from distutils.extension import Extension
from Cython.Distutils import build_ext
import subprocess
import numpy
import sys

# Code from https://github.com/rmcgibbo/npcuda-example/blob/master/cython/setup.py
compute_capability_args = [  # '-gencode=arch=compute_20,code=sm_20',
    '-gencode=arch=compute_30,code=sm_30',
    '-gencode=arch=compute_37,code=sm_37',
    '-gencode=arch=compute_52,code=sm_52',
    '-gencode=arch=compute_60,code=sm_60',
     '-gencode=arch=compute_61,code=sm_61',
    # '-gencode=arch=compute_70,code=sm_70', #untested
    '--ptxas-options=-v', '-c',
    '--default-stream=per-thread',
    '--compiler-options', "'-fPIC'"]


def find_in_path(name, path):
    "Find a file in a search path"
    # adapted fom http://code.activestate.com/recipes/52224-find-a-file-given-a-search-path/
    for dir in path.split(os.pathsep):
        binpath = pjoin(dir, name)
        if os.path.exists(binpath):
            return os.path.abspath(binpath)
    return None


def locate_cuda():
    """Locate the CUDA environment on the system

    Returns a dict with keys 'home', 'nvcc', 'include', and 'lib64'
    and values giving the absolute path to each directory.

    Starts by looking for the CUDAHOME env variable. If not found, everything
    is based on finding 'nvcc' in the PATH.
    """

    # first check if the CUDAHOME env variable is in use
    if 'CUDAHOME' in os.environ:
        home = os.environ['CUDAHOME']
        nvcc = pjoin(home, 'bin', 'nvcc')
    else:
        # otherwise, search the PATH for NVCC
        nvcc = find_in_path('nvcc', os.environ['PATH'])
        if nvcc is None:
            nvcc_is_in_path = find_in_path('/usr/local/cuda/bin', os.environ['PATH'])
            if not nvcc_is_in_path:
                print('WARNING: The nvcc binary could not be located in your $PATH. '
                      'Either add it to your path, or set $CUDAHOME')
                # raise EnvironmentError(
                #     'The CUDA  path could not be located in $PATH or $CUDAHOME. Either add it to your path, or set $CUDAHOME')
            home = '/usr/local/cuda'
            nvcc = pjoin(home, 'bin', 'nvcc')
        else:
            home = os.path.dirname(os.path.dirname(nvcc))

    cudaconfig = {'home': home, 'nvcc': nvcc,
                  'include': pjoin(home, 'include'),
                  'lib64': pjoin(home, 'lib64')}
    for k, v in cudaconfig.items():
        if not os.path.exists(v):
            raise EnvironmentError(
                'The CUDA  path could not be located in $PATH or $CUDAHOME. Either add it to your path, or set $CUDAHOME')

    return cudaconfig


CUDA = locate_cuda()

# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()


def customize_compiler_for_nvcc(self):
    """inject deep into distutils to customize how the dispatch
    to gcc/nvcc works.

    If you subclass UnixCCompiler, it's not trivial to get your subclass
    injected in, and still have the right customizations (i.e.
    distutils.sysconfig.customize_compiler) run on it. So instead of going
    the OO route, I have this. Note, it's kindof like a wierd functional
    subclassing going on."""

    # tell the compiler it can processes .cu
    self.src_extensions.append('.cu')

    # save references to the default compiler_so and _comple methods
    default_compiler_so = self.compiler_so
    super = self._compile

    # now redefine the _compile method. This gets executed for each
    # object but distutils doesn't have the ability to change compilers
    # based on source extension: we add it.
    def _compile(obj, src, ext, cc_args, extra_postargs, pp_opts):
        if os.path.splitext(src)[1] == '.cu':
            # use the cuda for .cu files
            self.set_executable('compiler_so', CUDA['nvcc'])
            # use only a subset of the extra_postargs, which are 1-1 translated
            # from the extra_compile_args in the Extension class
            postargs = extra_postargs['nvcc']
        else:
            postargs = extra_postargs['gcc']

        super(obj, src, ext, cc_args, postargs, pp_opts)
        # reset the default compiler_so, which we might have changed for cuda
        self.compiler_so = default_compiler_so

    # inject our redefined _compile method into the class
    self._compile = _compile


def include_headers(filename_list, sdist=False):
    """add hpp and h files to list if sdist is called"""
    c_extensions = ['.cu', ".c", ".C", ".cc", ".cpp", ".cxx", ".c++"]
    header_list = []
    for filename in filename_list:
        if sdist:
            header = list(os.path.splitext(filename))
            if header[1] in c_extensions:
                header[1] = '.hpp'
                header_list.append(''.join(header))
    if sdist:
        filename_list += ['tigre/Source/types_TIGRE.hpp','tigre/Source/errors.hpp']
    return filename_list + header_list


Ax_ext = Extension('_Ax',
                   sources=(include_headers(['tigre/Source/projection.cpp',
                                             'tigre/Source/Siddon_projection.cu',
                                             'tigre/Source/Siddon_projection_parallel.cu',
                                             'tigre/Source/ray_interpolated_projection.cu',
                                             'tigre/Source/ray_interpolated_projection_parallel.cu',
                                             'tigre/Source/_types.pxd',
                                             'tigre/Source/_Ax.pyx'], sdist=(sys.argv[1] == "sdist"))),
                   library_dirs=[CUDA['lib64']],
                   libraries=['cudart'],
                   language='c++',
                   runtime_library_dirs=[CUDA['lib64']],
                   # this syntax is specific to this build system
                   # we're only going to use certain compiler args with nvcc and not with gcc
                   # the implementation of this trick is in customize_compiler() below
                   extra_compile_args={'gcc': ['-O2'],
                                       'nvcc': compute_capability_args},
                   include_dirs=[numpy_include, CUDA['include'], 'Source'])

Atb_ext = Extension('_Atb',
                    sources=(include_headers(['tigre/Source/voxel_backprojection.cu',
                                              'tigre/Source/voxel_backprojection2.cu',
                                              'tigre/Source/voxel_backprojection_parallel.cu',
                                              'tigre/Source/_types.pxd',
                                              'tigre/Source/_Atb.pyx'], sdist=(sys.argv[1] == "sdist"))),
                    library_dirs=[CUDA['lib64']],
                    libraries=['cudart'],
                    language='c++',
                    runtime_library_dirs=[CUDA['lib64']],
                    # this syntax is specific to this build system
                    # we're only going to use certain compiler args with nvcc and not with gcc
                    # the implementation of this trick is in customize_compiler() below
                    extra_compile_args={'gcc': [],
                                        'nvcc': compute_capability_args},
                    include_dirs=[numpy_include, CUDA['include'], 'tigre/Source'])
tvdenoising_ext = Extension('_tvdenoising',
                            sources=(
                                include_headers(['tigre/Source/tvdenoising.cu',
                                                 'tigre/Source/_types.pxd',
                                                 'tigre/Source/_tvdenoising.pyx'], sdist=(sys.argv[1] == "sdist"))),
                            library_dirs=[CUDA['lib64']],
                            libraries=['cudart'],
                            language='c++',
                            runtime_library_dirs=[CUDA['lib64']],
                            # this syntax is specific to this build system
                            # we're only going to use certain compiler args with nvcc and not with gcc
                            # the implementation of this trick is in customize_compiler() below
                            extra_compile_args={'gcc': [],
                                                'nvcc': compute_capability_args},
                            include_dirs=[numpy_include, CUDA['include'], 'Source'])
minTV_ext = Extension('_minTV',
                      sources=(include_headers(['tigre/Source/POCS_TV.cu',
                                                'tigre/Source/_types.pxd',
                                                'tigre/Source/_minTV.pyx'], sdist=(sys.argv[1] == "sdist"))),
                      library_dirs=[CUDA['lib64']],
                      libraries=['cudart'],
                      language='c++',
                      runtime_library_dirs=[CUDA['lib64']],
                      # this syntax is specific to this build system
                      # we're only going to use certain compiler args with nvcc and not with gcc
                      # the implementation of this trick is in customize_compiler() below
                      extra_compile_args={'gcc': [],
                                          'nvcc': compute_capability_args},
                      include_dirs=[numpy_include, CUDA['include'], 'Source'])

AwminTV_ext = Extension('_AwminTV',
                        sources=(include_headers(['tigre/Source/POCS_TV2.cu',
                                                  # 'tigre/Source/_types.pxd',
                                                  'tigre/Source/_AwminTV.pyx'], sdist=(sys.argv[1] == "sdist"))),
                        library_dirs=[CUDA['lib64']],
                        libraries=['cudart'],
                        language='c++',
                        runtime_library_dirs=[CUDA['lib64']],
                        # this syntax is specific to this build system
                        # we're only going to use certain compiler args with nvcc and not with gcc
                        # the implementation of this trick is in customize_compiler() below
                        extra_compile_args={'gcc': [],
                                            'nvcc': compute_capability_args},
                        include_dirs=[numpy_include, CUDA['include'], 'Source'])


# run the customize_compiler
class custom_build_ext(build_ext):
    def build_extensions(self):
        customize_compiler_for_nvcc(self.compiler)
        build_ext.build_extensions(self)


setup(name='pytigre',
      version='0.1.7',
      author='Reuben Lindroos, Sam loescher',
      packages=find_packages(),
      scripts=['tigre/demos/launch.sh',
               'tests/runscript.sh'],
      include_package_data=True,
      ext_modules=[Ax_ext, Atb_ext, tvdenoising_ext, minTV_ext, AwminTV_ext],
      py_modules=['tigre.py'],
      # inject our custom trigger
      cmdclass={'build_ext': custom_build_ext},
      license_file='LICENSE.txt',
      license='BSD 3-Clause',
      # since the package has c code, the egg cannot be zipped
      zip_safe=False)
