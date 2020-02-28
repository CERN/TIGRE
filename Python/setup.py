from __future__ import absolute_import
import os
from os.path import join as pjoin
from setuptools import setup, find_packages, Extension
# from distutils.extension import Extension
from Cython.Distutils import build_ext
import subprocess
import numpy
import sys
import copy
import re

IS_WINDOWS = sys.platform == 'win32'


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
    ]


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

    Returns a dict with keys 'home', 'include', and 'lib64'
    and values giving the absolute path to each directory.

    Starts by looking for the CUDAHOME env variable. If not found, everything
    is based on finding 'nvcc' in the PATH.
    """

    '''Finds the CUDA install path.'''
    # Guess #1
    cuda_home = os.environ.get('CUDA_HOME') or os.environ.get('CUDA_PATH')
    if cuda_home is None:
        # Guess #2
        try:
            which = 'where' if IS_WINDOWS else 'which'
            nvcc = subprocess.check_output(
                [which, 'nvcc']).decode().rstrip('\r\n')
            cuda_home = os.path.dirname(os.path.dirname(nvcc))
        except Exception:
            # Guess #3
            if IS_WINDOWS:
                cuda_homes = glob.glob(
                    'C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v*.*')
                if len(cuda_homes) == 0:
                    cuda_home = ''
                else:
                    cuda_home = cuda_homes[0]
            else:
                cuda_home = '/usr/local/cuda'
            if not os.path.exists(cuda_home):
                cuda_home = None

    cudaconfig = {'home': cuda_home, 
                  'include': pjoin(cuda_home, 'include'),
                  'lib64': pjoin(cuda_home, 'lib/x64' if IS_WINDOWS else 'lib64')}
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


def _is_cuda_file(path):
    return os.path.splitext(path)[1] in ['.cu', '.cuh']


COMMON_MSVC_FLAGS = ['/MD', '/wd4819', '/EHsc']

COMMON_NVCC_FLAGS = [
    '-D__CUDA_NO_HALF_OPERATORS__',
    '-D__CUDA_NO_HALF_CONVERSIONS__',
    '-D__CUDA_NO_HALF2_OPERATORS__',
    '--expt-relaxed-constexpr'
]


def _join_cuda_home(*paths):
    return os.path.join(CUDA['home'], *paths)
	

class BuildExtension(build_ext, object):
    '''
    A custom :mod:`setuptools` build extension .

    This :class:`setuptools.build_ext` subclass takes care of passing the
    minimum required compiler flags (e.g. ``-std=c++11``) as well as mixed
    C++/CUDA compilation (and support for CUDA files in general).

    When using :class:`BuildExtension`, it is allowed to supply a dictionary
    for ``extra_compile_args`` (rather than the usual list) that maps from
    languages (``cxx`` or ``nvcc``) to a list of additional compiler flags to
    supply to the compiler. This makes it possible to supply different flags to
    the C++ and CUDA compiler during mixed compilation.
    '''

    @classmethod
    def with_options(cls, **options):
        '''
        Returns an alternative constructor that extends any original keyword
        arguments to the original constructor with the given options.
        '''
        def init_with_options(*args, **kwargs):
            kwargs = kwargs.copy()
            kwargs.update(options)
            return cls(*args, **kwargs)
        return init_with_options

    def __init__(self, *args, **kwargs):
        super(BuildExtension, self).__init__(*args, **kwargs)
        self.no_python_abi_suffix = kwargs.get("no_python_abi_suffix", False)

    def build_extensions(self):
        # Register .cu and .cuh as valid source extensions.
        self.compiler.src_extensions += ['.cu', '.cuh']
        # Save the original _compile method for later.
        if self.compiler.compiler_type == 'msvc':
            self.compiler._cpp_extensions += ['.cu', '.cuh']
            original_compile = self.compiler.compile
            original_spawn = self.compiler.spawn
        else:
            original_compile = self.compiler._compile

        def unix_wrap_compile(obj, src, ext, cc_args, extra_postargs, pp_opts):
            # Copy before we make any modifications.
            cflags = copy.deepcopy(extra_postargs)
            try:
                original_compiler = self.compiler.compiler_so
                if _is_cuda_file(src):
                    nvcc = _join_cuda_home('bin', 'nvcc')
                    if not isinstance(nvcc, list):
                        nvcc = [nvcc]
                    self.compiler.set_executable('compiler_so', nvcc)
                    if isinstance(cflags, dict):
                        cflags = cflags['nvcc']
                    cflags = COMMON_NVCC_FLAGS + ['--compiler-options',
                                                  "'-fPIC'"] + cflags + compute_capability_args
                elif isinstance(cflags, dict):
                    cflags = cflags['cxx']
                # NVCC does not allow multiple -std to be passed, so we avoid
                # overriding the option if the user explicitly passed it.
                if not any(flag.startswith('-std=') for flag in cflags):
                    cflags.append('-std=c++11')

                original_compile(obj, src, ext, cc_args, cflags, pp_opts)
            finally:
                # Put the original compiler back in place.
                self.compiler.set_executable('compiler_so', original_compiler)

        def win_wrap_compile(sources,
                             output_dir=None,
                             macros=None,
                             include_dirs=None,
                             debug=0,
                             extra_preargs=None,
                             extra_postargs=None,
                             depends=None):

            self.cflags = copy.deepcopy(extra_postargs)
            extra_postargs = None

            def spawn(cmd):
                # Using regex to match src, obj and include files
                src_regex = re.compile('/T(p|c)(.*)')
                src_list = [
                    m.group(2) for m in (src_regex.match(elem) for elem in cmd)
                    if m
                ]

                obj_regex = re.compile('/Fo(.*)')
                obj_list = [
                    m.group(1) for m in (obj_regex.match(elem) for elem in cmd)
                    if m
                ]

                include_regex = re.compile(r'((\-|\/)I.*)')
                include_list = [
                    m.group(1)
                    for m in (include_regex.match(elem) for elem in cmd) if m
                ]

                if len(src_list) >= 1 and len(obj_list) >= 1:
                    src = src_list[0]
                    obj = obj_list[0]
                    if _is_cuda_file(src):
                        nvcc = _join_cuda_home('bin', 'nvcc')
                        if isinstance(self.cflags, dict):
                            cflags = self.cflags['nvcc']
                        elif isinstance(self.cflags, list):
                            cflags = self.cflags
                        else:
                            cflags = []

                        cflags = COMMON_NVCC_FLAGS + cflags + compute_capability_args
                        for flag in COMMON_MSVC_FLAGS:
                            cflags = ['-Xcompiler', flag] + cflags
                        cmd = [nvcc, '-c', src, '-o', obj] + include_list + cflags
                    elif isinstance(self.cflags, dict):
                        cflags = COMMON_MSVC_FLAGS #+ self.cflags['cxx']
                        cmd += cflags
                    elif isinstance(self.cflags, list):
                        cflags = COMMON_MSVC_FLAGS + self.cflags
                        cmd += cflags

                return original_spawn(cmd)

            try:
                self.compiler.spawn = spawn
                return original_compile(sources, output_dir, macros,
                                        include_dirs, debug, extra_preargs,
                                        extra_postargs, depends)
            finally:
                self.compiler.spawn = original_spawn

        # Monkey-patch the _compile method.
        if self.compiler.compiler_type == 'msvc':
            self.compiler.compile = win_wrap_compile
        else:
            self.compiler._compile = unix_wrap_compile

        build_ext.build_extensions(self)

    def get_ext_filename(self, ext_name):
        # Get the original shared library name. For Python 3, this name will be
        # suffixed with "<SOABI>.so", where <SOABI> will be something like
        # cpython-37m-x86_64-linux-gnu. On Python 2, there is no such ABI name.
        # The final extension, .so, would be .lib/.dll on Windows of course.
        ext_filename = super(BuildExtension, self).get_ext_filename(ext_name)
        # If `no_python_abi_suffix` is `True`, we omit the Python 3 ABI
        # component. This makes building shared libraries with setuptools that
        # aren't Python modules nicer.
        if self.no_python_abi_suffix and sys.version_info >= (3, 0):
            # The parts will be e.g. ["my_extension", "cpython-37m-x86_64-linux-gnu", "so"].
            ext_filename_parts = ext_filename.split('.')
            # Omit the second to last element.
            without_abi = ext_filename_parts[:-2] + ext_filename_parts[-1:]
            ext_filename = '.'.join(without_abi)
        return ext_filename

    def _add_compile_flag(self, extension, flag):
        extension.extra_compile_args = copy.deepcopy(extension.extra_compile_args)
        if isinstance(extension.extra_compile_args, dict):
            for args in extension.extra_compile_args.values():
                args.append(flag)
        else:
            extension.extra_compile_args.append(flag)


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
                   runtime_library_dirs=[CUDA['lib64']] if not IS_WINDOWS else None,
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
                    runtime_library_dirs=[CUDA['lib64']] if not IS_WINDOWS else None,
                    include_dirs=[numpy_include, CUDA['include'], 'tigre/Source'])
tvdenoising_ext = Extension('_tvdenoising',
                            sources=(
                                include_headers(['tigre/Source/tvdenoising.cu',
                                                 'tigre/Source/_types.pxd',
                                                 'tigre/Source/_tvdenoising.pyx'], sdist=(sys.argv[1] == "sdist"))),
                            library_dirs=[CUDA['lib64']],
                            libraries=['cudart'],
                            language='c++',
                            runtime_library_dirs=[CUDA['lib64']] if not IS_WINDOWS else None,
                            include_dirs=[numpy_include, CUDA['include'], 'Source'])
minTV_ext = Extension('_minTV',
                      sources=(include_headers(['tigre/Source/POCS_TV.cu',
                                                'tigre/Source/_types.pxd',
                                                'tigre/Source/_minTV.pyx'], sdist=(sys.argv[1] == "sdist"))),
                      library_dirs=[CUDA['lib64']],
                      libraries=['cudart'],
                      language='c++',
                      runtime_library_dirs=[CUDA['lib64']] if not IS_WINDOWS else None,
                      include_dirs=[numpy_include, CUDA['include'], 'Source'])

AwminTV_ext = Extension('_AwminTV',
                        sources=(include_headers(['tigre/Source/POCS_TV2.cu',
                                                  # 'tigre/Source/_types.pxd',
                                                  'tigre/Source/_AwminTV.pyx'], sdist=(sys.argv[1] == "sdist"))),
                        library_dirs=[CUDA['lib64']],
                        libraries=['cudart'],
                        language='c++',
                        runtime_library_dirs=[CUDA['lib64']] if not IS_WINDOWS else None,
                        include_dirs=[numpy_include, CUDA['include'], 'Source'])


setup(name='pytigre',
      version='0.1.8',
      author='Reuben Lindroos, Sam loescher',
      packages=find_packages(),
      scripts=['tigre/demos/launch.sh',
               'tests/runscript.sh'],
      include_package_data=True,
      ext_modules=[Ax_ext, Atb_ext, tvdenoising_ext, minTV_ext, AwminTV_ext],
      py_modules=['tigre.py'],
      cmdclass={'build_ext': BuildExtension},
      license_file='LICENSE.txt',
      license='BSD 3-Clause',
      # since the package has c code, the egg cannot be zipped
      zip_safe=False)
