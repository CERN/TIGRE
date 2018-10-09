from __future__ import absolute_import
import  os
from os.path import join as pjoin
from setuptools import setup, find_packages
from distutils.extension import Extension
from Cython.Distutils import build_ext
from distutils.spawn import spawn, find_executable
import numpy
from sys import platform
PATH = os.environ.get('PATH')

# Code from https://github.com/rmcgibbo/npcuda-example/blob/master/cython/setup.py


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
    if platform == "win32":
        if 'CUDA_PATH' in os.environ:
            home = os.environ['CUDA_PATH']
            nvcc = pjoin(home, 'bin', 'nvcc.exe')
        else:
            raise EnvironmentError('CUDA_PATH could not be found in your environment variables.')
        home = os.path.dirname(os.path.dirname(nvcc))
        cudaconfig = {'home': home, 'nvcc': nvcc,
                      'include': pjoin(home, 'include'),
                      'x64': pjoin(home,'lib', 'x64')}
        for k, v in cudaconfig.iteritems():
            if not os.path.exists(v):
                raise EnvironmentError('The CUDA %s path could not be located in %s' % (k, v))


    #defaulting to linux for now.
    else:

    # first check if the CUDAHOME env variable is in use
        if 'CUDAHOME' in os.environ:
            home = os.environ['CUDAHOME']
            nvcc = pjoin(home, 'bin', 'nvcc')
        else:
            # otherwise, search the PATH for NVCC
            nvcc = find_in_path('nvcc', os.environ['PATH'])
            if nvcc is None:
                raise EnvironmentError('The nvcc binary could not be located in your $PATH. '
                                       'Either add it to your path, or set $CUDAHOME')
            home = os.path.dirname(os.path.dirname(nvcc))

        cudaconfig = {'home': home, 'nvcc': nvcc,
                      'include': pjoin(home, 'include'),
                      'lib64': pjoin(home, 'lib64')}
        for k, v in cudaconfig.iteritems():
            if not os.path.exists(v):
                raise EnvironmentError('The CUDA %s path could not be located in %s' % (k, v))

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

    self.set_executable('compiler_so', 'nvcc')
    self.set_executable('linker_so', CUDA['nvcc'])
    self._c_extensions.append('.cu')
    self._cpp_extensions.append('.cpp')
    # save references to the default compiler_so and _comple methods
    default_compiler_so = self.compiler_so
    super = self._compile

    # now redefine the _compile method. This gets executed for each
    # object but distutils doesn't have the ability to change compilers
    # based on source extension: we add it.

    def _compile(obj, src, ext, cc_args, extra_postargs, pp_opts):
        if os.path.splitext(src)[1] == '.cu':
            # use the cuda for .cu files
            self.set_executable('compiler_so', 'nvcc')
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

#NOTE: again, defaulting to linux.
lib_string = 'lib64'
extra_comp_args = {'gcc': [], 'nvcc':['-arch=sm_30', '--ptxas-options=-v', '-c']}
if platform == 'win32':
    lib_string = 'x64'
    extra_comp_args = []

Ax_ext = Extension('_Ax',
                   sources=(['tigre/Source/projection.cpp',
                             'tigre/Source/Siddon_projection.cu', 'tigre/Source/Siddon_projection_parallel.cu',
                             'tigre/Source/ray_interpolated_projection.cu', 'tigre/Source/ray_interpolated_projection_parallel.cu',
                             'tigre/Source/_types.pxd',
                             'tigre/Source/_Ax.pyx']),
                   library_dirs=[CUDA[lib_string]],
                   libraries=['cudart'],
                   language='c++',
                   #runtime_library_dirs=[CUDA[lib_string]],
                   # this syntax is specific to this build system
                   # we're only going to use certain compiler args with nvcc and not with gcc
                   # the implementation of this trick is in customize_compiler() below
                   extra_compile_args=extra_comp_args,
                                        #'nvcc': ['-arch=sm_30', '--ptxas-options=-v', '-c',
                                         #        '--compiler-options', '-fPIC', '-Wall', '-Wfatal-errors']},
                   include_dirs=[numpy_include, CUDA['include'], 'Source'])

Atb_ext = Extension('_Atb',
                    sources=(['tigre/Source/voxel_backprojection.cu', 'tigre/Source/voxel_backprojection2.cu',
                              'tigre/Source/voxel_backprojection_parallel.cu',
                              'tigre/Source/_types.pxd',
                              'tigre/Source/_Atb.pyx']),
                    library_dirs=[CUDA[lib_string]],
                    libraries=['cudart'],
                    language='c++',
                    #runtime_library_dirs=[CUDA[lib_string]],
                    # this syntax is specific to this build system
                    # we're only going to use certain compiler args with nvcc and not with gcc
                    # the implementation of this trick is in customize_compiler() below
                    extra_compile_args=extra_comp_args,
                                                  #'--compiler-options', "'-fPIC'"]},
                    include_dirs=[numpy_include, CUDA['include'], 'tigre/Source'])
tvdenoising_ext = Extension('_tvdenoising',
                    sources=(['tigre/Source/voxel_backprojection.cu', 'tigre/Source/tvdenoising.cu',
                              'tigre/Source/_types.pxd',
                              'tigre/Source/_tvdenoising.pyx']),
                    library_dirs=[CUDA[lib_string]],
                    libraries=['cudart'],
                    language='c++',
                    #runtime_library_dirs=[CUDA[lib_string]],
                    # this syntax is specific to this build system
                    # we're only going to use certain compiler args with nvcc and not with gcc
                    # the implementation of this trick is in customize_compiler() below
                    extra_compile_args=extra_comp_args,
                                                  #'--compiler-options', "'-fPIC'"
                    include_dirs=[numpy_include, CUDA['include'], 'Source'])

# run the customize_compiler
class custom_build_ext(build_ext):
    def build_extensions(self):

        customize_compiler_for_nvcc(self.compiler)
        self.compiler.spawn = self.spawn
        build_ext.build_extensions(self)

    def spawn(self, cmd, search_path=1, verbose=1, dry_run=0):
        """
        Perform any CUDA specific customizations before actually launching
        compile/link etc. commands.
        """
        print(cmd)

        if self.compiler.compiler_type == 'msvc':
            # There are several things we need to do to change the commands
            # issued by MSVCCompiler into one that works with nvcc. In the end,
            # it might have been easier to write our own CCompiler class for
            # nvcc, as we're only interested in creating a shared library to
            # load with ctypes, not in creating an importable Python extension.
            # - First, we replace the cl.exe or link.exe call with an nvcc
            #   call. In case we're running Anaconda, we search cl.exe in the
            #   original search path we captured further above -- Anaconda
            #   inserts a MSVC version into PATH that is too old for nvcc.
            cmd[:1] = ['nvcc', '--compiler-bindir',
                       os.path.dirname(find_executable("cl.exe", PATH))
                       or cmd[0]]

            # - Secondly, we fix a bunch of command line arguments.

            for idx, c in enumerate(cmd):
                # create .dll instead of .pyd files
                #if '.pyd' in c: cmd[idx] = c = c.replace('.pyd', '.dll')
                # replace /c by -c
                if c == '/c': cmd[idx] = '-c'
                # replace /DLL by --shared
                elif c == '/DLL': cmd[idx] = '--shared'
                # remove --compiler-options=-fPIC
                elif '-fPIC' in c: del cmd[idx]
                # replace /Tc... by ...
                elif c.startswith('/Tc'): cmd[idx] = c[3:]
                elif c.startswith('/Tp'):
                    cmd[idx] = c[3:]
                # replace /Fo... by -o ...
                elif c.startswith('/Fo'): cmd[idx:idx+1] = ['-o', c[3:]]
                # replace /LIBPATH:... by -L...
                elif c.startswith('/LIBPATH:'): cmd[idx] = '-L' + c[9:]
                # replace /OUT:... by -o ...
                elif c.startswith('/OUT:'): cmd[idx:idx+1] = ['-o', c[5:]]
                # remove /EXPORT:initlibcudamat or /EXPORT:initlibcudalearn
                elif c.startswith('/EXPORT:'): del cmd[idx]
                # replace cublas.lib by -lcublas
                elif c == 'cublas.lib': cmd[idx] = '-lcublas'
                #elif c.endswith('.obj'):cmd[idx] = c[:len(c)-4]+'.cpp'
            # - Finally, we pass on all arguments starting with a '/' to the
            #   compiler or linker, and have nvcc handle all other arguments
            if '--shared' in cmd:
                pass_on = '--linker-options='
                # we only need MSVCRT for a .dll, remove CMT if it sneaks in:
                cmd.append('/NODEFAULTLIB:libcmt.lib')
            else:
                pass_on = '--compiler-options='
            cmd = ([c for c in cmd if c[0] != '/'] +
                   [pass_on + ','.join(c for c in cmd if c[0] == '/')])
            # For the future: Apart from the wrongly set PATH by Anaconda, it
            # would suffice to run the following for compilation on Windows:
            # nvcc -c -O -o <file>.obj <file>.cu
            # And the following for linking:
            # nvcc --shared -o <file>.dll <file1>.obj <file2>.obj -lcublas
            # This could be done by a NVCCCompiler class for all platforms.

        print(cmd)
        spawn(cmd, search_path, verbose, dry_run)



setup(name='tigre',
      version = '0.0.2',
      author = 'Reuben Lindroos, Sam loescher',
      packages = find_packages(),
      include_package_data=True,
      ext_modules=[Ax_ext, Atb_ext,tvdenoising_ext],

      # inject our custom trigger
      cmdclass={'build_ext': custom_build_ext},

      # since the package has c code, the egg cannot be zipped
      zip_safe=False)
