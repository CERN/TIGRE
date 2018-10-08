# Copyright (c) 2017 Sony Corporation. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function

from setuptools import setup
from distutils.extension import Extension
import os
from os.path import dirname, realpath, join, isfile, splitext
from collections import namedtuple
import copy
import os
import shutil
import sys

root_dir = realpath(dirname(__file__))
a = dict()

__version__ = None
__short_version__ = None
__email__ = None
exec(open(os.path.join(root_dir, 'src', 'nnabla_ext',
                       'cuda', '_version.py')).read(), globals(), a)
if '__version__' in a:
    __version__ = a['__version__']
if '__short_version__' in a:
    __short_version__ = a['__short_version__']
if '__email__' in a:
    __email__ = a['__email__']
assert(__version__ is not None)
assert(__short_version__ is not None)
assert(__email__ is not None)

setup_requires = [
    'numpy>=1.12',
    'Cython>=0.24,<0.26',  # Requires python-dev.
]

whl_suffix = ''
if 'WHEEL_SUFFIX' in os.environ:
    whl_suffix += os.environ['WHEEL_SUFFIX']

install_requires = [
    'setuptools',
    'nnabla{}>={}'.format(whl_suffix, __short_version__),
]

LibInfo = namedtuple('LibInfo', ['file_name', 'path', 'name'])
ExtConfig = namedtuple('ExtConfig',
                       ['package_dir', 'packages', 'package_data',
                        'ext_modules', 'ext_opts'])


def get_libinfo():
    from six.moves.configparser import ConfigParser

    # Parse setup.cfg
    path_cfg = join(dirname(__file__), "setup.cfg")
    if not isfile(path_cfg):
        raise ValueError(
            "`setup.cfg` does not exist. Read installation document and install using CMake.")
    cfgp = ConfigParser()
    cfgp.read(path_cfg)

    # Read cpu lib info
    cpu_lib = LibInfo(None,
                      cfgp.get("cmake", "cpu_target_file"),
                      cfgp.get("cmake", "cpu_target_name"))
    print("CPU Library name:", cpu_lib.name)
    print("CPU Library file:", cpu_lib.path)

    # Read cuda lib info
    cuda_lib = LibInfo(cfgp.get("cmake", "cuda_target_file_name"),
                       cfgp.get("cmake", "cuda_target_file"),
                       cfgp.get("cmake", "cuda_target_name"))
    print("CUDA Library name:", cuda_lib.name)
    print("CUDA Library file name:", cuda_lib.file_name)
    print("CUDA Library file:", cuda_lib.path)
    return cpu_lib, cuda_lib


def get_cpu_extopts(lib):
    import numpy as np
    include_dir = realpath(join(dirname(__file__), '../include'))
    ext_opts = dict(
        include_dirs=[include_dir, np.get_include()],
        libraries=[lib.name],
        library_dirs=[dirname(lib.path)],
        language="c++",
        # The below definition breaks build. Use -Wcpp instead.
        # define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
    )
    if sys.platform != 'win32':
        ext_opts.update(dict(
            extra_compile_args=[
                '-std=c++11', '-Wno-sign-compare', '-Wno-unused-function', '-Wno-cpp'],
            runtime_library_dirs=['$ORIGIN/'],
        ))
    else:
        ext_opts.update(dict(extra_compile_args=['/W0']))
    return ext_opts


def cuda_config(root_dir, cuda_lib, ext_opts):
    # With CUDA
    src_dir = join(root_dir, 'src')
    path_cuda_pkg = join(src_dir, 'nnabla_ext', 'cuda')
    cuda_pkg = "nnabla_ext.cuda"
    package_dir = {cuda_pkg: path_cuda_pkg}
    packages = [cuda_pkg]

    cuda_lib_out = join(path_cuda_pkg, cuda_lib.file_name)
    shutil.copyfile(cuda_lib.path, cuda_lib_out)
    package_data = {cuda_pkg: [cuda_lib.file_name]}

    if sys.platform == 'win32':
        libdir = dirname(cuda_lib.path)
        libname, _ = splitext(cuda_lib.file_name)
        cuda_ext_lib_file_name = libname + '.lib'
        cuda_ext_lib_path = join(libdir, cuda_ext_lib_file_name)
        cuda_ext_lib_out = join(path_cuda_pkg, cuda_ext_lib_file_name)
        shutil.copyfile(cuda_ext_lib_path, cuda_ext_lib_out)
        package_data[cuda_pkg].append(cuda_ext_lib_file_name)

    cuda_ext_opts = copy.deepcopy(ext_opts)
    cuda_ext_opts['libraries'] += [cuda_lib.name]
    cuda_ext_opts['library_dirs'] += [dirname(cuda_lib.path)]
    ext_modules = [
        Extension(cuda_pkg + '.init',
                  [join(path_cuda_pkg, 'init.pyx')],
                  **cuda_ext_opts),
    ]
    return ExtConfig(package_dir, packages, package_data,
                     ext_modules, cuda_ext_opts)


def cudnn_config(root_dir, cuda_lib, cuda_ext_opts):
    src_dir = join(root_dir, 'src')
    path_cudnn_pkg = join(src_dir, 'nnabla_ext', 'cudnn')
    cudnn_pkg = 'nnabla_ext.cudnn'
    package_dir = {cudnn_pkg: path_cudnn_pkg}
    packages = [cudnn_pkg]
    ext_modules = [
        Extension(cudnn_pkg + '.init',
                  [join(path_cudnn_pkg, 'init.pyx')],
                  **cuda_ext_opts),
    ]
    return ExtConfig(package_dir, packages, {},
                     ext_modules, cuda_ext_opts)


def get_setup_config(root_dir):
    cpu_lib, cuda_lib = get_libinfo()

    packages = ['nnabla_ext']
    package_dir = {'nnabla_ext': join(root_dir, 'src', 'nnabla_ext')}
    package_data = {}
    ext_modules = []

    cuda_ext = cuda_config(root_dir, cuda_lib, get_cpu_extopts(cpu_lib))
    packages += cuda_ext.packages
    package_dir.update(cuda_ext.package_dir)
    package_data.update(cuda_ext.package_data)
    ext_modules += cuda_ext.ext_modules

    cudnn_ext = cudnn_config(root_dir, cuda_lib, cuda_ext.ext_opts)
    packages += cudnn_ext.packages
    package_dir.update(cudnn_ext.package_dir)
    package_data.update(cudnn_ext.package_data)
    ext_modules += cudnn_ext.ext_modules

    cuda_version = ''
    if 'WHL_NO_PREFIX' in os.environ and os.environ['WHL_NO_PREFIX'] == 'True':
        cuda_version = ''
    elif 'CUDA_VERSION_MAJOR' in os.environ:
        cuda_version = os.environ['CUDA_VERSION_MAJOR'] + \
            os.environ['CUDA_VERSION_MINOR']
    elif 'CUDAVER' in os.environ:
        cuda_version = os.environ['CUDAVER']

    if 'MULTI_GPU_SUFFIX' in os.environ:
        cuda_version += os.environ['MULTI_GPU_SUFFIX']

    pkg_name = 'nnabla_ext-cuda{}'.format(cuda_version)

    if 'WHEEL_SUFFIX' in os.environ:
        pkg_name += os.environ['WHEEL_SUFFIX']

    pkg_info = dict(
        name=pkg_name,
        description='A CUDA and cuDNN extension of NNabla',
        version=__version__,
        author_email=__email__,
        url="https://github.com/sony/nnabla-ext-cuda",
        license='Apache License 2.0',
        classifiers=[
                'Development Status :: 4 - Beta',
                'Intended Audience :: Developers',
                'Intended Audience :: Education',
                'Intended Audience :: Science/Research',
                'Topic :: Scientific/Engineering',
                'Topic :: Scientific/Engineering :: Artificial Intelligence',
                'License :: OSI Approved :: Apache Software License',
                'Programming Language :: C++',
                'Programming Language :: Python :: 2',
                'Programming Language :: Python :: 2.7',
                'Programming Language :: Python :: 3',
                'Programming Language :: Python :: 3.4',
                'Programming Language :: Python :: 3.5',
                'Programming Language :: Python :: 3.6',
                'Programming Language :: Python :: Implementation :: CPython',
                'Operating System :: Microsoft :: Windows',
                'Operating System :: POSIX :: Linux',
        ],
        keywords="deep learning artificial intelligence machine learning neural network cuda",
        python_requires='>=2.7,!=3.0.*,!=3.1.*,!=3.2.*,!=3.3.*',
    )
    return pkg_info, ExtConfig(package_dir, packages, package_data, ext_modules, {})


if __name__ == '__main__':
    from Cython.Build import cythonize

    pkg_info, cfg = get_setup_config(root_dir)

    # Cythonize
    ext_modules = cythonize(cfg.ext_modules, compiler_directives={
                            "embedsignature": True,
                            "c_string_type": 'str',
                            "c_string_encoding": "ascii"})

    # Setup
    setup(
        #setup_requires=setup_requires,
        #install_requires=install_requires,
        ext_modules=ext_modules,
        package_dir=cfg.package_dir,
        packages=cfg.packages,
        package_data=cfg.package_data,
        namespace_packages=['nnabla_ext'],
        **pkg_info)
