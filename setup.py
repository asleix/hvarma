import os
from setuptools import setup
from setuptools.command.install import install
import subprocess


def compile_and_copy_c_files():
    """Use the subprocess module to compile the C functions."""
    src_path = 'hvarma/ext_c/'
    subprocess.check_call('make', cwd=src_path, shell=True)

    assert os.path.exists('hvarma/ext_c/gradient.so'), \
        "Make sure shared object exists"


class CustomInstall(install):
    """Custom handler for the 'install' command."""
    def run(self):
        compile_and_copy_c_files()
        super().run()


setup (name = 'hvarma',
       version = '1.0',
       description = 'Horizontal-to-vertical ratio calculator',
       author = 'Aleix Segui',
       author_email = 'aleix.segui@estudiantat.upc.edu',
       url = 'https://docs.python.org/extending/building',
       long_description = 'Horizontal-to-vertical ratio calculator',
       packages=['hvarma'],
       install_requires=[
        'numpy>=1.19.0',
        'scipy>=1.5.1',
        'obspy>=1.2.2',
        'matplotlib>=3.3.0',
        'pytest'
        ],
        cmdclass={'install': CustomInstall},
        include_package_data=True,
)