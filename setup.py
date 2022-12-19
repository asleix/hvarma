import os
from setuptools import setup
from setuptools.command.install import install
import subprocess


def compile_and_copy_c_files():
    """Use the subprocess module to compile the C software."""
    src_path = 'hvarma/ext_c/'
    subprocess.check_call('make', cwd=src_path, shell=True)

    assert os.path.exists('hvarma/ext_c/gradient.so'), "Make sure shared object exists"

    # move into package directory
    import site, shutil
    package_dirs = site.getsitepackages()
    for pkg_dir in package_dirs:
        if 'hvarma' not in os.listdir(pkg_dir):
            os.mkdir(os.path.join(pkg_dir, 'hvarma'))
        ext_dir = os.path.join(pkg_dir, 'hvarma', 'c_ext')
        if not os.path.exists(ext_dir):
            os.mkdir(ext_dir)
        lib_path = os.path.join(ext_dir, 'gradient.so')
        if not os.path.exists(lib_path):
            shutil.copy('hvarma/ext_c/gradient.so', ext_dir)


class CustomInstall(install):
    """Custom handler for the 'install' command."""
    def run(self):
        compile_and_copy_c_files()
        super().run()


setup (name = 'hvarma',
       version = '1.0',
       description = 'This is a demo package',
       author = 'Aleix Segui',
       author_email = 'aleix.segui@estudiantat.upc.edu',
       url = 'https://docs.python.org/extending/building',
       long_description = '''
This is really just a demo package.
''',
       packages=['hvarma'],
       install_requires=[
        'numpy>=1.19.0',
        'scipy>=1.5.1',
        'obspy>=1.2.2',
        'matplotlib>=3.3.0',
        'pytest'
        ],
        cmdclass={'install': CustomInstall},
)