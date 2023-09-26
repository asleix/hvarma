from setuptools import setup, Extension

cmodule = Extension('gradient', sources=['hvarma/ext_c/gradient.c'])

setup (name = 'hvarma',
       version = '1.0',
       description = 'Horizontal-to-vertical ratio calculator',
       author = 'Aleix Segui',
       author_email = 'aleix.segui@student.ethz.ch',
       url = 'https://docs.python.org/extending/building',
       long_description = 'Horizontal-to-vertical ratio calculator',
       packages=['hvarma', 'hvarma/ext_c'],
       ext_modules = [cmodule],
       install_requires=[
        'numpy>=1.19.0',
        'scipy>=1.5.1',
        'obspy>=1.2.2',
        'matplotlib>=3.3.0',
        'pytest'
        ],
        include_package_data=True,
)