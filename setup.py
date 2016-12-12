from setuptools import setup

setup(
    name='rest',
    packages=['rest'],
    include_package_data=True,
    install_requires=[
        'flask', 'flask_restful', 'numpy', 'h5py', 'waitress'
    ],
    setup_requires=[
        'pytest-runner',
    ],
    tests_require=[
        'pytest',
    ],
)
