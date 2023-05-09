from setuptools import setup, find_packages


setup(name='bardic',
      version='0.4.1',
      description='A binomial RNA-DNA interaction caller.',
      url='http://github.com/dmitrymyl/BaRDIC',
      author='Dmitry Mylarshchikov',
      author_email='dmitrymyl@gmail.com',
      license='MIT',
      packages=find_packages(),
      install_requires=['numpy',
                        'scipy',
                        'pandas',
                        'bioframe>=0.3.0',
                        'tqdm',
                        'statsmodels',
                        'h5py'],
      python_requires='>=3.8',
      entry_points={'console_scripts':
                    ['bardic=bardic.__main__:main']},
      zip_safe=True)
