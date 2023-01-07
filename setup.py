from setuptools import setup


setup(name='bardic',
      version='0.3.0',
      description='A binomial RNA-DNA interaction caller.',
      url='http://github.com/dmitrymyl/BaRDIC',
      author='Dmitry Mylarshchikov',
      author_email='dmitrymyl@gmail.com',
      license='MIT',
      packages=['bardic'],
      install_requires=['numpy',
                        'scipy',
                        'pandas',
                        'bioframe>=0.3.0',
                        'tqdm',
                        'matplotlib',
                        'seaborn',
                        'statsmodels',
                        'h5py'],
      python_requires='>=3.8',
      entry_points={'console_scripts':
                    ['bardic=bardic.__main__:run']},
      zip_safe=True)
