from setuptools import setup


setup(name='bardic',
      version='0.1.0',
      description='A binomial RNA-DNA interaction caller.',
      url='http://github.com/dmitrymyl/BaRDIC',
      author='Dmitry Mylarshchikov',
      author_email='dmitrymyl@gmail.com',
      license='GPL-3.0',
      packages=['bardic'],
      install_requires=['numpy',
                        'scipy',
                        'pandas',
                        'bioframe>=0.3.0',
                        'ray',
                        'tqdm',
                        'matplotlib',
                        'seaborn',
                        'statsmodels'],
      python_requires='>=3.7',
      zip_safe=True)