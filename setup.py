from setuptools import setup

setup(name='exact_stats',
      version='0.1',
      description='Exact statistics for asssociation studies dealing with rare phenomena in which approximation models such as logistic regression are not valid.  Designed for genetic analysis such as GWAS and PheWAS, but applicable to any categorical case-control type studies.',
      url='https://github.com/phewas/exact_stats/',
      author='Joseph Bochenek',
      author_email='jbochenek@mykolab.ch',
      license='MIT',
      packages=['exact_stats'],
      zip_safe=True)
