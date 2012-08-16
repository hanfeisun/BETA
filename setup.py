#!/usr/bin/env python

import os
import sys
from distutils.core import setup, Extension

def main():
    if not float(sys.version[:3])>=2.5:
        sys.stderr.write("CRITICAL: Python version must be greater than or equal to 2.5! python 2.6.1 or newer is recommended!\n")
        sys.exit(1)
    setup(name="BETA-Package",
          version="0.1.5",
          description="BETA -- Binding and Expression Targets Analysis ",
          author='Su Wang',
          author_email='wangsu0623@gmail.com',
          package_dir={'BETA' : 'BETA'},
          packages=['BETA'],
          scripts=['bin/BETA'],

          classifiers=[
            'Development Status :: 4 - Beta',
            'Environment :: Console',
            'Environment :: Web Environment',
            'Intended Audience :: Developers',
            'License :: OSI Approved :: Artistic License',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: POSIX',
            'Programming Language :: Python',
            'Topic :: Database',
            ],
          )

if __name__ == '__main__':
    main()
