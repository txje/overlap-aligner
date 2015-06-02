from distutils.core import setup, Extension

extension_mod = Extension("ovlalign", ["ovlalignermodule.c"])

setup(name = "ovlalign", ext_modules=[extension_mod])
