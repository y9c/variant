from setuptools import setup, Extension

# Define the C extension
seqpy_extension = Extension(
    'variant.seqpy',
    sources=['variant/seqpy.c']
)

setup(
    ext_modules=[seqpy_extension],
)
