from distutils.command.build_ext import build_ext
from distutils.core import Distribution, Extension

# Define your C sources and dependencies
sources = ["variant/seqpy.c"]
include_dirs = []
extra_compile_args = []

# Define the Extension
extension = Extension(
    "seqpy",
    sources=sources,
    include_dirs=include_dirs,
    extra_compile_args=extra_compile_args,
)


def build():
    distribution = Distribution(
        {
            "name": "variant",
            "ext_modules": [extension],
            "cmdclass": {"build_ext": build_ext},
        }
    )

    distribution.package_dir = {"": "variant"}

    # Define the output directory for the compiled extensions
    output_dir = "variant"
    cmd = build_ext(distribution)
    cmd.build_lib = output_dir  # Direct output to variant directory
    cmd.inplace = 1  # Build extensions in place
    cmd.ensure_finalized()
    cmd.run()


if __name__ == "__main__":
    build()
