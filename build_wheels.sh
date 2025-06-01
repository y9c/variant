#!/bin/bash

# Script to build wheels locally using uv build
# This will help you test wheel building before pushing to CI

echo "Building wheels locally with uv build..."

# Set environment variables to force binary wheels only
export PIP_ONLY_BINARY=":all:"
export PIP_NO_BUILD_ISOLATION="false"

echo "Environment variables set:"
echo "PIP_ONLY_BINARY=$PIP_ONLY_BINARY"
echo "PIP_NO_BUILD_ISOLATION=$PIP_NO_BUILD_ISOLATION"

# Check if uv is available
if ! command -v uv &> /dev/null; then
    echo "Error: uv is not installed. Please install uv first."
    exit 1
fi

# Clean previous builds
echo "Cleaning previous builds..."
rm -rf dist/ build/ *.egg-info/

# Build both source distribution and wheel
echo "Building source distribution and wheel..."
uv build

# Test the wheel by installing it in a virtual environment
echo "Testing the built wheel..."
WHEEL_FILE=$(ls dist/*.whl | head -n 1)
if [ -f "$WHEEL_FILE" ]; then
    echo "Found wheel: $WHEEL_FILE"
    
    # Create a temporary virtual environment for testing
    TEMP_VENV=$(mktemp -d)
    echo "Creating test environment in $TEMP_VENV"
    
    uv venv "$TEMP_VENV"
    source "$TEMP_VENV/bin/activate"
    
    # Install the wheel with binary-only setting
    echo "Installing wheel (binary packages only)..."
    PIP_ONLY_BINARY=":all:" uv pip install "$WHEEL_FILE"
    
    # Test the C extension
    echo "Testing C extension..."
    python -c "import variant.seqpy; print('C extension test:', variant.seqpy.revcomp('ATCG'))"
    
    # Clean up
    deactivate
    rm -rf "$TEMP_VENV"
    
    echo "Build and test completed successfully!"
    echo "Built files:"
    ls -la dist/
else
    echo "Error: No wheel file found in dist/"
    exit 1
fi
