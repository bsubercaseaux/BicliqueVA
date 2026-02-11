#!/bin/bash
set -e

echo "========================================"
echo "    Installing BicliqueVA"
echo "========================================"

# Check for C++ compiler
if ! command -v g++ &> /dev/null; then
    echo "Error: g++ is not installed. Please install a C++ compiler."
    exit 1
fi

# Check if uv is installed
if ! command -v uv &> /dev/null; then
    echo "uv not found. Installing uv..."
    curl -LsSf https://astral.sh/uv/install.sh | sh
    
    # Add cargo bin to PATH for this session if it exists
    if [ -d "$HOME/.cargo/bin" ]; then
        export PATH="$HOME/.cargo/bin:$PATH"
    fi
    
    # Check again if uv is available
    if ! command -v uv &> /dev/null; then
        echo "Error: Failed to install uv or add it to PATH. Please install uv manually: https://github.com/astral-sh/uv"
        exit 1
    fi
else
    echo "Found uv: $(uv --version)"
fi

echo -e "\n[1/3] Compiling C++ components..."
cd src/bicliqueVA/partition_algs/cpp
make
if [ $? -ne 0 ]; then
    echo "Error: Compilation failed."
    exit 1
fi
cd - > /dev/null

echo -e "\n[2/3] Syncing dependencies..."
uv sync

echo -e "\n[3/3] Installing 'biva' command..."
uv tool install . --force

echo "========================================"
echo "    Installation Complete!"
echo "========================================"
echo "You can now run the 'biva' command from anywhere."
