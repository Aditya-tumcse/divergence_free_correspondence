#!/bin/bash
set -e
echo "Setting up RTK..."
curl -fsSL https://raw.githubusercontent.com/rtk-ai/rtk/master/install.sh | bash
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
export PATH="$HOME/.local/bin:$PATH"
export RTK_TELEMETRY_DISABLED=1
rtk init -g --auto-patch
echo "RTK setup complete ($(rtk --version))"