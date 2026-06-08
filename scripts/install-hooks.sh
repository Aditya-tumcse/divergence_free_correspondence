#!/usr/bin/env bash
# Configure git to use the project-tracked hooks from .githooks/
set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel)"

git config core.hooksPath .githooks
chmod +x "$REPO_ROOT"/.githooks/*

echo "Git hooks installed (core.hooksPath → .githooks/)"
echo "Active hooks:"
ls -1 "$REPO_ROOT/.githooks/" | grep -v '\.sample$' | sed 's/^/  /'
