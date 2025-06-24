#!/bin/bash
export PATH="/workspace/.pixi/envs/default/bin:/root/.pixi/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
export CONDA_PREFIX="/workspace/.pixi/envs/default"
export PIXI_PROJECT_VERSION="0.1.0"
export PIXI_PROJECT_ROOT="/workspace"
export PIXI_IN_SHELL="1"
export PIXI_PROJECT_MANIFEST="/workspace/pixi.toml"
export PIXI_PROJECT_NAME="pa_proj"
export PIXI_EXE="/usr/local/bin/pixi"
export CONDA_DEFAULT_ENV="pa_proj"
export PIXI_ENVIRONMENT_NAME="default"
export PIXI_ENVIRONMENT_PLATFORMS="osx-arm64"
export PIXI_PROMPT="(pa_proj) "
. "/workspace/.pixi/envs/default/etc/conda/activate.d/gdal-activate.sh"
. "/workspace/.pixi/envs/default/etc/conda/activate.d/libglib_activate.sh"
. "/workspace/.pixi/envs/default/etc/conda/activate.d/libxml2_activate.sh"
. "/workspace/.pixi/envs/default/etc/conda/activate.d/proj4-activate.sh"

# shellcheck shell=bash
pixi() {
    local first_arg="$1"
    local cmd="$PIXI_EXE $*"

    eval "$cmd"

    case "$first_arg" in
        add|a|remove|rm|install|i)
            eval "$($PIXI_EXE shell-hook --change-ps1 false)"
            hash -r
            ;;
    esac
}

export PS1="(pa_proj) ${PS1:-}"
exec "$@"
