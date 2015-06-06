#!/usr/bin/env bash
#
git submodule init
git submodule update
git pull origin master
git submodule update --recursive
