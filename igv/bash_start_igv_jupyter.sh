#!/bin/bash

jupyter-lab --no-browser &
python cors_server.py 8000 &
