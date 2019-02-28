#!/bin/bash

dir=`pwd`
cd $dir
root -l -b < x_run.C >& log.log
