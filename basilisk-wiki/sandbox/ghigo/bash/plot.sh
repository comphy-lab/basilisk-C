#!/bin/bash

function bplot()
{
    if test -f $1.c; then
	if test -f $1/plots; then
	    rm $1/plots
	fi
	make $1/plots
	display -update 1 $1/_plot$2.svg
    fi
}
