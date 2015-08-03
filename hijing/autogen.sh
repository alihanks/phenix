#!/bin/sh
#this is a generic file - you never need to change this one
srcdir=`dirname $0`
test -z "$srcdir" && srcdir=.

(cd $srcdir; aclocal -I ${OFFLINE_MAIN}/share;\
libtoolize --force; automake -a --add-missing; autoconf)

$srcdir/configure "$@"
