# Copyright (C) 2016-2020 Davidson Francis <davidsondfgl@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>

AR        = ar
INCLUDE   = $(CURDIR)/include
SRC       = $(CURDIR)/src
CFLAGS   +=  -O2
CFLAGS   +=  -I $(INCLUDE)
ARFLAGS   =  cr
LIB       =  libws.a
MCSS_DIR ?= /usr/bin/
MANPAGES  = $(CURDIR)/doc/man/man3
AFL_FUZZ ?= no

# Prefix
ifeq ($(PREFIX),)
    PREFIX := /usr/local
endif

# Detect machine type
MACHINE := $(shell uname -m)
ifeq ($(MACHINE), x86_64)
	LIBDIR = $(PREFIX)/lib64
else
	LIBDIR = $(PREFIX)/lib
endif

# Check if AFL fuzzing enabled
ifeq ($(AFL_FUZZ), yes)
    CC = afl-gcc
    CFLAGS := $(filter-out -O2,$(CFLAGS))
    CFLAGS += -DVERBOSE_MODE -DAFL_FUZZ -g
    $(info [+] AFL Fuzzing build enabled)
endif

# Source
C_SRC = $(wildcard $(SRC)/base64/*.c)    \
		$(wildcard $(SRC)/handshake/*.c) \
		$(wildcard $(SRC)/sha1/*.c)      \
		$(wildcard $(SRC)/*.c)

OBJ = $(C_SRC:.c=.o)

# Conflicts
.PHONY: doc
.PHONY: tests

# Paths
INCDIR = $(PREFIX)/include
BINDIR = $(PREFIX)/bin

# General objects
%.o: %.c
	$(CC99) $< $(CFLAGS) -c -o $@

# All
all: libws.a

# Library
libws.a: $(OBJ)
	$(AR) $(ARFLAGS) $(LIB) $^

# Clean
clean:
	@rm -f $(SRC)/base64/*.o
	@rm -f $(SRC)/handshake/*.o
	@rm -f $(SRC)/sha1/*.o
	@rm -f $(SRC)/*.o
	@rm -f $(LIB)

check:
