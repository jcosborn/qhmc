# makefile for bc library for Lua

include ../make.inc

# change these to reflect your Lua installation
LUA= ../$(LUADIR)
LUAINC= $(LUA)/src
LUALIB= $(LUA)/src
LUABIN= $(LUA)/src

# these will probably work if Lua has been installed globally
#LUA= /usr/local
#LUAINC= $(LUA)/include
#LUALIB= $(LUA)/lib
#LUABIN= $(LUA)/bin

# probably no need to change anything below here
CC= gcc
CFLAGS= $(INCS) $(WARN) -O2 $G
WARN= -std=c99 -pedantic -Wall -Wextra
INCS= -I$(LUAINC) -I.
#MAKESO= $(CC) -shared
#MAKESO= $(CC) -bundle -undefined dynamic_lookup

MYNAME= bc
MYLIB= l$(MYNAME)
#T= $(MYNAME).so
OBJS= $(MYLIB).o number.o
TEST= test.lua

LBC = liblbc.a

all: $(LBC)

$(LBC): $(OBJS)
	ar cru $@ $^

#all:	test

#test:	$T
#	$(LUABIN)/lua $(TEST)

o:	$(MYLIB).o

#so:	$T

#$T:	$(OBJS)
#	$(MAKESO) -o $@ $(OBJS)

realclean: clean
clean:
	rm -f $(OBJS) $T core core.* $(LBC)

doc:
	@echo "$(MYNAME) library:"
	@fgrep '/**' $(MYLIB).c | cut -f2 -d/ | tr -d '*' | sort | column

# eof
