# Makefile for extract_trajectory_parameters

SHELL = cmd
CC = gcc
CLINKER = $(CC)
CCFLAGS = -O3

.PHONY: all
all: extract_trajectory_parameters

.PHONY: extract_trajectory_parameters
extract_trajectory_parameters: extract_trajectory_parameters.exe

extract_trajectory_parameters.exe: extract_trajectory_parameters.obj
	@echo Linking $@ againgst $^ ...
	$(CLINKER) -o $@ $^ -static

extract_trajectory_parameters.obj: extract_trajectory_parameters.c
	@echo Compiling $@ ...
	$(CC) -o $@ -c $< $(CCFLAGS)

.PHONY: clean
clean: clean_tmp
	-del /q extract_trajectory_parameters.exe 1> NUL 2> NUL

.PHONY: clean_tmp
clean_tmp:
	-del /q extract_trajectory_parameters.obj 1> NUL 2> NUL

