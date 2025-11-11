# Makefile for extract_trajectory_parameters

SHELL = cmd
CC = gcc
CLINKER = $(CC)
CCFLAGS = -O2

.PHONY: all
all: extract_trajectory_parameters

.PHONY: extract_trajectory_parameters
extract_trajectory_parameters: extract_trajectory_parameters.exe

extract_trajectory_parameters.exe: extract_trajectory_parameters.obj
	@echo Linking $@ againgst $^ ...
	$(CLINKER) -o $@ $^ -static -s

extract_trajectory_parameters.obj: extract_trajectory_parameters.c
	@echo Compiling $@ ...
	$(CC) -o $@ -c $< $(CCFLAGS) -g0

.PHONY: veryclean
veryclean: clean
	-del /q extract_trajectory_parameters.exe 2> NUL

.PHONY: clean
clean:
	-del /q extract_trajectory_parameters.obj 2> NUL

